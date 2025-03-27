# This code is part of Qiskit.
#
# (C) Copyright IBM 2025.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from __future__ import annotations

__all__ = [
    "BasisConstructor",
    "BasisConstructorError",
    "GateCount",
    "LogFidelity",
]

import abc
import collections
import dataclasses
import heapq
import itertools
import logging
import math
import statistics
import typing

import rustworkx

from qiskit.circuit import ParameterVector, Parameter, QuantumCircuit, Gate
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGOpNode
from qiskit.transpiler import TransformationPass, TranspilerError

if typing.TYPE_CHECKING:
    from qiskit.circuit import EquivalenceLibrary, Instruction
    from qiskit.dagcircuit import DAGCircuit
    from qiskit.transpiler import Target, InstructionProperties

_LOGGER = logging.getLogger(__name__)


class BasisConstructorError(TranspilerError):
    pass


_PARAMETER_CACHE = {}


def _get_parameter_vector(num_params: int) -> ParameterVector:
    if (vector := _PARAMETER_CACHE.get(num_params)) is None:
        vector = _PARAMETER_CACHE[num_params] = ParameterVector(
            "_basis_constructor_internal_", num_params
        )
    return vector


class BasisConstructor(TransformationPass):
    def __init__(
        self,
        equivalences: EquivalenceLibrary,
        score: _ScoreComponent | typing.Iterable[_ScoreComponent],
        target: Target | None = None,
    ):
        """
        Args:
            equivalences: the equivalence library to use for translations.  A good choice is
                :data:`.DEFAULT_EQUIVALENCE_LIBRARY` (which is what the plugin form of this pass
                will use).
            score: a sequence of scoring components to use.  See :class:`.LogFidelity` and
                :class:`.GateCount` as the concrete options.
            target: the backend ``Target``.  Passing ``None`` suppresses translations entirely.
        """
        super().__init__()
        self._original_equivalences = equivalences
        self._equivalences_by_arity = None
        self._original_target = target
        self._target = None
        self._score = (score,) if isinstance(score, _ScoreComponent) else tuple(score)
        self._constructed = None

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        if self._original_target is None:
            return dag
        if self._target is None:
            infidelity_bin = next(
                (component.bin for component in self._score if isinstance(component, LogFidelity)),
                None,
            )
            self._target = _HomogenizedTarget.from_target(self._original_target, infidelity_bin)
        if self._equivalences_by_arity is None:
            self._equivalences_by_arity = _split_equivalences(self._original_equivalences)
        if self._constructed is None:
            self._constructed = {}
            for node_index in rustworkx.topological_sort(self._target.qargs_graph):
                self._constructed[node_index] = self._construct(node_index)
        return self._rewrite(dag, {q: i for i, q in enumerate(dag.qubits)})

    def _rewrite(self, dag, qubit_map):
        out = dag.copy_empty_like()
        for node in dag.topological_op_nodes():
            if node.is_control_flow():
                new_blocks = [
                    dag_to_circuit(
                        self._rewrite(
                            circuit_to_dag(block, copy_operations=False),
                            {
                                inner: qubit_map[outer]
                                for inner, outer in zip(block.qubits, node.qargs)
                            },
                        ),
                        copy_operations=False,
                    )
                    for block in node.op.blocks
                ]
                out.apply_operation_back(node.op.replace_blocks(new_blocks), node.qargs, node.cargs)
                continue
            if node.is_directive():
                out._apply_op_node_back(node)  # noqa: SLF001 (semi-public fast path)
                continue
            replacement = self._replacement_for(
                node.name, node.params, [qubit_map[qubit] for qubit in node.qargs]
            )
            # This is effectively `compose`, but restricted to the set of operations our
            # replacements will actually use.
            out.global_phase += replacement.global_phase
            for instruction in replacement.data:
                if getattr(node, "condition", None) is not None:
                    instruction = instruction.replace(
                        operation=instruction.operation.copy().c_if(*node.condition)
                    )
                instruction = instruction.replace(
                    qubits=tuple(
                        node.qargs[replacement.find_bit(qubit).index]
                        for qubit in instruction.qubits
                    ),
                    clbits=tuple(
                        node.cargs[replacement.find_bit(clbit).index]
                        for clbit in instruction.clbits
                    ),
                )
                out._apply_op_node_back(  # noqa: SLF001 (semi-public fast path)
                    DAGOpNode.from_instruction(instruction)
                )
        return out

    def _replacement_for(self, name, params, qargs) -> QuantumCircuit:
        homogenized = self._target.qargs_order(qargs) or self._construct_global(tuple(qargs))
        # TODO: improve failure messages - we could jump to a failure path where we find all the
        # failures and produce a complete accounting of everything that went wrong.
        if (
            homogenized is None
            or (sub_constructed := self._constructed.get(homogenized.index)) is None
        ):
            raise BasisConstructorError(f"no gates constructed for qubits {qargs}")
        if (replacement := sub_constructed.get((name, homogenized.qargs))) is None:
            raise BasisConstructorError(f"gate '{name}' not constructed for qubits {qargs}")
        if replacement.concrete_decomposition is None:
            concrete = replacement.abstract_decomposition.copy_empty_like()
            for instruction in replacement.abstract_decomposition.data:
                sub_qargs = [
                    qargs[
                        replacement.qargs[replacement.abstract_decomposition.find_bit(qubit).index]
                    ]
                    for qubit in instruction.qubits
                ]
                concrete.compose(
                    self._replacement_for(instruction.name, instruction.params, sub_qargs),
                    [
                        replacement.qargs[replacement.abstract_decomposition.find_bit(qubit).index]
                        for qubit in instruction.qubits
                    ],
                    instruction.clbits,
                    inplace=True,
                )
            _LOGGER.debug("concretizing '%s%s' as:\n%s", name, qargs, concrete)
            replacement.concrete_decomposition = concrete
        return replacement.concrete_decomposition.assign_parameters(
            params, inplace=False, strict=True, flat_input=True
        )

    def _construct_global(self, qargs: tuple[int, ...]) -> _HomogenizedQargs:
        node_index = self._target.add_explicit_node(qargs, [])
        self._constructed[node_index] = self._construct(node_index)
        return self._target.qargs_order(qargs)

    def _construct(self, node_index):
        entry = self._target.qargs_graph[node_index]
        base_circuit = QuantumCircuit(entry.num_qubits)

        def extract_sub_qargs(node_indices, cur_node, outer_from_current):
            for parent, _, parent_from_child in self._target.qargs_graph.in_edges(cur_node):
                sub_outer_from_current, sub_current_from_outer = (
                    {
                        parent: outer_from_current[child]
                        for child, parent in parent_from_child.items()
                    },
                    {
                        outer_from_current[child]: parent
                        for child, parent in parent_from_child.items()
                    },
                )
                node_indices[frozenset(sub_current_from_outer)] = (
                    parent,
                    sub_current_from_outer,
                )
                extract_sub_qargs(node_indices, parent, sub_outer_from_current)
            return node_indices

        node_index_from_homogenized_qargs = extract_sub_qargs(
            {}, node_index, {i: i for i in range(entry.num_qubits)}
        )

        def rule_satisfied(name, homogenized_qargs):
            if len(homogenized_qargs) == entry.num_qubits:
                return False
            if (
                lookup := node_index_from_homogenized_qargs.get(frozenset(homogenized_qargs))
            ) is None:
                return False
            parent, qarg_lookup = lookup
            local_qargs = tuple(qarg_lookup[q] for q in homogenized_qargs)
            return (name, local_qargs) in self._constructed.get(parent, {})

        def score_equivalence(equivalence, current_constructed):
            total_score = [0] * len(self._score)
            for instruction in equivalence.decomposition.data:
                num_qubits = len(instruction.qubits)
                qargs = tuple(
                    equivalence.decomposition.find_bit(q).index for q in instruction.qubits
                )
                if num_qubits == entry.num_qubits:
                    constructed = current_constructed
                    local_qargs = qargs
                else:
                    parent, qarg_lookup = node_index_from_homogenized_qargs[frozenset(qargs)]
                    constructed = self._constructed[parent]
                    local_qargs = tuple(qarg_lookup[q] for q in qargs)
                total_score = [
                    a + b
                    for a, b in zip(total_score, constructed[instruction.name, local_qargs].score)
                ]
            return _CandidateRule(
                tuple(total_score),
                equivalence.name,
                equivalence.qargs,
                equivalence.decomposition,
            )

        candidates = []
        for instruction in entry.instructions:
            qc = QuantumCircuit(instruction.example.num_qubits, instruction.example.num_clbits)
            qc.append(instruction.example, qc.qubits, qc.clbits)
            qc.assign_parameters(
                dict(
                    zip(
                        instruction.example.params,
                        _get_parameter_vector(len(instruction.example.params)),
                    )
                ),
                inplace=True,
                flat_input=True,
                strict=False,
            )
            heapq.heappush(
                candidates,
                _CandidateRule(
                    score=tuple(component.score_target(instruction) for component in self._score),
                    name=instruction.name,
                    qargs=instruction.qargs,
                    abstract_decomposition=qc,
                    # Rules that come from the `Target` are automatically concrete with themselves.
                    # Setting this forms the base case of the recursion in reconstruction passes.
                    concrete_decomposition=qc,
                ),
            )
        for instruction in self._target.instructions_global.get(entry.num_qubits, []):
            qc = QuantumCircuit(instruction.example.num_qubits, instruction.example.num_clbits)
            qc.append(instruction.example, qc.qubits, qc.clbits)
            qc.assign_parameters(
                dict(
                    zip(
                        instruction.example.params,
                        _get_parameter_vector(len(instruction.example.params)),
                    )
                ),
                inplace=True,
                flat_input=True,
                strict=False,
            )
            for permutation in itertools.permutations(range(entry.num_qubits)):
                heapq.heappush(
                    candidates,
                    _CandidateRule(
                        score=tuple(
                            component.score_target(instruction) for component in self._score
                        ),
                        name=instruction.name,
                        qargs=permutation,
                        abstract_decomposition=qc,
                        # Rules that come from the `Target` are automatically concrete with
                        # themselves.  Setting this forms the base case of the recursion in
                        # reconstruction passes.
                        concrete_decomposition=qc,
                    ),
                )

        constructed = {}
        used_by = collections.defaultdict(list)
        available_rules = []
        for equivalences in self._equivalences_by_arity[entry.num_qubits].values():
            for equivalence in equivalences:
                for permutation in itertools.permutations(range(entry.num_qubits)):
                    this_id = len(available_rules)
                    permuted_equivalence = dataclasses.replace(
                        dataclasses.replace(equivalence, qargs=permutation),
                        decomposition=base_circuit.compose(equivalence.decomposition, permutation),
                    )
                    rules_needed = 0
                    for (name, _), sub_qargs in permuted_equivalence.needs():
                        if rule_satisfied(name, sub_qargs):
                            continue
                        rules_needed += 1
                        used_by[name, sub_qargs].append(this_id)
                    if rules_needed == 0:
                        heapq.heappush(
                            candidates,
                            score_equivalence(permuted_equivalence, constructed),
                        )
                    else:
                        available_rules.append(_AvailableRule(permuted_equivalence, rules_needed))

        while candidates:
            best_rule = heapq.heappop(candidates)
            key = best_rule.key()
            if (previous := constructed.get(key)) is not None:
                _LOGGER.debug(
                    "skipping previously constructed %r as %s>=%s",
                    key,
                    best_rule.score,
                    previous.score,
                )
                continue
            _LOGGER.debug("constructing %r with score %s", key, best_rule.score)
            constructed[key] = best_rule
            for dependent in used_by[key]:
                rule = available_rules[dependent]
                rule.remaining -= 1
                if rule.remaining == 0 and rule.equivalence.key() not in constructed:
                    heapq.heappush(candidates, score_equivalence(rule.equivalence, constructed))
        return constructed


class _ScoreComponent(abc.ABC):
    __slots__ = ()

    @abc.abstractmethod
    def score_target(self, homogenized_instruction: _HomogenizedInstruction):
        raise NotImplementedError


class GateCount(_ScoreComponent):
    __slots__ = ("min_qubits",)

    def __init__(self, min_qubits=0):
        self.min_qubits = min_qubits

    def score_target(self, homogenized_instruction):
        return int(len(homogenized_instruction.qargs) >= self.min_qubits)


class LogFidelity(_ScoreComponent):
    __slots__ = ("bin",)

    def __init__(self, bin=0.0):
        if math.isnan(bin) or bin < 0.0:
            raise BasisConstructorError(f"bin must be non-negative, but was {bin}")
        self.bin = bin

    def score_target(self, homogenized_instruction):
        return homogenized_instruction.neg_log_fidelity


@dataclasses.dataclass
class _AvailableRule:
    equivalence: _EquivalenceRule
    remaining: int


@dataclasses.dataclass(order=True, eq=True)
class _CandidateRule:
    score: tuple[object, ...]
    name: str
    qargs: tuple[int, ...]
    abstract_decomposition: QuantumCircuit = dataclasses.field(compare=False)
    concrete_decomposition: QuantumCircuit | None = dataclasses.field(default=None, compare=False)

    def key(self):
        """Construction key of the rule.

        Several rules can have the same key, but only one will be chosen as the best."""
        return (self.name, self.qargs)


@dataclasses.dataclass
class _GraphEntry:
    num_qubits: int
    qargs: list[tuple[int, ...]]
    instructions: list[_HomogenizedInstruction]


@dataclasses.dataclass(frozen=True, eq=True)
class _HomogenizedInstruction:
    name: str
    qargs: tuple[int, ...]
    neg_log_fidelity: float
    example: Instruction = dataclasses.field(compare=False)


@dataclasses.dataclass(frozen=True, eq=True)
class _HomogenizedQargs:
    index: int
    qargs: tuple[int, ...]


@dataclasses.dataclass(frozen=True, eq=True)
class _QargsOrder:
    from_homogenized: tuple[int, ...]
    from_physical: dict[int, int]


class _HomogenizedTarget:
    # TODO: the homogenisation of the target is probably more broadly applicable than just here.
    # For example, the `BasisTranslator` effectively wants this same information with
    # `bin_width=None`, though specifics like the `qargs_graph` will differ.

    qargs_graph: rustworkx.PyDiGraph[_GraphEntry, None]
    """The graph stores all the unique instruction sets for qargs in its nodes, with the
    dependencies bewteen construction bases represented by directed edges.  As an example, qargs
    `(0, 1)` would depend on `(0,)` and `(1,)` having been constructed first."""
    node_index_from_qargs: dict[frozenset[int], int]
    natural_qargs_order: dict[frozenset[int], _QargsOrder]
    instructions_global: collections.defaultdict[list[_HomogenizedInstruction]]
    node_index_from_key: dict[object, int]
    """Mapping from each instruction-set key to the graph node index storing it."""

    @classmethod
    def empty(cls) -> typing.Self:
        self = cls()
        self.qargs_graph = rustworkx.PyDiGraph()
        self.natural_qargs_order = {}
        self.node_index_from_key = {}
        self.node_index_from_qargs = {}
        self.instructions_global = collections.defaultdict(list)
        return self

    @classmethod
    def from_target(cls, target: Target, bin_width: float | None) -> typing.Self:
        self = cls.empty()

        # Strictly there's no _need_ to bin based on the negative log fidelity itself - we could bin
        # on the error and avoid several unnecessary exponentials and logarithms - but it's pleasant
        # to have `_HomogenizedInstruction` always have its logical score field represent the same
        # value, and the cost of this doesn't seem to show up too much.
        bin_neg_log_fidelity = _neg_log_fidelity_bin(bin_width)

        instructions_per_arity_per_qargs = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )
        for name, properties in target.items():
            for qargs, instruction_properties in properties.items():
                # While we're looking through everything in the `Target`, we store the candidate
                # homogenised instructions with the binned negative log fidelity.  Once we've
                # completed the homogenisation, we can go back through and do the averaging.

                example = target.operation_from_name(name)
                # If the `example` isn't a `Gate`, then it's not going to be involved in any rule
                # rewriting, so we can just avoid binning on its fidelity.  This is effectively a
                # performance optimisation for the homogenisation, since `Measure`'s read-out error
                # is likely to be quite variable, and we don't want to introduce spurious
                # heterogeneity in the restricted target for fidelity differences that won't affect
                # any equivalence constructions.
                if isinstance(example, Gate):
                    neg_log_fidelity = bin_neg_log_fidelity(instruction_properties)
                else:
                    neg_log_fidelity = 0.0

                if qargs is None:
                    if not isinstance(num_qubits := getattr(example, "num_qubits", None), int):
                        # The `isinstance` check is to handle the case that we get a gate class
                        # where the `num_qubits` field is a property.
                        _LOGGER.info("ignoring %s due to unknown qubit count", example)
                        continue
                    self.instructions_global[num_qubits].append(
                        _HomogenizedInstruction(
                            name,
                            tuple(range(num_qubits)),
                            neg_log_fidelity,
                            example,
                        )
                    )
                else:
                    # The `BasisConstructor` algorithm handles gate direction by constructing the
                    # available bases for all permutations of multi-qubit qargs simultaneously
                    # (hence the `frozenset` key), but we still also need to know which concrete
                    # qargs we're working with.
                    instructions_per_arity_per_qargs[len(qargs)][frozenset(qargs)].append(
                        _HomogenizedInstruction(
                            name,
                            qargs,
                            neg_log_fidelity,
                            example,
                        )
                    )

        # Now we group all the property sets into "alike" ones.  Two groups of qubits are "alike" if
        # their set of interactions and associated fidelities (after rounding) are isomorphic, up to
        # renumbering of qubits.  We ensure that the zero case is always present to be the root
        # node, to simplify later logic.
        if 0 not in instructions_per_arity_per_qargs:
            _ = instructions_per_arity_per_qargs[0][frozenset()]
        for _, instructions_per_qargs in sorted(instructions_per_arity_per_qargs.items()):
            for qargs, instructions in instructions_per_qargs.items():
                self.add_explicit_node(qargs, instructions)

        # Now that the graph is fully constructed, we can go back and replace the negative log
        # fidelities derived from the binned errors with ones derived from the geometric mean of the
        # actual errors for all results that binned to that value.

        for entry in self.qargs_graph.nodes():

            def average(instruction):
                # Normal treatments of the geometric mean simply fail when encountering a zero.  We
                # can safely have it average out to a zero error; that's still a meaningful result.
                error = math.exp(
                    statistics.fmean(
                        _clipped_log(
                            _clip_error(
                                target[instruction.name][tuple(qargs[i] for i in instruction.qargs)]
                            )
                        )
                        for qargs in entry.qargs
                    )
                )
                if error <= 0.0:
                    return 0.0
                if error >= 1.0:
                    return math.inf
                return -math.log1p(-error)

            entry.instructions = [
                dataclasses.replace(instruction, neg_log_fidelity=average(instruction))
                for instruction in entry.instructions
            ]
        return self

    def add_explicit_node(self, qargs: tuple[int, ...], local_instructions: list) -> int:
        if not qargs:
            # We handle this specially so that the later loop-based logic can assume that the loop
            # always runs at least once.
            this_key = (self._sub_qargs_key(()), self._instruction_key(local_instructions, []))
            node_index = self.node_index_from_key[this_key] = self.qargs_graph.add_node(
                _GraphEntry(0, [()], local_instructions)
            )
            self.node_index_from_qargs[frozenset()] = node_index
            self.natural_qargs_order[frozenset()] = _QargsOrder((), {})
            return node_index

        # Obviously the scaling on this is absolute trash for multi-qubit instructions, but
        # since it's just 1 and 2 attempts for 1q and 2q respectively, and we don't really
        # deal with more than that, we can deal with it in the future.
        #
        # The loop is always entered at least once, so the `else` block can't `NameError`.
        qargs = tuple(sorted(qargs))
        for permutation in itertools.permutations(qargs):
            normalized_instructions = self._instruction_key(local_instructions, permutation)
            this_key = (self._sub_qargs_key(permutation), normalized_instructions)
            if (node_index := self.node_index_from_key.get(this_key)) is not None:
                self.qargs_graph[node_index].qargs.append(permutation)
                break
        else:
            # Pragmatically, humans like to define gates in terms of sorted indices where
            # possible, so this tends to make the resulting homogenised target easier to
            # read.  It also coincides with the first key attempted, which makes successful
            # dictionary lookups likely to be faster in the loop.
            permutation = qargs

            normalized_instructions = self._instruction_key(local_instructions, permutation)
            this_key = (self._sub_qargs_key(permutation), normalized_instructions)
            node_index = self.node_index_from_key[this_key] = self.qargs_graph.add_node(
                _GraphEntry(len(permutation), [permutation], normalized_instructions)
            )

            def add_predecessor_edges_to(child_index, sub_qargs_pairs):
                parent_qargs_key = frozenset(physical for physical, _ in sub_qargs_pairs)
                parent_index = self.node_index_from_qargs.get(parent_qargs_key)
                if parent_index is None:
                    for sub_qargs_pairs in itertools.combinations(
                        sub_qargs_pairs, len(sub_qargs_pairs) - 1
                    ):
                        add_predecessor_edges_to(child_index, sub_qargs_pairs)
                else:
                    parent_order = self.natural_qargs_order[parent_qargs_key]
                    parent_from_child = {
                        child: parent_order.from_physical[physical]
                        for physical, child in sub_qargs_pairs
                    }
                    self.qargs_graph.add_edge(parent_index, child_index, parent_from_child)

            # This time we include the 0 qargs case.
            physical_child = [(physical, child) for child, physical in enumerate(permutation)]
            for sub_qargs_pairs in itertools.combinations(physical_child, len(permutation) - 1):
                add_predecessor_edges_to(node_index, sub_qargs_pairs)

        self.node_index_from_qargs[frozenset(qargs)] = node_index
        self.natural_qargs_order[frozenset(qargs)] = _QargsOrder(
            permutation, {qubit: index for index, qubit in enumerate(permutation)}
        )
        return node_index

    def qargs_order(self, qargs: tuple[int, ...]) -> _HomogenizedQargs:
        lookup = frozenset(qargs)
        if (index := self.node_index_from_qargs.get(lookup, None)) is None:
            return None
        natural_order = self.natural_qargs_order[lookup]
        return _HomogenizedQargs(
            index, tuple(natural_order.from_physical[qubit] for qubit in qargs)
        )

    @staticmethod
    def _instruction_key(instructions, reverse_mapping: list[int]):
        forwards_mapping = {q: i for i, q in enumerate(reverse_mapping)}
        return frozenset(
            dataclasses.replace(entry, qargs=tuple(forwards_mapping[q] for q in entry.qargs))
            for entry in instructions
        )

    def _sub_qargs_key(self, permutation):
        def index_from_qargs(sub_qargs):
            if (node_index := self.node_index_from_qargs.get(frozenset(sub_qargs), None)) is None:
                # This is a subset of the qargs under consideration, and since we construct in
                # increasing instruction arity, this must refer to a node that can only have global
                # operations defined.
                return self.add_explicit_node(tuple(sub_qargs), [])
            return node_index

        # This part of the key always has the same elements, it's the order of them we
        # care about.  We assume that the zero-qargs case is always already present before this
        # function runs for the first time.
        return tuple(
            # TODO: the scaling here is needlessly awful for multi-q gates - there
            # should be a better way of building the predecessors into the key.
            index_from_qargs(sub_qargs)
            for k in range(len(permutation))
            for sub_qargs in itertools.combinations(permutation, k)
        )


def _neg_log_fidelity_bin(
    bin_width: float | None,
) -> typing.Callable[[InstructionProperties], float]:
    if bin_width is None:
        return lambda properties: 0.0

    def bin_neg_log_fidelity(properties):
        error = _clip_error(properties)
        if error == 0.0:
            return 0.0
        if error == 1.0:
            return math.inf
        if bin_width == 0.0:
            binned_error = error
        elif math.isinf(bin_width):
            binned_error = 0.0
        else:
            binned_error = math.exp(round(math.log(error) / bin_width) * bin_width)
        return -math.log1p(-binned_error) if binned_error < 1.0 else math.inf

    return bin_neg_log_fidelity


def _clip_error(properties: InstructionProperties | None):
    error = getattr(properties, "error", None) or 0.0
    if error <= 0.0:
        return 0.0
    if error >= 1.0:
        return 1.0
    return error


def _clipped_log(x):
    if x == 0.0:
        return -math.inf
    return math.log(x)


@dataclasses.dataclass
class _EquivalenceRule:
    name: str
    qargs: tuple[int, ...]
    decomposition: QuantumCircuit

    def key(self):
        """Construction key of the rule.

        Several rules can have the same key, but only one will be chosen as the best."""
        return (self.name, self.qargs)

    def needs(self):
        return {
            (
                (instruction.name, len(instruction.qubits)),
                tuple(self.decomposition.find_bit(q).index for q in instruction.qubits),
            )
            for instruction in self.decomposition.data
        }


def _split_equivalences(equivalences: EquivalenceLibrary):
    by_arity = collections.defaultdict(dict)

    def normalize_equivalence(key, equivalence) -> QuantumCircuit | None:
        if not all(isinstance(param, Parameter) for param in equivalence.params):
            # We rely on being able to bind all gate paramters, since we traverse the
            # `EquivalenceLibrary` differently to `BasisTranslator`, so for a first pass, it's
            # easier to just reject partially parametric rules.
            _LOGGER.debug(
                "skipping equivalence rule for %dq '%s' that is not fully parametric",
                key.num_qubits,
                key.name,
            )
            return None
        return equivalence.circuit.assign_parameters(
            dict(zip(equivalence.params, _get_parameter_vector(len(equivalence.params)))),
            inplace=False,
            strict=False,
            flat_input=True,
        )

    for key in equivalences.keys():
        by_arity[key.num_qubits][key.name] = [
            _EquivalenceRule(key.name, tuple(range(key.num_qubits)), circuit)
            # Private method so we can more easily access the base circuit, without needing to have
            # a fully paramteric gate object in hand, which isn't always easy to construct.
            for equivalence in equivalences._get_equivalences(key)  # noqa: SLF001
            if (circuit := normalize_equivalence(key, equivalence)) is not None
        ]
    return by_arity
