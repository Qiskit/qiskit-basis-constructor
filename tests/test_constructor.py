import cmath
import math
import unittest
import warnings

import ddt

from qiskit.circuit import EquivalenceLibrary, library as lib, QuantumCircuit, Parameter, Measure
from qiskit.circuit.classical import types
from qiskit.quantum_info import Operator
from qiskit.transpiler import Target, InstructionProperties
from qiskit.transpiler.passes import GatesInBasis

from qiskit_basis_constructor import (
    BasisConstructor,
    GateCount,
    LogFidelity,
    standard_equivalence_library,
)


@ddt.ddt
class TestConstructor(unittest.TestCase):
    def assertGatesInBasis(self, qc, target):  # noqa: N802 (`unittest` has different conventions)
        check = GatesInBasis(target=target)
        _ = check(qc)
        if not check.property_set["all_gates_in_basis"]:
            self.fail(f"Not all gates in basis for:\n{qc}")

    def test_noop(self):
        """Basically just an edge-case test for the 0q case."""
        qc = QuantumCircuit(global_phase=1.5)
        out = BasisConstructor(EquivalenceLibrary(), [GateCount()], Target(0))(qc)
        self.assertEqual(qc, out)

    def test_noop_bell(self):
        target = Target(2)
        target.add_instruction(lib.HGate(), {(0,): None, (1,): None})
        target.add_instruction(Measure(), {(0,): None, (1,): None})
        target.add_instruction(lib.CXGate(), {(0, 1): None})

        bell = QuantumCircuit(2, 2)
        bell.h(0)
        bell.cx(0, 1)
        bell.measure([0, 1], [0, 1])

        pass_ = BasisConstructor(EquivalenceLibrary(), target)
        out = pass_(bell)
        self.assertEqual(out, bell)

    @ddt.idata(
        (gate, target_direction, circuit_direction)
        for gate in ("cx", "cz", "ecr", "rzx", "rzz", "iswap")
        for target_direction in ((0, 1), (1, 0))
        for circuit_direction in ((0, 1), (1, 0))
    )
    @ddt.unpack
    def test_bell_conversion(self, gate_2q, target_direction, circuit_direction):
        """A simple Bell circuit should work, regardless of the direction of the available 2q gate
        in the target, and the direction that the 2q gate is used in the circuit."""
        target = Target(2)
        for gate_1q in (lib.RZGate(Parameter("a")), lib.SXGate()):
            target.add_instruction(gate_1q, {(0,): None, (1,): None})
        target.add_instruction(
            lib.get_standard_gate_name_mapping()[gate_2q], {target_direction: None}
        )

        bell = QuantumCircuit(2)
        bell.h(circuit_direction[0])
        bell.cx(*circuit_direction)

        pass_ = BasisConstructor(standard_equivalence_library(), [GateCount()], target)
        out = pass_(bell)
        self.assertGatesInBasis(out, target)
        self.assertEqual(Operator(out), Operator(bell))

    @ddt.idata(
        (gate, target_direction, circuit_direction)
        for gate in ("cx", "cz", "ecr", "rzx", "rzz", "iswap")
        for target_direction in ((2, 27), (27, 2))
        for circuit_direction in ((2, 27), (27, 2))
    )
    @ddt.unpack
    def test_bell_conversion_offset(self, gate_2q, target_direction, circuit_direction):
        """A simple Bell circuit should work, regardless of the direction of the available 2q gate
        in the target, and the direction that the 2q gate is used in the circuit.

        This also offsets the qubits from [0, 1] to make sure that homogenisation relabelling isn't
        creating spurious test passes."""
        target = Target(50)
        for instr_1q in (lib.RZGate(Parameter("a")), lib.SXGate(), Measure()):
            target.add_instruction(instr_1q, {(i,): None for i in target_direction})
        target.add_instruction(
            lib.get_standard_gate_name_mapping()[gate_2q], {target_direction: None}
        )

        bell = QuantumCircuit(50, 50)
        bell.h(circuit_direction[0])
        bell.cx(*circuit_direction)
        bell.measure(circuit_direction, circuit_direction)

        pass_ = BasisConstructor(standard_equivalence_library(), [GateCount()], target)
        out = pass_(bell)
        self.assertGatesInBasis(out, target)

        # Simple reduced-size simulator.
        qargs_map = {out.qubits[qarg]: i for i, qarg in enumerate(circuit_direction)}
        state = Operator([[1, 0j, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        state *= cmath.exp(1j * out.global_phase)
        for instruction in out.data:
            if instruction.name == "measure":
                continue
            state = state.compose(
                instruction.operation.to_matrix(), qargs=[qargs_map[q] for q in instruction.qubits]
            )
        bell_op = Operator([[1, 1, 0j, 0], [0, 0, 1, -1], [0, 0, 1, 1], [1, -1, 0, 0]])
        bell_op /= math.sqrt(2)

        self.assertEqual(state, bell_op)

    def test_bell_heterogeneous_target(self):
        expensive = InstructionProperties(error=1e-1)
        # The difference between "cheap" and "cheaper" here is just to make it reliable that certain
        # symmetric decompositions will get picked.
        cheap = InstructionProperties(error=2e-3)
        cheaper = InstructionProperties(error=1e-3)

        # The weights in this target are nonsensical, but the point is to test that specific
        # decompositions are being used and that they differ between different qubits.
        #
        # `h` is cheap on all qubits except qubit 1, so it should appear verbatim as a rule on
        # those.  On qubit 1, we make `s` and `sx` much cheaper, so that should be preferred
        # instead.
        target = Target(4)
        target.add_instruction(
            lib.HGate(), {(0,): cheap, (1,): expensive, (2,): cheap, (3,): cheap}
        )
        target.add_instruction(lib.SGate(), {(1,): cheaper, (2,): cheaper})
        target.add_instruction(lib.SdgGate(), {(1,): cheap, (2,): cheap})
        target.add_instruction(lib.SXGate(), {(1,): cheaper, (2,): cheaper})
        target.add_instruction(lib.SXdgGate(), {(1,): cheap, (2,): cheap})
        target.add_instruction(lib.XGate(), {(1,): cheaper, (2,): cheaper})
        target.add_instruction(lib.CXGate(), {(0, 1): cheap, (2, 3): expensive})
        target.add_instruction(lib.ECRGate(), {(1, 2): cheap})
        target.add_instruction(lib.CZGate(), {(1, 2): expensive, (3, 2): cheap})

        multi_bell = QuantumCircuit(4)
        for base in range(3):
            multi_bell.h(base)
            multi_bell.cx(base, base + 1)
            multi_bell.barrier()

        # We're building this circuit:
        #      ┌───┐      ░                                         ░               ░
        # q_0: ┤ H ├──■───░─────────────────────────────────────────░───────────────░─
        #      └───┘┌─┴─┐ ░  ┌───┐  ┌────┐┌───┐┌─────┐┌──────┐┌───┐ ░               ░
        # q_1: ─────┤ X ├─░──┤ S ├──┤ √X ├┤ S ├┤ Sdg ├┤0     ├┤ X ├─░───────────────░─
        #           └───┘ ░ ┌┴───┴─┐└────┘└───┘└─────┘│  Ecr │└───┘ ░ ┌───┐         ░
        # q_2: ───────────░─┤ √Xdg ├──────────────────┤1     ├──────░─┤ H ├─■───────░─
        #                 ░ └──────┘                  └──────┘      ░ ├───┤ │ ┌───┐ ░
        # q_3: ───────────░─────────────────────────────────────────░─┤ H ├─■─┤ H ├─░─
        #                 ░                                         ░ └───┘   └───┘ ░
        expected = QuantumCircuit(4)
        # Both `h(0)` and `cx(0, 1)` are already cheap and available.
        expected.h(0)
        expected.cx(0, 1)
        expected.barrier()
        # On qubits (1, 2), `ecr` is the only cheap 2q gate.  `h(1)` is expensive, while `s` and
        # `sx` are very cheap, so we expect this decomposition (plus the cx-from-ecr bit).
        expected.s(1)
        expected.sx(1)
        expected.s(1)
        expected.sdg(1)
        expected.sxdg(2)
        expected.ecr(1, 2)
        expected.x(1)
        expected.barrier()
        # Finally, `cz(3, 2)` is cheap, and the equivalences should happily switch it for free into
        # a `cz(2, 3)` (to make `cx(2, 3)`).  `h(2)` and `h(3)` are both cheap.
        expected.h(2)
        expected.h(3)
        expected.cz(3, 2)
        expected.h(3)
        expected.barrier()

        pass_ = BasisConstructor(
            standard_equivalence_library(), [LogFidelity(math.log(10)), GateCount()], target
        )
        self.assertEqual(pass_(multi_bell), expected)

    def test_parametric_circuit(self):
        a, b = Parameter("a"), Parameter("b")
        qc = QuantumCircuit(3)
        qc.rx(a, 0)
        qc.rzx(b, 0, 1)
        qc.ry(a + b, 1)
        qc.cry(a, 1, 2)

        target = Target(3)
        target.add_instruction(lib.RZGate(Parameter("z")), {(i,): None for i in range(3)})
        target.add_instruction(lib.RXGate(Parameter("x")), {(i,): None for i in range(3)})
        target.add_instruction(lib.RXXGate(Parameter("xx")), {(0, 1): None, (1, 2): None})

        out = BasisConstructor(standard_equivalence_library(), [GateCount(2), GateCount()], target)(
            qc
        )
        self.assertGatesInBasis(out, target)

        for parameters in [(0.25, 0.75), (-0.5, 0.25)]:
            self.assertEqual(
                Operator(out.assign_parameters(parameters)),
                Operator(qc.assign_parameters(parameters)),
            )

    def test_swap_to_two_2q(self):
        qc = QuantumCircuit(4)
        qc.swap(0, 1)
        qc.swap(1, 2)
        qc.swap(2, 3)

        # We don't care about the 1q gates, but the 2q links both have two available 2q gates, from
        # different equivalence classes that should let us effect the swap in only two 2q gates.
        target = Target(4)
        target.add_instruction(lib.RZGate(Parameter("a")), {(i,): None for i in range(4)})
        target.add_instruction(lib.RXGate(Parameter("a")), {(i,): None for i in range(4)})
        target.add_instruction(lib.CXGate(), {(0, 1): None, (2, 1): None})
        target.add_instruction(lib.CZGate(), {(2, 3): None})
        target.add_instruction(lib.DCXGate(), {(1, 2): None, (2, 3): None})
        target.add_instruction(lib.iSwapGate(), {(0, 1): None})

        out = BasisConstructor(standard_equivalence_library(), [GateCount(2), GateCount()], target)(
            qc
        )
        self.assertGatesInBasis(out, target)
        self.assertEqual(Operator(out), Operator(qc))
        num_2q_gates = len([gate for gate in out.data if len(gate.qubits) == 2])
        self.assertEqual(num_2q_gates, 6)  # Three swaps, with each swap taking two 2q gates.

    def test_conditional_operation(self):
        qc = QuantumCircuit(2, 2)
        expected = qc.copy_empty_like()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning, message=".*c_if.*")
            qc.cz(0, 1).c_if(qc.cregs[0], False)
            expected.h(1).c_if(expected.cregs[0], False)
            expected.cx(0, 1).c_if(expected.cregs[0], False)
            expected.h(1).c_if(expected.cregs[0], False)

        target = Target(2)
        target.add_instruction(lib.HGate(), {(0,): None, (1,): None})
        target.add_instruction(lib.CXGate(), {(0, 1): None})
        pass_ = BasisConstructor(standard_equivalence_library(), GateCount(), target)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning, message=".*c_if.*")
            out = pass_(qc)
        self.assertEqual(out, expected)

    def test_fidelity_binning(self):
        # The idea of this Target is that s.s is preferred on (0,) and z is preferred on (1,), but
        # if the fidelity binning is set to combine them, then z will be preferred on both.
        target = Target(2)
        target.add_instruction(
            lib.SGate(),
            {
                (0,): InstructionProperties(error=1e-3),
                (1,): InstructionProperties(error=6e-3),
            },
        )
        target.add_instruction(
            lib.ZGate(),
            {
                (0,): InstructionProperties(error=3e-3),
                (1,): InstructionProperties(error=6e-3),
            },
        )
        qc = QuantumCircuit(2)
        qc.z(0)
        qc.z(1)

        wide_expected = QuantumCircuit(2)
        wide_expected.z(0)
        wide_expected.z(1)

        narrow_expected = QuantumCircuit(2)
        narrow_expected.s(0)
        narrow_expected.s(0)
        narrow_expected.z(1)

        equivalences = standard_equivalence_library()
        wide_pass = BasisConstructor(equivalences, [LogFidelity(4 * math.log(10))], target)
        narrow_pass = BasisConstructor(equivalences, [LogFidelity(math.log(10) / 4)], target)

        self.assertEqual(wide_pass(qc), wide_expected)
        self.assertEqual(narrow_pass(qc), narrow_expected)

    def test_control_flow(self):
        target = Target(3)
        target.add_instruction(lib.HGate(), {(1,): None, (2,): None})
        target.add_instruction(lib.SGate(), {(1,): None})
        target.add_instruction(lib.XGate(), {(1,): None})
        target.add_instruction(lib.SXGate(), {(2,): None})
        target.add_instruction(lib.CZGate(), {(0, 1): None})
        target.add_instruction(lib.CXGate(), {(1, 2): None})

        qc = QuantumCircuit(3)
        a = qc.add_input("a", types.Uint(8))
        b = qc.add_var("b", False)
        with qc.if_test(b):
            qc.cx(0, 1)
        with qc.if_test(b) as else_:
            pass
        with else_:
            qc.cx(0, 1)
        with qc.while_loop(b):
            with qc.switch(a) as case:
                with case(0):
                    qc.cz(1, 2)
                with case(1):
                    qc.ecr(1, 2)

        expected = QuantumCircuit(3)
        expected.add_input(a)
        expected.add_var(b, False)
        with expected.if_test(b):
            expected.h(1)
            expected.cz(0, 1)
            expected.h(1)
        with expected.if_test(b) as else_:
            pass
        with else_:
            expected.h(1)
            expected.cz(0, 1)
            expected.h(1)
        with expected.while_loop(b):
            with expected.switch(a) as case:
                with case(0):
                    expected.h(2)
                    expected.cx(1, 2)
                    expected.h(2)
                with case(1):
                    expected.global_phase += -math.pi / 4
                    expected.s(1)
                    expected.sx(2)
                    expected.cx(1, 2)
                    expected.x(1)

        pass_ = BasisConstructor(standard_equivalence_library(), [GateCount()], target)
        self.assertEqual(pass_(qc), expected)
