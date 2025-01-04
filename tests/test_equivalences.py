import unittest

import ddt

from qiskit.circuit.library import get_standard_gate_name_mapping
from qiskit.quantum_info import Operator
from qiskit_basis_constructor import standard_equivalence_library

EQUIVALENCES = standard_equivalence_library()


@ddt.ddt
class TestEquivalenceLibrary(unittest.TestCase):
    @ddt.idata(EQUIVALENCES.keys())
    def test_equivalence_validity(self, key):
        # All the keys in our default library should be valid, standard Qiskit gates, so we can
        # reconstruct them from parameters from their base class.
        gate_base = get_standard_gate_name_mapping()[key.name]
        params = [0.5**i for i in range(len(gate_base.params))]
        gate = gate_base.base_class(*params)
        matrix = Operator(gate.to_matrix())
        rules = EQUIVALENCES.get_entry(gate)
        non_equals = [i for i, rule in enumerate(rules) if Operator(rule) != matrix]
        if non_equals:
            # Accessing the base rule rather than rebinding with a parametric gate makes it easier
            # to locate the problematic part of the library definition, since the `Parameter` names
            # will match exactly.
            all_rules = EQUIVALENCES._get_equivalences(key)  # noqa: SLF001
            bads = [all_rules[i].circuit for i in non_equals]
            self.fail(f"Failed rules for {key}:\n" + "\n".join(str(bad) for bad in bads))
