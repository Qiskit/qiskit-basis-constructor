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

import math
import typing

from qiskit.transpiler.preset_passmanagers.plugin import PassManagerStagePlugin
from qiskit.transpiler.passes import UnitarySynthesis, HighLevelSynthesis
from qiskit.transpiler import PassManager

from .constructor import BasisConstructor, GateCount, LogFidelity

if typing.TYPE_CHECKING:
    from qiskit.transpiler import PassManagerConfig


class BasisConstructorPlugin(PassManagerStagePlugin):
    def pass_manager(
        self, pass_manager_config: PassManagerConfig, optimization_level: int | None = None
    ) -> PassManager:
        # Delayed import is to avoid eager creation of the library during plugin scanning.
        from . import DEFAULT_EQUIVALENCE_LIBRARY

        equivalences = DEFAULT_EQUIVALENCE_LIBRARY

        optimization_level = optimization_level if optimization_level is not None else 2
        if optimization_level == 0:
            score = [GateCount(2), GateCount()]
        elif optimization_level == 1:
            # Use fidelity averaged over all gates.
            score = [LogFidelity(math.inf), GateCount(2), GateCount()]
        else:
            # Use nearest-nines fidelity averaging.
            score = [LogFidelity(math.log(10)), GateCount(2), GateCount()]

        return PassManager(
            [
                UnitarySynthesis(
                    approximation_degree=pass_manager_config.approximation_degree,
                    basis_gates=pass_manager_config.basis_gates,
                    coupling_map=pass_manager_config.coupling_map,
                    method=pass_manager_config.unitary_synthesis_method,
                    plugin_config=pass_manager_config.unitary_synthesis_plugin_config,
                    target=pass_manager_config.target,
                ),
                HighLevelSynthesis(
                    basis_gates=pass_manager_config.basis_gates,
                    coupling_map=pass_manager_config.coupling_map,
                    equivalence_library=equivalences,
                    hls_config=pass_manager_config.hls_config,
                    qubits_initially_zero=pass_manager_config.qubits_initially_zero,
                    target=pass_manager_config.target,
                    use_qubit_indices=True,
                ),
                BasisConstructor(
                    equivalences,
                    score,
                    target=pass_manager_config.target,
                ),
            ]
        )
