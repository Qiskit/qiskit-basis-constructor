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

__all__ = [
    "standard_equivalence_library",
    "BasisConstructor",
    "BasisConstructorError",
    "BasisConstructorPlugin",
    "GateCount",
    "LogFidelity",
    "DEFAULT_EQUIVALENCE_LIBRARY",
]
__version__ = "0.1.0"

from .constructor import BasisConstructor, BasisConstructorError, GateCount, LogFidelity
from .equivalence import standard_equivalence_library
from .plugin import BasisConstructorPlugin

_DEFAULT_EQUIVALENCE_LIBRARY = None


def __getattr__(name: str):
    if name == "DEFAULT_EQUIVALENCE_LIBRARY":
        global _DEFAULT_EQUIVALENCE_LIBRARY  # noqa: PLW0603
        if _DEFAULT_EQUIVALENCE_LIBRARY is None:
            _DEFAULT_EQUIVALENCE_LIBRARY = standard_equivalence_library()
        return _DEFAULT_EQUIVALENCE_LIBRARY
    raise AttributeError(f"module 'qiskit_basis_constructor' has no attribute '{name}'")
