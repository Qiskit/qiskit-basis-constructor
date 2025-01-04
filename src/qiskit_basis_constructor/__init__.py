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
