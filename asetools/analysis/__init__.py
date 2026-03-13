"""
Analysis package for ASEtools.

Provides analysis functionality for atomic structures including:
- VASP output analysis (convergence, energy, forces, parameters)
- Symmetry analysis (requires optional spglib dependency)
"""

# Import VASP analysis functions directly for backward compatibility
from .vasp import (
    check_energy_and_maxforce,
    check_outcar_convergence,
    classify_calculation_type,
    extract_comprehensive_metadata,
    extract_magnetic_moments,
    find_initial_structure,
    get_parameter_from_run,
)

# Lazy imports for optional dependencies
_SYMMETRY_AVAILABLE = None


def _check_symmetry_available():
    """Check if symmetry functionality is available."""
    global _SYMMETRY_AVAILABLE
    if _SYMMETRY_AVAILABLE is None:
        try:
            import spglib

            _SYMMETRY_AVAILABLE = True
        except ImportError:
            _SYMMETRY_AVAILABLE = False
    return _SYMMETRY_AVAILABLE


def symmetry_available() -> bool:
    """
    Check if symmetry analysis is available.

    Returns True if spglib is installed.
    """
    return _check_symmetry_available()


# Lazy load SymmetryAnalyzer
def __getattr__(name):
    if name == "SymmetryAnalyzer":
        from .symmetry import SymmetryAnalyzer

        return SymmetryAnalyzer
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    # Symmetry analysis
    "SymmetryAnalyzer",
    "check_energy_and_maxforce",
    # VASP analysis
    "check_outcar_convergence",
    "classify_calculation_type",
    "extract_comprehensive_metadata",
    "extract_magnetic_moments",
    "find_initial_structure",
    "get_parameter_from_run",
    "symmetry_available",
]
