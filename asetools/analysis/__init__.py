"""
Analysis package for ASEtools.

Provides analysis functionality for atomic structures including:
- Symmetry analysis (requires optional spglib dependency)

Modules are lazily imported when accessed.
"""

# Lazy imports to avoid loading optional dependencies at package import time
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
    if name == 'SymmetryAnalyzer':
        from .symmetry import SymmetryAnalyzer
        return SymmetryAnalyzer
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ['SymmetryAnalyzer', 'symmetry_available']
