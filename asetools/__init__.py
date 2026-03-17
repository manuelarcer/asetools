# ASEtools - Tools for computational materials science
"""
ASEtools: Python toolkit for computational materials science with VASP and ASE.

Subpackages:
    analysis        — VASP output analysis (convergence, energy, forces, symmetry)
    electronic      — DOS, PDOS, orbital projections, band center
    structure       — SurfaceAnalyzer, bond valence
    electrochemistry — Applied potential calculations
    pathways        — NEB, dimer method
    thermodynamics  — Ab initio thermodynamics
    workflow        — YAML-based multi-stage workflow runner
    database        — ASE database integration
    plotting        — PES and energy profile plotting
    parsers         — OUTCAR parser
    cli             — CLI scripts

Example usage::

    from asetools.analysis import check_outcar_convergence
    from asetools.electronic import DOS, extract_dos
    from asetools.structure import SurfaceAnalyzer
    from asetools.thermodynamics import ThermodynamicsCalculator
"""

import importlib as _importlib

__version__ = "0.2.0"

# Subpackages available for lazy loading
_SUBPACKAGES = {
    "analysis",
    "electronic",
    "structure",
    "electrochemistry",
    "pathways",
    "thermodynamics",
    "workflow",
    "database",
    "plotting",
    "parsers",
    "cli",
}

# Backward-compatible direct imports from thermodynamics
# (these were previously the only top-level exports)
from .thermodynamics.ab_initio import (  # noqa: E402
    AdsorbateSpecies,
    InterpolationModel,
    LatticeGasModel,
    SurfaceProperties,
    ThermodynamicsCalculator,
)


def __getattr__(name: str):
    """Lazy-load subpackages on attribute access."""
    if name in _SUBPACKAGES:
        return _importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "AdsorbateSpecies",
    "InterpolationModel",
    "LatticeGasModel",
    "SurfaceProperties",
    # Backward-compatible top-level classes
    "ThermodynamicsCalculator",
    # Metadata
    "__version__",
    # Subpackages
    "analysis",
    "cli",
    "database",
    "electrochemistry",
    "electronic",
    "parsers",
    "pathways",
    "plotting",
    "structure",
    "thermodynamics",
    "workflow",
]
