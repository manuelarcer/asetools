# ASEtools

**Python toolkit for computational materials science with VASP and ASE**

[![CI](https://github.com/manuelarcer/asetools/actions/workflows/ci.yml/badge.svg)](https://github.com/manuelarcer/asetools/actions/workflows/ci.yml)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

ASEtools provides analysis tools, workflow management, and CLI utilities for researchers working with VASP density functional theory calculations, built on the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/).

## Package Structure

```
asetools/
├── analysis/          # VASP OUTCAR parsing — convergence, energy, forces, metadata
│   ├── vasp.py        # Core analysis functions
│   └── symmetry.py    # Symmetry analysis (requires spglib)
├── electronic/        # Electronic structure
│   └── doscar.py      # DOS class — total/projected DOS, orbital projections, band centers
├── structure/         # Structure analysis
│   ├── adsorbate.py   # SurfaceAnalyzer for adsorbate placement
│   └── bond_valence.py# Bond valence sum calculations
├── electrochemistry/  # Potential-dependent energy corrections
│   └── appliedpotential.py
├── pathways/          # Reaction pathways
│   ├── neb.py         # NEB energy extraction
│   └── dimer.py       # Dimer method for saddle points
├── thermodynamics/    # Ab initio thermodynamics
│   └── ab_initio.py   # ThermodynamicsCalculator
├── workflow/          # YAML-based multi-stage workflow runner
│   ├── manager.py     # run_workflow, make_calculator
│   ├── calculatorsetuptools.py  # VASPConfigurationFromYAML, setup_initial_magmom
│   ├── constraints.py # ConstraintManager (Hookean + FixAtoms)
│   └── logger.py      # Workflow logging
├── database/          # ASE database integration
│   └── databases.py
├── parsers/           # Low-level VASP output parsers
│   └── vasp_outcar.py
├── plotting/          # PES and energy profile plotting
│   └── plots.py
└── cli/               # 17 command-line tools
```

## Installation

```bash
git clone https://github.com/manuelarcer/asetools.git
cd asetools
pip install -e .
```

**Dependencies:** `numpy`, `scipy`, `pandas`, `matplotlib`, `ase`, `pyyaml`

**Optional:** `vasp-interactive` (for ASE optimizers), `spglib` (for symmetry analysis)

## Quick Start

```python
from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce

converged, version = check_outcar_convergence('OUTCAR')
energy, max_force = check_energy_and_maxforce('OUTCAR')
print(f"Converged: {converged}, Energy: {energy:.4f} eV, Max Force: {max_force:.3f} eV/Å")
```

```python
from asetools.electronic.doscar import DOS

dos = DOS('DOSCAR')
d_center = dos.calculate_band_center([0], orbitals='all-d')
print(f"d-band center: {d_center:.3f} eV")
```

## CLI Tools

17 command-line tools installed via `pip install -e .`:

| Tool | Description |
|------|-------------|
| `getenergy` | Quick convergence/energy check from OUTCAR |
| `summaryfolders` | Batch analysis across subdirectories |
| `comparevaspparam` | Side-by-side VASP parameter comparison |
| `outcar-extract` | Extract a structure frame from OUTCAR |
| `inputparam-from-outcar` | Extract categorized VASP parameters from OUTCAR |
| `plotdos` | Plot DOS/PDOS from DOSCAR |
| `gibbsFE` | Gibbs free energy from vibrational analysis |
| `thermochem` | Thermochemistry corrections |
| `genpotentialdepcalc` | Generate potential-dependent calculation folders |
| `calculate-bader-charges` | Bader charge analysis |
| `reorder-atoms` | Reorder atoms by element/z-coordinate |
| `asegui` | Enhanced ASE GUI launcher |
| `view-outcars` | Batch structure visualization |
| `vasp2db` | Extract VASP results to pandas DataFrame pickle |
| `vasp2arc` | VASP to .arc format conversion |
| `vaspbackup` | Intelligent VASP file backup |
| `remove-slashes` | Fix corrupted VASP files |

```bash
getenergy OUTCAR
summaryfolders -m
reorder-atoms POSCAR --order Cu O H --z-order top-bottom
```

## Documentation

- [Package Overview & Quick Start](docs/index.md)
- [CLI Reference](docs/cli_reference.md) — all 17 tools with flags and usage
- [API Reference](docs/api_reference.md) — Python API for all subpackages
- [Ab Initio Thermodynamics Guide](docs/ab_initio_thermodynamics.md)
- [Constraints Quick Start](docs/constraints_quickstart.md)
- [VASP Calculation Guide](docs/vasp_calculation_guide.md)

## License

Open source. See repository for license details.
