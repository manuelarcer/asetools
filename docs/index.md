# ASEtools Documentation

**Python toolkit for computational materials science with VASP and ASE**

Version: 0.2.0

## Package Overview

ASEtools provides analysis tools, workflow management, and utilities for researchers working with VASP (Vienna Ab initio Simulation Package) density functional theory calculations, built on the Atomic Simulation Environment (ASE). It covers OUTCAR parsing, DOS analysis, thermochemistry corrections, structure manipulation, electrochemistry, and database integration.

## Package Structure

```
asetools/
├── analysis/          # VASP output analysis (convergence, energy, forces, metadata)
│   ├── vasp.py        # Core analysis functions
│   └── symmetry.py    # Symmetry analysis (requires spglib)
├── electronic/        # Electronic structure
│   └── doscar.py      # DOS class — total DOS, PDOS, orbital projections, band centers
├── structure/         # Structure analysis
│   ├── adsorbate.py   # SurfaceAnalyzer for adsorbate placement
│   └── bond_valence.py# BondValenceSum, BondValenceParameters
├── electrochemistry/  # Electrochemistry
│   └── appliedpotential.py  # Potential-dependent calculations with Fermi corrections
├── pathways/          # Reaction pathway methods
│   ├── neb.py         # NEB pathway energy extraction
│   └── dimer.py       # Dimer method for saddle point searches
├── thermodynamics/    # Ab initio thermodynamics
│   └── ab_initio.py   # ThermodynamicsCalculator and related classes
├── workflow/          # Workflow management
│   ├── manager.py     # YAML-based multi-stage workflow runner
│   ├── calculatorsetuptools.py  # VASPConfigurationFromYAML, setup_initial_magmom
│   ├── constraints.py # ConstraintManager (Hookean + FixAtoms)
│   ├── logger.py      # Workflow logging
│   └── sample_yaml/   # Example YAML configurations
├── database/          # Database integration
│   └── databases.py   # ASE-DB utilities with duplicate detection
├── plotting/          # Visualization
│   └── plots.py       # PES and energy profile plotting
├── parsers/           # File parsers
│   └── vasp_outcar.py # Extensible OUTCAR parser class hierarchy
├── cli/               # 17 command-line tools
└── data/              # Runtime data files
    ├── bvparm2020.cif
    └── surface_properties.json
```

## Installation

```bash
# Editable install
pip install -e .

# With development tools (pytest, ruff)
pip install -e ".[dev]"

# With symmetry analysis (spglib)
pip install -e ".[symmetry]"

# With vasp-interactive support
pip install -e ".[interactive]"
```

**Core requirements:** pandas, matplotlib, numpy, scipy, ase, pyyaml

**Python:** >= 3.8

## Quick Start

### Check a VASP calculation

```python
from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce

converged, vasp_version = check_outcar_convergence("OUTCAR")
energy, maxforce = check_energy_and_maxforce("OUTCAR")
print(f"Converged: {converged}, E={energy:.3f} eV, Fmax={maxforce:.4f} eV/Å")
```

### Analyze DOS

```python
from asetools.electronic.doscar import DOS

dos = DOS("DOSCAR")
print(f"Fermi energy: {dos.fermi_energy:.3f} eV, Atoms: {dos.natoms}")

# d-band center for atoms 0-4
center = dos.calculate_band_center(atoms=[0, 1, 2, 3, 4], orbitals="all-d")
```

### Apply Hookean constraints

```python
from asetools.workflow.constraints import ConstraintManager
from ase.io import read

atoms = read("CONTCAR")
cm = ConstraintManager(distance_factor=1.134)
cm.apply_from_json(atoms, "proton_mappings.json", k=20.0)
```

### Load YAML workflow configuration

```python
from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML, setup_initial_magmom
from ase.io import read

config = VASPConfigurationFromYAML("config.yaml", system="NCA")
atoms = read("POSCAR")
atoms = setup_initial_magmom(atoms, config.initial_magmom_data)
```

### Extract calculations to a DataFrame

```python
# Using ASE database
from asetools.database.databases import db_to_pandas
from ase.db import connect

db = connect("calculations.db")
df = db_to_pandas(db, columns=["name", "id", "energy", "free_energy", "magmom"])
```

### Command-line usage

```bash
getenergy OUTCAR                             # convergence + energy + max force
summaryfolders -m                            # batch summary with magnetic moments
plotdos --atoms 0-5 --orbitals d             # d-band PDOS for atoms 0-5
reorder-atoms POSCAR --order Cu O H          # reorder by element
thermochem OUTCAR --temp 298                 # Gibbs free energy corrections
```

## Documentation Index

| Document | Contents |
|----------|----------|
| [CLI Reference](cli_reference.md) | All 17 command-line tools with exact usage and flags |
| [API Reference](api_reference.md) | Python API for all subpackages |
| [Ab Initio Thermodynamics](ab_initio_thermodynamics.md) | Thermodynamics module deep-dive |
| [Constraints Quickstart](constraints_quickstart.md) | Hookean constraint quick start |
| [Constraints Migration](constraints_migration.md) | Migrating from manual to YAML constraints |
| [VASP Calculation Guide](vasp_calculation_guide.md) | Tips for spin-polarized and magnetic calculations |

## Development

```bash
pytest                  # run test suite
pytest -v               # verbose
ruff check asetools/    # lint
ruff format asetools/   # format
```

Test files live in `tests/`, test data in `tests/data/`. Configuration is in `pyproject.toml` under `[tool.pytest.ini_options]` and `[tool.ruff]`.
