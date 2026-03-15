# ASEtools Documentation

**Python toolkit for computational materials science with VASP and ASE**

Version: 0.1.0

## Package Overview

ASEtools provides analysis tools, workflow management, and utilities for researchers working with VASP density functional theory calculations, built on the Atomic Simulation Environment (ASE).

## Package Structure

```
asetools/
├── analysis/          # VASP output analysis (convergence, energy, forces, symmetry)
├── electronic/        # DOS, PDOS, orbital projections, band center calculations
├── structure/         # SurfaceAnalyzer, bond valence sum analysis
├── electrochemistry/  # Applied potential calculations with Fermi corrections
├── pathways/          # NEB and dimer method for transition states
├── thermodynamics/    # Ab initio thermodynamics (coverage, surface energy)
├── workflow/          # YAML-based multi-stage workflow runner
├── database/          # ASE database integration with deduplication
├── plotting/          # PES and energy profile plotting
├── parsers/           # OUTCAR parser
├── cli/               # 17 command-line tools
└── data/              # Runtime data (bvparm2020.cif, surface_properties.json)
```

## Installation

```bash
git clone https://github.com/manuelarcer/asetools.git
cd asetools
pip install -e ".[dev]"    # editable install with test deps
```

### Requirements

Core: `numpy`, `scipy`, `pandas`, `matplotlib`, `ase`, `pyyaml`

Optional:
- `vasp-interactive` — required for ASE optimizer workflows (FIRE, BFGS) and dimer method
- `spglib` — required for symmetry analysis
- `pytest` — for running tests

## Documentation Index

| Document | Contents |
|----------|----------|
| [CLI Reference](cli_reference.md) | All 17 command-line tools with usage and examples |
| [API Reference](api_reference.md) | Module-level Python API documentation |
| [Ab Initio Thermodynamics](ab_initio_thermodynamics.md) | Thermodynamics module deep-dive |
| [Constraint Management](constraints_quickstart.md) | Hookean constraint quick start |
| [Constraint Migration Guide](constraints_migration.md) | Migrating from manual to YAML constraints |
| [VASP Calculation Guide](vasp_calculation_guide.md) | Tips for spin-polarized and magnetic calculations |

## Quick Start

### Check a VASP calculation

```python
from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce

converged, version = check_outcar_convergence('OUTCAR')
energy, max_force = check_energy_and_maxforce('OUTCAR')
```

### Analyze DOS

```python
from asetools.electronic.doscar import DOS

dos = DOS('DOSCAR')
d_center = dos.calculate_band_center([0], orbitals='all-d')
```

### Run a workflow

```python
from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML
from ase.io import read

config = VASPConfigurationFromYAML('config.yaml')
atoms = read('POSCAR')
# run_workflow(atoms, config, 'optimization')
```

### Command-line usage

```bash
getenergy OUTCAR                        # quick convergence + energy
summaryfolders -m                       # batch analysis with magmoms
plotdos --atoms 0-5 --orbitals all-d    # d-band DOS plot
reorder_atoms POSCAR --order Cu O H     # reorder by element
thermochem OUTCAR --temp 298            # Gibbs free energy corrections
```

## Development

```bash
pytest                  # run test suite
ruff check asetools/    # lint
ruff format asetools/   # format
```

## License

Open source. See repository for details.
