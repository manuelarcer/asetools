# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ASEtools is a Python package for computational materials science, specifically designed for VASP (Vienna Ab initio Simulation Package) and ASE (Atomic Simulation Environment) integration. It provides analysis tools, workflow management, and utilities for DFT calculations.

## Development Commands

### Installation
```bash
pip install -e .  # Editable install for development
pip install -e ".[dev]"  # With test dependencies
```

### Testing
```bash
pytest  # Run all tests (configured in pyproject.toml)
pytest tests/test_analysis.py  # Run specific test file
pytest -v  # Verbose output
```

### Dependencies
Core dependencies are defined in `pyproject.toml`:
- **Required**: pandas, matplotlib, numpy, scipy, ase, pyyaml
- **Optional**:
  - `symmetry` - spglib for symmetry analysis
  - `interactive` - vasp-interactive for ASE optimizer support
  - `dev` - pytest for testing

## Code Architecture

### Package Structure
```
asetools/
├── analysis/          # VASP output analysis (convergence, energy, forces, symmetry)
│   ├── vasp.py        # Core VASP analysis functions
│   └── symmetry.py    # Symmetry analysis (requires spglib)
├── electronic/        # Electronic structure analysis
│   └── doscar.py      # DOS, PDOS, orbital projections, band center
├── structure/         # Structure analysis
│   ├── adsorbate.py   # SurfaceAnalyzer for adsorbate placement
│   └── bond_valence.py # Bond valence sum analysis
├── electrochemistry/  # Electrochemistry support
│   ├── appliedpotential.py  # Potential-dependent calculations
│   └── copy_appliedpotential.py
├── pathways/          # Reaction pathway analysis
│   ├── neb.py         # NEB calculations
│   └── dimer.py       # Dimer method
├── thermodynamics/    # Thermodynamics calculations
│   └── ab_initio.py   # Ab initio thermodynamics
├── workflow/          # Workflow management
│   ├── manager.py     # YAML-based multi-stage workflow runner
│   ├── calculatorsetuptools.py  # VASP calculator configuration
│   ├── constraints.py # Constraint management
│   ├── logger.py      # Workflow logging
│   └── sample_yaml/   # Example YAML configurations
├── database/          # Database tools
│   └── databases.py   # ASE database integration
├── plotting/          # Visualization
│   └── plots.py       # PES and energy profile plotting
├── cli/               # Command-line scripts
│   ├── getenergy.py   # Quick energy/convergence summary
│   ├── summaryfolders.py  # Batch analysis
│   ├── plotdos.py     # DOS plotting
│   ├── vasp2db.py     # VASP to database conversion
│   └── ...            # 17 CLI tools total
├── parsers/           # File parsers
│   └── vasp_outcar.py # OUTCAR parser
└── data/              # Runtime data files
    ├── bvparm2020.cif
    └── surface_properties.json
```

### Key Classes and Patterns
- **`SurfaceAnalyzer`**: Surface analysis and adsorbate placement
- **`BondValenceSum`**: Bond valence sum calculator
- **`VASPConfigurationFromYAML`**: Workflow configuration management
- **`ConstraintManager`**: Constraint application from JSON configs
- **`ThermodynamicsCalculator`**: Ab initio thermodynamics
- **`DOS`**: Electronic density of states analysis
- **Integration-first design**: All modules work with ASE Atoms objects
- **VASP-centric**: Functions expect standard VASP output files

### Testing Strategy
- Uses pytest with config in `pyproject.toml`
- Test files in `tests/`, test data in `tests/data/`
- Runtime data (bvparm2020.cif, surface_properties.json) stays in `asetools/data/`
- Integration tests using real VASP output files

### Documentation
- `docs/` directory contains guides:
  - `vasp_calculation_guide.md` - VASP calculation setup
  - `constraints_migration.md` - Constraint system migration
  - `constraints_quickstart.md` - Quick start for constraints
  - `ab_initio_thermodynamics.md` - Thermodynamics module guide

## Git Workflow
- Stage and commit changes (ask first) after making changes to the code
- Commit messages should be concise but descriptive (max 5 lines)
- When prompted to alter the code significantly, ask about creating a new branch
