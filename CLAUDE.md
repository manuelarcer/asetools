# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ASEtools is a Python package for computational materials science, specifically designed for VASP (Vienna Ab initio Simulation Package) and ASE (Atomic Simulation Environment) integration. It provides analysis tools, workflow management, and utilities for DFT calculations.

## Development Commands

### Installation
```bash
pip install -e .  # Editable install for development
```

### Testing
```bash
pytest  # Run all tests (configured in pytest.ini)
pytest asetools/tests/test_analysis.py  # Run specific test file
pytest -v  # Verbose output
```

### Dependencies
Core dependencies are defined in `setup.cfg`:
- pandas, matplotlib, numpy, scipy, ase
- Additional: pyyaml (for workflow management)

## Code Architecture

### Core Module Structure
- **`analysis.py`**: VASP output analysis (convergence, energy, forces, magnetic moments)
- **`doscar_analysis.py`**: Electronic structure analysis (DOS, PDOS, orbital projections)
- **`bond_valence.py`**: Bond valence sum analysis using Brown's equation with automatic valence determination
- **`adsorbate.py`**: Surface science toolkit with `SurfaceAnalyzer` class
- **`appliedpotential.py`**: Electrochemistry support with potential-dependent calculations
- **`databases.py`**: ASE database integration and DataFrame conversion
- **`manager/`**: YAML-based workflow management system for multi-stage calculations
- **`bin/`**: Command-line tools (14+ scripts for various VASP tasks)

### Key Classes and Patterns
- **`SurfaceAnalyzer`**: Main class for surface analysis and adsorbate placement
- **`BondValenceSum`**: Bond valence sum calculator with automatic valence optimization
- **`BondValenceParameters`**: Bond valence parameter database from Brown's accumulated table
- **`VASPConfigurationFromYAML`**: Workflow configuration management
- **Integration-first design**: All modules work with ASE Atoms objects
- **VASP-centric**: Functions expect standard VASP output files (OUTCAR, DOSCAR, POSCAR)

### Testing Strategy
- Uses pytest with real VASP output files in `asetools/data/`
- Integration tests rather than unit tests
- Test files: `test_analysis.py`, `test_database.py`, `test_doscar.py`, `test_freq.py`, `test_bond_valence.py`

### Command-Line Tools
Key scripts in `asetools/bin/`:
- **`getenergy.py`**: Quick energy/convergence summary
- **`summaryfolders.py`**: Batch analysis with DataFrame output
- **`asegui.py`**: Enhanced ASE GUI
- **`vaspbackup.py`**: VASP file backup system
- **`gibbsFE.py`**: Gibbs free energy corrections

## Common File Patterns

### VASP Integration
- Functions expect standard VASP filenames: `OUTCAR`, `DOSCAR`, `POSCAR`
- Supports both VASP 5 and VASP 6 formats
- Magnetic moment analysis throughout

### Configuration
- YAML-based workflow configuration in `manager/`
- Multi-stage calculation support
- System-specific parameter overrides

### Data Flow
1. Read VASP output files
2. Extract/analyze using core modules
3. Store in ASE database or export to pandas DataFrame
4. Visualize using `plots.py` module

## Development Notes

### Package Configuration
- Uses `setup.cfg` for metadata and dependencies
- Minimal `setup.py` for setuptools compatibility
- Scripts in `bin/` directory are NOT installed via pip - run directly from source
- `bin/` contains 15+ command-line tools that should be accessed via PATH or direct invocation

### Materials Science Focus
- Designed for surface science, electrochemistry, and catalysis
- Extensive DOS and electronic structure analysis
- Bond valence analysis for structural validation and valence state determination (with per-atom optimization)
- Reaction pathway analysis (NEB) support
- Applied potential calculations for electrochemistry

### Documentation
- Comprehensive README.md with API documentation and examples
- Inline docstrings in modules
- Command-line help available for all scripts

## Environment Notes
- You should run within python with /Users/juar/venv/workgeneric/bin/python

## Git Workflow
- Stage and commit changes (ask first) after making changes to the code. The commit message should be concise but descriptive (Max 5 lines)

## Memory for Code Interactions
- When prompted to alter the code in a significant way, ask if I want to create a new branch for the project and provide a name for it