# ASEtools

**A comprehensive Python toolkit for computational materials science with VASP and ASE**

[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![ASE](https://img.shields.io/badge/ASE-compatible-green.svg)](https://wiki.fysik.dtu.dk/ase/)
[![VASP](https://img.shields.io/badge/VASP-5%20%7C%206-red.svg)](https://www.vasp.at/)

## Overview

ASEtools is a specialized Python package designed for computational materials scientists working with density functional theory (DFT) calculations, particularly VASP (Vienna Ab initio Simulation Package). Built on top of the Atomic Simulation Environment (ASE), it provides a comprehensive suite of analysis tools, workflow management, and utilities for materials science research.

## Key Features

- **VASP Analysis Tools**: Convergence checking, energy/force extraction, magnetic moment analysis
- **Electronic Structure Analysis**: Comprehensive DOS analysis with orbital projections
- **Surface Science**: Automated adsorbate placement and surface analysis
- **Electrochemistry**: Applied potential calculations with Fermi shift corrections
- **Reaction Pathways**: NEB analysis and potential energy surface plotting
- **Workflow Management**: YAML-based configuration system for multi-stage calculations
- **Database Integration**: ASE database support with duplicate detection
- **Command-Line Tools**: 12+ specialized scripts for various VASP tasks

## Installation

### Requirements

```bash
pip install pandas numpy scipy matplotlib ase pyyaml
```

### Install ASEtools

```bash
git clone https://github.com/manuelarcer/asetools.git
cd asetools
pip install -e .
```

## Quick Start

### Basic VASP Analysis

```python
from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce

# Check convergence and extract properties
converged, version = check_outcar_convergence('OUTCAR')
energy, max_force = check_energy_and_maxforce('OUTCAR')

print(f"Converged: {converged}, Energy: {energy:.4f} eV, Max Force: {max_force:.3f} eV/Å")
```

### DOS Analysis

```python
from asetools.doscar_analysis import extract_dos, extract_pdos_perstate

# Extract density of states
dos_data = extract_dos('DOSCAR')

# Get projected DOS for specific atoms and orbitals
atoms_of_interest = [0, 1, 2]
states = ['s_states', 'p_states', 'd_states']
energy, dos_up, dos_down = extract_pdos_perstate(dos_data, atoms_of_interest, states)
```

### Surface Analysis and Adsorbate Placement

```python
from asetools.adsorbate import SurfaceAnalyzer
from ase.io import read

# Load surface structure
atoms = read('POSCAR')

# Analyze surface and place adsorbate
analyzer = SurfaceAnalyzer(atoms)
surface_sites = analyzer.find_surface_neighbors()
analyzer.add_adsorbate_to_mid(surface_sites[0], adsorbate='H', z_off=1.5)

print(f"Adsorbate distances: {analyzer.adsneighdistances}")
```

## Core Modules

### Analysis (`asetools.analysis`)
- **`check_outcar_convergence(outcar)`**: Validates VASP calculation convergence
- **`check_energy_and_maxforce(outcar, magmom=False)`**: Extracts energy, forces, and magnetic moments
- **`extract_magnetic_moments(outcar, atoms_list)`**: Gets magnetic moments for specific atoms
- **`get_parameter_from_run(outcar, parameter='ISIF')`**: Extracts VASP parameters

### DOS Analysis (`asetools.doscar_analysis`)
- **`extract_dos(doscarfile)`**: Complete DOSCAR parsing with spin support
- **`extract_pdos_perstate(data, atoms, states)`**: Projects DOS onto s/p/d orbitals
- **`extract_pdos_perorbital(data, atoms, orbitals)`**: Detailed orbital projections with crystalline field
- **`extract_fermi_e(doscarfile)`**: Simple Fermi energy extraction

### Surface Science (`asetools.adsorbate`)
- **`SurfaceAnalyzer`**: Comprehensive surface analysis toolkit
  - Surface atom identification
  - Triangular site detection
  - Automated adsorbate placement
  - Local environment analysis

### Applied Potential (`asetools.appliedpotential`)
- **`extract_corrected_energy_fermie(folders, calc_zero)`**: Potential-dependent energy corrections
- **`fit_data(X, Y, fit_type='polynomial')`**: Energy fitting with polynomial/spline methods
- **`get_energy_at_givenpotential(results, desiredU=0.)`**: Energy interpolation at specific potentials

### Workflow Management (`asetools.manager`)
- **`VASPConfigurationFromYAML`**: YAML-based configuration system
- **`run_workflow(atoms, config, workflow_name)`**: Multi-stage calculation execution
- **`make_calculator(config, overrides=None)`**: VASP calculator creation

### Database Tools (`asetools.databases`)
- **`add_config_to_db(db, outcar, idname=None)`**: Adds structures to ASE database
- **`check_if_exists_in_db(db, atoms)`**: Duplicate structure detection
- **`db_to_pandas(db)`**: Database to DataFrame conversion

### Reaction Pathways (`asetools.nebanalysis`, `asetools.plots`)
- **`extract_neb_data(folder_path, final)`**: NEB energy extraction
- **`plot_nebs(list_dfs)`**: Multi-pathway NEB plotting
- **`add_line_to_pes(ax, data)`**: Potential energy surface visualization

## Command-Line Tools

ASEtools includes 12+ specialized command-line utilities:

### Core Analysis Tools
- **`getenergy`**: Quick convergence and energy summary
- **`summaryfolders`**: Batch analysis of multiple calculations with DataFrame output
- **`comparevaspparam`**: Side-by-side VASP parameter comparison

### Structure and Visualization
- **`asegui`**: Enhanced ASE GUI with error recovery
- **`view_outcars`**: Batch structure visualization
- **`sortposcarbyelement`**: POSCAR sorting by element and coordinates

### Specialized Analysis
- **`gibbsFE`**: Gibbs free energy corrections from vibrational analysis
- **`genpotentialdepcalc`**: Generates potential-dependent calculation folders

### File Management
- **`vaspbackup`**: Intelligent VASP file backup system
- **`vasp_2_arc`**: Format conversion utilities
- **`remove_slashes`**: Corrupted file repair tool

### Usage Examples

```bash
# Quick energy summary
getenergy OUTCAR

# Batch analysis with magnetic moments
summaryfolders -m

# Sort POSCAR by elements
sortposcarbyelement -f POSCAR -axis z

# Backup VASP files
vaspbackup -p "OUTCAR*" -c
```

## Workflow Management

ASEtools provides a powerful YAML-based workflow system for complex multi-stage calculations:

```yaml
# config.yaml example
basic_config:
  xc: PBE
  encut: 400
  kpts: [4, 4, 4]

workflows:
  optimization:
    stages:
      - name: "rough_opt"
        ediff: 1e-4
        nsw: 100
      - name: "fine_opt"
        ediff: 1e-6
        nsw: 200
```

```python
from asetools.manager import VASPConfigurationFromYAML, run_workflow
from ase.io import read

config = VASPConfigurationFromYAML('config.yaml')
atoms = read('POSCAR')
run_workflow(atoms, config, 'optimization')
```

## Advanced Features

### Electrochemical Analysis
```python
from asetools.appliedpotential import extract_corrected_energy_fermie, fit_data

# Analyze potential-dependent calculations
folders = ['U_-0.5', 'U_0.0', 'U_0.5', 'U_1.0']
results = extract_corrected_energy_fermie(folders, calc_zero='U_0.0')

# Fit energy vs. potential
potentials = [-0.5, 0.0, 0.5, 1.0]
energies = [r['corrected_energy'] for r in results]
fit_results = fit_data(potentials, energies, fit_type='polynomial', order=2)
```

### Database Integration
```python
from asetools.databases import add_config_to_db, db_to_pandas
from ase.db import connect

db = connect('calculations.db')
add_config_to_db(db, 'OUTCAR', idname='surface_opt')
df = db_to_pandas(db)
print(df[['name', 'energy', 'magmom']])
```

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

## License

This project is open source. Please check the repository for license details.

## Citation

If you use ASEtools in your research, please consider citing this package and the underlying ASE framework.

## Support

For questions, issues, or feature requests, please visit the [GitHub repository](https://github.com/manuelarcer/asetools.git).
