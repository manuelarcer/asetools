# ASEtools

**A comprehensive Python toolkit for computational materials science with VASP and ASE**

[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![ASE](https://img.shields.io/badge/ASE-compatible-green.svg)](https://wiki.fysik.dtu.dk/ase/)
[![VASP](https://img.shields.io/badge/VASP-5%20%7C%206-red.svg)](https://www.vasp.at/)

## Overview

ASEtools is a specialized Python package designed for computational materials scientists working with density functional theory (DFT) calculations, particularly VASP (Vienna Ab initio Simulation Package). Built on top of the Atomic Simulation Environment (ASE), it provides a comprehensive suite of analysis tools, workflow management, and utilities for materials science research.

## Key Features

- **VASP Analysis Tools**: Convergence checking, energy/force extraction, magnetic moment analysis
- **Electronic Structure Analysis**: Comprehensive DOS analysis with orbital projections and band center calculations
- **Band Center Analysis**: d-band and p-band center calculations for catalysis and electronic structure studies
- **Bond Valence Analysis**: Bond valence sum calculations using Brown's equation for structural validation
- **Surface Science**: Automated adsorbate placement and surface analysis
- **Electrochemistry**: Applied potential calculations with Fermi shift corrections
- **Reaction Pathways**: NEB and dimer method support for transition state finding
- **Constraint Management**: Hookean spring constraints for O-H bonds and PCET reactions
- **Workflow Management**: YAML-based configuration with ASE optimizer support (FIRE, BFGS, LBFGS)
- **Database Integration**: ASE database support with duplicate detection
- **Command-Line Tools**: 15+ specialized scripts for various VASP tasks

## Installation

### Requirements

```bash
pip install pandas numpy scipy matplotlib ase pyyaml
```

**Optional dependencies:**
- `vasp-interactive`: Required for ASE optimizers (FIRE, BFGS) and dimer method
  ```bash
  pip install vasp-interactive
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
from asetools.doscar_analysis import DOS, calculate_band_center

# Modern DOS class approach (recommended)
dos = DOS('DOSCAR')

# Calculate d-band center for transition metal catalysis
d_center = dos.calculate_band_center([0], orbitals='all-d')
print(f"d-band center: {d_center:.3f} eV")

# Calculate t2g and eg orbital contributions
t2g_center = dos.calculate_band_center([0], orbitals='t2g')
eg_center = dos.calculate_band_center([0], orbitals='eg')

# Spin-resolved band center analysis
spin_resolved = dos.calculate_band_center([0], orbitals='all-d', spin_treatment='separate')
print(f"Spin-up d-center: {spin_resolved['up']:.3f} eV")
print(f"Spin-down d-center: {spin_resolved['down']:.3f} eV")

# Legacy approach
dos_data = extract_dos('DOSCAR')
atoms_of_interest = [0, 1, 2]
states = ['s_states', 'p_states', 'd_states']
energy, dos_up, dos_down = extract_pdos_perstate(dos_data, atoms_of_interest, states)
```

### Bond Valence Sum Analysis

```python
from asetools.bond_valence import BondValenceSum
from ase.io import read

atoms = read('POSCAR')

# Basic BVS calculation
bvs = BondValenceSum(atoms, valence_states={'Ti': 4, 'O': -2})
df = bvs.analyze_structure()
print(df[['element', 'expected_valence', 'calculated_bvs', 'deviation']])

# Automatic valence determination (per-element or per-atom)
bvs_auto = BondValenceSum(atoms, valence_states={'O': -2},
                         auto_determine_valence=True,
                         per_atom_valence=True)  # Allow mixed valences
df_auto = bvs_auto.analyze_structure()
bvs_auto.print_valence_optimization_summary()
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
- **`calculate_band_center(doscarfile, atoms, orbitals/states)`**: d-band and p-band center calculations
- **`DOS` class**: Modern object-oriented interface with comprehensive plotting and analysis methods

### Bond Valence Analysis (`asetools.bond_valence`)
- **`BondValenceParameters`**: Bond valence parameter database from Brown's accumulated table
- **`BondValenceSum`**: Bond valence sum calculator using Brown's equation
- **`get_parameters(element1, valence1, element2, valence2)`**: Parameter lookup with reliability ordering
- **`calculate_bvs()`**: Bond valence sum calculation with neighbor detection
- **`analyze_structure()`**: Comprehensive structural analysis with validation metrics
- **`auto_determine_valence`**: Automatically optimize metal valence states to minimize BVS deviation
- **`per_atom_valence`**: Enable per-atom valence optimization (allows mixed valences per element)
- **`get_valence_optimization_results()`**: Detailed optimization results for each optimized atom
- **`print_valence_optimization_summary()`**: Display optimization details and tried valences
- **`allowed_pairs`**: Specify element pairs to consider (e.g., only Ti-O bonds)
- **`exclude_same_element`**: Automatically exclude same-element pairs (default: True)
- **`get_allowed_pairs()`**: List allowed element pairs for BVS calculation
- **`print_allowed_pairs()`**: Display allowed pairs with statistics

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
- **ASE optimizer support**: FIRE, BFGS, LBFGS with VaspInteractive
- **Dimer method**: Transition state finding with MinModeAtoms

### Constraint Management (`asetools.constraints`)
- **`ConstraintManager`**: Hookean spring constraint system
- **`apply_stage_constraints(atoms, config)`**: YAML-based constraint application
- **Automatic O-H pair detection**: JSON configuration for PCET reactions

### Database Tools (`asetools.databases`)
- **`add_config_to_db(db, outcar, idname=None)`**: Adds structures to ASE database
- **`check_if_exists_in_db(db, atoms)`**: Duplicate structure detection
- **`db_to_pandas(db)`**: Database to DataFrame conversion

### Reaction Pathways (`asetools.nebanalysis`, `asetools.plots`)
- **`extract_neb_data(folder_path, final)`**: NEB energy extraction
- **`plot_nebs(list_dfs)`**: Multi-pathway NEB plotting
- **`add_line_to_pes(ax, data)`**: Potential energy surface visualization

## Command-Line Tools

ASEtools includes 15+ specialized command-line utilities:

### Core Analysis Tools
- **`getenergy`**: Quick convergence and energy summary
- **`summaryfolders`**: Batch analysis of multiple calculations with DataFrame output
- **`comparevaspparam`**: Side-by-side VASP parameter comparison

### Structure and Visualization
- **`asegui`**: Enhanced ASE GUI with error recovery
- **`view_outcars`**: Batch structure visualization
- **`sortposcarbyelement`**: POSCAR sorting by element and coordinates
- **`reorder_atoms`**: Flexible atom reordering by element and z-coordinate

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

# Reorder atoms by element and z-coordinate
reorder_atoms POSCAR --order Cu O H --z-order top-bottom

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

  # Optimization with Hookean constraints
  constrained_opt:
    stages:
      - name: OPT_WITH_CONSTRAINTS
        constraints:
          type: hookean
          config_file: oh_pairs.json  # JSON with O-H atom pairs
          spring_constant: 10.0
        steps:
          - name: optimization
            overrides: { nsw: 800 }
            optimizer: FIRE  # Requires vasp-interactive
            optimizer_kwargs: { fmax: 0.02, maxstep: 0.2 }
```

```python
from asetools.manager import VASPConfigurationFromYAML, run_workflow
from ase.io import read

config = VASPConfigurationFromYAML('config.yaml')
atoms = read('POSCAR')
run_workflow(atoms, config, 'optimization')
```

See `CONSTRAINTS_QUICKSTART.md` for detailed constraint examples.

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
