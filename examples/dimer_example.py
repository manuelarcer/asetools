#!/usr/bin/env python

"""
Example script demonstrating dimer method usage with asetools.

This example shows how to:
1. Set up a dimer calculation with VaspInteractive
2. Use MODECAR files for initial displacement
3. Apply ASE constraints (FixAtoms)
4. Run dimer optimization through the workflow manager
"""

import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms

# Example 1: Direct dimer calculation (requires VaspInteractive)
def example_direct_dimer():
    """Example of running dimer calculation directly."""

    print("=== Example 1: Direct Dimer Calculation ===")

    # Load your structure
    atoms = read('POSCAR')  # or CONTCAR, or any ASE-supported format

    # Apply constraints if needed (e.g., fix bottom layers of slab)
    # Fix atoms based on z-coordinate (bottom 2 layers)
    z_coords = atoms.get_positions()[:, 2]
    z_min = z_coords.min()
    z_threshold = z_min + 2.0  # Fix atoms within 2 Å of bottom

    mask = z_coords < z_threshold
    constraint = FixAtoms(mask=mask)
    atoms.set_constraint(constraint)

    print(f"Applied FixAtoms constraint to {np.sum(mask)} atoms")

    # Set up VaspInteractive calculator
    from vasp_interactive import VaspInteractive

    calc_kwargs = {
        'encut': 300,
        'kspacing': 1.0,
        'ediff': 1e-6,
        'nsw': 2000,
        'ibrion': -1,  # Let ASE handle optimization
        'isif': 0,
        'prec': 'Normal',
        'xc': 'PBE',
    }

    # Use VaspInteractive with context manager
    with VaspInteractive(**calc_kwargs) as calc:
        atoms.calc = calc

        # Set up and run dimer calculation
        from asetools.dimer import setup_dimer_atoms
        from ase.mep import MinModeTranslate

        # MODECAR will be used if present, otherwise random displacement
        d_atoms = setup_dimer_atoms(atoms)

        # Run dimer optimization
        dimer_opt = MinModeTranslate(d_atoms, logfile='dimer.log', trajectory='dimer.traj')
        dimer_opt.run(fmax=0.001)

        print("Dimer optimization completed")


# Example 2: Using the workflow manager system
def example_workflow_dimer():
    """Example using the YAML workflow system."""

    print("=== Example 2: Workflow Manager Dimer ===")

    # Create a YAML configuration file
    yaml_content = """
basic:
  encut: 300
  ediff: !!float 1e-6
  nsw: 2000
  ibrion: -1
  prec: Normal
  xc: PBE

systems:
  default:
    kspacing: 1.0
    isif: 0

workflows:
  dimer_search:
    stages:
      - name: DIMER_OPT
        steps:
          - name: dimer_calculation
            overrides: {}
            optimizer: DIMER
            optimizer_kwargs:
              fmax: 0.001
              initial_eigenmode_method: 'displacement'
              displacement_method: 'vector'
              logfile: 'workflow_dimer.log'
              trajectory: 'workflow_dimer.traj'

globals:
  VASP_PP_PATH: "/path/to/your/pseudopotentials"
"""

    # Save YAML configuration
    with open('dimer_config.yaml', 'w') as f:
        f.write(yaml_content)

    # Load structure and apply constraints
    atoms = read('POSCAR')

    # Apply constraints (same as example 1)
    z_coords = atoms.get_positions()[:, 2]
    mask = z_coords < (z_coords.min() + 2.0)
    atoms.set_constraint(FixAtoms(mask=mask))

    # Run workflow
    from asetools.manager.calculatorsetuptools import VASPConfigurationFromYAML
    from asetools.manager.manager import run_workflow

    cfg = VASPConfigurationFromYAML('dimer_config.yaml', 'default')
    run_workflow(atoms, cfg, 'dimer_search')

    print("Workflow dimer optimization completed")


# Example 3: Working with MODECAR files
def example_modecar_usage():
    """Example of MODECAR file handling."""

    print("=== Example 3: MODECAR File Handling ===")

    from asetools.dimer import read_modecar, write_modecar, generate_displacement_vector

    atoms = read('POSCAR')

    # Generate displacement vector manually
    displacement_vector = generate_displacement_vector(atoms, method='random')

    # Save to MODECAR
    write_modecar(displacement_vector, atoms, 'MODECAR')
    print("Wrote MODECAR file")

    # Read back MODECAR
    if os.path.exists('MODECAR'):
        disp_read = read_modecar('MODECAR')
        print(f"Read MODECAR: shape={disp_read.shape}, magnitude={np.linalg.norm(disp_read):.6f}")

        # MODECAR will be automatically used in dimer calculations
        print("MODECAR will be automatically detected in dimer calculations")


# Example 4: Advanced constraint setup
def example_advanced_constraints():
    """Example of advanced constraint setups for surfaces."""

    print("=== Example 4: Advanced Constraint Setup ===")

    atoms = read('POSCAR')

    # Different constraint strategies:

    # 1. Fix by layer (z-coordinate)
    positions = atoms.get_positions()
    z_coords = positions[:, 2]
    z_sorted = np.sort(np.unique(z_coords))

    # Fix bottom 2 layers
    n_layers_fix = 2
    if len(z_sorted) >= n_layers_fix:
        z_threshold = z_sorted[n_layers_fix - 1] + 0.1
        mask_layers = z_coords <= z_threshold
        constraint_layers = FixAtoms(mask=mask_layers)
        print(f"Layer constraint: fixing {np.sum(mask_layers)} atoms in bottom {n_layers_fix} layers")

    # 2. Fix by element (e.g., substrate atoms)
    symbols = atoms.get_chemical_symbols()
    substrate_elements = ['Al', 'Ni']  # Example substrate
    mask_elements = np.array([s in substrate_elements for s in symbols])
    constraint_elements = FixAtoms(mask=mask_elements)
    print(f"Element constraint: fixing {np.sum(mask_elements)} substrate atoms")

    # 3. Fix by distance from a point (e.g., center of surface)
    center = np.mean(positions, axis=0)
    distances = np.linalg.norm(positions - center, axis=1)
    mask_distance = distances > 5.0  # Fix atoms more than 5 Å from center
    constraint_distance = FixAtoms(mask=mask_distance)
    print(f"Distance constraint: fixing {np.sum(mask_distance)} atoms far from center")

    # Apply one of the constraints
    atoms.set_constraint(constraint_layers)

    # For dimer calculations, you can also use the mask parameter for additional control
    dimer_mask = [False] * len(atoms)
    # Only allow specific atoms to participate in dimer search
    adsorbate_indices = [i for i, s in enumerate(symbols) if s in ['C', 'O', 'H']]
    for idx in adsorbate_indices:
        dimer_mask[idx] = True

    print(f"Dimer mask: {np.sum(dimer_mask)} atoms will participate in dimer search")

    return atoms, dimer_mask


# Example 5: Analysis of dimer results
def example_dimer_analysis():
    """Example of analyzing dimer calculation results."""

    print("=== Example 5: Dimer Result Analysis ===")

    from asetools.dimer import check_dimer_convergence, extract_saddle_point_info
    from ase.io import Trajectory

    # Analyze completed dimer calculation
    if os.path.exists('dimer.traj'):
        traj = Trajectory('dimer.traj', 'r')

        print(f"Dimer trajectory contains {len(traj)} images")

        # Get final structure
        final_atoms = traj[-1]

        # If you have the MinModeAtoms object (d_atoms), you can check convergence
        # conv_info = check_dimer_convergence(d_atoms)
        # saddle_info = extract_saddle_point_info(d_atoms)

        print(f"Final energy: {final_atoms.get_potential_energy():.6f} eV")

        # Analyze trajectory
        energies = [atoms.get_potential_energy() for atoms in traj]
        print(f"Energy change: {energies[-1] - energies[0]:.6f} eV")

        # Check for convergence
        forces = final_atoms.get_forces()
        max_force = np.max(np.linalg.norm(forces, axis=1))
        print(f"Final max force: {max_force:.4f} eV/Å")

        traj.close()

    # Check if final MODECAR was saved
    if os.path.exists('MODECAR_final'):
        print("Final eigenvector saved to MODECAR_final")


def main():
    """Run examples based on available files."""

    print("Dimer Method Examples with asetools\\n")

    if not os.path.exists('POSCAR'):
        print("No POSCAR file found. Creating example files...")
        # Create a simple example structure
        from ase.build import fcc111, add_adsorbate

        slab = fcc111('Al', size=(2, 2, 2), vacuum=7.0)
        add_adsorbate(slab, 'O', height=2.0, position='hcp')

        write('POSCAR', slab, format='vasp')
        print("Created example POSCAR file with Al(111) + O adsorbate")

    # Run examples that don't require VaspInteractive
    try:
        example_modecar_usage()
        example_advanced_constraints()
        example_dimer_analysis()

        print("\\n=== Summary ===")
        print("Examples completed. To run full dimer calculations:")
        print("1. Set up VaspInteractive environment")
        print("2. Configure VASP pseudopotentials")
        print("3. Use example_direct_dimer() or example_workflow_dimer()")
        print("4. Monitor convergence in log files and trajectories")

    except Exception as e:
        print(f"Example failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()