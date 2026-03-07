#!/usr/bin/env python
"""
Example usage of vasp2db tool and how to work with the generated database.

This script demonstrates:
1. Creating a database with vasp2db
2. Loading and exploring the database
3. Accessing structures and metadata
4. Filtering and analysis
"""

import pickle
import pandas as pd
from ase.io import write

# ============================================================
# Part 1: Creating the database (command-line)
# ============================================================
print("=" * 60)
print("Part 1: Creating the database")
print("=" * 60)
print("""
Run vasp2db from the command line:

# Basic usage
vasp2db --paths CuCu_2 NiCu_1 NiCu_2 --output database.pkl

# With all folders and relative energies
vasp2db --paths */ --output all_calcs.pkl -r -v

# Skip problematic folders
vasp2db --paths calc* --output db.pkl --skip-errors --verbose

# Custom initial structure pattern
vasp2db --paths folder* --output db.pkl --initial-pattern "POSCAR_init"
""")

# ============================================================
# Part 2: Loading and exploring the database
# ============================================================
print("\n" + "=" * 60)
print("Part 2: Loading and exploring the database")
print("=" * 60)

# Example database path - adjust to your actual database
database_path = 'NH2_var_sites_Cu7.pkl'

print(f"\nLoading database from: {database_path}")
print("Code:")
print("""
import pickle
import pandas as pd

with open('database.pkl', 'rb') as f:
    df = pickle.load(f)
""")

# Load the database
try:
    with open(database_path, 'rb') as f:
        df = pickle.load(f)

    print(f"\nSuccessfully loaded database with {len(df)} calculations")

    # Show available columns
    print(f"\nAvailable columns ({len(df.columns)}):")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i:2d}. {col}")

    # ============================================================
    # Part 3: Basic data exploration
    # ============================================================
    print("\n" + "=" * 60)
    print("Part 3: Basic data exploration")
    print("=" * 60)

    # Summary table
    print("\nSummary table:")
    summary_cols = ['Path', 'CalcType', 'Formula', 'Energy', 'MaxForce', 'Converged', 'Rel.E']
    available_cols = [col for col in summary_cols if col in df.columns]
    print(df[available_cols].to_string(index=False))

    # Energy statistics
    if 'Energy' in df.columns:
        print(f"\nEnergy statistics:")
        print(f"  Mean: {df['Energy'].mean():.3f} eV")
        print(f"  Std:  {df['Energy'].std():.3f} eV")
        print(f"  Min:  {df['Energy'].min():.3f} eV")
        print(f"  Max:  {df['Energy'].max():.3f} eV")

    # Convergence statistics
    if 'Converged' in df.columns:
        conv_pct = (df['Converged'].sum() / len(df)) * 100
        print(f"\nConvergence: {df['Converged'].sum()}/{len(df)} ({conv_pct:.1f}%)")

    # ============================================================
    # Part 4: Accessing structures
    # ============================================================
    print("\n" + "=" * 60)
    print("Part 4: Accessing ASE Atoms structures")
    print("=" * 60)

    print("\nCode to access structures:")
    print("""
# Access initial structure of first calculation
initial_atoms = df.loc[0, 'InitialStructure']

# Access final structure of first calculation
final_atoms = df.loc[0, 'FinalStructure']

# Number of atoms
natoms = len(final_atoms)

# Write structure to file
from ase.io import write
write('structure.vasp', final_atoms, format='vasp', direct=True)
    """)

    if 'FinalStructure' in df.columns and df.loc[0, 'FinalStructure'] is not None:
        final_atoms = df.loc[0, 'FinalStructure']
        print(f"\nExample: First calculation")
        print(f"  Number of atoms: {len(final_atoms)}")
        print(f"  Chemical formula: {final_atoms.get_chemical_formula()}")
        print(f"  Cell volume: {final_atoms.get_volume():.2f} Å³")

    # ============================================================
    # Part 5: Filtering and analysis
    # ============================================================
    print("\n" + "=" * 60)
    print("Part 5: Filtering and analysis")
    print("=" * 60)

    # Filter converged calculations
    if 'Converged' in df.columns:
        converged_df = df[df['Converged'] == True]
        print(f"\nConverged calculations: {len(converged_df)}/{len(df)}")

    # Find lowest energy configuration
    if 'Energy' in df.columns:
        min_idx = df['Energy'].idxmin()
        min_energy_path = df.loc[min_idx, 'Path']
        min_energy = df.loc[min_idx, 'Energy']
        print(f"\nLowest energy configuration:")
        print(f"  Path: {min_energy_path}")
        print(f"  Energy: {min_energy:.3f} eV")

        if 'Rel.E' in df.columns and df.loc[min_idx, 'Rel.E'] is not None:
            print(f"  Rel.E: {df.loc[min_idx, 'Rel.E']:.3f} eV")

    # Filter by calculation type
    if 'CalcType' in df.columns:
        calc_types = df['CalcType'].unique()
        print(f"\nCalculation types found: {', '.join(calc_types)}")

        for calc_type in calc_types:
            subset = df[df['CalcType'] == calc_type]
            print(f"  {calc_type}: {len(subset)} calculations")

    # Filter by maximum force
    if 'MaxForce' in df.columns:
        fmax_threshold = 0.02
        well_converged = df[df['MaxForce'] < fmax_threshold]
        print(f"\nWell converged (fmax < {fmax_threshold}): {len(well_converged)}/{len(df)}")

    # ============================================================
    # Part 6: Advanced usage - VASP parameters
    # ============================================================
    print("\n" + "=" * 60)
    print("Part 6: Accessing VASP parameters and metadata")
    print("=" * 60)

    print("\nCode examples:")
    print("""
# Check VASP parameters
print(df[['Path', 'ENCUT', 'EDIFFG', 'IBRION', 'NSW']].head())

# Access INCAR content
incar_content = df.loc[0, 'INCAR_full']

# Access POTCAR information
potcar_info = df.loc[0, 'POTCAR_info']
print("Pseudopotentials used:")
for pp in potcar_info:
    print(f"  {pp}")
    """)

    # Show VASP parameters table
    vasp_param_cols = ['Path', 'ENCUT', 'EDIFFG', 'IBRION', 'NSW', 'ISIF']
    available_param_cols = [col for col in vasp_param_cols if col in df.columns]
    if available_param_cols:
        print("\nVASP parameters (first 3 calculations):")
        print(df[available_param_cols].head(3).to_string(index=False))

    # Show POTCAR info example
    if 'POTCAR_info' in df.columns and df.loc[0, 'POTCAR_info']:
        print(f"\nPseudopotentials in first calculation:")
        for pp in df.loc[0, 'POTCAR_info'][:5]:  # Show first 5
            print(f"  {pp}")
        if len(df.loc[0, 'POTCAR_info']) > 5:
            print(f"  ... and {len(df.loc[0, 'POTCAR_info']) - 5} more")

    # ============================================================
    # Part 7: Exporting results
    # ============================================================
    print("\n" + "=" * 60)
    print("Part 7: Exporting results")
    print("=" * 60)

    print("\nCode examples:")
    print("""
# Export summary to CSV (without structures)
summary_cols = ['Path', 'CalcType', 'Energy', 'MaxForce', 'Converged', 'Rel.E']
df[summary_cols].to_csv('summary.csv', index=False)

# Export lowest energy structure
min_idx = df['Energy'].idxmin()
best_structure = df.loc[min_idx, 'FinalStructure']
write('best_structure.vasp', best_structure, format='vasp', direct=True)

# Create a filtered database (only converged)
converged_df = df[df['Converged'] == True]
with open('converged_only.pkl', 'wb') as f:
    pickle.dump(converged_df, f)
    """)

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("""
The vasp2db tool provides:
- Comprehensive VASP parameter extraction
- Full structure storage (initial + final) as ASE Atoms objects
- INCAR and POTCAR metadata for reproducibility
- Automatic calculation type classification
- Relative energy calculations
- Easy filtering and analysis with pandas

For more information:
- Run: vasp2db --help
- See: asetools/bin/vasp2db.py for implementation
    """)

except FileNotFoundError:
    print(f"\nDatabase file not found: {database_path}")
    print("Please create a database first using vasp2db command.")
    print("\nExample:")
    print("  vasp2db --paths folder1 folder2 --output database.pkl -r -v")

except Exception as e:
    print(f"\nError loading database: {e}")
