#!/usr/bin/env python
"""
vasp2db: Extract VASP calculation data to pandas DataFrame (pickle format)

Scans VASP calculation directories, extracts comprehensive metadata (structures,
parameters, INCAR, POTCAR info), and saves to a pickle file containing a pandas
DataFrame with ASE Atoms objects.

Usage:
    vasp2db --paths folder1 folder2 folder3 --output database.pkl
    vasp2db --paths */ --output all_calcs.pkl -r -v
"""

import argparse
import os
import sys
import pickle
import pandas as pd
from pathlib import Path
from ase.io import read

# Add parent directory to path to import asetools
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from asetools.analysis import (
    extract_comprehensive_metadata,
    find_initial_structure,
    classify_calculation_type
)


def process_single_calc(path, initial_pattern, verbose):
    """
    Process a single calculation directory.

    Args:
        path: Path to calculation directory (relative or absolute)
        initial_pattern: Glob pattern for initial structure files
        verbose: Print progress messages

    Returns:
        dict: Record with all extracted data

    Raises:
        FileNotFoundError: If OUTCAR not found
        Exception: For other errors during processing
    """
    # Convert to absolute path
    abs_path = os.path.abspath(path)

    # Check OUTCAR exists
    outcar_path = os.path.join(abs_path, 'OUTCAR')
    if not os.path.exists(outcar_path):
        raise FileNotFoundError(f"OUTCAR not found in {path}")

    # Extract comprehensive metadata
    if verbose:
        print(f"  Extracting metadata from OUTCAR...")
    metadata = extract_comprehensive_metadata(outcar_path)

    # Find and load initial structure
    initial_atoms = None
    try:
        initial_struct_path = find_initial_structure(abs_path, initial_pattern)
        initial_atoms = read(initial_struct_path)
        if verbose:
            print(f"  Initial: {os.path.basename(initial_struct_path)}")
    except FileNotFoundError as e:
        if verbose:
            print(f"  WARNING: {e}")
    except Exception as e:
        if verbose:
            print(f"  WARNING: Could not read initial structure: {e}")

    # Load final structure from OUTCAR
    final_atoms = None
    try:
        final_atoms = read(outcar_path, format='vasp-out', index=-1)
        if verbose and metadata.get('Energy') is not None:
            print(f"  Final: OUTCAR (E={metadata['Energy']:.3f} eV)")
        elif verbose:
            print(f"  Final: OUTCAR")
    except Exception as e:
        if verbose:
            print(f"  WARNING: Could not read final structure: {e}")

    # Build record
    record = {
        'Path': path,
        'AbsPath': abs_path,
        'InitialStructure': initial_atoms,
        'FinalStructure': final_atoms,
        **metadata  # Unpack all metadata fields
    }

    return record


def build_database(paths, initial_pattern='*.vasp', verbose=False, skip_errors=False):
    """
    Build pandas DataFrame from VASP calculation directories.

    Args:
        paths: List of directory paths (relative or absolute)
        initial_pattern: Glob pattern for initial structures
        verbose: Print progress messages
        skip_errors: Continue on errors vs stop

    Returns:
        pd.DataFrame: DataFrame with all extracted data
    """
    records = []
    skipped = []

    for i, path in enumerate(paths):
        if verbose:
            print(f"[{i+1}/{len(paths)}] Processing: {path}")

        try:
            record = process_single_calc(path, initial_pattern, verbose)
            records.append(record)
        except Exception as e:
            error_msg = f"{path}: {e}"
            if skip_errors:
                print(f"  ERROR: {error_msg}")
                skipped.append(error_msg)
                continue
            else:
                raise RuntimeError(f"Error processing {path}: {e}") from e

    # Print summary of skipped folders
    if skipped:
        print(f"\n{'='*60}")
        print(f"Skipped {len(skipped)} folder(s) due to errors:")
        for err in skipped:
            print(f"  - {err}")
        print(f"{'='*60}\n")

    df = pd.DataFrame(records)
    return df


def add_relative_energies(df):
    """
    Add Rel.E column with energies relative to minimum.

    Args:
        df: DataFrame with 'Energy' column

    Returns:
        pd.DataFrame: DataFrame with added 'Rel.E' column
    """
    energies = df['Energy'].values

    # Filter out None or invalid values
    valid_energies = [e for e in energies if e is not None and isinstance(e, (int, float))]

    if valid_energies:
        min_energy = min(valid_energies)
        df['Rel.E'] = df['Energy'].apply(
            lambda x: x - min_energy if (x is not None and isinstance(x, (int, float))) else None
        )
    else:
        df['Rel.E'] = None

    return df


def save_database(df, output_path, verbose):
    """
    Save DataFrame to pickle file.

    Args:
        df: DataFrame to save
        output_path: Output pickle file path
        verbose: Print file size information
    """
    with open(output_path, 'wb') as f:
        pickle.dump(df, f, protocol=pickle.HIGHEST_PROTOCOL)

    if verbose:
        file_size = os.path.getsize(output_path) / (1024**2)  # MB
        print(f"\nSaved database to: {output_path}")
        print(f"File size: {file_size:.2f} MB")


def print_summary(df):
    """
    Print summary statistics of database.

    Args:
        df: DataFrame to summarize
    """
    print(f"\n{'='*60}")
    print(f"Database Summary")
    print(f"{'='*60}")
    print(f"Total calculations: {len(df)}")

    # Calculation types
    if 'CalcType' in df.columns:
        print(f"\nCalculation types:")
        calc_types = df['CalcType'].value_counts()
        for calc_type, count in calc_types.items():
            print(f"  {calc_type}: {count}")

    # Convergence status
    if 'Converged' in df.columns:
        print(f"\nConvergence status:")
        conv_counts = df['Converged'].value_counts()
        for status, count in conv_counts.items():
            print(f"  {status}: {count}")

    # Energy range
    if 'Energy' in df.columns:
        valid_energies = df['Energy'].dropna()
        if len(valid_energies) > 0:
            print(f"\nEnergy range:")
            print(f"  Min: {valid_energies.min():.3f} eV")
            print(f"  Max: {valid_energies.max():.3f} eV")
            print(f"  Range: {valid_energies.max() - valid_energies.min():.3f} eV")

    # Relative energies if present
    if 'Rel.E' in df.columns:
        valid_rel_e = df['Rel.E'].dropna()
        if len(valid_rel_e) > 0:
            print(f"\nRelative energy range:")
            print(f"  Max: {valid_rel_e.max():.3f} eV")

    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Extract VASP calculation data to pandas DataFrame (pickle format)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  vasp2db --paths CuCu_2 NiCu_1 NiCu_2 --output database.pkl
  vasp2db --paths */ --output all_calcs.pkl --relative-energy -v
  vasp2db --paths calc1 calc2 --output db.pkl --skip-errors

Loading the database:
  import pandas as pd
  import pickle

  with open('database.pkl', 'rb') as f:
      df = pickle.load(f)

  # Access structures
  initial = df.loc[0, 'InitialStructure']  # ASE Atoms object
  final = df.loc[0, 'FinalStructure']

  # View summary
  print(df[['Path', 'CalcType', 'Energy', 'Converged', 'Rel.E']])
        """
    )

    parser.add_argument(
        '--paths',
        nargs='+',
        required=True,
        help='List of calculation directories (relative or absolute paths)'
    )

    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output pickle file path (e.g., database.pkl)'
    )

    parser.add_argument(
        '--initial-pattern',
        default='*.vasp',
        help='Glob pattern for initial structure files (default: *.vasp)'
    )

    parser.add_argument(
        '--relative-energy', '-r',
        action='store_true',
        help='Calculate relative energies (subtract minimum)'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output showing progress'
    )

    parser.add_argument(
        '--skip-errors',
        action='store_true',
        help='Skip folders with errors instead of stopping'
    )

    args = parser.parse_args()

    # Build database
    print(f"Scanning {len(args.paths)} folder(s)...")
    df = build_database(
        args.paths,
        args.initial_pattern,
        args.verbose,
        args.skip_errors
    )

    if len(df) == 0:
        print("ERROR: No calculations successfully processed!")
        sys.exit(1)

    # Calculate relative energies if requested
    if args.relative_energy:
        if args.verbose:
            print("\nCalculating relative energies...")
        df = add_relative_energies(df)

    # Save to pickle
    save_database(df, args.output, args.verbose)

    # Print summary
    print_summary(df)

    print(f"Successfully processed {len(df)} calculation(s).")


if __name__ == "__main__":
    main()
