#!/usr/bin/env python

"""
Reorder atoms in a configuration file by chemical symbol and optionally by z-coordinate.

This script reads an atomic structure file (POSCAR, CIF, XYZ, etc.) and reorders
the atoms according to user-specified chemical symbol ordering and/or z-coordinate.
"""

import argparse
import sys
from ase.io import read, write
from ase import Atoms
from ase.constraints import FixAtoms


def get_fixed_indices(atoms):
    """Extract indices of fixed atoms from constraints."""
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            return constraint.index
    return []


def reorder_atoms(atoms, element_order=None, z_order=None):
    """
    Reorder atoms by chemical symbol and optionally by z-coordinate.

    Parameters
    ----------
    atoms : ase.Atoms
        Input atomic structure
    element_order : list of str, optional
        Order of chemical symbols (e.g., ['Cu', 'O', 'H']).
        If None, uses alphabetical order.
    z_order : str, optional
        Secondary ordering by z-coordinate: 'top-bottom' (descending z)
        or 'bottom-top' (ascending z). If None, no z-ordering is applied.

    Returns
    -------
    ase.Atoms
        Reordered atomic structure with preserved constraints
    """
    fixed_indices = get_fixed_indices(atoms)

    # Determine element ordering
    if element_order is None:
        # Default: alphabetical order
        element_order = sorted(set(atoms.get_chemical_symbols()))
    else:
        # Validate that all elements in structure are in element_order
        unique_symbols = set(atoms.get_chemical_symbols())
        ordered_symbols = set(element_order)
        missing = unique_symbols - ordered_symbols
        if missing:
            print(f"Warning: Elements {missing} in structure but not in element_order. "
                  f"They will be appended at the end.")
            # Append missing elements in alphabetical order
            element_order = element_order + sorted(list(missing))

    # Group atoms by element
    sorted_atoms_list = []
    for element in element_order:
        atoms_of_element = [atom for atom in atoms if atom.symbol == element]

        # Apply z-ordering within each element group if requested
        if z_order is not None and len(atoms_of_element) > 0:
            if z_order == 'top-bottom':
                # Sort by z descending (highest z first)
                atoms_of_element = sorted(atoms_of_element,
                                         key=lambda atom: atom.position[2],
                                         reverse=True)
            elif z_order == 'bottom-top':
                # Sort by z ascending (lowest z first)
                atoms_of_element = sorted(atoms_of_element,
                                         key=lambda atom: atom.position[2],
                                         reverse=False)

        sorted_atoms_list.extend(atoms_of_element)

    # Create new Atoms object with sorted atoms
    new_atoms = Atoms(cell=atoms.cell, pbc=atoms.pbc)
    for atom in sorted_atoms_list:
        new_atoms.append(atom)

    # Preserve constraints by mapping old indices to new indices
    if len(fixed_indices) > 0:
        # Create mapping from old index to new index
        index_mapping = {}
        for old_idx, atom in enumerate(atoms):
            for new_idx, new_atom in enumerate(new_atoms):
                # Match by position and symbol
                if (atom.symbol == new_atom.symbol and
                    all(abs(atom.position - new_atom.position) < 1e-6)):
                    index_mapping[old_idx] = new_idx
                    break

        new_fixed_indices = [index_mapping[old_idx] for old_idx in fixed_indices
                            if old_idx in index_mapping]

        if new_fixed_indices:
            new_atoms.set_constraint(FixAtoms(indices=new_fixed_indices))

    return new_atoms


def main():
    parser = argparse.ArgumentParser(
        description='Reorder atoms in a configuration file by chemical symbol and z-coordinate.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Reorder by alphabetical order of elements
  %(prog)s POSCAR

  # Reorder with custom element order
  %(prog)s POSCAR --order Cu O H

  # Reorder alphabetically, then by z-coordinate (top to bottom)
  %(prog)s POSCAR --z-order top-bottom

  # Custom element order + z-ordering within each element
  %(prog)s POSCAR --order Cu O H --z-order bottom-top

  # Specify output file name
  %(prog)s POSCAR --order Cu O --output reordered_POSCAR

  # Read CIF, write VASP format
  %(prog)s structure.cif --format vasp --output POSCAR
        """
    )

    parser.add_argument(
        'input_file',
        help='Input structure file (POSCAR, CIF, XYZ, etc.)'
    )

    parser.add_argument(
        '--order', '-o',
        nargs='+',
        metavar='ELEMENT',
        help='Order of chemical symbols (e.g., Cu O H). Default: alphabetical'
    )

    parser.add_argument(
        '--z-order', '-z',
        choices=['top-bottom', 'bottom-top'],
        help='Secondary ordering by z-coordinate within each element group. '
             'top-bottom: highest z first; bottom-top: lowest z first'
    )

    parser.add_argument(
        '--output', '-out',
        default='reordered_POSCAR',
        help='Output file name (default: reordered_POSCAR)'
    )

    parser.add_argument(
        '--format', '-f',
        default='vasp',
        help='Output format (default: vasp). See ASE documentation for supported formats.'
    )

    args = parser.parse_args()

    # Read input structure
    try:
        atoms = read(args.input_file)
    except Exception as e:
        print(f"Error reading file '{args.input_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # Reorder atoms
    reordered_atoms = reorder_atoms(atoms,
                                    element_order=args.order,
                                    z_order=args.z_order)

    # Write output
    try:
        write(args.output, reordered_atoms, format=args.format)
        print(f"Reordered structure written to '{args.output}'")

        # Print summary
        print(f"\nAtom ordering summary:")
        element_order = args.order if args.order else sorted(set(atoms.get_chemical_symbols()))
        for element in element_order:
            count = sum(1 for atom in reordered_atoms if atom.symbol == element)
            if count > 0:
                print(f"  {element}: {count} atoms")

        if args.z_order:
            print(f"\nZ-ordering: {args.z_order}")

    except Exception as e:
        print(f"Error writing file '{args.output}': {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
