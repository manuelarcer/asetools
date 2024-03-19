#!/usr/bin/env python

from ase.io import read, write
from ase import Atoms
from ase.constraints import FixAtoms

def get_fixed_indices(atoms):
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            return constraint.index
    return []

def sort_atoms(filename, element_order):
    atoms = read(filename)
    fixed_indices = get_fixed_indices(atoms)

    sorted_atoms = Atoms(cell=atoms.cell, pbc=atoms.pbc)
    index_map = {}  # Maps old indices to new indices
    new_fixed_indices = []

    # Sort the atoms and track index changes
    for element in element_order:
        for i, atom in enumerate(atoms):
            if atom.symbol == element:
                sorted_atoms.append(atom)
                index_map[i] = len(sorted_atoms) - 1

    # Update the fixed indices based on the sorted atoms
    for old_index in fixed_indices:
        new_fixed_indices.append(index_map[old_index])

    if new_fixed_indices:
        sorted_atoms.set_constraint(FixAtoms(indices=new_fixed_indices))

    return sorted_atoms

def main():
    import sys
    if len(sys.argv) < 3:
        print("Usage: script.py <filename> <element1> <element2> ...")
        sys.exit(1)

    filename = sys.argv[1]
    element_order = sys.argv[2:]

    sorted_atoms = sort_atoms(filename, element_order)

    # Save the sorted structure with preserved constraints
    write('sorted_POSCAR', sorted_atoms, format='vasp')

if __name__ == "__main__":
    main()

