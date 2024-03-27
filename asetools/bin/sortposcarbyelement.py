#!/usr/bin/env python

import sys
from ase.io import read, write
from ase import Atoms
from ase.constraints import FixAtoms

def get_fixed_indices(atoms):
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            return constraint.index
    return []

def sort_atoms(filename, element_order, sort_axis='z', ascending=True):
    atoms = read(filename)
    fixed_indices = get_fixed_indices(atoms)

    # Sorting by element type first
    sorted_atoms_by_type = []
    for element in element_order:
        sorted_atoms_by_type.extend([atom for atom in atoms if atom.symbol == element])

    # Sorting within each type by specified axis
    axis_indices = {'x': 0, 'y': 1, 'z': 2}
    axis_index = axis_indices.get(sort_axis, 2)  # Default to z-axis if unspecified
    sorted_atoms_by_type_and_axis = []
    for element in element_order:
        atoms_of_element = [atom for atom in sorted_atoms_by_type if atom.symbol == element]
        sorted_atoms_of_element = sorted(atoms_of_element, key=lambda atom: atom.position[axis_index], reverse=not ascending)
        sorted_atoms_by_type_and_axis.extend(sorted_atoms_of_element)

    # Create new Atoms object with sorted atoms
    new_atoms = Atoms(cell=atoms.cell, pbc=atoms.pbc)
    for atom in sorted_atoms_by_type_and_axis:
        new_atoms.append(atom)

    # Update the fixed indices based on the new sorted order
    new_fixed_indices = []
    for old_index in fixed_indices:
        for new_index, atom in enumerate(new_atoms):
            if atom.index == old_index:
                new_fixed_indices.append(new_index)
                break

    if new_fixed_indices:
        new_atoms.set_constraint(FixAtoms(indices=new_fixed_indices))

    return new_atoms

def main():
    if len(sys.argv) < 4:
        print("Usage: script.py <filename> <element1> <element2> ... [-axis x|y|z] [-desc]")
        sys.exit(1)

    filename = sys.argv[1]
    element_order = sys.argv[2:]

    # Detect if sorting by axis is specified
    sort_axis = 'z'  # Default axis
    ascending = True  # Default order
    if '-axis' in element_order:
        axis_index = element_order.index('-axis')
        sort_axis = element_order[axis_index + 1].lower()
        element_order = element_order[:axis_index]  # Remove axis options from element_order list
        if '-desc' in element_order:
            ascending = False
            element_order.remove('-desc')

    sorted_atoms = sort_atoms(filename, element_order, sort_axis, ascending)

    # Save the sorted structure with preserved constraints
    write('sorted_POSCAR', sorted_atoms, format='vasp')

if __name__ == "__main__":
    main()



