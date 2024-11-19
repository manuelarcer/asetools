#!/usr/bin/env python

from ase.io import read, write
import sys

# User input is the VASP (POSCAR/CONTCAR) file
input = sys.argv[1]

name = input.split('.')[0]

atoms = read(input)
write(f'{name}.arc', atoms, format='dmol-arc')