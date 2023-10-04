#!/usr/bin/env python

import os
import shutil
import argparse
import numpy as np
from ase.io import read

VALANCES = {'Cu': 11, 'Zn': 12, 'C': 4, 'O': 6, 'H': 1}
EPSILON = 0.001

def get_sum_electrons(poscar):
    atoms = read(poscar)
    symbols = atoms.get_chemical_symbols()
    elecpersymb = [VALANCES[symb] for symb in symbols]
    return sum(elecpersymb)

def modify_python_file(pythonfile, ne):
    with open(pythonfile, 'r') as file:
        content = file.readlines()
        for k, line in enumerate(content):
            if "nelect =" in line:
                parts = line.split('nelect =')
                content[k] = f"{parts[0]}nelect = {ne:.2f} # Modified value\n"
                break

    with open(pythonfile, 'w') as file:
        file.writelines(content)

def main():
    parser = argparse.ArgumentParser(description='Preparation of folders for Potential-dependent calculations')
    parser.add_argument('poscar', help='Path to POSCAR file')
    parser.add_argument('pythonfile', help='Script used for the calculations')
    parser.add_argument('--lower', type=float, default=-1, help='Lower limit for the number of electrons')
    parser.add_argument('--higher', type=float, default=1, help='Higher limit for the number of electrons')
    parser.add_argument('--step', type=float, default=0.2, help='Step value for electron count')

    args = parser.parse_args()

    numelec = get_sum_electrons(args.poscar)
    nelec_min = numelec + args.lower
    nelec_max = numelec + args.higher

    print(f'Num_elec (Neutral) = {numelec} ; lower = {nelec_min} ; higher = {nelec_max}')

    for i, ne in enumerate(np.arange(nelec_min, nelec_max + EPSILON, args.step)):
        name = f'{i:02d}_nelec_{ne:.2f}'
        os.mkdir(name)
        shutil.copyfile(args.pythonfile, os.path.join(name, args.pythonfile))
        shutil.copyfile(args.poscar, os.path.join(name, args.poscar))

        os.chdir(name)
        modify_python_file(args.pythonfile, ne)
        os.chdir('..')

if __name__ == "__main__":
    main()
