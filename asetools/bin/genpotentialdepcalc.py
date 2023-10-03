#!/usr/bin/env python

from ase.io import read
import os, shutil, sys
import numpy as np

def get_sum_electrons(poscar):
    valences = {'Cu':11, 'Zn':12, 'C':4, 'O':6, 'H':1}
    atoms = read(poscar)
    symbols = atoms.get_chemical_symbols()
    elecpersymb = [valences[symb] for symb in symbols ]
    return sum(elecpersymb)

print()
print('Preparation of folders for Potential-dependent calculations')
print('USAGE:   genporentialdepcalc.py   POSCAR   pythonfile.py   lower_n_e(OPT)   higher_n_e(OPT)  step(OPT)')
print('pythonfile.py, is the script used for the calculations, all ase-related lines for the calculations are here')
print('lower_n_e, represent the lower limit for the number of electrons added to the neutral system')
print()

## Parameters ########
poscar = sys.argv[1]
pythonfile = sys.argv[2]

if len(sys.argv) == 3:
    lower = -1
    higher = 1
    step = 0.2
elif len(sys.argv) == 4:
    lower = float( sys.argv[3] )
    higher = -lower
    step = 0.2
elif len(sys.argv) == 5:
    lower = float( sys.argv[3] )
    higher = float( sys.argv[4] )
    step = 0.2
elif len(sys.argv) == 6:
    lower = float( sys.argv[3] )
    higher = float( sys.argv[4] )
    step = float( sys.argv[5] )
else:
    print('Please verify the input parameters')
######################

numelec = get_sum_electrons(poscar)
nelec_min = numelec + lower
nelec_max = numelec + higher

print(f'Num_elec (Neutral) = {numelec} ; lower = {nelec_min} ; higher = {nelec_max}')

for i, ne in enumerate( np.arange(nelec_min, nelec_max + 0.001, step) ):
    
    name = f'{i:02d}_nelec_{ne:.2f}'
    os.mkdir(name)
    shutil.copyfile(pythonfile, name+'/'+pythonfile)

    os.chdir(name)
    with open(pythonfile, 'r') as file:
        content = file.readlines()
        for k, line in enumerate(content):
            if "nelect =" in line:
                parts = line.split('nelect =')
                content[k] = f"{parts[0]}nelect = {ne:.2f} # Modified value\n"
                break

    with open(pythonfile, 'w') as file:
        file.writelines(content)

    os.chdir('..')