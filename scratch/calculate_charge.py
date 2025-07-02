#!/usr/bin/env python

import sys, os
from ase import Atoms
from ase.io import read
import numpy as np

# This script calculate the charges of all atoms doing this operation:
# CH = CHG_bader - ZVAL

atoms = read('CONTCAR',format='vasp')
symbols = atoms.get_chemical_symbols()
coord = atoms.get_positions()
# Number of different atoms:
ndat = len(np.unique(symbols))

if not os.path.exists('ACF.dat'):
        print( 'ACF.dat file does not exist, calculating bader charge...' )
        os.system('chgsum.pl AECCAR0 AECCAR2')
        os.system('bader -vac off CHGCAR -ref CHGCAR_sum')

print('Extracting bader charge from ACF.dat file and calculating final charge...')
# Determine the Zval of the different types of atoms
outcar = open('OUTCAR', 'r')
outlines = outcar.readlines()
elem = []
for i, line in enumerate(outlines):
        if 'POTCAR:' in line:
                if len(elem) < ndat:
                        elem.append(line.split()[2].split('_')[0])
                #for k in range(ndat):
                        # Double split here in case the element potential is given as Ti_pv, for example.
                #        print(outlines[i+1+k].split())
                #        elem.append(outlines[i+1+k].split()[2].split('_')[0])
        if 'Ionic Valenz' in line:
                e_line = outlines[i+1]
                words = e_line.split()
                e_val = [float(i) for i in words[2:]]


# The dictionary containing the Z valance of all different elements
ZVAL = {}
for i, el in enumerate(elem):
        ZVAL[el] = e_val[i]

# Now saving the Zval of all atoms in model:
#init_chg = np.zeros(len(symbols))
#for i, el in enumerate(symbols):
#       init_chg[i] = ZVAL[el]

# Determine charge from ACF file of all atoms
acf = open('ACF.dat','r')
acflines = acf.readlines()
chg = np.zeros(len(symbols))
for i, line in enumerate(acflines[2:len(symbols)+2]):
        chg[i] = line.split()[4]

print(ZVAL)
# FINALLY.
# Determine and write the final charge of all atoms in model
print('Writing resulting_charges.txt ...')
results = open('resulting_charges.txt','w')
results.write('#\tSymbol\tX\tY\tZ\tBader_Chg\tFinal_Chg\n')
for i, el in enumerate(symbols):
        results.write(str(i)+'\t'+el+'\t'+str(coord[i,0])+'\t'+
                str(coord[i,1])+'\t'+str(coord[i,2])+'\t'+
                str(chg[i])+'\t'+str(ZVAL[el]-chg[i])+'\n')
results.close()
print('Done!!')

