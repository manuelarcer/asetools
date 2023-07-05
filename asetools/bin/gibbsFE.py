#!/usr/bin/env python

# This script can extract the relevant information to calculate the corrections
# to the energy for GFE
# Harmonic limit approximation. Just like it is implemented in ASE
# https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html

import math, sys
from ase.io import read, write
import numpy as np
from ase.thermochemistry import HarmonicThermo

def determineshift(shiftarray, k=1 ):
        newshift = []
        for shiftat in shiftarray:
                newshift.append( [ k * sh for sh in shiftat ] )
        return newshift


print()
print('This should be a folder with a finite differences calculations (Vib Analysis)')
print('USAGE:   gibbsFE.py   OUTCAR_FILE   WRITE_VIB   TEMP')
print('WRITE_VIB is either y or n \nTEMP in K (298 K is the default)')
print()

outcar = sys.argv[1]
writevib = sys.argv[2]

if len(sys.argv) == 4:
    T = float( sys.argv[3] )
else:
    T = 298         # K

kB = 0.0000861733   # eV/K

output = open(outcar, 'r')
lines = output.readlines()
vib = {}
zpe = 0 ; cpT = 0 ; S = 0
vasp6 = False
finished = False
for line in lines:
    if 'General timing and accounting informations' in line:
        finished = True
    if 'vasp.6' in line:
        vasp6 = True

if finished:
    for line in lines:
        if ' f  = ' in line:
            index = line.split()[0]
            vib[index] = {}
            vibenergy = float( line.split()[-2] ) / 1000  # this is now in eV
            freq = float( line.split()[-4] )
            vib[index]['e'] = vibenergy
            vib[index]['freq'] = freq
            zpe += vibenergy        # For ZPE energy, we need to divide by 2 at the end
            cpT += vibenergy / ( math.exp( vibenergy / kB / T ) - 1)
            S += kB * ( vibenergy / ( kB*T * ( math.exp(vibenergy/kB/T)-1 )  ) - math.log( 1 - math.exp(-vibenergy/kB/T
) ) )
        elif 'f/i=' in line:
            print('there is an imaginary frequency')
            index = line.split()[0]
            vib[index+'i'] = {}
            vibenergy = float( line.split()[-2] ) / 1000  # this is now in eV
            freq = float( line.split()[-4] )
            vib[index+'i']['e'] = vibenergy
            vib[index+'i']['freq'] = freq
else:
    print('Finite differences did not complete!! Please check')
    sys.exit()

print()
#results = open('GFEcorrections.log', 'a')
print ('-----------------------')
print ('{:>3} {:>12} {:>6}'.format('#', 'Freq[cm-1]', 'E[eV]'))
print ('-----------------------')
for i in vib:
    #results.write('{:.0} {:.1f} {:.3f}'.format(vib[i], vib['freq'], vib['e']))
    print ('{:>3} {:12.1f} {:6.3f}'.format(i, vib[i]['freq'], vib[i]['e']))
print ('-----------------------')
#print('\n')
print('ZPE, eV = {:.3f}'.format( zpe / 2. ))
print('S, eV/K = {:.6f}'.format( S ))
print('Temperature is = {:.1f}'.format( T ) )
print('Thermal correction (0->T), eV = {:.3f}'.format( cpT ) )
print('Entropy correction (-S*T), eV = {:.3f}'.format( -S * T ) )
print('\n')
print('Three together {:.3f} {:.3f} {:.3f}'.format(zpe / 2., cpT, -S * T))
print('\n')

print('*************************************************')
print('---------------  Harmonic Limit  ----------------')
print('*************************************************')
print()
energies = []
for key in vib.keys():
      energies.append( vib[key]['e'] )
harm_lim = HarmonicThermo(energies, potentialenergy=0.0)
# The following line already gives the details of the S, CpT and G energy calculations
harm_lim.get_helmholtz_energy(T, verbose=True)

### WRITE VIB section

if writevib == 'y' or writevib == 'Y' or writevib == 'yes':
	print('Writing vibrational frequencies to traj files...')
	if vasp6:
		try:
			atoms = read('POSCAR', format='vasp')
		except:
			print('Could not read structure file')
	else:
		atoms = read(outcar, format='vasp-out')
	nat = len(atoms)
	for i, line in enumerate(lines):
		count = 0
		if 'f  =' in line or 'f/i=' in line:
			name1 = line.split()[0]
			if 'f/i=' in line:
				name2 = 'IMG'
			else:
				name2 = 'f'
			shift = []
			while count < nat:
				#shift = np.array( [ float(sh) for sh in lines[i+count+2].split()[3:] ] )
				shift.append( [ float(sh) for sh in lines[i+count+2].split()[3:] ] )
				count += 1
			vib = []
			for k in np.arange(-1,1,0.2):
				atomscopy = atoms.copy()
				newshift = determineshift( shift, k )
				atomscopy.translate(newshift)
				vib.append(atomscopy)
			write('vib_'+name1+'_'+name2+'.traj', vib , format='traj')