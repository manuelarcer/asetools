#!/usr/bin/env python

import math
import sys
import argparse
from ase.io import read, write
import numpy as np
from ase.thermochemistry import HarmonicThermo, IdealGasThermo

def extract_vib_info(lines):
    vib = {}
    vasp6 = any('vasp.6' in line for line in lines)
    for line in lines:
        if ' f  = ' in line or 'f/i=' in line:
            index = line.split()[0]
            vibenergy = float(line.split()[-2]) / 1000
            freq = float(line.split()[-4])
            if ' f  = ' in line:
                vib[index] = {'e': vibenergy, 'freq': freq}
            elif 'f/i=' in line:
                vib[index] = {'e': vibenergy, 'freq': -freq}
    return vib, vasp6

def compute_corrections(vib, T):
    kB = 0.0000861733
    zpe = sum([v['e'] for v in vib.values() if v['freq'] > 0]) / 2.0
    cpT = sum([v['e'] / (math.exp(v['e'] / kB / T) - 1) for v in vib.values() if v['freq'] >= 0])
    S = sum([kB * (v['e'] / (kB * T * (math.exp(v['e'] / kB / T) - 1)) - math.log(1 - math.exp(-v['e'] / kB / T))) for v in vib.values() if v['freq'] >= 100])
    return zpe, cpT, S

def write_vib_files(vib, vasp6, lines, outcar):
    if vasp6:
        try:
            atoms = read('POSCAR', format='vasp')
        except:
            print('Could not read structure file')
            return
    else:
        atoms = read(outcar, format='vasp-out')

    nat = len(atoms)
    for i, line in enumerate(lines):
        if 'f  =' in line or 'f/i=' in line:
            name1 = line.split()[0]
            name2 = 'IMG' if 'f/i=' in line else 'f'
            shift = [list(map(float, lines[i+count+2].split()[3:])) for count in range(nat)]
            vib = []
            for k in np.arange(-1, 1, 0.2):
                atomscopy = atoms.copy()
                newshift = determineshift(shift, k)
                atomscopy.translate(newshift)
                vib.append(atomscopy)
                write(f'vib_{name1}_{name2}.traj', vib, format='traj')

def determineshift(shiftarray, k=1):
    return [[k * sh for sh in shiftat] for shiftat in shiftarray]

def is_close(a, b, threshold=0.001):
    return abs(a - b) < threshold

def remove_repeated_energies(energies, atoms, geom):  # It happens that for gas-phase I found repeated vib
    # energies, is the list of vibrational energies
    # atoms, is the atoms object
    # geom, is the geometry of the molecule : 'linear' or 'nonlinear'
    if geom == 'linear':
        numvib = 3 * len(atoms) - 5
    elif geom == 'nonlinear':
        numvib = 3 * len(atoms) - 6
    
    test = len(energies) == numvib
    if test:
        print('Number of energies as expected')
        return energies
    else:
        print()
        print('WARNING: The number of energies is different from the expected')
        print('Removing similar vibrational energies')
        diff = len(energies) - numvib
        invertenergies = [energies[i] for i in range(len(energies)-1, -1, -1)]
        print(invertenergies)
        newenergies = []
        if len(invertenergies) > numvib:
            for i, e in enumerate(invertenergies):
                if i == 0:
                    newenergies.append(e)
                if i > 0 and not is_close(e, invertenergies[i-1]):
                    newenergies.append(e)
        return newenergies

def main():
    parser = argparse.ArgumentParser(description='Calculate corrections to the energy for GFE.')
    parser.add_argument('outcar', type=str, help='OUTCAR file path')
    parser.add_argument('--writevib', type=str, choices=['y', 'n'], default='n', help='Whether to write vibrations (y/n)')
    parser.add_argument('--temp', type=float, default=298, help='Temperature in K')
    parser.add_argument('--gas', type=bool, default=False, help='Ideal Gas-phase Thermo')
    parser.add_argument('--symnum', type=int, default=2, help='Symmetry Number')
    parser.add_argument('--geom', type=str, choices=['monoatomic', 'linear', 'nonlinear'], default='linear', help='Geometry of molecule')
    parser.add_argument('--pressure', type=float, default=1.0, help='Gas pressure (bar)')
    args = parser.parse_args()

    Pa = 100000.    # 1 bar in Pa

    lines = open(args.outcar, 'r').readlines()
    vib, vasp6 = extract_vib_info(lines)
    zpe, cpT, S = compute_corrections(vib, args.temp)

    ## Display results
    print ('-----------------------')
    print ('{:>3} {:>12} {:>6}'.format('#', 'Freq[cm-1]', 'E[eV]'))
    print ('-----------------------')
    for key, value in vib.items():
        print ('{:>3} {:12.1f} {:6.3f}'.format(key, value['freq'], value['e']))
    print ('-----------------------')

    e_dft = read(args.outcar, format='vasp-out', index=0).get_potential_energy()
    print('E, eV = {:.3f}'.format(e_dft))
    print('ZPE, eV = {:.3f}'.format(zpe))
    print('S, eV/K = {:.6f}'.format(S))
    print('Temperature is = {:.1f}'.format(args.temp))
    print('Thermal correction (0->T), eV = {:.3f}'.format(cpT))
    print('Entropy correction (-S*T), eV = {:.3f}'.format(-S * args.temp))
    print('\n')
    print('All values together: E_tot   E_ZPE   CpT   -S*T')
    print('{:.3f} {:.3f} {:.3f} {:.3f}'.format(e_dft, zpe, cpT, -S * args.temp))
    print('\n')

    print('*************************************************')
    print('---------------  Harmonic Limit  ----------------')
    print('*************************************************')
    print()
    energies = [value['e'] for value in vib.values()]
    harm_lim = HarmonicThermo(energies, potentialenergy=0.0)
    # The following line already gives the details of the S, CpT and G energy calculations
    harm_lim.get_helmholtz_energy(args.temp, verbose=True)


    ##### Ideal Gas Thermo
    if args.gas:
        energies = [value['e'] for value in vib.values() if value['freq'] >= 0]
        print(energies)
        atoms = read(args.outcar, format='vasp-out', index=0)
        energies = remove_repeated_energies(energies, atoms, args.geom)
        thermo = IdealGasThermo(vib_energies=energies,
                        atoms=atoms,
                        #potentialenergy=atoms.get_potential_energy(),
                        geometry=args.geom,      # monoatomic, linear, nonlinear
                        symmetrynumber=args.symnum,   # CO2: 2, H2O: 2, CO: 1, H2: 2
                        spin=0)     # Different for radicals or unpaired electrons            
        G = thermo.get_gibbs_energy(temperature=args.temp, pressure=args.pressure * Pa)
    
    if args.writevib == 'y':
        write_vib_files(vib, vasp6, lines, args.outcar)

if __name__ == "__main__":
    main()
