# File containing the functions to analyze OUTCARs and VASP runs

from ase.io import read
import pandas as pd
import numpy as np
    
def check_outcar_convergence(outcar, verbose=False):
    try:
        out = open(outcar, 'r')
    except:
        print('check_outcar_convergence --> OUTCAR file is damaged')
        return False
    lines = out.readlines()
    ibrion = None
    nsw = None
    opt = False
    convergence = False
    vasp = ''
    for line in lines:
        if 'vasp.6' in line:
            vasp = 'vasp6'
        if 'vasp.5' in line:
            vasp = 'vasp5'
        if 'IBRION' in line:
            ibrion = int( line.split()[2] )
            if ibrion == 1 or ibrion == 2 or ibrion == 3:
                opt = True
            elif verbose:
                print(f'IBRION --> {ibrion}, not an OPTIMIZATION job')
                opt = False
            else:    
                opt = False
        elif 'NSW' in line:
            nsw = int( line.split()[2] )
            if nsw > 0:
                opt = True
        elif 'reached required accuracy' in line and opt:
            convergence = True
        elif not opt and 'total amount of memory' in line:
            convergence = True
    if verbose:
        print(f'IBRION --> {ibrion}, NSW --> {nsw}')    
    
    if ibrion == 1 or ibrion == 2 or ibrion == 3:
        if nsw > 0 and convergence and verbose:
            print('Optimization Job --> CONVERGED')
        elif nsw > 0 and not convergence and verbose:
            print('Optimization Job --> *NOT* converged')
            
    return convergence, vasp

def check_energy_and_maxforce(outcar, magmom=False, verbose=False):

    convergence, vasp = check_outcar_convergence(outcar, verbose=verbose)
    try:
        atoms = read(outcar, format='vasp-out', index=-1)
        energy = atoms.get_potential_energy()
        vecforces = atoms.get_forces()
        forces = [np.linalg.norm(f) for f in vecforces]
        maxforce = max( forces )
        if magmom:
            mm = atoms.get_magnetic_moment() 
            return energy, maxforce, mm
        else:
            return energy, maxforce
    except:
        print('Missing or damaged OUTCAR file')
        return 9999.99, 9.99

def extract_magnetic_moments(outcar, listatoms, verbose=False):
    # listatoms: is a list with the indexes of atoms of interest
    #convergence, vasp = check_outcar_convergence(outcar, verbose=verbose)
    try:
        atoms = read(outcar, format='vasp-out', index=-1)
        #energy = atoms.get_potential_energy()
        mm = atoms.get_magnetic_moment()
        return [atoms.get_magnetic_moments()[i] for i in listatoms]
    except:
        print('Missing or damaged OUTCAR file')
        return []


    