# TODO
# File containing the functions to analyze OUTCARs and VASP runs

from ase.io import read
import pandas as pd
import numpy as np
    
def check_outcar_convergence(outcar, verbose=True):
    out = open(outcar, 'r')
    lines = out.readlines()
    ibrion = None
    nsw = None
    opt = False
    convergence = False
    for line in lines:
        if 'IBRION' in line:
            ibrion = int( line.split()[2] )
            if ibrion == 1 or ibrion == 2:
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
        elif 'reached required accuracy - stopping' in line and opt:
            convergence = True
    if verbose:
        print(f'IBRION --> {ibrion}, NSW --> {nsw}')    
    
    if ibrion == 1 or ibrion == 2:
        if nsw > 0 and convergence and verbose:
            print('Optimization Job --> CONVERGED')
        elif nsw > 0 and not convergence and verbose:
            print('Optimization Job --> *NOT* converged')
            
    return convergence

def check_energy_and_maxforce(outcar, magmom=False, verbose=True):
    ## TODO need to add Total MagMom functionally for output

    check_outcar_convergence(outcar, verbose=verbose)
    try:
        atoms = read(outcar, format='vasp-out', index=-1)
    except:
        print('Missing or damaged OUTCAR file')
        pass

    energy = atoms.get_potential_energy()
    vecforces = atoms.get_forces()
    forces = [np.linalg.norm(f) for f in vecforces]
    maxforce = max( forces )
    return energy, maxforce
    