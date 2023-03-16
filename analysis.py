# TODO
# File containing the functions to analyze OUTCARs and VASP runs

from ase.io import read
import pandas as pd
    
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
        if nsw > 0 and convergence:
            print('Optimization Job --> CONVERGED')
        elif nsw > 0 and not convergence and verbose:
            print('Optimization Job --> *NOT* converged')
            
    return convergence

def check_status_calc(outcar):
    
    check_outcar_convergence(outcar, verbose=True)
    try:
        atoms = read(outcar, format='vasp-out', index=-1)
    except:
        print('Missing or damaged OUTCAR file')
        pass
    
    return max( atoms.get_forces() )
    