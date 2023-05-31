#!/usr/bin/env python

from ase.io import read
import pandas as pd
import numpy as np
import glob, os
from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence

def main():
    folders = glob.glob('*/')

    dic = {'Config': [], 'Converged':[], 'MaxForce': [], 'Energy':[], 'MagMom':[]}
    for f in sorted(folders):
        if os.path.exists(f+'OUTCAR'):
            converged = check_outcar_convergence(f+'OUTCAR', verbose=False)
            energy, maxforce = check_energy_and_maxforce(f+'OUTCAR', magmom=True, verbose=False)
            dic['Config'].append(f)
            dic['Converged'].append(converged)
            dic['MaxForce'].append(maxforce)
            dic['Energy'].append(energy)

        else:
            print('No OUTCAR in ', f)
    
    df = pd.DataFrame.from_dict(dic) 
    print(df)

if __name__ == "__main__":
    main()
