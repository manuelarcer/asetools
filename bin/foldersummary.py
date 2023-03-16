#!/usr/bin/env python

from ase.io import read
import pandas as pd
import numpy as np
import glob, os
from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence


folders = glob.glob('*/')

dic = {'Config': [], 'Converged':[], 'MaxForce': [], 'Energy':[]}
for f in folders:
    if os.path.exists(f+'OUTCAR'):
        converged = check_outcar_convergence(f+'OUTCAR')
        energy, maxforce = check_energy_and_maxforce(f+'OUTCAR', magmom=False, verbose=False)
        dic['Config'].append(f)
        dic['Converged'].append(converged)
        dic['MaxForce'].append(maxforce)
        dic['Energy'].append(energy)

    else:
        print('No OUTCAR in ', f)
    
df = pd.DataFrame.from_dict(dic) 
print(df)