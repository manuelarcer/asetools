#!/usr/bin/env python

import argparse
from ase.io import read
import pandas as pd
import numpy as np
import glob, os
from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence

def main():

    ## Block from chatgpt
    parser = argparse.ArgumentParser(description='Process folders.')
    parser.add_argument('-m', '--magmom', action='store_true',
                        help='extract and present magnetic moments')
    args = parser.parse_args()
    ##

    folders = glob.glob('*/')
    if args.magmom:
        dic = {'Config': [], 'Converged':[], 'MaxForce': [], 'Energy':[], 'MagMom':[]}
    else:
        dic = {'Config': [], 'Converged':[], 'MaxForce': [], 'Energy':[]}
    for f in sorted(folders):
        #foundoutcar = False
        #if os.path.exists(f+'OUTCAR'):
        #    outcar = f+'OUTCAR'
        #    foundoutcar = True
        #elif os.path.exists(f+'OUTCAR.gz'):
        #    outcar = f+'OUTCAR.gz'
        #    foundoutcar = True
        if os.path.exists(f+'OUTCAR'):
            converged = check_outcar_convergence(f+'OUTCAR', verbose=False)
            if args.magmom:
                energy, maxforce, magmom = check_energy_and_maxforce(f+'OUTCAR', magmom=args.magmom, verbose=False)
                dic['MagMom'].append(round(magmom,3))
            else:
                energy, maxforce = check_energy_and_maxforce(f+'OUTCAR', magmom=False, verbose=False)
            dic['Config'].append(f)
            dic['Converged'].append(converged)
            dic['MaxForce'].append(round(maxforce,3))
            dic['Energy'].append(round(energy,3))
        else:
            print('No OUTCAR in ', f)
    
    dic['Rel.E'] = []
    for e in dic['Energy']:
        dic['Rel.E'].append(e - min(dic['Energy']))

    df = pd.DataFrame.from_dict(dic) 
    print(df)

if __name__ == "__main__":
    main()
