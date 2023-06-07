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
    dic = {'Config': [], 'Converged':[], 'MaxForce': [], 'Energy':[], 'MagMom':[]}
    for f in sorted(folders):
        outcar = False
        if os.path.exists(f+'OUTCAR'):
            outcar = f+'OUTCAR'
        elif os.path.exists(f+'OUTCAR.gz'):
            outcar = f+'OUTCAR.gz'
        if not outcar:
            converged = check_outcar_convergence(outcar, verbose=False)
            if args.magmom:
                energy, maxforce, magmom = check_energy_and_maxforce(outcar, magmom=args.magmom, verbose=False)
            else:
                energy, maxforce = check_energy_and_maxforce(outcar, magmom=False, verbose=False)
                magmom = 'NA'
            dic['Config'].append(f)
            dic['Converged'].append(converged)
            dic['MaxForce'].append(round(maxforce,3))
            dic['Energy'].append(round(energy,3))
            if type(magmom) == str:
                dic['MagMom'].append(magmom)
            else:
                dic['MagMom'].append(round(magmom,3))

        else:
            print('No OUTCAR in ', f)
    
    dic['Rel.E'] = []
    for e in dic['Energy']:
        dic['Rel.E'].append(e - min(dic['Energy']))

    df = pd.DataFrame.from_dict(dic) 
    print(df)

if __name__ == "__main__":
    main()
