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
    parser.add_argument('-f', '--fast', action='store_true',
                        help='fast mode (reading out file)')
    args = parser.parse_args()
    ##

    folders = glob.glob('*/')
    if args.magmom:
        dic = {'Config': [], 'Converged':[], 'MaxForce': [], 'Energy':[], 'MagMom':[]}
    else:
        dic = {'Config': [], 'Converged':[], 'MaxForce': [], 'Energy':[]}

    # Define alternative filenames to look for when fast mode is enabled
    alternative_filenames = ['vasp.out', 'out.txt']

    for f in sorted(folders):
        
        ########## Block for fast mode ########
        if args.fast:
            for alt in alternative_filenames:
                if os.path.exists(f + alt):
                    vaspout = open(f + alt, 'r')
                    break
            if not vaspout:
                print('No OUT file in', f)
                continue
            
            incar = open(f + 'INCAR', 'r')
            incarlines = incar.readlines()
            ibrion = None
            for line in incarlines:
                if 'IBRION' in line:
                    ibrion = int(line.split()[2])
                    break
            lines = vaspout.readlines()
            if ibrion == 1 or ibrion == 2 or ibrion == 3:
                for line in lines:
                    if 'reached required accuracy' in line:
                        converged = True
                        break
                    else:
                        converged = False
            else:
                for line in lines:
                    if 'E0=' in line:
                        converged = True
                        break
                    else:
                        converged = False
        #########################################

        else:
            if os.path.exists(f + 'OUTCAR'):
                try:
                    converged = check_outcar_convergence(f + 'OUTCAR', verbose=False)
                    if args.magmom:
                        energy, maxforce, magmom = check_energy_and_maxforce(f + 'OUTCAR', magmom=args.magmom, verbose=False)
                        dic['MagMom'].append(round(magmom, 3))
                    else:
                        energy, maxforce = check_energy_and_maxforce(f + 'OUTCAR', magmom=False, verbose=False)

                    dic['Config'].append(f)
                    dic['Converged'].append(converged)
                    dic['MaxForce'].append(round(maxforce, 3))
                    dic['Energy'].append(round(energy, 3))
                except ValueError as e:
                    print(f'Error processing {f}: {e}. OUTCAR may be incomplete or damaged.')
            else:
                print('No OUTCAR in', f)

    dic['Rel.E'] = []
    for e in dic['Energy']:
        dic['Rel.E'].append(e - min(dic['Energy'], default=0))  # default=0 to handle empty Energy list

    df = pd.DataFrame.from_dict(dic)
    print(df)

if __name__ == "__main__":
    main()
