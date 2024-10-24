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
        ispyatoms = False
        if f+'log.info' in glob.glob(f+'log.info'):
            ispyatoms = True
            finishedpyatoms = False
            # Read log.info to check if the global processes have finished
            logfile = open(f+'log.info', 'r')
            lines = logfile.readlines()
            optsteps = []; e = []; fmax = []
            for line in lines:
                if 'pyatoms.jobs.job] Procedure' in line:
                    optsteps.append(int(line.split()[-5]))
                if 'success: energy' in line:               # Check if the job has finished, this may not be the best
                    finishedpyatoms = True
                if 'energy' in line:
                    e.append( float(line.split()[-3].split(',')[0]) )   # example: energy -481.276399, force 0.013.
                    fmax.append( float(line.split()[-1][:-1]) )
            if not finishedpyatoms:
                print(f, 'Pyatoms Job Not Finished')
                continue
            else:
                print(f, 'Pyatoms Job Finished')
                dic['Config'].append(f)
                dic['Converged'].append(True)
                dic['MaxForce'].append( round(fmax[-1], 3) )
                dic['Energy'].append( round(e[-1], 3) )

        ########## Block for fast mode ########
        elif args.fast:
            converged = fast_mode_check(f, alternative_filenames)
            if converged:
                dic['Config'].append(f)
                dic['Converged'].append(converged)
                dic['MaxForce'].append('N/A')
                dic['Energy'].append('N/A')
            else:
                print(f, 'not converged')
            continue
        #########################################

        elif not ispyatoms:
            if os.path.exists(f + 'OUTCAR'):
                try:
                    converged = check_outcar_convergence(f + 'OUTCAR', verbose=False)
                    if not converged:
                        print(f, 'not converged')
                        energy, maxforce = 'N/A', 'N/A'
                        continue
                    else:
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
        if e == 'N/A':
            dic['Rel.E'].append('N/A')
        else:
            dic['Rel.E'].append(e - min(dic['Energy'], default=0))  # default=0 to handle empty Energy list

    df = pd.DataFrame.from_dict(dic)

    # Add these lines to display the entire DataFrame
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)

    print(df)

def fast_mode_check(f, alternative_filenames):
    for alt in alternative_filenames:
        foundout = False
        if os.path.exists(f + alt):
            vaspout = open(f + alt, 'r')
            foundout = True
            break
    if not foundout:
        print(f'No OUT file ({alternative_filenames}) in', f)
        return False
    incar = open(f + 'INCAR', 'r')
    incarlines = incar.readlines()
    ibrion = None
    for line in incarlines:
        if 'IBRION' in line:
            ibrion = int(line.split()[2])
            break
    lines = vaspout.readlines()[-5:]
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
    return converged


if __name__ == "__main__":
    main()
