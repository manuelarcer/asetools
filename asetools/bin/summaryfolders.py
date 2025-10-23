#!/usr/bin/env python

import argparse
from ase.io import read
import pandas as pd
import numpy as np
import glob, os
from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence, get_parameter_from_run

def is_summary_up_to_date(magmom_requested=False):
    """Check if summary.log exists and is newer than all folders, and if magmom flag status matches"""
    if not os.path.exists('summary.log'):
        return False

    # Check if magnetic moments flag status has changed
    with open('summary.log', 'r') as f:
        log_content = f.read()
        log_has_magmom = 'MagMom' in log_content

    # If magmom is requested but log doesn't have it, or vice versa, regenerate
    if magmom_requested != log_has_magmom:
        return False

    summary_mtime = os.path.getmtime('summary.log')
    folders = glob.glob('*/')

    for folder in folders:
        # Check if any OUTCAR in folders is newer than summary.log
        outcar_path = os.path.join(folder, 'OUTCAR')
        if os.path.exists(outcar_path):
            if os.path.getmtime(outcar_path) > summary_mtime:
                return False
        # Also check for log.info files (pyatoms)
        loginfo_path = os.path.join(folder, 'log.info')
        if os.path.exists(loginfo_path):
            if os.path.getmtime(loginfo_path) > summary_mtime:
                return False

    return True

def main():
    ## Block from chatgpt
    parser = argparse.ArgumentParser(description='Process folders.')
    parser.add_argument('-m', '--magmom', action='store_true',
                        help='extract and present magnetic moments')
    parser.add_argument('-f', '--fast', action='store_true',
                        help='fast mode (reading out file)')
    args = parser.parse_args()
    ##

    # Check if summary.log is up to date
    if is_summary_up_to_date(magmom_requested=args.magmom):
        print("Summary is up to date. Reading from summary.log:")
        print()
        with open('summary.log', 'r') as f:
            print(f.read())
        return

    folders = glob.glob('*/')
    if args.magmom:
        dic = {'Config': [], 'ISIF':[], 'Converged':[], 'ENCUT': [], 'Target-fmax': [], 'MaxForce': [], 'Energy':[], 'MagMom':[]}
    else:
        dic = {'Config': [], 'ISIF':[], 'Converged':[], 'ENCUT': [], 'Target-fmax': [], 'MaxForce': [], 'Energy':[]}

    # Define alternative filenames to look for when fast mode is enabled
    alternative_filenames = ['vasp.out', 'out.txt']
    not_converged = []
    for f in sorted(folders):
        ispyatoms = False
        if f+'log.info' in glob.glob(f+'log.info'):
            ispyatoms = True
            finishedpyatoms = False
            finishedprocess = False
            # Read log.info to check if the global processes have finished
            logfile = open(f+'log.info', 'r')
            lines = logfile.readlines()
            optsteps = []; e = []; fmax = []
            for line in lines:
                if 'pyatoms.jobs.job]' in line and 'Procedure step' in line:
                    optsteps.append(int(line.split()[-5]))
                #else:
                #    optsteps = 1
                # I think if multistep Opt there will be multiple 'VaspGeomOptProcedure' and 'VaspGeomOptJob' in log.info
                # for a single step procedure, there will be only one 'VaspGeomOptJob'
                if 'VaspGeomOptProcedure' in line and 'completed successfully' in line:      # this work for GeomOpt only
                    finishedprocess = True
                elif 'VaspGeomOptJob' in line and 'completed successfully' in line:           # this work for GeomOpt only
                    finishedprocess = True

                if finishedprocess:
                    if 'Closing global processes' in line:
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
                dic['Converged'].append((True, 'pyatoms'))
                dic['ENCUT'].append('N/A')
                dic['Target-fmax'].append('N/A')
                dic['MaxForce'].append( round(fmax[-1], 3) )
                dic['Energy'].append( round(e[-1], 3) )

        ########## Block for fast mode ########
        elif args.fast:
            converged = fast_mode_check(f, alternative_filenames)
            if converged:
                dic['Config'].append(f)
                dic['Converged'].append(converged)
                dic['ENCUT'].append('N/A')
                dic['Target-fmax'].append('N/A')
                dic['MaxForce'].append('N/A')
                dic['Energy'].append('N/A')
            else:
                print(f, 'not converged')
                not_converged.append(f)
            continue
        #########################################

        elif not ispyatoms:
            if os.path.exists(f + 'OUTCAR'):
                try:
                    converged = check_outcar_convergence(f + 'OUTCAR', verbose=False)
                    if not converged[0]:
                        print(f, 'not converged')
                        not_converged.append(f)
                    
                    if args.magmom:
                        energy, maxforce, magmom = check_energy_and_maxforce(f + 'OUTCAR', magmom=args.magmom, verbose=False)
                        dic['MagMom'].append(round(magmom, 3))
                    else:
                        energy, maxforce = check_energy_and_maxforce(f + 'OUTCAR', magmom=False, verbose=False)
                    
                    isif, _ = get_parameter_from_run(f + 'OUTCAR', check_converg=False, parameter='ISIF')
                    encut, _ = get_parameter_from_run(f + 'OUTCAR', check_converg=False, parameter='ENCUT')
                    ediffg, _ = get_parameter_from_run(f + 'OUTCAR', check_converg=False, parameter='EDIFFG')
                    
                    # Convert EDIFFG to positive value for target fmax (VASP uses negative for forces)
                    target_fmax = abs(ediffg) if ediffg is not None else 'N/A'

                    dic['Config'].append(f)
                    dic['ISIF'].append(isif)
                    dic['Converged'].append(converged[0])
                    dic['ENCUT'].append(encut)
                    dic['Target-fmax'].append(target_fmax)
                    dic['MaxForce'].append(round(maxforce, 3))
                    dic['Energy'].append(round(energy, 3))
                    

                except ValueError as e:
                    print(f'Error processing {f}: {e}. OUTCAR may be incomplete or damaged.')
                    not_converged.append(f)
            else:
                print('No OUTCAR in', f)
                not_converged.append(f)

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
    print()
    not_converged = [f.split('/')[0] for f in not_converged]
    print('Not converged:')
    print(' '.join(not_converged))
    
    # Write output to summary.log
    with open('summary.log', 'w') as f:
        f.write(str(df))
        f.write('\n\n')
        f.write('Not converged:\n')
        f.write(' '.join(not_converged))
        f.write('\n')

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
