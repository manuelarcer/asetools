#!/usr/bin/env python

import os
import argparse
from ase.io import read
from ase.visualize import view
from asetools.analysis import check_outcar_convergence

def main():
    parser = argparse.ArgumentParser(description='Visualize the last configuration in OUTCAR files.')
    parser.add_argument('--converged', action='store_true', help='Only consider converged calculations')
    parser.add_argument('--folders', '-f', nargs='+', help='List specific folders to process. If not provided, all folders in the current directory will be processed.')
    parser.add_argument('--pyatoms', action='store_false', help='Assume that Pyatoms was not used for this run')
    args = parser.parse_args()
    
    if args.folders:
        # User has specified folders to process
        folders_to_process = args.folders
        valid_folders = []
        for folder in folders_to_process:
            if os.path.isdir(folder):   
                valid_folders.append(folder)
            else:
                print(f'WARNING: {folder} is not a valid directory, skipping...')
        if not valid_folders:
            print('No valid folders provided. Exiting')
            return
    else:
        # Process all folders in the current directory
        valid_folders = [f for f in sorted( os.listdir('.') ) if os.path.isdir(f)]

    listatoms = []
    for f in valid_folders:
        if args.pyatoms:        # Pyatoms was used
            # List all step_0* with OUTCARs within and use the last instance with OUTCAR in it
            ste_dirs = [d for d in os.listdir(f) if d.startswith('step_0') and os.path.isdir(os.path.join(f, d))]
            if not ste_dirs:
                print(f'No step_0* directories found in {f}, skipping...')
                continue
            ste_dirs.sort()
            for d in ste_dirs:
                if os.path.isfile(os.path.join(f, d, 'OUTCAR')):
                    last_ste_dir = d
            outcar_path = os.path.join(f, last_ste_dir, 'OUTCAR')
        else:               # Pyatoms was not used
            outcar_path = os.path.join(f, 'OUTCAR')
        
        if os.path.isfile(outcar_path):
            converged, _ = check_outcar_convergence(outcar_path, verbose=False)
            if not args.converged or converged:
                try:
                    atoms = read(outcar_path, format='vasp-out', index=-1)
                    listatoms.append(atoms)
                except Exception as e:
                    print(f'Error reading {outcar_path}: {e}')
        else:
            print(f'OUTCAR not found in {f}, skipping...')
            
    if listatoms:
        view(listatoms)
    else:
        print("No suitable OUTCAR files found.")

if __name__ == "__main__":
    main()