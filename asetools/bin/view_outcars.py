#!/usr/bin/env python

import os
import argparse
from ase.io import read
from ase.visualize import view
from asetools.analysis import check_outcar_convergence

def main():
    parser = argparse.ArgumentParser(description='Visualize the last configuration in OUTCAR files.')
    parser.add_argument('--converged', action='store_true', help='Only consider converged calculations')
    args = parser.parse_args()
    
    listatoms = []
    for f in sorted( os.listdir('.') ):
        outcar_path = os.path.join(f, 'OUTCAR')
        if os.path.isfile(outcar_path):
            converged, _ = check_outcar_convergence(outcar_path, verbose=False)
            if not args.converged or converged:
                atoms = read(outcar_path, format='vasp-out', index=-1)
                listatoms.append(atoms)
    if listatoms:
        view(listatoms)
    else:
        print("No suitable OUTCAR files found.")

if __name__ == "__main__":
    main()

