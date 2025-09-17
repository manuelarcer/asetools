#!/usr/bin/env python

from ase.io import read, write
import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python outcar_extract.py OUTCAR [index] [output_file] [format]")
        print("       python outcar_extract.py OUTCAR                    # last frame to final.vasp")
        print("       python outcar_extract.py OUTCAR 5                 # 5th frame to final.vasp")
        print("       python outcar_extract.py OUTCAR 5 freq/CONTCAR.vasp  # 5th frame to specific file")
        print("       python outcar_extract.py OUTCAR -1 final.extxyz extxyz  # with format")
        sys.exit(1)

    input_file = sys.argv[1]

    # Parse optional arguments
    index = -1
    output_file = 'final.vasp'
    format_type = None

    if len(sys.argv) > 2:
        index = int(sys.argv[2])
    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    if len(sys.argv) > 4:
        format_type = sys.argv[4]

    if not os.path.exists(input_file):
        print(f'ERROR: Input file not found: {input_file}')
        sys.exit(1)

    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read and write structure
    atoms = read(input_file, index=index)

    if format_type:
        write(output_file, atoms, format=format_type)
    else:
        write(output_file, atoms)

    print(f'Wrote frame {index} from {input_file} -> {output_file}')

if __name__ == "__main__":
    main()