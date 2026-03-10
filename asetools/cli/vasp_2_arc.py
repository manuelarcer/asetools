#!/usr/bin/env python

import sys
import argparse
from ase.io import read, write

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert VASP output files to ARC format.')
    parser.add_argument('input_file', help='Path to the VASP output file (e.g., OUTCAR)')
    parser.add_argument('output_file', nargs='?', help='Path to the output ARC file (optional)')
    args = parser.parse_args()

    if not args.output_file:
        args.output_file = args.input_file.split('.')[0] + '.arc'

    return args

def main():
    args = parse_arguments()
    
    input_file = args.input_file
    output_file = args.output_file
    
    try:
        atoms = read(input_file)
        write(output_file, atoms, format='dmol-arc')
        print(f'Successfully converted {input_file} to {output_file}')
    except Exception as e:
        print(f'Error during conversion: {e}', file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()