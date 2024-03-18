#!/usr/bin/env python3

# This script creates a backup folder with user provided name and copies the CONTCAR, OUTCAR, and vasprun.xml files to the backup folder.

import os
import shutil
import argparse

def main():
    parser = argparse.ArgumentParser(description='Backup VASP files.')
    parser.add_argument('backupname', type=str, help='name of backup folder')
    args = parser.parse_args()

    if os.path.exists(args.backupname):
        print('Backup folder already exists')
        return

    try:
        os.makedirs(args.backupname)
        print(f'Created backup folder: {args.backupname}')

        for filename in ['POSCAR', 'CONTCAR', 'OUTCAR', 'vasprun.xml', 'vasp.out', 'INCAR', 'KPOINTS']:
            shutil.copy(filename, args.backupname)
            print(f'Copied {filename} to {args.backupname}')

        print('Backup completed successfully.')
    except FileNotFoundError as e:
        print(f'Error: {e}')
    except Exception as e:
        print(f'An unexpected error occurred: {e}')

if __name__ == "__main__":
    main()

