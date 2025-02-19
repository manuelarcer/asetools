#!/usr/bin/env python3

# This script creates a backup folder with user provided name and copies the CONTCAR, OUTCAR, and vasprun.xml files to the backup folder.

import os
import shutil
import argparse
import glob
import gzip

def compress_file(original_file):
    with open(original_file, 'rb') as f_in:
        with gzip.open(original_file + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(original_file)  # Remove the original file after compression

def main():
    parser = argparse.ArgumentParser(description='Backup VASP files.')
    parser.add_argument('backupname', type=str, help='name of backup folder')
    args = parser.parse_args()

    wildcard_patterns = [
        '*.vasp', 
        '*.traj', 
        '*.cif', 
        '*.xyz', 
        '*.json', 
        '*.py',
        '*.out' 
        '*.sh', 
        '*.txt', 
        '*.log', 
        'OUTCAR*',
        'OSZICAR*',
        'INCAR*',
        'KPOINTS*',
        'CONTCAR*',
        'POSCAR*',
    ]

    if os.path.exists(args.backupname):
        print('Backup folder already exists')
        return

    try:
        os.makedirs(args.backupname)
        print(f'Created backup folder: {args.backupname}')

        file_list = []
        if not os.path.exists('vasprun.xml'):
            print('vasprun.xml not found in the current directory')
        else:
            file_list.append('vasprun.xml')

        # Add all files matching wildcard patterns
        for pattern in wildcard_patterns:
            file_list.extend(glob.glob(pattern))

        for filename in file_list:
            # Make sure CONTCAR stays in the current directory
            if 'CONTCAR' in filename:
                shutil.copy(filename, os.path.join(args.backupname, 'CONTCAR'))
                print(f'Copied {filename} to {args.backupname}')
            else:
                shutil.move(filename, args.backupname)
                print(f'Moved {filename} to {args.backupname}')

        # Compress specific files within the backup directory
        for file in glob.glob(os.path.join(args.backupname, 'OUTCAR*')):
            compress_file(file)
        compress_file(os.path.join(args.backupname, 'vasprun.xml'))


        print('Backup and compression completed successfully.')
    except FileNotFoundError as e:
        print(f'Error: {e}')
    except Exception as e:
        print(f'An unexpected error occurred: {e}')

if __name__ == "__main__":
    main()

