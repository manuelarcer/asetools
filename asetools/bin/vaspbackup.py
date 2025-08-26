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

    # Define which patterns should be copied (kept in place or copied) vs moved
    actions = {
        'copy': [
            'CONTCAR*',      # keep a copy of CONTCAR in the backup
            'POSCAR*',
            'KPOINTS*',
            'INCAR*',
            '*.cif',
            '*.xyz',
            '*.json',
            '*.py',
            '*.traj',
            '*.out',
            '*.txt',
            '*.log',
            '*.err',
            '*.info',
            '*.arc',         # requested patterns to copy
            '*.vasp',
            '*.in',
            '*.sh',
        ],
        'move': [
            'OUTCAR*',
            'OSZICAR*',
        ],
    }

    if os.path.exists(args.backupname):
        print('Backup folder already exists')
        return

    try:
        os.makedirs(args.backupname)
        print(f'Created backup folder: {args.backupname}')

        # Collect files to copy and move according to actions dict
        files_to_copy = []
        files_to_move = []

        # vasprun.xml is handled specially: include it in move list if present
        if not os.path.exists('vasprun.xml'):
            print('vasprun.xml not found in the current directory')
        else:
            files_to_move.append('vasprun.xml')

        for pattern in actions['copy']:
            files_to_copy.extend(glob.glob(pattern))

        for pattern in actions['move']:
            files_to_move.extend(glob.glob(pattern))

        # Deduplicate while preserving order
        def dedupe(seq):
            seen = set()
            out = []
            for x in seq:
                if x not in seen:
                    seen.add(x)
                    out.append(x)
            return out

        files_to_copy = dedupe(files_to_copy)
        files_to_move = dedupe(files_to_move)

        # Perform copy operations
        for filename in files_to_copy:
            # ensure destination filename for CONTCAR remains 'CONTCAR'
            if 'CONTCAR' in filename:
                shutil.copy(filename, os.path.join(args.backupname, 'CONTCAR'))
                print(f'Copied {filename} to {args.backupname}')
            else:
                shutil.copy(filename, args.backupname)
                print(f'Copied {filename} to {args.backupname}')

        # Perform move operations
        for filename in files_to_move:
            # skip if file was already copied
            if filename in files_to_copy:
                continue
            shutil.move(filename, args.backupname)
            print(f'Moved {filename} to {args.backupname}')

        # Compress specific files within the backup directory
        for file in glob.glob(os.path.join(args.backupname, 'OUTCAR*')):
            compress_file(file)
        # Compress vasprun.xml only if it was moved into the backup
        if os.path.exists(os.path.join(args.backupname, 'vasprun.xml')):
            compress_file(os.path.join(args.backupname, 'vasprun.xml'))


        print('Backup and compression completed successfully.')
    except FileNotFoundError as e:
        print(f'Error: {e}')
    except Exception as e:
        print(f'An unexpected error occurred: {e}')

if __name__ == "__main__":
    main()

