#!/usr/bin/env python

import sys
from ase.io import read
from ase.visualize import view
from ase.gui.gui import GUI
import io

def preprocess_file_content(filename):
    """
    Remove all forward slashes from the file content and overwrite the original file.
    :param filename: str, the path to the file.
    """
    try:
        # Read the original content of the file
        with open(filename, 'r') as file:
            content = file.read()

        # Remove all forward slashes from the content
        modified_content = content.replace('/', '')

        # Overwrite the original file with the modified content
        with open(filename, 'w') as file:
            file.write(modified_content)

        print(f"File {filename} has been successfully cleaned of slashes.")
        file.close()

    except Exception as e:
        print(f"Error processing file {filename}: {e}")

def read_and_visualize_structure(filename):
    """
    Read the structure from a file and visualize it using ASE GUI.
    If a KeyError occurs due to element symbols containing slashes, fix the file content and retry.
    :param filename: str, the path to the structure file.
    """
    # Determine type of file, e.g., POSCAR, CONTCAR, OUTCAR, etc.
    is_poscar = 'POSCAR' in filename or 'CONTCAR' in filename or '.vasp' in filename
    is_outcar = 'OUTCAR' in filename

    try:
        if is_outcar:
            atoms = read(filename, format='vasp-out', index=':')
        elif is_poscar:
            atoms = read(filename, format='vasp')

        view(atoms)

    except KeyError as e:
        # Check if the error is due to a '/'
        if "/" in str(e):
            print(f"Detected problematic '/' in file: {filename} (CONTCAR can cause the issue in OUTCAR). Attempting to fix...")
            try:
                # Attempt to fix the file content
                if is_outcar:
                    file2fix = 'CONTCAR'
                    preprocess_file_content(file2fix)
                    atoms = read(filename, format='vasp-out', index=':')
                elif is_poscar:
                    file2fix = filename
                    preprocess_file_content(file2fix)
                    atoms = read(filename)
                view(atoms)
            except Exception as fix_e:
                print(f"Failed to fix and open file due to: {fix_e}")
        else:
            print(f"Unhandled KeyError: {e}")

    except Exception as e:
        print(f"Error while reading the file: {e}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python asegui.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    read_and_visualize_structure(filename)

if __name__ == "__main__":
    main()