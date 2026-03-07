#!/usr/bin/env python

import sys

def remove_slashes_from_file(input_file, output_file=None):
    """
    Read the content of the input file, remove all forward slashes '/',
    and write the modified content to the output file or overwrite the input file if
    no output file is specified.

    :param input_file: str, the path to the file to be modified.
    :param output_file: str or None, the path to the file where the modified content will be written.
                        If None, the input file will be overwritten.
    """
    try:
        # Read the content of the input file
        with open(input_file, 'r') as file:
            content = file.read()
        
        # Remove all forward slashes
        modified_content = content.replace('/', '')
        
        # If no output file is specified, overwrite the input file
        if output_file is None:
            output_file = input_file
        
        # Write the modified content to the output file
        with open(output_file, 'w') as file:
            file.write(modified_content)
        
        print(f"Successfully removed '/' from {input_file} and saved to {output_file}")

    except Exception as e:
        print(f"Failed to process the file due to: {e}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python remove_slashes.py <input_file> [output_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    remove_slashes_from_file(input_file, output_file)

if __name__ == "__main__":
    main()