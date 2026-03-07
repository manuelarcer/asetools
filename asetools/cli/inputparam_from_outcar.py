#!/usr/bin/env python

import re
import argparse

# List of PARAM to check from OUTCAR
list_param = {'StartParam': ['PREC', 'ISTART', 'ICHARG', 'ISPIN'],
              'Electronic': ['ENCUT', 'NELM', 'EDIFF', 'LREAL', 'LMAXMIX', 'VOSKOWN'],
              'IonicRelax': ['EDIFFG', 'NSW', 'IBRION', 'NFREE', 'ISIF', 'ISYM', 'LCORR', 'POTIM'],
              'System': ['POMASS', 'ZVAL', 'RWIGS', 'VCA', 'NELECT', 'NUPDOWN'],
              'DOS': ['EMIN', 'EMAX', 'EFERMI', 'ISMEAR', 'IALGO', 'SIGMA'], 
              'ElectronicRelax': ['IALGO'],
              'Write': ['LWAVE', 'LCHARG', 'LVTOT', 'LVHAR', 'LELF', 'LORBIT'],
              'DipoleCorr': ['LMONO', 'LDIPOL', 'IDIPOL', 'EPSILON'],
              'ExCorr':  ['GGA', 'VOSKOWN', 'LHFCALC', 'LHFONE', 'AEXX'],
              'LinearResp': ['LEPSILON']
            }

def extract_parameters(outcar_text, param_dict):
    found_params = {key: {} for key in param_dict.keys()}
    
    # Loop through each category in param_dict
    for category, params in param_dict.items():
        for param in params:
            # Find all instances of the parameter and its value in the OUTCAR text
            regex_pattern = rf"^\s*({param})\s*=\s*([^\s;]+)"  # Match parameter and value until space or semicolon
            matches = re.findall(regex_pattern, outcar_text, re.MULTILINE)
            for param_name, param_value in matches:
                found_params[category][param_name] = param_value.strip()
    
    return found_params

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract parameters from an OUTCAR file.")
    parser.add_argument('outcar_file', type=str, help='Path to the OUTCAR file to analyze')
    
    # Parse arguments
    args = parser.parse_args()

    # Read the OUTCAR file
    with open(args.outcar_file, 'r') as file:
        outcar_text = file.read()
    
    # Extract parameters
    extracted_params = extract_parameters(outcar_text, list_param)

    # Print extracted parameters
    for category, params in extracted_params.items():
        print(f"{category}:")
        for param, value in params.items():
            print(f"  {param} = {value}")
        print()

if __name__ == '__main__':
    main()