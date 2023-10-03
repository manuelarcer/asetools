import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def extract_dos(doscarfile):
    # doscarfile: Is the path to the DOSCAR file
    
    doscar = open(doscarfile, 'r')
    lines = doscar.readlines()
    nat = int( lines[0].split()[0] )        # Number of atoms in DOSCAR
    typedos = int( lines[0].split()[2] )
    if typedos == 0:
        partialdos = False
    elif typedos == 1:
        partialdos = True
    
    ##repline = lines[5]

    data = { 'energy':[], 'DOSup':[], 'DOSdown':[] }
    for i in range(nat):
        # Define the data dictionary with the order of states as they appear in the DOSCAR
        data['at-'+str(i)] = {'energy':[], 's+':[], 's-':[], 'py+':[], 'py-':[],
            'pz+':[], 'pz-':[], 'px+':[], 'px-':[],
            'dxy+':[], 'dxy-':[], 'dyz+':[], 'dyz-':[],
            'dz2+':[], 'dz2-':[], 'dxz+':[], 'dxz-':[],
            'dx2+':[], 'dx2-':[] }

    goDOS = False; gopDOS = False
    count = -1
    for i, line in enumerate(lines):
        if goDOS:
            data['energy'].append( float(line.split()[0]) - fermie)
            data['DOSup'].append( float(line.split()[1]) )
            data['DOSdown'].append( -float(line.split()[2]) )
        if gopDOS and line != repline:
            for k, key in enumerate( data['at-'+str(count)] ):
                if k == 0:
                    data['at-'+str(count)]['energy'].append( float(line.split()[k]) - fermie )
                else:
                    data['at-'+str(count)][key].append( float(line.split()[k]) )
        if i == 5:
            repline = line  # Line with 5 values
            fermie = float( line.split()[3] )
            goDOS = True
        if i > 5 and line == repline and partialdos:
            goDOS = False
            gopDOS = True
            count += 1

    return data

def extract_fermi_e(doscarfile):
    try:
        with open(doscarfile, 'r') as file:
            lines = file.readlines()
            fermie = float( lines[5].split()[3] )
        return fermie
    except FileNotFoundError:
        print(f"'{doscarfile}' not found.")
    except IndexError:
        print("Unexpected file format. Couldn't extract Fermi energy.")
    except ValueError:
        print("Error converting Fermi energy to float.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    return None


def extract_pdos_perstate(data, atoms, states):
    # data: First, extract 'data' with extract_dos()
    # atoms: list atoms of interest, e.g., [0, 10]
    # states: list of states of interest, e.g., ['s_states', 'p_states', 'd_states']

    dicstates = {
            's_states':{'s+': 1, 's-': 2},
            'p_states':{'py+': 3, 'py-': 4,
            'pz+': 5, 'pz-': 6,
            'px+': 7, 'px-': 8},
            'd_states':{
            'dxy+': 9, 'dxy-': 10,
            'dyz+': 11, 'dyz-': 12,
            'dz2+': 13, 'dz2-': 14,
            'dxz+': 15, 'dxz-': 16,
            'dx2+': 17, 'dx2-': 18
            }
        }
    sum_plus = np.zeros(len(data['at-0']['py+']))
    sum_minus = np.zeros(len(data['at-0']['py+']))
    e = data['at-0']['energy']
    for at in atoms:
        for sss in states:
            for key in dicstates[sss]:     ## Loop through py+, py-, pz+, pz-, etc 
                if key[-1] == '+':
                    sum_plus += data['at-'+str(at)][key]
                elif key[-1] == '-':
                    sum_minus -= data['at-'+str(at)][key]
    return e, sum_plus, sum_minus

def extract_pdos_perorbital(data, atoms, orbitals):
    # data: First, extract 'data' with extract_dos()
    # atoms: list atoms of interest, e.g., [0, 10]
    # orbitals: list of orbitals of interest, e.g., ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+'] 
    #         or 'all-d' to consider all d-states
    #         or 't2g' for ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-']
    #         or 'eg' for ['dz2+', 'dz2-', 'dx2+', 'dx2-']

    if orbitals == 'all-s':
        orbitals = ['s+', 's-']
    elif orbitals == 'all-p':
        orbitals = ['py+', 'py-', 'pz+', 'pz-', 'px+', 'px-']
    elif orbitals == 'all-d':
        orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-', 'dz2+', 'dz2-', 'dx2+', 'dx2-']
    elif orbitals == 'all':
        orbitals = [
            's+', 's-', 'py+', 'py-', 'pz+', 'pz-', 'px+', 'px-',
            'dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-', 'dz2+', 'dz2-', 'dx2+', 'dx2-'
            ]
    elif orbitals == 't2g':
        orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-']
    elif orbitals == 'eg':
        orbitals = ['dz2+', 'dz2-', 'dx2+', 'dx2-']
    
    sum_plus = np.zeros(len(data['at-0']['py+']))
    sum_minus = np.zeros(len(data['at-0']['py+']))
    e = data['at-0']['energy']
    for at in atoms:
        for sss in orbitals: 
            if sss[-1] == '+':
                sum_plus += data['at-'+str(at)][sss]
            elif sss[-1] == '-':
                sum_minus -= data['at-'+str(at)][sss]
    return e, sum_plus, sum_minus