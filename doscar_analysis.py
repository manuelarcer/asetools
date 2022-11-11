import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def extract_total_dos(doscarfile):
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
            data['DOSdown'].append( float(line.split()[2]) )
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
        if i > 5 and line == repline:
            goDOS = False
            gopDOS = True
            count += 1

    dfdos = pd.DataFrame()
    dfdos['energy'] = data['energy']
    dfdos['DOSup'] = data['DOSup']
    dfdos['DOSdown'] = data['DOSdown']

    return dfdos