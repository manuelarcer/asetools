#!/usr/bin/env python

from ase.io import read
import pandas as pd
import matplotlib.pyplot as plt
import os

def extract_neb_data(folder_path, final):
    ## Extracts the energy of each image in the NEB calculation
    ## folder_path: path to the folder containing the NEB calculation
    ## final: the final image number, e.g. 5 for 00 to 05
    ## Returns a pandas DataFrame with the energy of each image

    frange = range(0, final+1)
    res = {'i':[], 'E':[]}

    for f in frange:
        res['i'].append(f)
        fname = os.path.join(folder_path, '{:02d}'.format(f) + '/OUTCAR')
        atoms = read(fname, format='vasp-out', index='-1')
        e = atoms.get_potential_energy()
        res['E'].append(e)

    nebdf = pd.DataFrame.from_dict(res, orient='columns')
    return nebdf

def plot_nebs(list_dfs=[], font='large'):
    # list_dfs: list of pandas DataFrames with the NEB data
    # font: font size for the labels and ticks
    
    if font == 'large':
        f_label = 18
        f_tick = 16
    elif font == 'medium':
        f_label = 16
        f_tick = 14
    elif font == 'small':
        f_label = 14
        f_tick = 12

    # Determine the NEB_df with the most images
    max_images = 0
    r_min_y = 9999 ; r_max_y = -9999
    for df in list_dfs:
        if len(df) > max_images:
            max_images = len(df)
        
        rel_e = df['E'] - df['E'].iloc[0]       # Relative energy wrt the first image
        min_y = min(rel_e)
        max_y = max(rel_e)
        r_min_y = min(r_min_y, min_y)           # Find the min and max relative energies
        r_max_y = max(r_max_y, max_y)
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    ax.axhline(y=0, color='gray', linestyle='-.', linewidth=0.5)
    for df in list_dfs:
        npoints = len(df['i'])
        minE = min(df['E'])
        ax.plot(range(npoints), df['E'] - df['E'].iloc[0], '-o', linewidth=1.5, markersize=6)
    #xlim(0, max_images)
    ax.set_xlim(-0.5, max_images-0.5)
    ax.set_ylim(r_min_y - 0.1, r_max_y + 0.1)         
    
    # Draw gray horizontal lines at Zero

    plt.ylabel(f'E$_i$ - E$_0$ (eV)', fontsize=f_label)
    plt.xlabel('Image', fontsize=f_label)
    plt.tick_params(axis='both', which='major', labelsize=f_tick)

    plt.tight_layout()
    plt.show()