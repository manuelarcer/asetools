# Module for analysis of potential-dependent calculations

from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce
from asetools.doscar_analysis import extract_fermi_e
from scipy import odr
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import os

U_SHE = 4.43  # U(SHE) constant

def get_num_elect(outcar):
    with open(outcar, 'r') as file:
        for line in file:
            if 'NELECT' in line:
                return float(line.split()[-1])  # Float since it may be fractional

def extract_fermi_shift(folder):
    try:
        with open(os.path.join(folder, 'vasp.out'), 'r') as file:
            for line in file.readlines()[-40:]:
                if 'FERMI_SHIFT' in line:
                    return float(line.split('=')[-1])
    except:
        print('WARNING ! Error while reading "vasp.out", no FERMI_SHIFT stored')
        return None

def extract_corrected_energy_fermie(listfolders, calc_zero):
    results = {'nelect': [], 'e': [], 'fe': [], 'U': []}
    ref_nelect = get_num_elect(os.path.join(calc_zero, 'OUTCAR'))

    for folder in listfolders:
        if check_outcar_convergence(os.path.join(folder, 'OUTCAR'), verbose=False):
            print(f'GOOD news, calculation at {folder} CONVERGED')
            atoms = read(os.path.join(folder, 'OUTCAR'), format='vasp-out', index=-1)
            energy_original = atoms.get_potential_energy()
            nelect = get_num_elect(os.path.join(folder, 'OUTCAR'))
            
            fermie_original = extract_fermi_e(os.path.join(folder, 'DOSCAR'))
            fermi_shift = extract_fermi_shift(folder)
            
            fermie_corr = fermie_original + (fermi_shift if fermi_shift else 0)
            energy_corr = energy_original - (nelect - ref_nelect) * (fermi_shift if fermi_shift else 0)

            results['e'].append(energy_corr)
            results['fe'].append(fermie_corr)
            results['nelect'].append(nelect - ref_nelect)
            results['U'].append(-fermie_corr - U_SHE)
        else:
            print(f'Oh NOOOO!, something wrong with calculation at {folder}')
    
    return results

# defining a function with a fix constant value. This should be the value at 0 added electrons
def custom_polynomial(beta, x, fixed_constant=0):
    y = fixed_constant
    for i in range(1, len(beta) + 1):
        y += beta[i-1] * (x ** i)
    return y

def fit_polynomial(results, order=3, energy_ref=0, plot=False, ploterrors=False):
    poly_model = odr.Model(lambda beta, x: custom_polynomial(beta, x, energy_ref))
    data = odr.Data(results['nelect'], results['e'])
    odr_obj = odr.ODR(data, poly_model, beta0=[1.0]*(order))
    output = odr_obj.run()
    if sum([plot, ploterrors]) == 2:
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
        plot_fit(results, output, energy_ref, ax1)
        plot_errors(results, output, energy_ref, ax2)
        plt.tight_layout()
        plt.show()
    elif sum([plot, ploterrors]) == 1:
        fig, ax = plt.subplots(figsize=(5,5))
        if plot:
            plot_fit(results, output, energy_ref, ax)
        elif ploterrors:
            plot_errors(results, output, energy_ref, ax)
        plt.tight_layout()
        plt.show()
    return output

def plot_errors(results, polyfit, energy_ref, ax):
    # Input, results from the 'extract_corrected_energy_fermie' and polyfit from 'fit_polynomial'
    # energy_ref is the energy of the neutral system
    # ax is the axis pyplot object
     
    parameters = np.append(polyfit.beta[::-1], energy_ref)
    poly = np.poly1d( parameters )
    poly_y = poly( results['nelect'] )
    errors = poly_y - results['e']
    ax.bar( results['nelect'], errors )
    lower = min(results['nelect'])
    higher = max(results['nelect'])
    shift = (higher - lower) / 10
    ax.plot( [lower - shift, higher + shift], [0, 0], '-k', linewidth=1.5 )
    ax.set_xlabel(r'$\Delta$N$_{elec}$')
    ax.set_ylabel(r'E$_{fit}$ - E$_{DFT}$, eV')
    return ax

def plot_fit(results, polyfit, energy_ref, ax):
    
    parameters = np.append(polyfit.beta[::-1], energy_ref)
    poly = np.poly1d(parameters)
    x = np.linspace(min(results['nelect']), max(results['nelect']), 100)
    ax.plot(results['nelect'], results['e'], 'ok', label="original data")
    ax.plot(x, poly(x), '-k', label="fit")
    ax.legend()
    ax.set_xlabel(r'$\Delta$N$_{elec}$')
    ax.set_ylabel(r'E, eV')
    return ax

def print_results(results):
    # Input is the dictionary from "extract_corrected_energy_fermie"
    for i, nelec in enumerate(results['nelect']):
        print(nelec, results['e'][i], results['fe'][i], results['U'][i])
