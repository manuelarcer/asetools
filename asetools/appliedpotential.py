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

def fit_polynomial(results, order=3, fixed_constant=0, plot=False):
    poly_model = odr.Model(lambda beta, x: custom_polynomial(beta, x, fixed_constant))
    data = odr.Data(results['nelect'], results['e'])
    odr_obj = odr.ODR(data, poly_model, beta0=[1.0]*(order))
    output = odr_obj.run()
    
    if plot:
        parameters = np.append(output.beta[::-1], fixed_constant)
        poly = np.poly1d(parameters)
        x = np.linspace(min(results['nelect']), max(results['nelect']), 100)
        plt.plot(results['nelect'], results['e'], 'ok', label="original data")
        plt.plot(x, poly(x), '-k', label="fit")
        plt.legend()
        plt.xlabel(r'$\Delta$N$_{elec}$')
        plt.ylabel(r'E, eV')
        plt.tight_layout()
        plt.show()
    
    return output
