# Module for analysis of potential-dependent calculations

from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce
from asetools.doscar_analysis import extract_fermi_e
from scipy import odr
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt

def get_num_elect(outcar):
    with open(outcar, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if 'NELECT' in line:
                nelect = float( line.split()[-1] )  # Float since it may be fractional
                break
    return nelect

def extract_corrected_energy_fermie(listfolders, calc_zero):
    # The quadratic fit relies on calculation of the derivatives of the energy with respect to the num of elec
    # these calculations should be done around the neutral system (prob better if the difference in num elec is small)
    # Input should be the folders of the calculations:
    #   with decreased num of elec (calc), neutral system (calc_zero), and the one with more electrons (calc_plus)

    results = {'nelect': [], 'e': [], 'fe': [], 'U': []}

    ref_nelect = get_num_elect(calc_zero+'/OUTCAR')

    for folder in listfolders:
        convergence = check_outcar_convergence(folder+'/OUTCAR', verbose=False)
        
        if convergence:
            print(f'GOOD news, calculation at {folder} CONVERGED')
            atoms = read(folder+'/OUTCAR', format='vasp-out', index=-1)
            energy_original = atoms.get_potential_energy()
            nelect = get_num_elect(folder+'/OUTCAR')

            #energy, maxforce = check_energy_and_maxforce(folder+'/OUTCAR', magmom=False, verbose=False)
            
            fermie_original = extract_fermi_e(folder+'/DOSCAR')
            fermi_shift = None
            try:
                with open(folder+'/vasp.out', 'r') as file:
                    lines = file.readlines()
                    for line in lines[-40:]:
                        if 'FERMI_SHIFT' in line:
                            fermi_shift = float( line.split('=')[-1] )
            
            except:
                print('WARNING ! Error while reading "vasp.out", no FERMI_SHIFT stored')
                fermi_shift = None

            fermie_corr = fermie_original + fermi_shift
            energy_corr = energy_original - ( nelect - ref_nelect) * fermi_shift

            results['e'].append(energy_corr)
            results['fe'].append(fermie_corr)
            results['nelect'].append(nelect - ref_nelect)
            results['U'].append( -fermie_corr - 4.43 )   # U(SHE) = 4.43 V
            
        else:
            print(f'Oh NOOOO!, something wrong with calculation at {folder}')
    
    return results

def custom_polynomial(beta, x, fixed_constant=0):
    """
    Custom polynomial function with a fixed constant term.
    """
    # The first coefficient is fixed
    y = fixed_constant
    for i in range(1, len(beta) + 1):
        y += beta[i-1] * (x ** i)
    return y

def fit_polynomial(results, order=3, fixed_constant=0, plot=False):
    # Define the custom model
    poly_model = odr.Model(lambda beta, x: custom_polynomial(beta, x, fixed_constant))

    data = odr.Data(results['nelect'], results['e'])
    odr_obj = odr.ODR(data, poly_model, beta0=[1.0]*(order))  # Initial guesses for coefficients
    output = odr_obj.run()  # running ODR fitting
    if plot:
        parameters = np.append(output.beta[::-1], fixed_constant)
        poly = np.poly1d(parameters)
        x = np.arange( min(results['nelect']), max( results['nelect'] ), 0.01 )
        poly_y = poly(x)
        plt.plot(results['nelect'], results['e'], 'ok', label="original data")
        plt.plot(x, poly_y, '-k', label="fit")
        plt.legend()
        plt.xlabel(r'$\Delta$N$_{elec}$')
        plt.ylabel(r'E, eV')
        plt.tight_layout()
        plt.show()
    return output

def first_and_second_derivatives(results):
    # input is the "results" dictionary from the "extract_corrected_energy_fermie" function
    # The shape is like this: "{'nelect': [], 'e': [], 'fe': []}"     

    dE_dNe = ( results['e'][2] - results['e'][0] ) / ( results['nelect'][2] - results['nelect'][0] )
    d2E_dNE2 = ( results['e'][2] - 2 * results['e'][1] + results['e'][0] ) / ( results['nelect'][2] - results['nelect'][1] )**2

    return dE_dNe, d2E_dNE2

