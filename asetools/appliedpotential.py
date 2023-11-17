# Module for analysis of potential-dependent calculations

from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce
from asetools.doscar_analysis import extract_fermi_e
from scipy.interpolate import UnivariateSpline
from scipy import odr
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import os

U_SHE = 4.43  # U(SHE) constant
VALANCES = {'Cu': 11, 'Zn': 12, 'C': 4, 'O': 6, 'H': 1}

def get_sum_electrons(poscar):
    atoms = read(poscar)
    symbols = atoms.get_chemical_symbols()
    elecpersymb = [VALANCES[symb] for symb in symbols]
    return sum(elecpersymb)

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
    
def correct_energy_fermishift(folder):
    fshift = extract_fermi_shift(folder)
    e, maxf = check_energy_and_maxforce(folder+'OUTCAR', magmom=False, verbose=False)
    numelec_neutral = get_sum_electrons(folder+'POSCAR')
    numelec = get_num_elect(folder+'OUTCAR')
    charge = numelec - numelec_neutral
    return e + charge * fshift

def extract_corrected_energy_fermie(listfolders, calc_zero):
    results = {'nelect': [], 'e': [], 'fe': [], 'U': []}
    ref_nelect = get_num_elect(os.path.join(calc_zero, 'OUTCAR'))

    for folder in listfolders:
        if check_outcar_convergence(os.path.join(folder, 'OUTCAR'), verbose=False):
            #print(f'GOOD news, calculation at {folder} CONVERGED')
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
    
    # Convert lists to numpy arrays
    for key in results:
        results[key] = np.array(results[key])

    return results

# defining a function with a fix constant value. This should be the value at 0 added electrons
def custom_polynomial(beta, x, fixed_constant=0):
    y = fixed_constant
    for i in range(1, len(beta) + 1):
        y += beta[i-1] * np.power(x, i)
    return y

def fitenergy_polynomial(results, order=3, energy_ref=0, plot=False, ploterrors=False):
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

def fit_data(X, Y, fit_type='polynomial', order=2, ref_value=None, plot=False, ploterrors=False):
    
    if fit_type == 'polynomial':
        data = odr.Data(X, Y)
        
        if ref_value != None:
            poly_model = odr.Model(lambda beta, x: custom_polynomial(beta, x, fixed_constant=ref_value))
            odr_obj = odr.ODR(data, poly_model, beta0=[1.0]*(order))
            output = odr_obj.run()
            Y_fit = custom_polynomial(output.beta, X, ref_value)
        elif ref_value == None:
            poly_model = odr.polynomial(order)
            odr_obj = odr.ODR(data, poly_model)
            output = odr_obj.run()
            Y_fit = np.polyval(output.beta, X)    
        
        fit_result = output

    elif fit_type == 'spline':
        spline = UnivariateSpline(X, Y, k=order)
        Y_fit = spline(X)
        fit_result = spline

    else:
        raise ValueError("Invalid fit_type. Choose either 'polynomial' or 'spline'.")
    
    if ref_value is not None:
        parameters = np.append(fit_result.beta[::-1], ref_value)
    else:
        parameters = fit_result.beta[::-1]
    poly = np.poly1d(parameters)
    # Total Sum of Squares
    SST = np.sum((Y - np.mean(Y))**2)
    # Residual Sum of Squares
    SSR = np.sum((Y - poly(X))**2)
    # R-squared
    R_squared = 1 - (SSR / SST)
    print(f'Sum of squared-residuals : {SSR:.5f}')
    print(f'R-squared : {R_squared:.5f}')

    ### Plotting section
    if sum([plot, ploterrors]) == 2:
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
        plot_fit(X, Y, fit_result, ref_value, ax1)
        plot_errors(X, Y, fit_result, ref_value, ax2)
        plt.tight_layout()
        plt.show()
    elif sum([plot, ploterrors]) == 1:
        fig, ax = plt.subplots(figsize=(5,5))
        if plot:
            plot_fit(X, Y, fit_result, ref_value, ax)
        elif ploterrors:
            plot_errors(X, Y, fit_result, ref_value, ax)
        
        plt.tight_layout()
        plt.show()
    
    return fit_result

def interpolate_new_x(fit_result, new_X, fit_type, ref_value=None):
    if fit_type == 'polynomial':
        if ref_value is not None:
            # Using custom polynomial if a reference value was provided
            return custom_polynomial(fit_result.beta, new_X, ref_value)
        else:
            # Using np.polyval for standard polynomial evaluation
            return np.polyval(fit_result.beta, new_X)
    elif fit_type == 'spline':
        # Using the spline object directly for interpolation
        return fit_result(new_X)
    else:
        raise ValueError("Invalid fit_type. Choose either 'polynomial' or 'spline'.")

def get_energy_at_givenpotential(results, fit_type='polynomial', e_ref=None, order=2, desiredU=0.):
    # e_ref, is the energy of the neutral system (reference)

    # Fit Nelec to U
    if fit_type == 'polynomial':  # Check if it's an ODR output (polynomial fit)
        potential_fit = fit_data(results['U'], results['nelect'], fit_type=fit_type, ref_value=None, order=order )
        parameters = potential_fit.beta[::-1]    # reversing the order of items
        poly = np.poly1d( parameters )
        nelec = poly( desiredU )
    elif fit_type == 'spline':  # Check if it's a spline
        # There are problems if spline tries to fit data with x in decreasing order
        if results['U'][0] > results['U'][-1]:
            X = results['U'][::-1]
            Y = results['nelect'][::-1]
        else:
            X = results['U']
            Y = results['nelect']
        potential_fit = fit_data(X, Y, fit_type=fit_type, ref_value=None, order=order )
        nelec = potential_fit( desiredU )
    else:
        raise ValueError("Unknown fit type. The fit_result object type is not recognized.")

    # Fit E to Nelec
    potential_fit = fit_data(results['nelect'], results['e'], fit_type=fit_type, ref_value=e_ref, order=order)
    if isinstance(potential_fit, odr.Output):
        parameters = np.append(potential_fit.beta[::-1], e_ref)
        poly = np.poly1d( parameters )
        e_pred = poly( nelec )
    elif isinstance(potential_fit, UnivariateSpline):
        e_pred = potential_fit( nelec )
    else:
        raise ValueError("Unknown fit type. The fit_result object type is not recognized.")
    
    #print(f'at U = {desiredU:0.3f}; Nelect = {nelec:0.2f}; Energy = {e_pred:.3f}')
    return e_pred

def plot_errors(X, Y, fit_result, energy_ref, ax):
    # Input, results from the 'extract_corrected_energy_fermie' and polyfit from 'fit_polynomial'
    # energy_ref is the energy of the neutral system
    # ax is the axis pyplot object
    if isinstance(fit_result, odr.Output):  # Check if it's an ODR output (polynomial fit)
        if energy_ref is not None:
            parameters = np.append(fit_result.beta[::-1], energy_ref)
        else:
            parameters = fit_result.beta[::-1]
        #parameters = np.append(fit_result.beta[::-1], energy_ref)
        poly = np.poly1d( parameters )
        poly_y = poly( X )
        errors = poly_y - Y
    elif isinstance(fit_result, UnivariateSpline):  # Check if it's a spline
        errors = fit_result(X) - Y
    else:
        raise ValueError("Unknown fit type. The fit_result object type is not recognized.")

    width = ( max(X) - min(X) ) / len(X)
    ax.bar( X, errors, width=width )
    lower = min(X)
    higher = max(X)
    shift = (higher - lower) / 10
    ax.plot( [lower - shift, higher + shift], [0, 0], '-k', linewidth=1.5 )
    ax.set_xlabel(r'X')
    ax.set_ylabel(r'Y$_{fit}$ - Y$_{DFT}$, eV')
    return ax

def plot_fit(X, Y, fit_result, energy_ref, ax):
    x = np.linspace(min(X), max(X), 100)
    if isinstance(fit_result, odr.Output):  # Check if it's an ODR output (polynomial fit)
        #parameters = np.append(fit_result.beta[::-1], energy_ref)
        if energy_ref is not None:
            parameters = np.append(fit_result.beta[::-1], energy_ref)
        else:
            parameters = fit_result.beta[::-1]
        poly = np.poly1d(parameters)
        ax.plot(x, poly(x), '-k', label="fit")
    elif isinstance(fit_result, UnivariateSpline):  # Check if it's a spline
        ax.plot(x, fit_result(x), '-k', label="spline fit")
    else:
        raise ValueError("Unknown fit type. The fit_result object type is not recognized.")

    ax.plot(X, Y, 'ok', label="original data")
    ax.legend()
    ax.set_xlabel(r'X')
    ax.set_ylabel(r'Y')
    return ax

def print_results(results):
    # Input is the dictionary from "extract_corrected_energy_fermie"
    for i, nelec in enumerate(results['nelect']):
        print(nelec, results['e'][i], results['fe'][i], results['U'][i])