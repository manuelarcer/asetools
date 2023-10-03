# Module for analysis of potential-dependent calculations

from asetools.analysis import check_outcar_convergence, check_energy_and_maxforce
from asetools.doscar_analysis import extract_fermi_e
from ase.io import read
import os

def get_num_elect(outcar):
    with open(outcar, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if 'NELECT' in line:
                nelect = float( line.split()[-1] )  # Float since it may be fractional
                break
    return nelect

def extract_corrected_energy_fermie(calc_minus, calc_zero, calc_plus):
    # The quadratic fit relies on calculation of the derivatives of the energy with respect to the num of elec
    # these calculations should be done around the neutral system (prob better if the difference in num elec is small)
    # Input should be the folders of the calculations:
    #   with decreased num of elec (calc), neutral system (calc_zero), and the one with more electrons (calc_plus)

    results = {'nelect': [], 'e': [], 'fe': []}
    for folder in [calc_minus, calc_zero, calc_plus]:
        convergence = check_outcar_convergence(folder+'/OUTCAR', verbose=False)
        
        if convergence:
            print(f'GOOD news, calculation at {folder} CONVERGED')
            atoms = read(folder+'/OUTCAR', format='vasp-out', index=-1)
            nelect = get_num_elect(folder+'/OUTCAR')
            energy_original = atoms.get_potential_energy()
            results['nelect'].append(nelect)

            #energy, maxforce = check_energy_and_maxforce(folder+'/OUTCAR', magmom=False, verbose=False)
            
            fermie_original = extract_fermi_e(folder+'/DOSCAR')

            try:
                with open(folder+'/vasp.out', 'r') as file:
                    lines = file.readlines()
                    for line in lines[-20]:
                        if 'FERMI_SHIFT' in line:
                            fermi_shift = float( line.split('=')[-1] )
                            break
            
            except:
                print('WARNING ! Error while reading "vasp.out", no FERMI_SHIFT stored')
                fermi_shift = None

            fermie_corr = fermie_original + fermi_shift
            energy_corr = energy_original + nelect * fermi_shift

            results['e'].append(energy_corr)
            results['fe'].append(fermie_corr)
            

        else:
            print(f'Oh NOOOO!, something wrong with calculation at {folder}')
    
    return results

        


