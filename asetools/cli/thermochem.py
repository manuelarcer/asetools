#!/usr/bin/env python

import math
import sys
import argparse
from ase.io import read, write
import numpy as np
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
import logging
import datetime
from io import StringIO
from contextlib import redirect_stdout

def setup_logging():
    """Configure logging to write to both file and console"""
    logfile = f"FREQ_ANALYSIS_{datetime.datetime.now():%Y%m%d_%H%M}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",  # Simple format for frequency analysis output
        handlers=[
            logging.FileHandler(logfile),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def extract_vib_info(lines):
    """Extract vibrational frequencies and eigenvectors from OUTCAR.

    Returns:
        vib: dict with mode info including energy, frequency, and eigenvectors
        vasp6: bool indicating if VASP6 format
    """
    vib = {}
    vasp6 = any('vasp.6' in line for line in lines)

    for i, line in enumerate(lines):
        if ' f  = ' in line or 'f/i=' in line:
            index = line.split()[0]
            vibenergy = float(line.split()[-2]) / 1000
            freq = float(line.split()[-4])

            # Parse eigenvectors (displacements)
            eigenvectors = []
            j = i + 2  # Start from 2 lines after the frequency line
            while j < len(lines):
                parts = lines[j].split()
                # Each atom line has format: x y z dx dy dz
                if len(parts) == 6:
                    try:
                        dx, dy, dz = float(parts[3]), float(parts[4]), float(parts[5])
                        eigenvectors.append([dx, dy, dz])
                        j += 1
                    except (ValueError, IndexError):
                        break
                else:
                    break

            if ' f  = ' in line:
                vib[index] = {'e': vibenergy, 'freq': freq, 'eigenvectors': eigenvectors}
            elif 'f/i=' in line:
                vib[index] = {'e': vibenergy, 'freq': -freq, 'eigenvectors': eigenvectors}

    return vib, vasp6

def compute_corrections(vib, T):
    kB = 0.0000861733
    zpe = sum([v['e'] for v in vib.values() if v['freq'] > 0]) / 2.0
    cpT = sum([v['e'] / (math.exp(v['e'] / kB / T) - 1) for v in vib.values() if v['freq'] >= 0])
    S = sum([kB * (v['e'] / (kB * T * (math.exp(v['e'] / kB / T) - 1)) - math.log(1 - math.exp(-v['e'] / kB / T))) for v in vib.values() if v['freq'] >= 100])
    return zpe, cpT, S

def write_vib_files(vib, vasp6, lines, outcar):
    if vasp6:
        try:
            atoms = read('POSCAR', format='vasp')
        except:
            print('Could not read structure file')
            return
    else:
        atoms = read(outcar, format='vasp-out')

    nat = len(atoms)
    for i, line in enumerate(lines):
        if 'f  =' in line or 'f/i=' in line:
            name1 = line.split()[0]
            name2 = 'IMG' if 'f/i=' in line else 'f'
            shift = [list(map(float, lines[i+count+2].split()[3:])) for count in range(nat)]
            vib = []
            for k in np.arange(-1, 1, 0.2):
                atomscopy = atoms.copy()
                newshift = determineshift(shift, k)
                atomscopy.translate(newshift)
                vib.append(atomscopy)
                write(f'vib_{name1}_{name2}.traj', vib, format='traj')

def determineshift(shiftarray, k=1):
    return [[k * sh for sh in shiftat] for shiftat in shiftarray]

def parse_atom_indices(atom_string):
    """Parse atom indices from string like '0,1,2' or '0-4,7,9-11'.

    Args:
        atom_string: String with comma-separated indices or ranges

    Returns:
        set of atom indices (0-indexed)
    """
    if atom_string is None:
        return None

    indices = set()
    parts = atom_string.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = part.split('-')
            indices.update(range(int(start), int(end) + 1))
        else:
            indices.add(int(part))
    return indices


def calculate_mode_character(eigenvectors, selected_atoms):
    """Calculate what fraction of a mode's displacement comes from selected atoms.

    Args:
        eigenvectors: List of [dx, dy, dz] for each atom
        selected_atoms: Set of atom indices (0-indexed)

    Returns:
        character: Fraction of total displacement from selected atoms (0-1)
    """
    if not selected_atoms or not eigenvectors:
        return 1.0  # If no selection, include all modes

    total_disp = 0.0
    selected_disp = 0.0

    for atom_idx, (dx, dy, dz) in enumerate(eigenvectors):
        displacement = np.sqrt(dx**2 + dy**2 + dz**2)
        total_disp += displacement
        if atom_idx in selected_atoms:
            selected_disp += displacement

    if total_disp == 0:
        return 0.0

    return selected_disp / total_disp


def is_close(a, b, threshold=0.001):
    return abs(a - b) < threshold

def filter_modes_by_character(vib, selected_atoms, threshold, logger):
    """Filter vibrational modes based on character from selected atoms.

    Args:
        vib: Dict of vibrational mode info
        selected_atoms: Set of atom indices to consider
        threshold: Minimum character fraction to keep mode (0-1)
        logger: Logger for output

    Returns:
        filtered_vib: Dict of modes that pass the threshold
        mode_analysis: Dict with character info for all modes
    """
    if selected_atoms is None:
        return vib, {}

    filtered_vib = {}
    mode_analysis = {}

    logger.info('')
    logger.info('*************************************************')
    logger.info('----------  Mode Selection Analysis  ------------')
    logger.info('*************************************************')
    logger.info(f'Selected atoms: {sorted(selected_atoms)}')
    logger.info(f'Character threshold: {threshold:.2f}')
    logger.info('')
    logger.info('{:>3} {:>12} {:>10} {:>8}'.format('#', 'Freq[cm-1]', 'Character', 'Keep?'))
    logger.info('-' * 40)

    for key, value in vib.items():
        character = calculate_mode_character(value['eigenvectors'], selected_atoms)
        mode_analysis[key] = character
        keep = character >= threshold

        status = 'YES' if keep else 'NO'
        logger.info('{:>3} {:12.1f} {:10.3f} {:>8}'.format(
            key, value['freq'], character, status))

        if keep:
            filtered_vib[key] = value

    logger.info('-' * 40)
    logger.info(f'Kept {len(filtered_vib)}/{len(vib)} modes')
    logger.info('')

    return filtered_vib, mode_analysis


def remove_repeated_energies(energies, atoms, geom, logger):  # It happens that for gas-phase I found repeated vib
    # energies, is the list of vibrational energies
    # atoms, is the atoms object
    # geom, is the geometry of the molecule : 'linear' or 'nonlinear'
    # logger, is the logger object for output
    if geom == 'linear':
        numvib = 3 * len(atoms) - 5
    elif geom == 'nonlinear':
        numvib = 3 * len(atoms) - 6

    test = len(energies) == numvib
    if test:
        logger.info('Number of energies as expected')
        return energies
    else:
        logger.info('')
        logger.info('WARNING: The number of energies is different from the expected')
        logger.info('Removing similar vibrational energies')
        diff = len(energies) - numvib
        invertenergies = [energies[i] for i in range(len(energies)-1, -1, -1)]
        logger.info(invertenergies)
        newenergies = []
        if len(invertenergies) > numvib:
            for i, e in enumerate(invertenergies):
                if i == 0:
                    newenergies.append(e)
                if i > 0 and not is_close(e, invertenergies[i-1]):
                    newenergies.append(e)
        return newenergies

def main():
    parser = argparse.ArgumentParser(
        description='Calculate corrections to the energy for GFE.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard usage (all atoms)
  thermochem.py OUTCAR --temp 298

  # Select only adsorbate atoms (indices 0-4)
  thermochem.py OUTCAR --select-atoms 0-4 --threshold 0.5

  # Select specific atoms with custom threshold
  thermochem.py OUTCAR --select-atoms 0,1,2,3,4 --threshold 0.7
        """)
    parser.add_argument('outcar', type=str, help='OUTCAR file path')
    parser.add_argument('--writevib', type=str, choices=['y', 'n'], default='n', help='Whether to write vibrations (y/n)')
    parser.add_argument('--temp', type=float, default=298, help='Temperature in K')
    parser.add_argument('--gas', type=bool, default=False, help='Ideal Gas-phase Thermo')
    parser.add_argument('--symnum', type=int, default=2, help='Symmetry Number')
    parser.add_argument('--geom', type=str, choices=['monoatomic', 'linear', 'nonlinear'], default='linear', help='Geometry of molecule')
    parser.add_argument('--pressure', type=float, default=1.0, help='Gas pressure (bar)')
    parser.add_argument('--select-atoms', type=str, default=None,
                        help='Atom indices to consider (e.g., "0-4" or "0,1,2,3,4"). Only modes dominated by these atoms will be used.')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Minimum character threshold (0-1) for mode selection. Only modes with >= this fraction of displacement from selected atoms are kept. Default: 0.5')
    args = parser.parse_args()

    # Setup logging
    logger = setup_logging()

    Pa = 100000.    # 1 bar in Pa

    lines = open(args.outcar, 'r').readlines()
    vib_all, vasp6 = extract_vib_info(lines)

    # Parse selected atoms and filter modes if requested
    selected_atoms = parse_atom_indices(args.select_atoms)
    if selected_atoms is not None:
        vib, mode_analysis = filter_modes_by_character(vib_all, selected_atoms, args.threshold, logger)
    else:
        vib = vib_all
        mode_analysis = {}

    zpe, cpT, S = compute_corrections(vib, args.temp)

    ## Display results
    logger.info('-----------------------')
    logger.info('{:>3} {:>12} {:>6}'.format('#', 'Freq[cm-1]', 'E[eV]'))
    logger.info('-----------------------')
    for key, value in vib.items():
        logger.info('{:>3} {:12.1f} {:6.3f}'.format(key, value['freq'], value['e']))
    logger.info('-----------------------')

    e_dft = read(args.outcar, format='vasp-out', index=0).get_potential_energy()
    logger.info('E, eV = {:.3f}'.format(e_dft))
    logger.info('ZPE, eV = {:.3f}'.format(zpe))
    logger.info('S, eV/K = {:.6f}'.format(S))
    logger.info('Temperature is = {:.1f}'.format(args.temp))
    logger.info('Thermal correction (0->T), eV = {:.3f}'.format(cpT))
    logger.info('Entropy correction (-S*T), eV = {:.3f}'.format(-S * args.temp))
    logger.info('\n')
    logger.info('All values together: E_tot   E_ZPE   CpT   -S*T')
    logger.info('{:.3f} {:.3f} {:.3f} {:.3f}'.format(e_dft, zpe, cpT, -S * args.temp))
    logger.info('\n')

    logger.info('*************************************************')
    logger.info('---------------  Harmonic Limit  ----------------')
    logger.info('*************************************************')
    logger.info('')
    energies = [value['e'] for value in vib.values() if value['freq'] >= 0]
    harm_lim = HarmonicThermo(energies, potentialenergy=0.0)
    # The following line already gives the details of the S, CpT and G energy calculations
    # Capture the verbose output from ASE and log it
    f = StringIO()
    with redirect_stdout(f):
        harm_lim.get_helmholtz_energy(args.temp, verbose=True)
    helmholtz_output = f.getvalue()
    # Print to console and log to file
    print(helmholtz_output, end='')
    for line in helmholtz_output.splitlines():
        logger.info(line)


    ##### Ideal Gas Thermo
    if args.gas:
        energies = [value['e'] for value in vib.values() if value['freq'] >= 0]
        logger.info(energies)
        atoms = read(args.outcar, format='vasp-out', index=0)
        energies = remove_repeated_energies(energies, atoms, args.geom, logger)
        thermo = IdealGasThermo(vib_energies=energies,
                        atoms=atoms,
                        #potentialenergy=atoms.get_potential_energy(),
                        geometry=args.geom,      # monoatomic, linear, nonlinear
                        symmetrynumber=args.symnum,   # CO2: 2, H2O: 2, CO: 1, H2: 2
                        spin=0)     # Different for radicals or unpaired electrons
        G = thermo.get_gibbs_energy(temperature=args.temp, pressure=args.pressure * Pa)

    if args.writevib == 'y':
        write_vib_files(vib, vasp6, lines, args.outcar)

if __name__ == "__main__":
    main()
