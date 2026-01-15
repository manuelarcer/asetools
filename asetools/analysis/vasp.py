# File containing the functions to analyze OUTCARs and VASP runs

from ase.io import read
import pandas as pd
import numpy as np
import re
    
def check_outcar_convergence(outcar, verbose=False):
    try:
        out = open(outcar, 'r')
    except:
        if verbose:
            print('check_outcar_convergence --> OUTCAR file not found or damaged')
        return False, ''
    lines = out.readlines()
    ibrion = None
    nsw = None
    opt = False
    convergence = False
    vasp = ''
    has_general_timing = False

    for line in lines:
        if 'vasp.6' in line:
            vasp = 'vasp6'
        if 'vasp.5' in line:
            vasp = 'vasp5'
        if 'IBRION' in line:
            ibrion = int( line.split()[2] )
            if ibrion == 1 or ibrion == 2 or ibrion == 3:
                opt = True
            elif verbose:
                print(f'Not an OPTIMIZATION job')
                opt = False
            else:
                opt = False
        elif 'NSW' in line:
            nsw = int( line.split()[2] )
            if nsw > 0:
                opt = True
        elif 'reached required accuracy' in line and opt:
            convergence = True
        elif 'General timing and accounting informations for this job:' in line:
            has_general_timing = True
            if not opt:
                convergence = True

    # Special handling for ASE optimizers: IBRION=-1 with NSW>0
    # ASE controls optimization externally, VASP just runs single-point calculations
    # Check if VASP completed normally (has General timing section)
    if ibrion == -1 and nsw > 0 and has_general_timing:
        convergence = True
        if verbose:
            print('ASE optimizer job (IBRION=-1, NSW>0) --> Assuming converged (VASP completed normally)')

    if verbose:
        print(f'IBRION --> {ibrion}, NSW --> {nsw}')

    if ibrion == 1 or ibrion == 2 or ibrion == 3:
        if nsw > 0 and convergence and verbose:
            print('Optimization Job --> CONVERGED')
        elif nsw > 0 and not convergence and verbose:
            print('Optimization Job --> *NOT* converged')

    return convergence, vasp

def check_energy_and_maxforce(outcar, magmom=False, verbose=False):

    convergence, vasp = check_outcar_convergence(outcar, verbose=verbose)
    try:
        atoms = read(outcar, format='vasp-out', index=-1)
        energy = atoms.get_potential_energy()
        vecforces = atoms.get_forces()
        forces = [np.linalg.norm(f) for f in vecforces]
        maxforce = max( forces )
        if magmom:
            mm = atoms.get_magnetic_moment() 
            return energy, maxforce, mm
        else:
            return energy, maxforce
    except:
        print('Missing or damaged OUTCAR file')
        return 9999.99, 9.99

def extract_magnetic_moments(outcar, listatoms, verbose=False):
    # listatoms: is a list with the indexes of atoms of interest
    convergence, vasp = check_outcar_convergence(outcar, verbose=verbose)
    try:
        atoms = read(outcar, format='vasp-out', index=-1)
        #energy = atoms.get_potential_energy()
        mm = atoms.get_magnetic_moment()
        return [round(atoms.get_magnetic_moments()[i],2) for i in listatoms]
    except:
        print('Missing or damaged OUTCAR file')
        return []

def get_parameter_from_run(outcar, check_converg=True, parameter='ISIF'):
    # First check convergence
    convergence = 'Convergence-not-Checked'
    if check_converg:
        convergence, vasp = check_outcar_convergence(outcar, verbose=False)

    pattern = rf'^\s*{re.escape(parameter)}\s*=\s*(\S+)'
    with open(outcar, 'r') as f:
        for line in f:
            m = re.match(pattern, line)
            if m:
                val = m.group(1)
                # try int → float → fallback to string
                for caster in (int, float):
                    try:
                        return caster(val), convergence
                    except ValueError:
                        continue
                return val, convergence  # non-numeric string
    if check_converg:
        raise ValueError(f"Parameter '{parameter}' not found in {outcar!r}")
    return None, convergence


def classify_calculation_type(outcar_path, calc_dir):
    """
    Classify VASP calculation type based on OUTCAR parameters and auxiliary files.

    Args:
        outcar_path: Path to OUTCAR file
        calc_dir: Directory containing calculation files

    Returns:
        str: Calculation type - one of:
            'dimer', 'neb', 'ase-optimizer', 'finite-diff', 'md',
            'single-point', 'cell-relax', 'optimization', 'unknown'
    """
    import os
    import glob

    # Check for special files first (highest priority)
    # DIMER calculation
    if os.path.exists(os.path.join(calc_dir, 'DIMER.log')):
        return 'dimer'

    # ASE optimizer logs
    ase_optimizer_logs = ['BFGS.log', 'FIRE.log', 'LBFGS.log', 'GPMIN.log',
                          'MDMIN.log', 'QUASINEWTON.log']
    for log_name in ase_optimizer_logs:
        if os.path.exists(os.path.join(calc_dir, log_name)):
            return 'ase-optimizer'

    # Check OUTCAR for IMAGES tag (NEB calculation)
    try:
        with open(outcar_path, 'r') as f:
            for line in f:
                if 'IMAGES' in line and '=' in line:
                    # Extract IMAGES value to confirm it's not just a comment
                    parts = line.split('=')
                    if len(parts) > 1:
                        try:
                            images = int(parts[1].split()[0])
                            if images > 0:
                                return 'neb'
                        except (ValueError, IndexError):
                            pass
    except (OSError, IOError):
        pass

    # Extract VASP parameters
    try:
        ibrion, _ = get_parameter_from_run(outcar_path, check_converg=False, parameter='IBRION')
        nsw, _ = get_parameter_from_run(outcar_path, check_converg=False, parameter='NSW')
        isif, _ = get_parameter_from_run(outcar_path, check_converg=False, parameter='ISIF')
        nfree, _ = get_parameter_from_run(outcar_path, check_converg=False, parameter='NFREE')
    except Exception:
        return 'unknown'

    # Convert None to defaults
    if ibrion is None:
        ibrion = -1
    if nsw is None:
        nsw = 0
    if isif is None:
        isif = 0
    if nfree is None:
        nfree = 0

    # Apply classification rules (priority order)
    # Frequency/phonon calculation
    if nfree > 0:
        return 'finite-diff'

    # Molecular dynamics
    if ibrion == 0:
        return 'md'

    # Single-point calculation
    if nsw <= 1:
        return 'single-point'

    # Cell relaxation (ISIF allows cell shape/volume changes)
    if ibrion in [1, 2, 3] and isif in [3, 4, 5, 6, 7]:
        return 'cell-relax'

    # Geometry optimization (ISIF only allows atomic positions)
    if ibrion in [1, 2, 3, 5, 6]:
        return 'optimization'

    return 'unknown'


def find_initial_structure(calc_dir, pattern='*.vasp'):
    """
    Find initial structure file in calculation directory using pattern matching.

    Args:
        calc_dir: Directory to search for initial structure
        pattern: Glob pattern for initial structure file (default: '*.vasp')

    Returns:
        str: Absolute path to initial structure file

    Raises:
        FileNotFoundError: If no structure file matching pattern or POSCAR fallback found
    """
    import os
    import glob

    # Try pattern matching first
    search_pattern = os.path.join(calc_dir, pattern)
    matches = sorted(glob.glob(search_pattern))

    if matches:
        return os.path.abspath(matches[0])

    # Fallback to POSCAR
    poscar_path = os.path.join(calc_dir, 'POSCAR')
    if os.path.exists(poscar_path):
        return os.path.abspath(poscar_path)

    # No structure file found
    raise FileNotFoundError(
        f"No initial structure found in {calc_dir} "
        f"(searched for pattern '{pattern}' and fallback 'POSCAR')"
    )


def extract_comprehensive_metadata(outcar_path, incar_path=None, potcar_path=None):
    """
    Extract comprehensive metadata from VASP calculation files.

    Extracts VASP parameters, energies, forces, convergence status, INCAR content,
    and POTCAR information from OUTCAR and related files.

    Args:
        outcar_path: Path to OUTCAR file
        incar_path: Path to INCAR file (default: auto-detect from OUTCAR directory)
        potcar_path: Path to POTCAR file (optional, TITEL extracted from OUTCAR)

    Returns:
        dict: Dictionary containing:
            - CalcType (str): Calculation type classification
            - Formula (str): Chemical formula
            - ENCUT, KSPACING, EDIFF, EDIFFG (float): VASP parameters
            - GGA (str): Exchange-correlation functional
            - IBRION, ISPIN, NSW, ISIF, NFREE (int): VASP parameters
            - Energy (float): Final energy in eV
            - TotMagMom (float): Total magnetic moment
            - MaxForce (float): Maximum force in eV/Å
            - Converged (bool): Convergence status
            - VASPVersion (str): VASP version (vasp5/vasp6)
            - INCAR_full (str or None): Complete INCAR file content
            - POTCAR_info (list): List of POTCAR TITEL lines
    """
    import os

    calc_dir = os.path.dirname(outcar_path)
    metadata = {}

    # Extract basic VASP parameters
    param_names = ['ENCUT', 'KSPACING', 'EDIFF', 'EDIFFG', 'GGA',
                   'IBRION', 'ISPIN', 'NSW', 'ISIF', 'NFREE']

    for param in param_names:
        try:
            value, _ = get_parameter_from_run(outcar_path, check_converg=False, parameter=param)
            metadata[param] = value
        except Exception:
            metadata[param] = None

    # Extract energy, forces, and magnetic moment
    try:
        energy, maxforce, magmom = check_energy_and_maxforce(outcar_path, magmom=True, verbose=False)
        metadata['Energy'] = energy
        metadata['MaxForce'] = maxforce
        metadata['TotMagMom'] = magmom
    except Exception:
        metadata['Energy'] = None
        metadata['MaxForce'] = None
        metadata['TotMagMom'] = None

    # Check convergence and get VASP version
    try:
        converged, vasp_version = check_outcar_convergence(outcar_path, verbose=False)
        metadata['Converged'] = converged
        metadata['VASPVersion'] = vasp_version
    except Exception:
        metadata['Converged'] = None
        metadata['VASPVersion'] = None

    # Extract chemical formula from final structure
    try:
        atoms = read(outcar_path, format='vasp-out', index=-1)
        metadata['Formula'] = atoms.get_chemical_formula()
    except Exception:
        metadata['Formula'] = None

    # Classify calculation type
    try:
        metadata['CalcType'] = classify_calculation_type(outcar_path, calc_dir)
    except Exception:
        metadata['CalcType'] = 'unknown'

    # Read full INCAR content
    if incar_path is None:
        incar_path = os.path.join(calc_dir, 'INCAR')

    try:
        with open(incar_path, 'r') as f:
            metadata['INCAR_full'] = f.read()
    except (OSError, IOError):
        metadata['INCAR_full'] = None

    # Extract POTCAR TITEL lines from OUTCAR
    potcar_info = []
    try:
        with open(outcar_path, 'r') as f:
            for line in f:
                if 'TITEL' in line and '=' in line:
                    # Extract TITEL line format: "   TITEL  = PAW_PBE Zn 06Sep2000"
                    parts = line.split('=', 1)
                    if len(parts) > 1:
                        titel = parts[1].strip()
                        potcar_info.append(titel)
    except (OSError, IOError):
        pass

    metadata['POTCAR_info'] = potcar_info

    return metadata




    