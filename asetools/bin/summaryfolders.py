#!/usr/bin/env python

import argparse
from ase.io import read
import pandas as pd
import numpy as np
import glob, os
from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence, get_parameter_from_run


# Parameters to extract for pyatoms mode
PYATOMS_PARAMS = ['GGA', 'ENCUT', 'EDIFF', 'EDIFFG', 'ISPIN', 'ISMEAR', 'SIGMA', 'IVDW']


def get_pyatoms_steps(folder):
    """Get sorted list of step folders in a pyatoms calculation folder.

    Args:
        folder: Path to pyatoms calculation folder (e.g., 'H2O/')

    Returns:
        list: Sorted list of step folder paths (e.g., ['H2O/step_00/', 'H2O/step_01/'])
    """
    step_folders = sorted(glob.glob(os.path.join(folder, 'step_*/')))
    return step_folders


def extract_step_parameters(step_folder):
    """Extract VASP parameters from a pyatoms step OUTCAR.

    Args:
        step_folder: Path to step folder containing OUTCAR

    Returns:
        dict: Dictionary with parameter names as keys and values
    """
    outcar_path = os.path.join(step_folder, 'OUTCAR')
    params = {}

    if not os.path.exists(outcar_path):
        return {p: 'N/A' for p in PYATOMS_PARAMS}

    for param in PYATOMS_PARAMS:
        try:
            value, _ = get_parameter_from_run(outcar_path, check_converg=False, parameter=param)
            params[param] = value if value is not None else 'N/A'
        except Exception:
            params[param] = 'N/A'

    return params


def extract_step_energy(step_folder):
    """Extract final energy from a pyatoms step OUTCAR.

    Args:
        step_folder: Path to step folder containing OUTCAR

    Returns:
        float or str: Energy value or 'N/A' if not found
    """
    outcar_path = os.path.join(step_folder, 'OUTCAR')

    if not os.path.exists(outcar_path):
        return 'N/A'

    try:
        energy, _ = check_energy_and_maxforce(outcar_path, magmom=False, verbose=False)
        return round(energy, 4)
    except Exception:
        return 'N/A'


def check_pyatoms_finished(folder):
    """Check if a pyatoms calculation has finished successfully.

    Args:
        folder: Path to pyatoms calculation folder

    Returns:
        tuple: (finished: bool, final_energy: float or None, final_fmax: float or None)
    """
    loginfo_path = os.path.join(folder, 'log.info')

    if not os.path.exists(loginfo_path):
        return False, None, None

    finished = False
    final_energy = None
    final_fmax = None

    with open(loginfo_path, 'r') as f:
        lines = f.readlines()

    process_completed = False
    for line in lines:
        # Check for procedure or job completion
        if ('VaspGeomOptProcedure' in line or 'VaspGeomOptJob' in line) and 'completed successfully' in line:
            process_completed = True

        # Check for global processes closing (final confirmation)
        if process_completed and 'Closing global processes' in line:
            finished = True

        # Extract energy and force
        if 'energy' in line and 'force' in line:
            try:
                final_energy = float(line.split()[-3].split(',')[0])
                final_fmax = float(line.split()[-1].rstrip('.'))
            except (ValueError, IndexError):
                pass

    return finished, final_energy, final_fmax


def run_pyatoms_mode():
    """Run summary in pyatoms mode - analyze multi-step pyatoms calculations."""
    folders = sorted(glob.glob('*/'))

    # Filter to only folders containing pyatoms calculations (have log.info)
    pyatoms_folders = [f for f in folders if os.path.exists(os.path.join(f, 'log.info'))]

    if not pyatoms_folders:
        print("No pyatoms calculations found (no folders with log.info)")
        return

    # Get step count from first folder to establish expected structure
    first_steps = get_pyatoms_steps(pyatoms_folders[0])
    n_steps = len(first_steps)

    if n_steps == 0:
        print(f"No step folders found in {pyatoms_folders[0]}")
        return

    # Verify all folders have the same number of steps
    for folder in pyatoms_folders:
        folder_steps = get_pyatoms_steps(folder)
        if len(folder_steps) != n_steps:
            print(f"Error: {folder} has {len(folder_steps)} steps, expected {n_steps}")
            print("All pyatoms calculations must have the same number of steps.")
            return

    # === Print Parameters Table ===
    print("=" * 80)
    print("STEP PARAMETERS (from first calculation)")
    print("=" * 80)

    params_data = {'Step': []}
    for param in PYATOMS_PARAMS:
        params_data[param] = []

    for i, step_folder in enumerate(first_steps):
        params = extract_step_parameters(step_folder)
        params_data['Step'].append(i)
        for param in PYATOMS_PARAMS:
            params_data[param].append(params.get(param, 'N/A'))

    params_df = pd.DataFrame(params_data)
    print(params_df.to_string(index=False))
    print()

    # === Print Energy Summary Table ===
    print("=" * 80)
    print("ENERGY SUMMARY")
    print("=" * 80)

    energy_data = {'Config': []}
    for i in range(n_steps):
        energy_data[f'E_step{i}'] = []
    energy_data['Converged'] = []

    not_converged = []

    for folder in pyatoms_folders:
        folder_name = folder.rstrip('/')
        step_folders = get_pyatoms_steps(folder)

        # Check if calculation finished
        finished, _, _ = check_pyatoms_finished(folder)

        if not finished:
            print(f"{folder_name}: Not finished")
            not_converged.append(folder_name)
            continue

        energy_data['Config'].append(folder_name)
        energy_data['Converged'].append(True)

        # Extract energy from each step
        for i, step_folder in enumerate(step_folders):
            energy = extract_step_energy(step_folder)
            energy_data[f'E_step{i}'].append(energy)

    # Calculate relative energies based on final step
    if energy_data['Config']:
        final_step_col = f'E_step{n_steps - 1}'
        final_energies = [e for e in energy_data[final_step_col] if e != 'N/A']

        if final_energies:
            min_energy = min(final_energies)
            energy_data['Rel.E'] = []
            for e in energy_data[final_step_col]:
                if e == 'N/A':
                    energy_data['Rel.E'].append('N/A')
                else:
                    energy_data['Rel.E'].append(round(e - min_energy, 4))
        else:
            energy_data['Rel.E'] = ['N/A'] * len(energy_data['Config'])

    energy_df = pd.DataFrame(energy_data)
    print(energy_df.to_string(index=True, max_rows=None, max_cols=None, line_width=1000))
    print()

    print("Not converged/finished:")
    print(' '.join(not_converged) if not_converged else "(none)")

    # Write to summary_pyatoms.log
    with open('summary_pyatoms.log', 'w') as f:
        f.write("STEP PARAMETERS (from first calculation)\n")
        f.write("=" * 80 + "\n")
        f.write(params_df.to_string(index=False))
        f.write("\n\n")
        f.write("ENERGY SUMMARY\n")
        f.write("=" * 80 + "\n")
        f.write(energy_df.to_string(index=True, max_rows=None, max_cols=None, line_width=1000))
        f.write("\n\n")
        f.write("Not converged/finished:\n")
        f.write(' '.join(not_converged) if not_converged else "(none)")
        f.write("\n")

    print(f"\nOutput saved to summary_pyatoms.log")

def is_summary_up_to_date(magmom_requested=False):
    """Check if summary.log exists and is newer than all folders, and if magmom flag status matches"""
    if not os.path.exists('summary.log'):
        return False

    # Check if magnetic moments flag status has changed
    with open('summary.log', 'r') as f:
        log_content = f.read()
        log_has_magmom = 'MagMom' in log_content

    # Only regenerate if magmom is requested but log doesn't have it
    # If log has magmom but not requested, still use existing log (saves time)
    if magmom_requested and not log_has_magmom:
        return False

    summary_mtime = os.path.getmtime('summary.log')
    folders = glob.glob('*/')

    for folder in folders:
        # Check if any OUTCAR in folders is newer than summary.log
        outcar_path = os.path.join(folder, 'OUTCAR')
        if os.path.exists(outcar_path):
            if os.path.getmtime(outcar_path) > summary_mtime:
                return False
        # Also check for log.info files (pyatoms)
        loginfo_path = os.path.join(folder, 'log.info')
        if os.path.exists(loginfo_path):
            if os.path.getmtime(loginfo_path) > summary_mtime:
                return False

    return True

def main():
    ## Block from chatgpt
    parser = argparse.ArgumentParser(description='Process folders.')
    parser.add_argument('-m', '--magmom', action='store_true',
                        help='extract and present magnetic moments')
    parser.add_argument('-f', '--fast', action='store_true',
                        help='fast mode (reading out file)')
    parser.add_argument('-p', '--pyatoms', action='store_true',
                        help='pyatoms mode: analyze multi-step pyatoms calculations')
    args = parser.parse_args()
    ##

    # Handle pyatoms mode separately
    if args.pyatoms:
        run_pyatoms_mode()
        return

    # Check if summary.log is up to date
    if is_summary_up_to_date(magmom_requested=args.magmom):
        print("Summary is up to date. Reading from summary.log:")
        print()
        with open('summary.log', 'r') as f:
            print(f.read())
        return

    folders = glob.glob('*/')
    if args.magmom:
        dic = {'Config': [], 'ISIF':[], 'Converged':[], 'ENCUT': [], 'Target-fmax': [], 'MaxForce': [], 'Energy':[], 'MagMom':[]}
    else:
        dic = {'Config': [], 'ISIF':[], 'Converged':[], 'ENCUT': [], 'Target-fmax': [], 'MaxForce': [], 'Energy':[]}

    # Define alternative filenames to look for when fast mode is enabled
    alternative_filenames = ['vasp.out', 'out.txt']
    not_converged = []
    for f in sorted(folders):
        ispyatoms = False
        if f+'log.info' in glob.glob(f+'log.info'):
            ispyatoms = True
            finishedpyatoms = False
            finishedprocess = False
            # Read log.info to check if the global processes have finished
            logfile = open(f+'log.info', 'r')
            lines = logfile.readlines()
            optsteps = []; e = []; fmax = []
            for line in lines:
                if 'pyatoms.jobs.job]' in line and 'Procedure step' in line:
                    optsteps.append(int(line.split()[-5]))
                #else:
                #    optsteps = 1
                # I think if multistep Opt there will be multiple 'VaspGeomOptProcedure' and 'VaspGeomOptJob' in log.info
                # for a single step procedure, there will be only one 'VaspGeomOptJob'
                if 'VaspGeomOptProcedure' in line and 'completed successfully' in line:      # this work for GeomOpt only
                    finishedprocess = True
                elif 'VaspGeomOptJob' in line and 'completed successfully' in line:           # this work for GeomOpt only
                    finishedprocess = True

                if finishedprocess:
                    if 'Closing global processes' in line:
                        finishedpyatoms = True
                
                if 'energy' in line:
                    e.append( float(line.split()[-3].split(',')[0]) )   # example: energy -481.276399, force 0.013.
                    fmax.append( float(line.split()[-1][:-1]) )
            if not finishedpyatoms:
                print(f, 'Pyatoms Job Not Finished')
                continue
            else:
                print(f, 'Pyatoms Job Finished')
                dic['Config'].append(f)
                dic['Converged'].append((True, 'pyatoms'))
                dic['ENCUT'].append('N/A')
                dic['Target-fmax'].append('N/A')
                dic['MaxForce'].append( round(fmax[-1], 3) )
                dic['Energy'].append( round(e[-1], 3) )

        ########## Block for fast mode ########
        elif args.fast:
            converged = fast_mode_check(f, alternative_filenames)
            if converged is None:
                # Timeout or error reading files (likely OneDrive not downloaded)
                print(f, 'skipped (files not accessible)')
                not_converged.append(f)
            elif converged:
                dic['Config'].append(f)
                dic['Converged'].append(converged)
                dic['ENCUT'].append('N/A')
                dic['Target-fmax'].append('N/A')
                dic['MaxForce'].append('N/A')
                dic['Energy'].append('N/A')
            else:
                print(f, 'not converged')
                not_converged.append(f)
            continue
        #########################################

        elif not ispyatoms:
            if os.path.exists(f + 'OUTCAR'):
                try:
                    # First check if an ASE optimizer was used (check for optimizer log files)
                    ase_converged, ase_fmax, optimizer_type = check_ase_optimizer_convergence(f)

                    if ase_converged is not None:
                        # ASE optimizer was used - use its convergence status
                        # (ASE optimizers use constraint-adjusted forces, not raw VASP forces)
                        converged_status = ase_converged
                        final_fmax = ase_fmax
                        if not converged_status:
                            print(f, f'not converged ({optimizer_type}: fmax={ase_fmax:.4f})')
                            not_converged.append(f)
                    else:
                        # No ASE optimizer - use standard OUTCAR convergence check
                        converged = check_outcar_convergence(f + 'OUTCAR', verbose=False)
                        converged_status = converged[0]
                        energy, maxforce_outcar = check_energy_and_maxforce(f + 'OUTCAR', magmom=False, verbose=False)
                        final_fmax = maxforce_outcar
                        if not converged_status:
                            print(f, 'not converged')
                            not_converged.append(f)

                    # Get energy and forces
                    if args.magmom:
                        energy, maxforce, magmom = check_energy_and_maxforce(f + 'OUTCAR', magmom=args.magmom, verbose=False)
                        dic['MagMom'].append(round(magmom, 3))
                    else:
                        energy, maxforce = check_energy_and_maxforce(f + 'OUTCAR', magmom=False, verbose=False)

                    isif, _ = get_parameter_from_run(f + 'OUTCAR', check_converg=False, parameter='ISIF')
                    encut, _ = get_parameter_from_run(f + 'OUTCAR', check_converg=False, parameter='ENCUT')
                    ediffg, _ = get_parameter_from_run(f + 'OUTCAR', check_converg=False, parameter='EDIFFG')

                    # Convert EDIFFG to positive value for target fmax (VASP uses negative for forces)
                    target_fmax = abs(ediffg) if ediffg is not None else 'N/A'

                    dic['Config'].append(f)
                    dic['ISIF'].append(isif)
                    dic['Converged'].append(converged_status)
                    dic['ENCUT'].append(encut)
                    dic['Target-fmax'].append(target_fmax)
                    # Use ASE fmax if available, otherwise use OUTCAR maxforce
                    dic['MaxForce'].append(round(final_fmax if ase_converged is not None else maxforce, 3))
                    dic['Energy'].append(round(energy, 3))


                except ValueError as e:
                    print(f'Error processing {f}: {e}. OUTCAR may be incomplete or damaged.')
                    not_converged.append(f)
            else:
                print('No OUTCAR in', f)
                not_converged.append(f)

    dic['Rel.E'] = []
    for e in dic['Energy']:
        if e == 'N/A':
            dic['Rel.E'].append('N/A')
        else:
            dic['Rel.E'].append(e - min(dic['Energy'], default=0))  # default=0 to handle empty Energy list

    df = pd.DataFrame.from_dict(dic)

    # Use to_string() to ensure all columns display in one line
    print(df.to_string(index=True, max_rows=None, max_cols=None, line_width=1000))
    print()
    not_converged = [f.split('/')[0] for f in not_converged]
    print('Not converged:')
    print(' '.join(not_converged))
    
    # Write output to summary.log
    with open('summary.log', 'w') as f:
        f.write(df.to_string(index=True, max_rows=None, max_cols=None, line_width=1000))
        f.write('\n\n')
        f.write('Not converged:\n')
        f.write(' '.join(not_converged))
        f.write('\n')

def check_ase_optimizer_convergence(folder, fmax_threshold=0.02):
    """
    Check if an ASE optimizer converged by reading optimizer log files.

    ASE optimizers (BFGS, FIRE, LBFGS, etc.) write their own log files that contain
    the actual convergence information based on constraint-adjusted forces.

    Args:
        folder: Path to folder containing optimizer log
        fmax_threshold: Maximum force threshold for convergence (default: 0.02 eV/Å)

    Returns:
        tuple: (converged: bool or None, final_fmax: float or None, optimizer_type: str or None)
               Returns (None, None, None) if no ASE optimizer log found
    """
    # Common ASE optimizer log files
    optimizer_logs = ['BFGS.log', 'FIRE.log', 'LBFGS.log', 'GPMIN.log',
                      'MDMIN.log', 'QUASINEWTON.log', 'DIMER.log']

    for log_name in optimizer_logs:
        log_path = os.path.join(folder, log_name)
        if os.path.exists(log_path):
            try:
                with open(log_path, 'r') as f:
                    lines = f.readlines()

                # Parse the last optimization step
                # Format: "BFGS:  step  time  energy  fmax"
                optimizer_type = log_name.replace('.log', '')
                final_fmax = None

                for line in reversed(lines):
                    line = line.strip()
                    if line.startswith(optimizer_type + ':'):
                        parts = line.split()
                        try:
                            # Last column is typically fmax
                            final_fmax = float(parts[-1])
                            break
                        except (ValueError, IndexError):
                            continue

                if final_fmax is not None:
                    converged = final_fmax <= fmax_threshold
                    return (converged, final_fmax, optimizer_type)

            except (OSError, IOError) as e:
                continue

    # No ASE optimizer log found
    return (None, None, None)


def fast_mode_check(f, alternative_filenames):
    """Check convergence in fast mode with error handling for OneDrive files.

    Args:
        f: Folder path
        alternative_filenames: List of alternative output filenames to check

    Returns:
        Convergence status (True/False) or None if file read fails/times out
    """
    try:
        # Find and open output file
        foundout = False
        vaspout = None
        for alt in alternative_filenames:
            filepath = f + alt
            if os.path.exists(filepath):
                try:
                    # Try to open with timeout protection
                    vaspout = open(filepath, 'r')
                    # Test read to check if file is accessible (will fail for OneDrive placeholders)
                    _ = vaspout.readline()
                    vaspout.seek(0)  # Reset to beginning
                    foundout = True
                    break
                except (OSError, TimeoutError, IOError) as e:
                    if vaspout:
                        vaspout.close()
                    continue

        if not foundout:
            print(f'No accessible OUT file ({alternative_filenames}) in {f}')
            return None

        # Read INCAR to get IBRION
        ibrion = None
        incar_path = f + 'INCAR'
        if os.path.exists(incar_path):
            try:
                with open(incar_path, 'r') as incar:
                    # Test read first to check accessibility
                    incarlines = incar.readlines()
                    for line in incarlines:
                        if 'IBRION' in line:
                            ibrion = int(line.split()[2])
                            break
            except (OSError, TimeoutError, IOError) as e:
                print(f'Warning: Cannot read INCAR in {f} (file may not be downloaded)')
                if vaspout:
                    vaspout.close()
                return None

        # Read output file (last 5 lines to check convergence)
        try:
            lines = vaspout.readlines()[-5:]
            vaspout.close()
        except (OSError, TimeoutError, IOError) as e:
            print(f'Warning: Cannot read output file in {f} (file may not be downloaded)')
            if vaspout:
                vaspout.close()
            return None

        # Check convergence based on IBRION
        if ibrion in [1, 2, 3]:
            for line in lines:
                if 'reached required accuracy' in line:
                    return True
            return False
        else:
            for line in lines:
                if 'E0=' in line:
                    return True
            return False

    except TimeoutError as e:
        print(f'Timeout reading files in {f} (OneDrive files may not be downloaded)')
        return None
    except Exception as e:
        print(f'Error reading files in {f}: {e}')
        return None


if __name__ == "__main__":
    main()
