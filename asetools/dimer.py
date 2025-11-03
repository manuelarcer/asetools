#!/usr/bin/env python

"""
Dimer method implementation for ASE with VaspInteractive support.

This module provides utilities for dimer calculations including:
- MODECAR file I/O (compatible with Henkelman's VTST tools)
- Integration with ASE's dimer method (MinModeAtoms, DimerControl)
- Support for VaspInteractive calculator
- ASE constraints compatibility
"""

import os
import numpy as np
from ase.io import read, write
from ase.mep import DimerControl, MinModeAtoms
import logging

logger = logging.getLogger(__name__)

def read_modecar(filename='MODECAR'):
    """
    Read displacement vector from MODECAR file.

    MODECAR file format is similar to POSCAR but contains displacement vectors
    instead of positions. Each line represents the displacement vector for
    the corresponding atom in x, y, z directions.

    Parameters:
    -----------
    filename : str
        Path to MODECAR file (default: 'MODECAR')

    Returns:
    --------
    numpy.ndarray
        Displacement vectors as (N, 3) array where N is number of atoms
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"MODECAR file not found: {filename}")

    try:
        # Read MODECAR using ASE's POSCAR reader, but interpret as displacement vectors
        modecar_atoms = read(filename, format='vasp')
        displacement_vector = modecar_atoms.get_positions()

        # Normalize the displacement vector to unit magnitude
        magnitude = np.linalg.norm(displacement_vector)
        if magnitude > 1e-10:
            displacement_vector = displacement_vector / magnitude
        else:
            logger.warning("MODECAR displacement vector has zero magnitude, using random")
            displacement_vector = np.random.random(displacement_vector.shape) - 0.5
            displacement_vector = displacement_vector / np.linalg.norm(displacement_vector)

        logger.info(f"Successfully read MODECAR from {filename}")
        logger.info(f"Displacement vector magnitude: {np.linalg.norm(displacement_vector):.6f}")

        return displacement_vector

    except Exception as e:
        logger.error(f"Failed to read MODECAR file {filename}: {e}")
        raise


def write_modecar(displacement_vector, atoms, filename='MODECAR'):
    """
    Write displacement vector to MODECAR file.

    Parameters:
    -----------
    displacement_vector : numpy.ndarray
        Displacement vectors as (N, 3) array
    atoms : ase.Atoms
        Atoms object for structure information (cell, symbols, etc.)
    filename : str
        Output filename (default: 'MODECAR')
    """
    try:
        # Create a copy of atoms with displacement vectors as positions
        modecar_atoms = atoms.copy()
        modecar_atoms.set_positions(displacement_vector)

        # Write in VASP format
        write(filename, modecar_atoms, format='vasp')
        logger.info(f"Successfully wrote MODECAR to {filename}")

    except Exception as e:
        logger.error(f"Failed to write MODECAR file {filename}: {e}")
        raise


def generate_displacement_vector(atoms, method='random', direction=None, magnitude=0.01):
    """
    Generate initial displacement vector for dimer calculation.

    Parameters:
    -----------
    atoms : ase.Atoms
        Atoms object
    method : str
        Method for generating displacement ('random', 'direction', 'modecar')
    direction : list or numpy.ndarray, optional
        Specific direction vector for 'direction' method
    magnitude : float
        Magnitude of displacement (default: 0.01)

    Returns:
    --------
    numpy.ndarray
        Displacement vector as (N, 3) array
    """
    n_atoms = len(atoms)

    if method == 'modecar':
        # Try to read from MODECAR file
        try:
            return read_modecar('MODECAR')
        except FileNotFoundError:
            logger.warning("MODECAR file not found, falling back to random displacement")
            method = 'random'

    if method == 'random':
        # Generate random displacement vector
        displacement_vector = np.random.random((n_atoms, 3)) - 0.5
        displacement_vector = displacement_vector / np.linalg.norm(displacement_vector) * magnitude
        logger.info(f"Generated random displacement vector with magnitude {magnitude}")

    elif method == 'direction':
        if direction is None:
            raise ValueError("Direction vector must be provided for 'direction' method")

        direction = np.array(direction)
        if direction.shape != (n_atoms, 3):
            raise ValueError(f"Direction shape {direction.shape} doesn't match atoms shape ({n_atoms}, 3)")

        displacement_vector = direction / np.linalg.norm(direction) * magnitude
        logger.info(f"Generated directional displacement vector with magnitude {magnitude}")

    else:
        raise ValueError(f"Unknown displacement method: {method}")

    return displacement_vector


def setup_dimer_atoms(atoms, displacement_vector=None, dimer_control_kwargs=None):
    """
    Set up MinModeAtoms object for dimer calculation.

    Parameters:
    -----------
    atoms : ase.Atoms
        Atoms object with calculator attached
    displacement_vector : numpy.ndarray, optional
        Initial displacement vector. If None, will try MODECAR or random
    dimer_control_kwargs : dict, optional
        Parameters for DimerControl initialization

    Returns:
    --------
    ase.mep.MinModeAtoms
        MinModeAtoms object ready for dimer optimization
    """
    if dimer_control_kwargs is None:
        dimer_control_kwargs = {}

    # Set default dimer control parameters
    default_control_kwargs = {
        'initial_eigenmode_method': 'displacement',
        'displacement_method': 'vector',
        'logfile': 'dimer.log',
        'mask': None,  # None means all atoms participate
    }
    default_control_kwargs.update(dimer_control_kwargs)

    # Create DimerControl object
    d_control = DimerControl(**default_control_kwargs)

    # Create MinModeAtoms object
    d_atoms = MinModeAtoms(atoms, d_control)

    # Set up initial displacement vector
    if displacement_vector is None:
        # Try to read from MODECAR, fall back to random
        try:
            displacement_vector = read_modecar('MODECAR')
            logger.info("Using displacement vector from MODECAR file")
        except FileNotFoundError:
            logger.info("No MODECAR found, generating random displacement vector")
            displacement_vector = generate_displacement_vector(atoms, method='random')

    # Apply displacement vector
    d_atoms.displace(displacement_vector=displacement_vector)

    logger.info("Successfully set up MinModeAtoms for dimer calculation")
    logger.info(f"Displacement vector magnitude: {np.linalg.norm(displacement_vector):.6f}")

    return d_atoms


def check_dimer_convergence(d_atoms, force_threshold=0.01, eigenvalue_threshold=-0.01):
    """
    Check convergence criteria for dimer calculation.

    Parameters:
    -----------
    d_atoms : ase.mep.MinModeAtoms
        MinModeAtoms object from dimer calculation
    force_threshold : float
        Maximum force threshold for convergence (eV/Å)
    eigenvalue_threshold : float
        Maximum eigenvalue threshold (should be negative for saddle point)

    Returns:
    --------
    dict
        Convergence information including status and current values
    """
    try:
        forces = d_atoms.get_forces()
        max_force = np.max(np.linalg.norm(forces, axis=1))

        # Get eigenvalue if available
        eigenvalue = getattr(d_atoms, 'eigenvalue', None)

        force_converged = max_force < force_threshold
        eigenvalue_converged = eigenvalue is not None and eigenvalue < eigenvalue_threshold

        converged = force_converged and (eigenvalue is None or eigenvalue_converged)

        convergence_info = {
            'converged': converged,
            'max_force': max_force,
            'force_threshold': force_threshold,
            'force_converged': force_converged,
            'eigenvalue': eigenvalue,
            'eigenvalue_threshold': eigenvalue_threshold,
            'eigenvalue_converged': eigenvalue_converged,
        }

        logger.info(f"Dimer convergence check: max_force={max_force:.4f}, eigenvalue={eigenvalue}")

        return convergence_info

    except Exception as e:
        logger.error(f"Error checking dimer convergence: {e}")
        return {'converged': False, 'error': str(e)}


def extract_saddle_point_info(d_atoms):
    """
    Extract information about the saddle point from converged dimer calculation.

    Parameters:
    -----------
    d_atoms : ase.mep.MinModeAtoms
        Converged MinModeAtoms object

    Returns:
    --------
    dict
        Saddle point information including energy, eigenvalue, eigenvector
    """
    try:
        atoms = d_atoms.atoms
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        max_force = np.max(np.linalg.norm(forces, axis=1))

        # Get dimer-specific information
        eigenvalue = getattr(d_atoms, 'eigenvalue', None)
        eigenvector = getattr(d_atoms, 'eigenvector', None)

        saddle_info = {
            'energy': energy,
            'max_force': max_force,
            'eigenvalue': eigenvalue,
            'eigenvector': eigenvector,
            'positions': atoms.get_positions().copy(),
            'cell': atoms.get_cell().copy(),
            'symbols': atoms.get_chemical_symbols(),
        }

        logger.info(f"Extracted saddle point: E={energy:.6f} eV, eigenvalue={eigenvalue}")

        return saddle_info

    except Exception as e:
        logger.error(f"Error extracting saddle point info: {e}")
        return {'error': str(e)}


def save_dimer_trajectory(d_atoms, filename='dimer_trajectory.traj'):
    """
    Save dimer calculation trajectory.

    Parameters:
    -----------
    d_atoms : ase.mep.MinModeAtoms
        MinModeAtoms object
    filename : str
        Output trajectory filename
    """
    try:
        from ase.io import Trajectory

        # Check if trajectory is already being saved
        if hasattr(d_atoms, 'trajectory') and d_atoms.trajectory is not None:
            logger.info("Trajectory already being saved by optimizer")
            return

        # Create new trajectory
        traj = Trajectory(filename, 'w')
        traj.write(d_atoms.atoms)
        traj.close()

        logger.info(f"Saved dimer trajectory to {filename}")

    except Exception as e:
        logger.error(f"Error saving dimer trajectory: {e}")


# Utility functions for integration with workflow manager

def validate_dimer_kwargs(optimizer_kwargs):
    """
    Validate and process dimer-specific optimizer kwargs.

    Parameters:
    -----------
    optimizer_kwargs : dict
        Raw optimizer kwargs from YAML configuration

    Returns:
    --------
    tuple
        (dimer_control_kwargs, optimizer_init_kwargs, optimizer_run_kwargs)
    """
    # Separate dimer-specific parameters
    dimer_control_params = {
        'initial_eigenmode_method', 'displacement_method', 'mask',
        'logfile', 'displacement_vector'
    }

    # Parameters that go to optimizer initialization
    optimizer_init_params = {'logfile', 'trajectory'}

    # Parameters that go to optimizer run method
    optimizer_run_params = {'fmax', 'steps'}

    dimer_control_kwargs = {}
    optimizer_init_kwargs = {}
    optimizer_run_kwargs = {}

    for key, value in optimizer_kwargs.items():
        if key in dimer_control_params:
            if key != 'displacement_vector':  # Handle separately
                dimer_control_kwargs[key] = value
        elif key in optimizer_init_params:
            optimizer_init_kwargs[key] = value
        elif key in optimizer_run_params:
            optimizer_run_kwargs[key] = value
        else:
            # Default to run kwargs for unknown parameters
            optimizer_run_kwargs[key] = value

    # Set defaults
    if 'fmax' not in optimizer_run_kwargs:
        optimizer_run_kwargs['fmax'] = 0.01

    if 'logfile' not in optimizer_init_kwargs:
        optimizer_init_kwargs['logfile'] = 'dimer_optimization.log'

    return dimer_control_kwargs, optimizer_init_kwargs, optimizer_run_kwargs