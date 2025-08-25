#!/usr/bin/env python

from ase.io import read, write
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from pathlib import Path
from ase.neb import NEB
from ase.mep import idpp_interpolate
from ase.data import atomic_numbers, covalent_radii

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

def redistribute_images_evenly(images, mic=True):
    """
    Redistribute images along the path to be more evenly spaced.
    Uses cumulative distance along the path to achieve equal spacing.
    
    Parameters:
    -----------
    images : list of ase.Atoms
        List of ASE Atoms objects representing NEB images
    mic : bool
        Use minimum image convention for periodic boundary conditions
        
    Returns:
    --------
    list of ase.Atoms
        Redistributed images with even spacing along the path
    """
    n_images = len(images)
    if n_images < 3:
        return images
    
    # Calculate cumulative distances along the path
    distances = [0.0]
    for i in range(1, n_images):
        # Calculate distance between consecutive images
        pos1 = images[i-1].get_positions()
        pos2 = images[i].get_positions()
        
        if mic:
            # Use minimum image convention
            cell = images[i].get_cell()
            diff = pos2 - pos1
            diff = diff - np.round(diff @ np.linalg.inv(cell)) @ cell
        else:
            diff = pos2 - pos1
        
        dist = np.sqrt(np.sum(diff**2))
        distances.append(distances[-1] + dist)
    
    # Total path length
    total_length = distances[-1]
    
    # Create evenly spaced target distances
    target_distances = np.linspace(0, total_length, n_images)
    
    # Redistribute images
    redistributed_images = []
    redistributed_images.append(images[0].copy())  # Keep initial image
    
    for i in range(1, n_images - 1):
        target_dist = target_distances[i]
        
        # Find which segment this target distance falls in
        for j in range(len(distances) - 1):
            if distances[j] <= target_dist <= distances[j + 1]:
                # Interpolate between images j and j+1
                if distances[j + 1] - distances[j] > 1e-10:
                    alpha = (target_dist - distances[j]) / (distances[j + 1] - distances[j])
                else:
                    alpha = 0.0
                
                pos1 = images[j].get_positions()
                pos2 = images[j + 1].get_positions()
                
                if mic:
                    cell = images[j].get_cell()
                    diff = pos2 - pos1
                    diff = diff - np.round(diff @ np.linalg.inv(cell)) @ cell
                    new_pos = pos1 + alpha * diff
                else:
                    new_pos = pos1 + alpha * (pos2 - pos1)
                
                new_image = images[j].copy()
                new_image.set_positions(new_pos)
                redistributed_images.append(new_image)
                break
    
    redistributed_images.append(images[-1].copy())  # Keep final image
    return redistributed_images

def interpolate_neb_images(initial_atoms, final_atoms, n_images, method='idpp', mic=True):
    """
    Generate intermediate images between initial and final structures.
    
    Parameters:
    -----------
    initial_atoms : ase.Atoms
        Initial structure
    final_atoms : ase.Atoms
        Final structure  
    n_images : int
        Total number of images including endpoints
    method : str
        Interpolation method ('idpp' or 'linear')
    mic : bool
        Use minimum image convention for periodic systems
        
    Returns:
    --------
    list of ase.Atoms
        List of interpolated images
    """
    if n_images < 2:
        raise ValueError("Need at least 2 images (initial and final)")
    
    # Validate that structures are compatible
    if len(initial_atoms) != len(final_atoms):
        raise ValueError("Initial and final structures must have same number of atoms")
    
    if initial_atoms.get_chemical_symbols() != final_atoms.get_chemical_symbols():
        raise ValueError("Initial and final structures must have same atomic species")
    
    # Create image list
    images = [initial_atoms.copy() for _ in range(n_images)]
    images[-1] = final_atoms.copy()
    
    if method == 'idpp':
        try:
            # Use IDPP interpolation with robust parameters
            idpp_interpolate(images, mic=mic, fmax=0.05, steps=200)
            print(f"IDPP interpolation completed for {n_images} images")
        except Exception as e:
            print(f"IDPP interpolation failed: {e}")
            print("Falling back to linear interpolation")
            method = 'linear'
    
    if method == 'linear':
        # Linear interpolation as fallback
        pos_initial = initial_atoms.get_positions()
        pos_final = final_atoms.get_positions()
        
        if mic and initial_atoms.get_pbc().any():
            # Handle periodic boundary conditions
            cell = initial_atoms.get_cell()
            diff = pos_final - pos_initial
            diff = diff - np.round(diff @ np.linalg.inv(cell)) @ cell
            pos_final = pos_initial + diff
        
        for i in range(1, n_images - 1):
            alpha = i / (n_images - 1)
            new_pos = pos_initial + alpha * (pos_final - pos_initial)
            images[i].set_positions(new_pos)
        
        print(f"Linear interpolation completed for {n_images} images")
    
    return images

def check_atomic_distances(atoms, shrink_factor=0.8, mic=True):
    """
    Check for atoms that are too close to each other based on covalent radii.
    
    Parameters:
    -----------
    atoms : ase.Atoms
        ASE Atoms object to check
    shrink_factor : float
        Factor to multiply sum of covalent radii (default: 0.8)
    mic : bool
        Use minimum image convention for periodic systems
        
    Returns:
    --------
    list of tuples
        List of (i, j, distance, threshold) for atom pairs that are too close
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    n_atoms = len(atoms)
    
    close_pairs = []
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # Get atomic numbers and covalent radii
            z_i = atomic_numbers[symbols[i]]
            z_j = atomic_numbers[symbols[j]]
            
            # Sum of covalent radii as threshold
            threshold = shrink_factor * (covalent_radii[z_i] + covalent_radii[z_j])
            
            # Calculate distance
            pos_i = positions[i]
            pos_j = positions[j]
            
            if mic and atoms.get_pbc().any():
                # Use minimum image convention
                cell = atoms.get_cell()
                diff = pos_j - pos_i
                diff = diff - np.round(diff @ np.linalg.inv(cell)) @ cell
                distance = np.linalg.norm(diff)
            else:
                distance = np.linalg.norm(pos_j - pos_i)
            
            if distance < threshold:
                close_pairs.append((i, j, distance, threshold))
    
    return close_pairs

def check_neb_images_sanity(output_dir, shrink_factor=0.8, mic=True):
    """
    Check all NEB images for atoms that are too close to each other.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing NEB image folders (00, 01, 02, etc.)
    shrink_factor : float
        Factor to multiply sum of covalent radii (default: 0.8)
    mic : bool
        Use minimum image convention for periodic systems
        
    Returns:
    --------
    dict
        Dictionary with image names as keys and lists of close pairs as values
    """
    output_path = Path(output_dir)
    problems = {}
    
    # Find all image directories (00, 01, 02, etc.)
    image_dirs = sorted([d for d in output_path.iterdir() 
                        if d.is_dir() and d.name.isdigit() and len(d.name) == 2])
    
    for image_dir in image_dirs:
        poscar_file = image_dir / "POSCAR"
        if poscar_file.exists():
            try:
                atoms = read(str(poscar_file))
                close_pairs = check_atomic_distances(atoms, shrink_factor, mic)
                
                if close_pairs:
                    problems[image_dir.name] = close_pairs
                    
            except Exception as e:
                print(f"Warning: Could not check {poscar_file}: {e}")
    
    return problems

def setup_neb_calculation(initial_file='IS/CONTCAR', final_file='FS/CONTCAR', output_dir='.', n_images=5, use_even_spacing=True, mic=True, check_distances=True, shrink_factor=0.8):
    """
    Set up a complete NEB calculation folder structure.
    
    Parameters:
    -----------
    initial_file : str
        Path to initial structure file (default: 'IS/CONTCAR')
    final_file : str  
        Path to final structure file (default: 'FS/CONTCAR')
    output_dir : str
        Directory where NEB folders will be created (default: '.')
    n_images : int
        Number of intermediate images (default: 5, creates 7 total folders 00-06)
    use_even_spacing : bool
        Apply even spacing redistribution after interpolation (default: True)
    mic : bool
        Use minimum image convention for periodic systems (default: True)
    check_distances : bool
        Perform sanity check for atoms that are too close (default: True)
    shrink_factor : float
        Factor to multiply sum of covalent radii for distance check (default: 0.8)
        
    Returns:
    --------
    bool
        True if setup successful, False otherwise
    """
    
    output_path = Path(output_dir)
    if not output_path.exists():
        try:
            output_path.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        except Exception as e:
            print(f"Error creating output directory {output_dir}: {e}")
            return False
    
    # Check if initial and final files exist
    initial_path = Path(initial_file)
    final_path = Path(final_file)
    
    if not initial_path.exists():
        print(f"Error: Initial structure file not found at {initial_file}")
        return False
    
    if not final_path.exists():
        print(f"Error: Final structure file not found at {final_file}")
        return False
    
    # Read initial and final structures
    try:
        initial_atoms = read(str(initial_path))
        final_atoms = read(str(final_path))
        print(f"Successfully read structures from {initial_file} and {final_file}")
    except Exception as e:
        print(f"Error reading structure files: {e}")
        return False
    
    # Total number of images including endpoints
    total_images = n_images + 2
    
    # Generate interpolated images
    try:
        images = interpolate_neb_images(initial_atoms, final_atoms, total_images, method='idpp', mic=mic)
    except Exception as e:
        print(f"Error during interpolation: {e}")
        return False
    
    # Apply even spacing if requested
    if use_even_spacing:
        try:
            images = redistribute_images_evenly(images, mic=mic)
            print(f"Images redistributed for even spacing along reaction path")
        except Exception as e:
            print(f"Error during redistribution: {e}")
            print("Continuing with IDPP spacing...")
    
    # Create directories and POSCAR files
    for i, atoms in enumerate(images):
        image_dir = output_path / f"{i:02d}"
        try:
            image_dir.mkdir(exist_ok=True)
            poscar_path = image_dir / "POSCAR"
            write(str(poscar_path), atoms, format='vasp')
            print(f"Created {poscar_path}")
        except Exception as e:
            print(f"Error creating image {i:02d}: {e}")
            return False
    
    # Perform sanity check for atomic distances if requested
    if check_distances:
        print("\nPerforming sanity check for atomic distances...")
        problems = check_neb_images_sanity(output_dir, shrink_factor, mic)
        
        if problems:
            print("WARNING: Found atoms that may be too close to each other:")
            for image, close_pairs in problems.items():
                print(f"\nImage {image}:")
                for i, j, distance, threshold in close_pairs:
                    print(f"  Atoms {i}-{j}: distance={distance:.3f} Å < threshold={threshold:.3f} Å")
            
            print(f"\nNote: Threshold = {shrink_factor} × (r_i + r_j) where r_i, r_j are covalent radii")
            print("Consider adjusting structures or using different interpolation settings")
        else:
            print(f"✓ All atomic distances OK (threshold = {shrink_factor} × covalent radii)")
    
    spacing_method = "even spacing" if use_even_spacing else "IDPP spacing"
    print(f"\nSuccessfully set up NEB calculation with {total_images} images using {spacing_method}")
    print(f"Created folders: {', '.join([f'{i:02d}' for i in range(total_images)])}")
    print(f"Each folder contains a POSCAR file for the corresponding image")
    
    return True