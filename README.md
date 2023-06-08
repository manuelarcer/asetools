# ASETOOLS 

## Description:
This package ...

## Installation

As simple as 

`pip install git+https://github.com/manuelarcer/asetools.git`

or

`pip install git+git@github.com:manuelarcer/asetools.git`

### Requirements

- pandas
- numpy
- ase
- matplotlib

## Package's modules

### Module `analysis`

The `analysis.py` module contains functions designed to analyze OUTCAR files produced by the VASP (Vienna Ab initio Simulation Package) software. These functions can be used to check the convergence of a calculation, as well as to extract information about the energy and forces calculated in the simulation. The module makes use of the Atomic Simulation Environment (ASE) to read information from OUTCAR files.

Here's a brief summary of the two functions in the module:

`check_outcar_convergence(outcar, verbose=False)`: This function checks if a calculation represented by the given OUTCAR file has converged. It returns a tuple containing a boolean value that indicates whether the calculation has converged and a string indicating the VASP version used for the calculation.

`check_energy_and_maxforce(outcar, magmom=False, verbose=False)`: This function returns the final potential energy and maximum force from the calculation represented by the OUTCAR file. If the `magmom` parameter is set to True, it also returns the magnetic moment of the system.

Here's an example of how you can use these functions:

```python
from my_package.analysis import check_outcar_convergence, check_energy_and_maxforce

outcar_file = "/path/to/your/OUTCAR"

# Check if the calculation has converged
convergence, vasp_version = check_outcar_convergence(outcar_file)
print(f"Convergence: {convergence}")
print(f"VASP version: {vasp_version}")

# Get the final potential energy and maximum force
energy, max_force = check_energy_and_maxforce(outcar_file)
print(f"Final potential energy: {energy} eV")
print(f"Maximum force: {max_force} eV/A")

# If you also want to get the magnetic moment
energy, max_force, mag_moment = check_energy_and_maxforce(outcar_file, magmom=True)
print(f"Magnetic moment: {mag_moment} Bohr magnetons")
```


This Python code reads an OUTCAR file, checks if the associated VASP calculation has converged, and extracts the final potential energy and maximum force. If you want to obtain the magnetic moment as well, you can set the `magmom` parameter to True when calling the `check_energy_and_maxforce` function. Replace "/path/to/your/OUTCAR" with the actual path to your OUTCAR file.

### Module `doscar_analysis`

The `doscar_analysis` module provides functionality to analyze the Density of States (DOS) from VASP calculations. It allows users to extract the DOS from DOSCAR files and also to extract projected DOS (pDOS) per state or per orbital. Here is a description of the functions contained in the module:

`extract_dos(doscarfile)`: This function extracts the DOS from a given DOSCAR file. It returns a dictionary containing the energy, DOSup and DOSdown.

`extract_pdos_perstate(data, atoms, states)`: This function extracts the pDOS per state for specific atoms. The user can specify a list of atoms and states (s, p, d) of interest.

`extract_pdos_perorbital(data, atoms, orbitals)`: This function extracts the pDOS per orbital for specific atoms. The user can specify a list of atoms and orbitals of interest.

Here is an example usage of the module for the README file:

```python
from doscar_analysis import extract_dos, extract_pdos_perstate, extract_pdos_perorbital

# Specify the path to your DOSCAR file
doscarfile = "path/to/your/DOSCAR"

# Extract the DOS
dos_data = extract_dos(doscarfile)

# Print the extracted data
print(dos_data)

# Specify atoms and states of interest
atoms = [0, 1, 2]  # Indices of atoms of interest
states = ['s_states', 'p_states', 'd_states']  # States of interest

# Extract the PDOS per state
e, sum_plus, sum_minus = extract_pdos_perstate(dos_data, atoms, states)

# Print the extracted data
print(f"Energy: {e}")
print(f"Sum Plus: {sum_plus}")
print(f"Sum Minus: {sum_minus}")

# Specify atoms and orbitals of interest
atoms = [0, 1, 2]  # Indices of atoms of interest
orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+']  # Orbitals of interest

# Extract the PDOS per orbital
e, sum_plus, sum_minus = extract_pdos_perorbital(dos_data, atoms, orbitals)

# Print the extracted data
print(f"Energy: {e}")
print(f"Sum Plus: {sum_plus}")
print(f"Sum Minus: {sum_minus}")
```

### Module `adsorbate`

The SurfaceAnalyzer class provided in your script is a powerful tool designed to analyze the surface of a molecular structure. It takes an ASE Atoms object as input and provides several methods to analyze and manipulate the surface of the structure.

Here's a brief summary of what each method does:

`__init__(self, atoms)`: This is the constructor of the SurfaceAnalyzer class. It initializes the instance with a given set of atoms and calculates the surface atoms and their neighbors.

`find_surface_atoms(self)`: This method finds and returns the indices of the atoms that are on the surface of the structure.

`find_surface_neighbors(self, thr=3.0)`: This method identifies neighboring atoms on the surface. The neighbor threshold can be adjusted with the thr parameter.

`midpoint_three_atoms(self, three)`: This method calculates the midpoint of a triangle formed by three atoms, given their indices.

`add_adsorbate_to_mid(self, three, adsorbate='H', z_off=1., thr=4.4)`: This method adds an adsorbate atom at the midpoint of a triangle formed by three surface atoms. The type of atom added, its z offset, and the threshold for finding neighbors around the adsorbate can be adjusted.

`get_cluster_around_adsorbate(self, thr=4.4)`: This method identifies a cluster of atoms around the adsorbate. The threshold for inclusion in the cluster can be adjusted.

`adsorbate_to_cluster_neighbors(self, thr=4.4, prec=3)`: This method calculates the distances from the adsorbate atom to its neighbors in the cluster.

```python
from ase.build import molecule
from adsorbate import SurfaceAnalyzer

# Load a molecule structure using ASE
atoms = molecule("CH3CH2OH")

# Initialize a SurfaceAnalyzer object
analyzer = SurfaceAnalyzer(atoms)

# Find the indices of surface atoms
print(analyzer.surface_indices)

# Add an adsorbate to the surface
analyzer.add_adsorbate_to_mid(analyzer.surface_indices[:3])

# Print distances from the adsorbate to its neighbors
print(analyzer.adsneighdistances)
```

In this example, we load a molecule structure using the ASE build module. We then initialize a `SurfaceAnalyzer` object with the atoms of this molecule. We print the indices of the surface atoms. Finally, we add an adsorbate to the surface at the midpoint of the first three surface atoms and print the distances from the adsorbate to its neighbors.

Please replace "my_module" with the actual name of your module. This example assumes that the SurfaceAnalyzer class is imported from a module called "my_module".

## Package's executables

### `getenergy.py`

The `getenergy.py` script is a tool that enables users to analyze the convergence, energy, and maximum force of an OUTCAR file (a common output file in DFT calculations, especially with VASP software). This is particularly useful for users conducting materials simulations and those who require a simplified method of checking their OUTCAR files' status.

The script uses `asetools.analysis` library functions, specifically `check_energy_and_maxforce` and `check_outcar_convergence`. It takes an OUTCAR file as an argument and determines whether the calculation has converged by analyzing the OUTCAR file. If the convergence is successful, the energy and maximum force of the system are returned.

If no argument is provided, the script defaults to analyzing a file named "OUTCAR" in the current directory.

### `summaryfolders.py`

This script is used for analysis and summarization of results from multiple OUTCAR files (output files from DFT calculations, commonly generated with the VASP software), specifically across different configuration folders. It determines convergence status, calculates the energy, maximum force, and optionally, the magnetic moments, then outputs the data as a dataframe.

The script utilizes `glob`, `pandas`, `numpy`, and `asetools.analysis` libraries to process and analyze multiple OUTCAR files. It performs a directory-wide search for directories containing OUTCAR files and then carries out the following analyses:

1. Check if the calculations in the OUTCAR files have converged.
2. Extracts the energy, maximum force, and magnetic moments (optional) from the OUTCAR files.
3. Calculates the relative energy with respect to the minimum energy across all the configurations.
4. Prints out a `DataFrame` that includes all the extracted and calculated information.
5. If the `-m` or `--magmom` flag is included, the script will also calculate the magnetic moments.

The script assumes that each directory contains one OUTCAR file.

## Extra content

Extra content can be found in the folder 'extra', within package's folder.

### 1. vasp_outcar_parser.py

This is a file to fix one of the bugs while reading an OUTCAR file generated by VASP 6. This has been apparent elsewhere (https://gitlab.com/ase/ase/-/issues/1119). I had the same issue. Copy this file into ase/io/vasp_parser/ to overwrite the existing one.

The fix was submitted here 
https://gitlab.com/ase/ase/-/commit/9e9234a785c97ca3d49a0f3fca4fec2bdc52eb2b

And the full file can be found in this other page
https://gitlab.com/ase/ase/-/blob/9e9234a785c97ca3d49a0f3fca4fec2bdc52eb2b/ase/io/vasp_parsers/vasp_outcar_parsers.py


  
