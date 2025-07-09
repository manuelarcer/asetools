"""
Bond Valence Sum (BVS) calculations using Brown's equation.

This module provides functionality to calculate bond valence sums for atomic structures
using the bond valence parameters from Brown's accumulated table (bvparm2020.cif).

The bond valence equation is: bond_valence = exp((R0-R)/B)
where R is the bond length, R0 and B are bond valence parameters.

Classes:
    BondValenceParameters: Handles parsing and lookup of bond valence parameters
    BondValenceSum: Calculates bond valence sums for atomic structures

Author: ASEtools development team
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
from ase import Atoms
from ase.neighborlist import NeighborList
import os
import re
from collections import defaultdict


class BondValenceParameters:
    """
    Class to handle bond valence parameters from Brown's accumulated table.
    
    Loads and provides lookup functionality for bond valence parameters
    used in Brown's bond valence equation: bond_valence = exp((R0-R)/B)
    """
    
    def __init__(self, parameter_file: Optional[str] = None):
        """
        Initialize BondValenceParameters.
        
        Args:
            parameter_file: Path to bond valence parameter file (default: bvparm2020.cif)
        """
        if parameter_file is None:
            # Default to the included parameter file
            current_dir = os.path.dirname(os.path.abspath(__file__))
            parameter_file = os.path.join(current_dir, 'data', 'bvparm2020.cif')
        
        self.parameter_file = parameter_file
        self.parameters = {}
        self.references = {}
        self._load_parameters()
    
    def _load_parameters(self):
        """Load bond valence parameters from the CIF file."""
        try:
            with open(self.parameter_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"Bond valence parameter file not found: {self.parameter_file}")
        
        # Parse references first
        self._parse_references(lines)
        
        # Parse parameter data
        self._parse_parameters(lines)
    
    def _parse_references(self, lines: List[str]):
        """Parse the reference section of the CIF file."""
        in_ref_section = False
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Check if we're entering the reference section
            if line.startswith('loop_') and i + 1 < len(lines) and '_valence_ref_id' in lines[i + 1]:
                in_ref_section = True
                continue
            
            if in_ref_section:
                # Check if we're entering the parameter section
                if line.startswith('loop_') and i + 1 < len(lines) and '_valence_param_' in lines[i + 1]:
                    break
                
                if line and not line.startswith('_'):
                    # Parse reference line: "a 'Brown and Altermatt, (1985)...'"
                    parts = line.split("'", 1)
                    if len(parts) >= 2:
                        ref_id = parts[0].strip()
                        ref_text = parts[1].rstrip("'")
                        self.references[ref_id] = ref_text
    
    def _parse_parameters(self, lines: List[str]):
        """Parse the parameter data section of the CIF file."""
        in_param_section = False
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Check if we're entering the parameter section
            if line.startswith('loop_') and i + 1 < len(lines) and '_valence_param_' in lines[i + 1]:
                in_param_section = True
                continue
            
            if in_param_section and line and not line.startswith('_'):
                # Parse parameter line format:
                # Element1 Valence1 Element2 Valence2 R0 B Reference Details
                parts = line.split()
                if len(parts) >= 7:
                    element1 = parts[0]
                    try:
                        valence1 = int(parts[1])
                    except ValueError:
                        continue
                    
                    element2 = parts[2]
                    try:
                        valence2 = int(parts[3])
                    except ValueError:
                        continue
                    
                    try:
                        R0 = float(parts[4])
                        B = float(parts[5])
                    except ValueError:
                        continue
                    
                    reference = parts[6]
                    details = ' '.join(parts[7:]) if len(parts) > 7 else ''
                    
                    # Create key for parameter lookup
                    key = (element1, valence1, element2, valence2)
                    
                    # Store parameters in order of reliability (first entry is most reliable)
                    if key not in self.parameters:
                        self.parameters[key] = []
                    
                    self.parameters[key].append({
                        'R0': R0,
                        'B': B,
                        'reference': reference,
                        'details': details,
                        'reference_text': self.references.get(reference, 'Unknown reference')
                    })
    
    def get_parameters(self, element1: str, valence1: int, element2: str, valence2: int,
                      most_reliable: bool = True, exclude_unchecked: bool = True) -> Union[Dict, List[Dict]]:
        """
        Get bond valence parameters for a given element pair.
        
        Args:
            element1: First element symbol
            valence1: Valence state of first element
            element2: Second element symbol
            valence2: Valence state of second element
            most_reliable: If True, return only the most reliable parameter set
            exclude_unchecked: If True, exclude parameters marked as 'unchecked'
            
        Returns:
            Dictionary or list of dictionaries with R0, B, reference, details
            
        Raises:
            ValueError: If no parameters found for the given element pair
        """
        # Try both orientations of the element pair
        key1 = (element1, valence1, element2, valence2)
        key2 = (element2, valence2, element1, valence1)
        
        params = None
        if key1 in self.parameters:
            params = self.parameters[key1]
        elif key2 in self.parameters:
            params = self.parameters[key2]
        
        if params is None:
            raise ValueError(f"No bond valence parameters found for {element1}({valence1})-{element2}({valence2})")
        
        # Filter out unchecked parameters if requested
        if exclude_unchecked:
            params = [p for p in params if 'unchecked' not in p['details'].lower()]
        
        if not params:
            raise ValueError(f"No reliable bond valence parameters found for {element1}({valence1})-{element2}({valence2})")
        
        if most_reliable:
            return params[0]  # Most reliable is first in the list
        else:
            return params
    
    def list_available_pairs(self) -> List[Tuple[str, int, str, int]]:
        """
        List all available element pairs with their valence states.
        
        Returns:
            List of tuples (element1, valence1, element2, valence2)
        """
        return list(self.parameters.keys())
    
    def get_elements(self) -> List[str]:
        """
        Get list of all elements with bond valence parameters.
        
        Returns:
            List of element symbols
        """
        elements = set()
        for element1, _, element2, _ in self.parameters.keys():
            elements.add(element1)
            elements.add(element2)
        return sorted(list(elements))
    
    def get_valence_states(self, element: str) -> List[int]:
        """
        Get available valence states for a given element.
        
        Args:
            element: Element symbol
            
        Returns:
            List of valence states
        """
        valences = set()
        for element1, valence1, element2, valence2 in self.parameters.keys():
            if element1 == element:
                valences.add(valence1)
            if element2 == element:
                valences.add(valence2)
        return sorted(list(valences))


class BondValenceSum:
    """
    Class to calculate bond valence sums for atomic structures.
    
    Uses Brown's bond valence equation: bond_valence = exp((R0-R)/B)
    where R is the bond length, R0 and B are bond valence parameters.
    """
    
    def __init__(self, atoms: Atoms, valence_states: Optional[Dict[str, int]] = None,
                 custom_parameters: Optional[Dict[str, Dict[str, float]]] = None,
                 distance_cutoff: float = 3.5):
        """
        Initialize BondValenceSum calculator.
        
        Args:
            atoms: ASE Atoms object
            valence_states: Dict mapping element symbols to valence states
            custom_parameters: Dict with custom R0 and B values for specific elements
            distance_cutoff: Maximum distance for neighbor detection (Angstrom)
        """
        self.atoms = atoms
        self.distance_cutoff = distance_cutoff
        
        # Load bond valence parameters
        self.bv_params = BondValenceParameters()
        
        # Set default valence states (most common oxidation states)
        self.default_valences = {
            'H': 1, 'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': -3, 'O': -2, 'F': -1,
            'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': -2, 'Cl': -1,
            'K': 1, 'Ca': 2, 'Sc': 3, 'Ti': 4, 'V': 5, 'Cr': 3, 'Mn': 2,
            'Fe': 3, 'Co': 2, 'Ni': 2, 'Cu': 2, 'Zn': 2, 'Ga': 3, 'Ge': 4,
            'As': 3, 'Se': -2, 'Br': -1, 'Rb': 1, 'Sr': 2, 'Y': 3, 'Zr': 4,
            'Nb': 5, 'Mo': 6, 'Tc': 7, 'Ru': 4, 'Rh': 3, 'Pd': 2, 'Ag': 1,
            'Cd': 2, 'In': 3, 'Sn': 4, 'Sb': 3, 'Te': -2, 'I': -1, 'Cs': 1,
            'Ba': 2, 'La': 3, 'Ce': 3, 'Pr': 3, 'Nd': 3, 'Pm': 3, 'Sm': 3,
            'Eu': 3, 'Gd': 3, 'Tb': 3, 'Dy': 3, 'Ho': 3, 'Er': 3, 'Tm': 3,
            'Yb': 3, 'Lu': 3, 'Hf': 4, 'Ta': 5, 'W': 6, 'Re': 7, 'Os': 4,
            'Ir': 3, 'Pt': 2, 'Au': 1, 'Hg': 2, 'Tl': 3, 'Pb': 2, 'Bi': 3,
            'Po': 4, 'At': -1, 'Rn': 0, 'Fr': 1, 'Ra': 2, 'Ac': 3, 'Th': 4,
            'Pa': 5, 'U': 6, 'Np': 5, 'Pu': 4, 'Am': 3, 'Cm': 3, 'Bk': 3,
            'Cf': 3, 'Es': 3, 'Fm': 3, 'Md': 3, 'No': 3, 'Lr': 3
        }
        
        # Use provided valence states or defaults
        if valence_states is not None:
            self.valence_states = valence_states.copy()
        else:
            self.valence_states = {}
            
        # Apply custom parameters if provided
        self.custom_parameters = custom_parameters or {}
        
        # Initialize analysis results
        self.bond_valence_sums = None
        self.bond_data = None
        self.neighbors = None
        
    def _get_valence_state(self, element: str) -> int:
        """Get valence state for an element."""
        if element in self.valence_states:
            return self.valence_states[element]
        elif element in self.default_valences:
            return self.default_valences[element]
        else:
            raise ValueError(f"No valence state specified for element {element}")
    
    def _find_neighbors(self):
        """Find neighbors for each atom using distance cutoff with proper minimum image convention."""
        # Use full distance_cutoff for NeighborList to ensure we don't miss neighbors
        # due to periodic boundary conditions
        cutoffs = [self.distance_cutoff] * len(self.atoms)
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(self.atoms)
        
        self.neighbors = {}
        for i in range(len(self.atoms)):
            indices, offsets = nl.get_neighbors(i)
            self.neighbors[i] = []
            
            for j, offset in zip(indices, offsets):
                # Calculate distance with proper minimum image convention
                distance = self.atoms.get_distance(i, j, mic=True, vector=False)
                if distance <= self.distance_cutoff:
                    self.neighbors[i].append({
                        'index': j,
                        'distance': distance,
                        'offset': offset
                    })
    
    def _calculate_bond_valence(self, element1: str, valence1: int, element2: str, valence2: int,
                               distance: float) -> float:
        """
        Calculate bond valence for a specific bond.
        
        Args:
            element1: First element symbol
            valence1: Valence state of first element
            element2: Second element symbol
            valence2: Valence state of second element
            distance: Bond length in Angstrom
            
        Returns:
            Bond valence value
        """
        # Check for custom parameters
        custom_key = f"{element1}-{element2}"
        if custom_key in self.custom_parameters:
            params = self.custom_parameters[custom_key]
            R0 = params['R0']
            B = params['B']
        else:
            # Get parameters from database
            try:
                params = self.bv_params.get_parameters(element1, valence1, element2, valence2)
                R0 = params['R0']
                B = params['B']
            except ValueError:
                # If no parameters found, return 0 (no contribution)
                return 0.0
        
        # Calculate bond valence using Brown's equation
        bond_valence = np.exp((R0 - distance) / B)
        return bond_valence
    
    def calculate_bvs(self) -> Dict[int, float]:
        """
        Calculate bond valence sums for all atoms.
        
        Returns:
            Dictionary mapping atom indices to bond valence sums
        """
        if self.neighbors is None:
            self._find_neighbors()
        
        self.bond_valence_sums = {}
        self.bond_data = {}
        
        for i in range(len(self.atoms)):
            element1 = self.atoms[i].symbol
            valence1 = self._get_valence_state(element1)
            
            bvs = 0.0
            bonds = []
            
            for neighbor in self.neighbors[i]:
                j = neighbor['index']
                distance = neighbor['distance']
                
                element2 = self.atoms[j].symbol
                valence2 = self._get_valence_state(element2)
                
                # Calculate bond valence
                bond_valence = self._calculate_bond_valence(element1, valence1, element2, valence2, distance)
                
                if bond_valence > 0:
                    bvs += bond_valence
                    bonds.append({
                        'neighbor_index': j,
                        'neighbor_element': element2,
                        'neighbor_valence': valence2,
                        'distance': distance,
                        'bond_valence': bond_valence
                    })
            
            self.bond_valence_sums[i] = bvs
            self.bond_data[i] = bonds
        
        return self.bond_valence_sums
    
    def analyze_structure(self) -> pd.DataFrame:
        """
        Analyze the structure and return a summary DataFrame.
        
        Returns:
            DataFrame with analysis results for each atom
        """
        if self.bond_valence_sums is None:
            self.calculate_bvs()
        
        analysis_data = []
        
        for i in range(len(self.atoms)):
            element = self.atoms[i].symbol
            expected_valence = abs(self._get_valence_state(element))
            calculated_bvs = self.bond_valence_sums[i]
            deviation = calculated_bvs - expected_valence
            coordination = len(self.bond_data[i])
            
            analysis_data.append({
                'atom_index': i,
                'element': element,
                'expected_valence': expected_valence,
                'calculated_bvs': calculated_bvs,
                'deviation': deviation,
                'coordination_number': coordination,
                'relative_deviation': deviation / expected_valence if expected_valence > 0 else 0
            })
        
        return pd.DataFrame(analysis_data)
    
    def get_bond_details(self, atom_index: int) -> List[Dict]:
        """
        Get detailed bond information for a specific atom.
        
        Args:
            atom_index: Index of the atom
            
        Returns:
            List of dictionaries with bond details
        """
        if self.bond_data is None:
            self.calculate_bvs()
        
        return self.bond_data.get(atom_index, [])
    
    def print_summary(self):
        """Print a summary of the bond valence analysis."""
        df = self.analyze_structure()
        
        print("Bond Valence Sum Analysis Summary")
        print("=" * 50)
        print(f"Total atoms: {len(self.atoms)}")
        print(f"Distance cutoff: {self.distance_cutoff:.2f} Å")
        print()
        
        # Group by element
        for element in sorted(df['element'].unique()):
            element_data = df[df['element'] == element]
            print(f"{element} atoms ({len(element_data)}):")
            print(f"  Expected valence: {element_data['expected_valence'].iloc[0]:.2f}")
            print(f"  Average BVS: {element_data['calculated_bvs'].mean():.3f} ± {element_data['calculated_bvs'].std():.3f}")
            print(f"  Average deviation: {element_data['deviation'].mean():.3f} ± {element_data['deviation'].std():.3f}")
            print(f"  Average coordination: {element_data['coordination_number'].mean():.1f}")
            print()