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
                 distance_cutoff: float = 3.5,
                 allowed_pairs: Optional[List[Tuple[str, str]]] = None,
                 exclude_same_element: bool = True,
                 auto_determine_valence: bool = False,
                 per_atom_valence: bool = False):
        """
        Initialize BondValenceSum calculator.
        
        Args:
            atoms: ASE Atoms object
            valence_states: Dict mapping element symbols to valence states
            custom_parameters: Dict with custom R0 and B values for specific elements
            distance_cutoff: Maximum distance for neighbor detection (Angstrom)
            allowed_pairs: List of element pairs to consider, e.g., [('Ti', 'O'), ('Fe', 'O')]
                          If None, will use default meaningful pairs based on elements present
            exclude_same_element: If True, exclude same-element pairs (e.g., Ti-Ti, O-O)
            auto_determine_valence: If True, automatically determine best valence states for metals
                                   by minimizing BVS deviation from expected values
            per_atom_valence: If True, optimize valence for each atom individually (allows mixed valences)
        """
        self.atoms = atoms
        self.distance_cutoff = distance_cutoff
        self.exclude_same_element = exclude_same_element
        self.auto_determine_valence = auto_determine_valence
        self.per_atom_valence = per_atom_valence
        
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
        
        # Set up allowed element pairs
        if allowed_pairs is not None:
            # Convert to set of tuples for fast lookup, including both orientations
            self.allowed_pairs = set()
            for pair in allowed_pairs:
                self.allowed_pairs.add(tuple(sorted(pair)))
        else:
            # Generate default allowed pairs based on elements present
            self.allowed_pairs = self._get_default_allowed_pairs()
        
        # Initialize analysis results
        self.bond_valence_sums = None
        self.bond_data = None
        self.neighbors = None
        
        # For auto-determined valence states
        self.optimized_valences = None
        self.valence_optimization_results = None
        self.per_atom_optimized_valences = None  # For per-atom optimization
        
    def _get_valence_state(self, element: str, atom_index: Optional[int] = None) -> int:
        """Get valence state for an element, with optional per-atom lookup."""
        # For per-atom optimization, check if we have a specific valence for this atom
        if self.per_atom_valence and atom_index is not None:
            per_atom_key = f"{element}_{atom_index}"
            if per_atom_key in self.valence_states:
                return self.valence_states[per_atom_key]
            # Also check the per_atom_optimized_valences dict
            if hasattr(self, 'per_atom_optimized_valences') and self.per_atom_optimized_valences:
                if atom_index in self.per_atom_optimized_valences:
                    return self.per_atom_optimized_valences[atom_index]
        
        # Standard element-wide valence lookup
        if element in self.valence_states:
            return self.valence_states[element]
        elif element in self.default_valences:
            return self.default_valences[element]
        else:
            raise ValueError(f"No valence state specified for element {element}")
    
    def _get_default_allowed_pairs(self) -> set:
        """
        Generate default allowed element pairs based on elements present in the structure.
        
        Returns sensible pairs like metal-oxygen, metal-halogen, etc. and excludes
        same-element pairs and non-meaningful pairs like metal-metal or hydrogen-metal.
        """
        elements = set(self.atoms.get_chemical_symbols())
        allowed_pairs = set()
        
        # Define element categories
        metals = {'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 
                 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
                 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce',
                 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi'}
        
        anions = {'O', 'F', 'Cl', 'Br', 'I', 'S', 'Se', 'Te', 'N', 'P', 'As'}
        
        # Generate meaningful pairs
        for elem1 in elements:
            for elem2 in elements:
                pair = tuple(sorted([elem1, elem2]))
                
                # Skip same-element pairs if requested
                if self.exclude_same_element and elem1 == elem2:
                    continue
                
                # Include metal-anion pairs (most common for BVS analysis)
                if (elem1 in metals and elem2 in anions) or (elem2 in metals and elem1 in anions):
                    allowed_pairs.add(pair)
                
                # Include metal-metal pairs only if specifically common (e.g., mixed oxides)
                # but this is rare, so we'll be conservative
                
                # Skip hydrogen-metal pairs (usually not relevant for BVS)
                if ('H' in pair and any(elem in metals for elem in pair)):
                    continue
                
                # Add other meaningful pairs if they have bond valence parameters
                try:
                    val1 = self._get_valence_state(elem1)
                    val2 = self._get_valence_state(elem2)
                    self.bv_params.get_parameters(elem1, val1, elem2, val2)
                    # If we get here, parameters exist - consider adding the pair
                    # But be selective to avoid unwanted pairs
                    if not ('H' in pair and any(elem in metals for elem in pair)):
                        allowed_pairs.add(pair)
                except (ValueError, KeyError):
                    # No parameters available, skip this pair
                    continue
        
        return allowed_pairs
    
    def _is_pair_allowed(self, element1: str, element2: str) -> bool:
        """
        Check if a specific element pair is allowed for bond valence calculation.
        
        Args:
            element1: First element
            element2: Second element
            
        Returns:
            True if the pair should be considered for BVS calculation
        """
        pair = tuple(sorted([element1, element2]))
        return pair in self.allowed_pairs
    
    def _get_possible_valences(self, element: str) -> List[int]:
        """
        Get possible valence states for an element from the bond valence parameter database.
        
        Args:
            element: Element symbol
            
        Returns:
            List of possible valence states
        """
        return self.bv_params.get_valence_states(element)
    
    def _is_metal_element(self, element: str) -> bool:
        """
        Determine if an element is a metal (for auto-valence determination).
        
        Args:
            element: Element symbol
            
        Returns:
            True if element is considered a metal
        """
        metals = {'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 
                 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
                 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce',
                 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi'}
        return element in metals
    
    def _determine_best_valence(self, atom_index: int, element: str) -> Tuple[int, float, Dict]:
        """
        Determine the best valence state for a metal atom by minimizing BVS deviation.
        
        Args:
            atom_index: Index of the atom
            element: Element symbol
            
        Returns:
            Tuple of (best_valence, best_deviation, optimization_results)
        """
        if self.neighbors is None:
            self._find_neighbors()
        
        possible_valences = self._get_possible_valences(element)
        if not possible_valences:
            # Fall back to default valence if no options available
            default_val = self._get_valence_state(element)
            return default_val, 0.0, {'tried_valences': [default_val], 'deviations': [0.0]}
        
        optimization_results = {'tried_valences': [], 'deviations': [], 'bvs_values': []}
        best_valence = None
        best_deviation = float('inf')
        
        for valence in possible_valences:
            # Temporarily set this valence
            original_valence = self.valence_states.get(element, self.default_valences.get(element, valence))
            self.valence_states[element] = valence
            
            # Calculate BVS for this atom with this valence
            bvs = 0.0
            for neighbor in self.neighbors[atom_index]:
                j = neighbor['index']
                distance = neighbor['distance']
                
                element2 = self.atoms[j].symbol
                valence2 = self._get_valence_state(element2)
                
                # Check if this element pair is allowed
                if not self._is_pair_allowed(element, element2):
                    continue
                
                # Calculate bond valence
                bond_valence = self._calculate_bond_valence(element, valence, element2, valence2, distance)
                if bond_valence > 0:
                    bvs += bond_valence
            
            # Calculate deviation from expected BVS (which equals the valence magnitude)
            expected_bvs = abs(valence)
            deviation = abs(bvs - expected_bvs)
            
            optimization_results['tried_valences'].append(valence)
            optimization_results['deviations'].append(deviation)
            optimization_results['bvs_values'].append(bvs)
            
            if deviation < best_deviation:
                best_deviation = deviation
                best_valence = valence
            
            # Restore original valence
            if element in self.valence_states:
                self.valence_states[element] = original_valence
        
        return best_valence, best_deviation, optimization_results
    
    def _auto_determine_valences_per_atom(self):
        """
        Automatically determine the best valence state for each individual atom.
        Allows different atoms of the same element to have different valences.
        """
        if self.neighbors is None:
            self._find_neighbors()
        
        self.per_atom_optimized_valences = {}
        self.valence_optimization_results = {}
        
        symbols = self.atoms.get_chemical_symbols()
        
        # Group atoms by element for reporting
        elements = set(symbols)
        element_results = {}
        
        for element in elements:
            if self._is_metal_element(element):
                atom_indices = [i for i, sym in enumerate(symbols) if sym == element]
                element_results[element] = {}
                
                for atom_idx in atom_indices:
                    best_valence, best_deviation, opt_results = self._determine_best_valence(atom_idx, element)
                    
                    # Store per-atom optimized valence
                    self.per_atom_optimized_valences[atom_idx] = best_valence
                    
                    # Store results for reporting
                    element_results[element][atom_idx] = {
                        'best_valence': best_valence,
                        'best_deviation': best_deviation,
                        'optimization_results': opt_results
                    }
                    
                    # Update valence_states for this specific atom
                    # We'll use a special key format for per-atom valences
                    self.valence_states[f"{element}_{atom_idx}"] = best_valence
        
        self.valence_optimization_results = element_results

    def _auto_determine_valences(self):
        """
        Automatically determine the best valence states for all metal atoms.
        """
        if self.per_atom_valence:
            return self._auto_determine_valences_per_atom()
        
        if self.neighbors is None:
            self._find_neighbors()
        
        self.optimized_valences = {}
        self.valence_optimization_results = {}
        
        # Get unique elements in the structure
        elements = set(self.atoms.get_chemical_symbols())
        
        for element in elements:
            if self._is_metal_element(element):
                # Get all atom indices for this element
                atom_indices = [i for i, sym in enumerate(self.atoms.get_chemical_symbols()) if sym == element]
                
                # Determine best valence for each atom of this element
                element_results = {}
                valence_votes = {}
                
                for atom_idx in atom_indices:
                    best_valence, best_deviation, opt_results = self._determine_best_valence(atom_idx, element)
                    element_results[atom_idx] = {
                        'best_valence': best_valence,
                        'best_deviation': best_deviation,
                        'optimization_results': opt_results
                    }
                    
                    # Vote for the best valence across all atoms of this element
                    if best_valence not in valence_votes:
                        valence_votes[best_valence] = 0
                    valence_votes[best_valence] += 1
                
                # Choose the valence with the lowest average deviation across ALL atoms
                if valence_votes:
                    # Get all possible valences tried for this element
                    all_valences = set()
                    for result in element_results.values():
                        all_valences.update(result['optimization_results']['tried_valences'])
                    
                    # Calculate average deviation for each valence across ALL atoms
                    valence_avg_deviations = {}
                    for valence in all_valences:
                        total_deviation = 0
                        count = 0
                        
                        for result in element_results.values():
                            opt_results = result['optimization_results']
                            if valence in opt_results['tried_valences']:
                                valence_idx = opt_results['tried_valences'].index(valence)
                                deviation = opt_results['deviations'][valence_idx]
                                total_deviation += deviation
                                count += 1
                        
                        if count > 0:
                            valence_avg_deviations[valence] = total_deviation / count
                    
                    # Find valence with lowest average deviation across all atoms
                    best_element_valence = min(valence_avg_deviations.keys(), 
                                             key=lambda v: valence_avg_deviations[v])
                    
                    self.optimized_valences[element] = best_element_valence
                    self.valence_optimization_results[element] = element_results
                    
                    # Update the valence_states with the optimized value
                    self.valence_states[element] = best_element_valence
    
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
        
        # Auto-determine valences if requested
        if self.auto_determine_valence:
            self._auto_determine_valences()
        
        self.bond_valence_sums = {}
        self.bond_data = {}
        
        for i in range(len(self.atoms)):
            element1 = self.atoms[i].symbol
            valence1 = self._get_valence_state(element1, i)
            
            bvs = 0.0
            bonds = []
            
            for neighbor in self.neighbors[i]:
                j = neighbor['index']
                distance = neighbor['distance']
                
                element2 = self.atoms[j].symbol
                valence2 = self._get_valence_state(element2, j)
                
                # Check if this element pair is allowed
                if not self._is_pair_allowed(element1, element2):
                    continue
                
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
            current_valence = self._get_valence_state(element, i)
            expected_valence = abs(current_valence)
            calculated_bvs = self.bond_valence_sums[i]
            deviation = calculated_bvs - expected_valence
            coordination = len(self.bond_data[i])
            
            # Base analysis data
            atom_data = {
                'atom_index': i,
                'element': element,
                'used_valence': current_valence,
                'expected_valence': expected_valence,
                'calculated_bvs': calculated_bvs,
                'deviation': deviation,
                'coordination_number': coordination,
                'relative_deviation': deviation / expected_valence if expected_valence > 0 else 0
            }
            
            # Add valence optimization results if available
            if self.auto_determine_valence and self.valence_optimization_results:
                if element in self.valence_optimization_results:
                    if i in self.valence_optimization_results[element]:
                        opt_results = self.valence_optimization_results[element][i]
                        atom_data.update({
                            'valence_optimized': True,
                            'optimization_deviation': opt_results['best_deviation'],
                            'tried_valences': str(opt_results['optimization_results']['tried_valences']),
                            'tried_deviations': str([f"{d:.3f}" for d in opt_results['optimization_results']['deviations']])
                        })
                    else:
                        atom_data.update({
                            'valence_optimized': False,
                            'optimization_deviation': None,
                            'tried_valences': None,
                            'tried_deviations': None
                        })
                else:
                    atom_data.update({
                        'valence_optimized': False,
                        'optimization_deviation': None,
                        'tried_valences': None,
                        'tried_deviations': None
                    })
            
            analysis_data.append(atom_data)
        
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
    
    def get_allowed_pairs(self) -> List[Tuple[str, str]]:
        """
        Get the list of element pairs that are considered in BVS calculations.
        
        Returns:
            List of allowed element pairs as tuples
        """
        return sorted(list(self.allowed_pairs))
    
    def print_allowed_pairs(self):
        """Print the allowed element pairs for bond valence calculations."""
        pairs = self.get_allowed_pairs()
        print("Allowed element pairs for bond valence calculations:")
        for i, pair in enumerate(pairs):
            print(f"  {i+1:2d}. {pair[0]}-{pair[1]}")
        print(f"\nTotal: {len(pairs)} pairs")
        if self.exclude_same_element:
            print("Note: Same-element pairs are excluded")
    
    def get_valence_optimization_results(self) -> Dict:
        """
        Get detailed results from valence optimization.
        
        Returns:
            Dictionary with optimization results for each element
        """
        if not self.auto_determine_valence:
            return {}
        
        if self.valence_optimization_results is None:
            self.calculate_bvs()  # This will trigger optimization
        
        return self.valence_optimization_results or {}
    
    def print_valence_optimization_summary(self):
        """Print a summary of the valence optimization results."""
        if not self.auto_determine_valence:
            print("Valence optimization was not performed.")
            return
        
        opt_results = self.get_valence_optimization_results()
        if not opt_results:
            print("No valence optimization results available.")
            return
        
        print("Valence Optimization Summary")
        print("=" * 40)
        
        for element, element_results in opt_results.items():
            print(f"\n{element} atoms:")
            if element in self.optimized_valences:
                print(f"  Optimized valence: {self.optimized_valences[element]:+d}")
            
            for atom_idx, result in element_results.items():
                opt_data = result['optimization_results']
                best_valence = result['best_valence']
                best_deviation = result['best_deviation']
                
                print(f"  Atom {atom_idx}:")
                print(f"    Best valence: {best_valence:+d} (deviation: {best_deviation:.3f})")
                print(f"    Tried valences: {opt_data['tried_valences']}")
                print(f"    Deviations: {[f'{d:.3f}' for d in opt_data['deviations']]}")
                print(f"    BVS values: {[f'{bvs:.3f}' for bvs in opt_data['bvs_values']]}")
    
    def print_summary(self):
        """Print a summary of the bond valence analysis."""
        df = self.analyze_structure()
        
        print("Bond Valence Sum Analysis Summary")
        print("=" * 50)
        print(f"Total atoms: {len(self.atoms)}")
        print(f"Distance cutoff: {self.distance_cutoff:.2f} Å")
        if self.auto_determine_valence:
            print("Valence optimization: ENABLED")
        else:
            print("Valence optimization: DISABLED")
        print()
        
        # Group by element
        for element in sorted(df['element'].unique()):
            element_data = df[df['element'] == element]
            print(f"{element} atoms ({len(element_data)}):")
            
            # Show valence information
            used_valence = element_data['used_valence'].iloc[0]
            expected_valence = element_data['expected_valence'].iloc[0]
            print(f"  Used valence: {used_valence:+d} (expected: {expected_valence:.2f})")
            
            # Show optimization info if available
            if self.auto_determine_valence and 'valence_optimized' in element_data.columns:
                optimized_count = element_data['valence_optimized'].sum()
                if optimized_count > 0:
                    avg_opt_deviation = element_data[element_data['valence_optimized']]['optimization_deviation'].mean()
                    print(f"  Optimized atoms: {optimized_count}/{len(element_data)} (avg deviation: {avg_opt_deviation:.3f})")
            
            print(f"  Average BVS: {element_data['calculated_bvs'].mean():.3f} ± {element_data['calculated_bvs'].std():.3f}")
            print(f"  Average deviation: {element_data['deviation'].mean():.3f} ± {element_data['deviation'].std():.3f}")
            print(f"  Average coordination: {element_data['coordination_number'].mean():.1f}")
            print()
        
        # Print detailed valence optimization results if available
        if self.auto_determine_valence:
            print("-" * 50)
            self.print_valence_optimization_summary()