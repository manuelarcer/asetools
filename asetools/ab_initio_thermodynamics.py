"""
Ab Initio Thermodynamics Module

This module provides functionality for calculating equilibrium adsorbate coverage and 
surface energies under varying temperature and pressure conditions using ab initio 
thermodynamics principles.

The module supports both simple coverage-dependent interpolation models and lattice gas 
models with adsorbate-adsorbate interactions, based on approaches from:
- Zhu, Meng, and Gao, J. Phys. Chem. C 121 (2017) 5629-5634
- Standard Langmuir isotherm and competitive adsorption theory

Classes:
    AdsorbateSpecies: Handles adsorbate properties (entropy, dissociation)
    SurfaceProperties: Manages surface energies and areas for different facets
    InterpolationModel: Simple coverage-Eads interpolation from CSV data
    LatticeGasModel: Interaction-based coverage dependence
    ThermodynamicsCalculator: Main calculation engine

Author: ASEtools development team
"""

import numpy as np
import pandas as pd
import scipy.optimize as opt
import math
import random
import os
import json
from typing import Dict, List, Tuple, Optional, Union
import logging

# Constants
EV_TO_KJ_MOL = 96.485  # eV to kJ/mol conversion
P_STANDARD = 100000    # Standard pressure in Pa
R_GAS = 0.008314463    # Gas constant in kJ/(mol·K)

logger = logging.getLogger(__name__)


class AdsorbateSpecies:
    """
    Handles adsorbate species properties including entropy functions and dissociation behavior.
    
    Entropy is calculated using the empirical relationship: S = a*T^b - R*ln(P/P0)
    where a and b are species-specific parameters.
    """
    
    # Default entropy parameters from notebooks
    DEFAULT_ENTROPY_PARAMS = {
        'CO': {'a': 88.4037336109, 'b': 0.141268928},
        'O': {'a': 87.9953566884, 'b': 0.1480182342},
        'H': {'a': 87.9953566884, 'b': 0.1480182342},  # Placeholder, user should provide
        'N': {'a': 87.9953566884, 'b': 0.1480182342},  # Placeholder, user should provide
    }
    
    # Default dissociation behavior (True for H2, O2, N2)
    DEFAULT_DISSOCIATION = {
        'CO': False,   # Molecular adsorption
        'O': True,     # From O2 dissociation
        'H': True,     # From H2 dissociation  
        'N': True,     # From N2 dissociation
    }
    
    def __init__(self, name: str, entropy_params: Optional[Dict[str, float]] = None, 
                 dissociative: Optional[bool] = None):
        """
        Initialize adsorbate species.
        
        Args:
            name: Adsorbate species name (e.g., 'CO', 'O', 'H')
            entropy_params: Dict with 'a' and 'b' parameters for S = a*T^b
            dissociative: Whether adsorption is dissociative (None = use default)
        """
        self.name = name
        
        if entropy_params is None:
            self.entropy_params = self.DEFAULT_ENTROPY_PARAMS.get(name, 
                {'a': 88.0, 'b': 0.14})  # Generic default
        else:
            self.entropy_params = entropy_params
            
        if dissociative is None:
            self.dissociative = self.DEFAULT_DISSOCIATION.get(name, False)
        else:
            self.dissociative = dissociative
    
    def entropy(self, temperature: float, pressure: float) -> float:
        """
        Calculate entropy of gas phase adsorbate at given T,P.
        
        Args:
            temperature: Temperature in K
            pressure: Pressure in Pa
            
        Returns:
            Entropy in kJ/(mol·K)
        """
        a = self.entropy_params['a']
        b = self.entropy_params['b'] 
        
        S = a * temperature**b / 1000 - R_GAS * np.log(pressure / P_STANDARD)
        return S


class SurfaceProperties:
    """
    Manages surface energies and areas for different metal facets.
    """
    
    def __init__(self, surface_data_file: Optional[str] = None):
        """
        Initialize surface properties.
        
        Args:
            surface_data_file: Path to JSON file with surface data (None = use default)
        """
        if surface_data_file is None:
            # Use default database file
            current_dir = os.path.dirname(os.path.abspath(__file__))
            surface_data_file = os.path.join(current_dir, 'data', 'surface_properties.json')
        
        with open(surface_data_file, 'r') as f:
            data = json.load(f)
        self.surface_energies = data['surface_energies']
        self.surface_areas = data['surface_areas']
    
    def get_surface_energy(self, metal: str, facet: str) -> float:
        """Get surface energy in eV/Å²."""
        return self.surface_energies[metal][facet]
    
    def get_surface_area(self, metal: str, facet: str) -> float:
        """Get surface area per atom in Å²."""
        return self.surface_areas[metal][facet]
    
    def get_available_metals(self) -> List[str]:
        """Get list of available metals."""
        return list(self.surface_energies.keys())
    
    def get_available_facets(self, metal: str) -> List[str]:
        """Get list of available facets for a metal."""
        return list(self.surface_energies[metal].keys())


class InterpolationModel:
    """
    Simple coverage-dependent adsorption energy model using interpolation.
    Accepts data from CSV file or directly as pandas DataFrame.
    Data must contain columns: ['Catalyst', 'facet', 'Ads', 'cov', 'Eads']
    """
    
    def __init__(self, data_source: Union[str, pd.DataFrame]):
        """
        Initialize interpolation model from CSV file or DataFrame.
        
        Args:
            data_source: Either path to CSV file or pandas DataFrame with adsorption energy data
        """
        if isinstance(data_source, str):
            self.data = pd.read_csv(data_source)
        elif isinstance(data_source, pd.DataFrame):
            self.data = data_source.copy()
        else:
            raise ValueError("data_source must be either a file path (str) or pandas DataFrame")
        
        self._validate_data()
    
    def _validate_data(self):
        """Validate CSV data format."""
        required_columns = ['Catalyst', 'facet', 'Ads', 'cov', 'Eads']
        if not all(col in self.data.columns for col in required_columns):
            raise ValueError(f"CSV must contain columns: {required_columns}")
        
        # Convert facet to string to handle both string and numeric facet names
        self.data['facet'] = self.data['facet'].astype(str)
    
    def get_adsorption_energy(self, catalyst: str, facet: str, adsorbate: str, 
                            coverage: float) -> float:
        """
        Get adsorption energy for given conditions using interpolation.
        
        Args:
            catalyst: Metal catalyst name
            facet: Surface facet 
            adsorbate: Adsorbate species name
            coverage: Coverage (0-1)
            
        Returns:
            Adsorption energy in eV
        """
        # Filter data for specific conditions
        mask = ((self.data['Catalyst'] == catalyst) & 
                (self.data['facet'] == facet) & 
                (self.data['Ads'] == adsorbate))
        subset = self.data[mask]
        
        if len(subset) == 0:
            raise ValueError(f"No data found for {catalyst}/{facet}/{adsorbate}")
        
        # Interpolate
        coverages = subset['cov'].values
        energies = subset['Eads'].values
        
        # Sort by coverage
        sort_idx = np.argsort(coverages)
        coverages = coverages[sort_idx]
        energies = energies[sort_idx]
        
        return np.interp(coverage, coverages, energies)


class LatticeGasModel:
    """
    Lattice gas model with adsorbate-adsorbate interactions.
    Implements the approach from the Pd notebook with interaction parameters.
    """
    
    def __init__(self, interaction_params: Dict):
        """
        Initialize lattice gas model.
        
        Args:
            interaction_params: Dictionary with interaction parameters
                Format: {
                    'metal': {
                        'facet': {
                            'coordination': z_value,
                            'zero_coverage_energies': {'species': energy_eV},
                            'self_interactions': {'species': w_value_eV},
                            'cross_interactions': {('species1', 'species2'): w_value_eV}
                        }
                    }
                }
        """
        self.params = interaction_params
    
    def get_adsorption_energy(self, catalyst: str, facet: str, adsorbate: str,
                            coverages: Dict[str, float]) -> float:
        """
        Calculate coverage-dependent adsorption energy using lattice gas model.
        
        Args:
            catalyst: Metal catalyst
            facet: Surface facet
            adsorbate: Target adsorbate species
            coverages: Dict of coverage values for all species
            
        Returns:
            Adsorption energy in eV
        """
        facet_params = self.params[catalyst][facet]
        z = facet_params['coordination']  # Coordination number
        E_zero = facet_params['zero_coverage_energies'][adsorbate]
        
        # Calculate interaction energy
        interaction_energy = 0.0
        
        # Self-interaction
        if adsorbate in facet_params['self_interactions']:
            w_self = facet_params['self_interactions'][adsorbate]
            interaction_energy += w_self * coverages[adsorbate]
        
        # Cross-interactions with other species
        for other_species, coverage in coverages.items():
            if other_species != adsorbate:
                key1 = (adsorbate, other_species)
                key2 = (other_species, adsorbate)
                if key1 in facet_params['cross_interactions']:
                    w_cross = facet_params['cross_interactions'][key1]
                    interaction_energy += w_cross * coverage
                elif key2 in facet_params['cross_interactions']:
                    w_cross = facet_params['cross_interactions'][key2]
                    interaction_energy += w_cross * coverage
        
        # For dissociative adsorption (O), multiply by 2 as in notebook
        if adsorbate in ['O', 'H', 'N']:  # Dissociative species
            return 2 * E_zero - 2 * z * interaction_energy
        else:  # Molecular adsorption
            return E_zero - z * interaction_energy


class ThermodynamicsCalculator:
    """
    Main calculator for ab initio thermodynamic calculations.
    """
    
    def __init__(self, metal: str, adsorbates: List[str], 
                 surface_properties: Optional[SurfaceProperties] = None):
        """
        Initialize thermodynamics calculator.
        
        Args:
            metal: Metal catalyst name
            adsorbates: List of adsorbate species
            surface_properties: SurfaceProperties object (None = use default)
        """
        self.metal = metal
        self.adsorbates = adsorbates
        self.surface_props = surface_properties or SurfaceProperties()
        self.species = {name: AdsorbateSpecies(name) for name in adsorbates}
        self.model = None  # Will be set by load_* methods
        
        # Validate metal exists
        if metal not in self.surface_props.get_available_metals():
            raise ValueError(f"Metal '{metal}' not in surface database")
    
    def load_interpolation_model(self, data_source: Union[str, pd.DataFrame]):
        """
        Load interpolation model from CSV file or DataFrame.
        
        Args:
            data_source: Either path to CSV file or pandas DataFrame with adsorption energy data
        """
        self.model = InterpolationModel(data_source)
        self.model_type = 'interpolation'
    
    def load_lattice_gas_model(self, interaction_params: Dict):
        """Load lattice gas model with interaction parameters.""" 
        self.model = LatticeGasModel(interaction_params)
        self.model_type = 'lattice_gas'
    
    def _ads_energy(self, adsorbate: str, facet: str, coverages: Union[Dict[str, float], float]) -> float:
        """
        Calculate adsorption energy for given conditions.
        
        Args:
            adsorbate: Adsorbate species name
            facet: Surface facet
            coverages: Coverage values (dict for lattice gas, float for interpolation)
            
        Returns:
            Adsorption energy in kJ/mol
        """
        if self.model is None:
            raise ValueError("No adsorption model loaded. Use load_interpolation_model() or load_lattice_gas_model()")
        
        if self.model_type == 'interpolation':
            # For single adsorbate interpolation
            if isinstance(coverages, dict):
                coverage = coverages.get(adsorbate, 0.0)
            else:
                coverage = coverages
            eads_ev = self.model.get_adsorption_energy(self.metal, facet, adsorbate, coverage)
        
        elif self.model_type == 'lattice_gas':
            # For multi-adsorbate lattice gas
            if not isinstance(coverages, dict):
                raise ValueError("Lattice gas model requires coverage dict")
            eads_ev = self.model.get_adsorption_energy(self.metal, facet, adsorbate, coverages)
        
        return eads_ev * EV_TO_KJ_MOL
    
    def _system_equilibrium_single(self, coverage: float, adsorbate: str, facet: str,
                                 temperature: float, pressure: float) -> float:
        """
        Equilibrium equation for single adsorbate system.
        
        Returns:
            Objective function value (should be minimized to 0)
        """
        species = self.species[adsorbate]
        
        # Adsorption energy at this coverage
        eads = self._ads_energy(adsorbate, facet, coverage)
        
        # Entropy
        entropy = species.entropy(temperature, pressure)
        
        # Equilibrium equations
        if species.dissociative:
            # Dissociative: A2(g) + 2* ⇌ 2A*
            left = coverage**2
            right = (pressure / P_STANDARD) * np.exp(-(eads - temperature * (-entropy)) / (R_GAS * temperature)) * (1 - coverage)**2
        else:
            # Molecular: A(g) + * ⇌ A*
            left = coverage
            right = (pressure / P_STANDARD) * np.exp(-(eads - temperature * (-entropy)) / (R_GAS * temperature)) * (1 - coverage)
        
        return abs(left - right)
    
    def _system_equilibrium_multi(self, coverages_list: List[float], facet: str,
                                temperature: float, pressures: Dict[str, float]) -> float:
        """
        Equilibrium equation for multi-adsorbate competitive adsorption.
        
        Args:
            coverages_list: List of coverage values for each adsorbate
            facet: Surface facet
            temperature: Temperature in K
            pressures: Dict of pressures for each adsorbate
            
        Returns:
            Objective function value
        """
        # Convert list to dict
        coverages = {adsorbate: cov for adsorbate, cov in zip(self.adsorbates, coverages_list)}
        
        total_objective = 0.0
        
        for adsorbate, coverage in coverages.items():
            species = self.species[adsorbate]
            pressure = pressures[adsorbate]
            
            # Adsorption energy with interactions
            eads = self._ads_energy(adsorbate, facet, coverages)
            
            # Entropy
            entropy = species.entropy(temperature, pressure)
            
            # Available sites
            total_coverage = sum(coverages.values())
            available_sites = 1 - total_coverage
            
            # Equilibrium equations
            if species.dissociative:
                left = coverage**2
                right = (pressure / P_STANDARD) * np.exp(-(eads - temperature * (-entropy)) / (R_GAS * temperature)) * available_sites**2
            else:
                left = coverage
                right = (pressure / P_STANDARD) * np.exp(-(eads - temperature * (-entropy)) / (R_GAS * temperature)) * available_sites
            
            total_objective += abs(left - right)
        
        return total_objective
    
    def _constraint_total_coverage(self, coverages_list: List[float]) -> float:
        """Constraint: total coverage ≤ 1."""
        return 1.0 - sum(coverages_list)
    
    def calculate_equilibrium_coverage(self, facet: str, temperature: float, 
                                     pressures: Union[float, Dict[str, float]], 
                                     max_iterations: int = 100) -> Dict[str, float]:
        """
        Calculate equilibrium coverage for given conditions.
        
        Args:
            facet: Surface facet
            temperature: Temperature in K
            pressures: Pressure (Pa) - float for single adsorbate, dict for multiple
            max_iterations: Maximum optimization attempts
            
        Returns:
            Dict of equilibrium coverages for each adsorbate
        """
        if len(self.adsorbates) == 1:
            # Single adsorbate system
            adsorbate = self.adsorbates[0]
            pressure = pressures if isinstance(pressures, (int, float)) else pressures[adsorbate]
            
            # Solve equilibrium
            best_solution = None
            best_obj = float('inf')
            
            for _ in range(max_iterations):
                initial_guess = random.random()
                
                try:
                    solution = opt.minimize(
                        self._system_equilibrium_single,
                        initial_guess,
                        args=(adsorbate, facet, temperature, pressure),
                        method='SLSQP',
                        bounds=[(0, 1)],
                        options={'eps': 1e-12, 'ftol': 1e-15, 'maxiter': 5000}
                    )
                    
                    if solution.fun < best_obj and solution.success:
                        best_solution = solution
                        best_obj = solution.fun
                        
                    # Check for high coverage convergence
                    if solution.x[0] >= 0.999:
                        return {adsorbate: 1.0}
                    
                    # Good enough convergence
                    if solution.fun <= 1e-5 and not math.isnan(solution.fun):
                        return {adsorbate: round(solution.x[0], 4)}
                        
                except:
                    continue
            
            # Use best solution found
            if best_solution is not None and best_obj <= 1e-4:
                return {adsorbate: round(best_solution.x[0], 4)}
            else:
                logger.warning(f"Convergence not achieved for T={temperature}, facet={facet}. Setting coverage to 1.")
                return {adsorbate: 1.0}
        
        else:
            # Multi-adsorbate system
            if not isinstance(pressures, dict):
                raise ValueError("Multi-adsorbate system requires pressure dict")
            
            constraint = {'type': 'ineq', 'fun': self._constraint_total_coverage}
            
            best_solution = None
            best_obj = float('inf')
            
            for _ in range(max_iterations):
                # Random initial guess
                initial_coverages = [random.random() for _ in self.adsorbates]
                # Normalize to satisfy constraint
                total = sum(initial_coverages)
                if total > 1:
                    initial_coverages = [c/total * 0.9 for c in initial_coverages]
                
                try:
                    solution = opt.minimize(
                        self._system_equilibrium_multi,
                        initial_coverages,
                        args=(facet, temperature, pressures),
                        method='SLSQP',
                        bounds=[(0, 1) for _ in self.adsorbates],
                        constraints=constraint,
                        options={'eps': 1e-8, 'ftol': 1e-12, 'maxiter': 5000}
                    )
                    
                    if solution.fun < best_obj and solution.success:
                        best_solution = solution
                        best_obj = solution.fun
                    
                    # Good convergence
                    if solution.fun <= 1e-4 and not math.isnan(solution.fun):
                        return {adsorbate: round(cov, 4) for adsorbate, cov in 
                               zip(self.adsorbates, solution.x)}
                        
                except:
                    continue
            
            # Use best solution
            if best_solution is not None:
                return {adsorbate: round(cov, 4) for adsorbate, cov in 
                       zip(self.adsorbates, best_solution.x)}
            else:
                logger.warning(f"Convergence failed for T={temperature}, facet={facet}")
                return {adsorbate: 0.0 for adsorbate in self.adsorbates}
    
    def calculate_surface_energy(self, facet: str, coverages: Dict[str, float],
                               temperature: float, pressures: Union[float, Dict[str, float]]) -> float:
        """
        Calculate surface energy including adsorbate contributions.
        
        Args:
            facet: Surface facet
            coverages: Coverage dict for each adsorbate
            temperature: Temperature in K
            pressures: Pressure(s) in Pa
            
        Returns:
            Surface energy in eV/Å²
        """
        # Clean surface energy
        gamma_clean = self.surface_props.get_surface_energy(self.metal, facet)
        
        # Area per atom
        area_per_atom = self.surface_props.get_surface_area(self.metal, facet)
        
        # Adsorbate contribution
        adsorbate_contribution = 0.0
        
        for adsorbate, coverage in coverages.items():
            if coverage > 0:
                # Get pressure for this adsorbate
                if isinstance(pressures, dict):
                    pressure = pressures[adsorbate]
                else:
                    pressure = pressures
                
                # Adsorption energy
                eads = self._ads_energy(adsorbate, facet, coverages)  # kJ/mol
                eads_ev = eads / EV_TO_KJ_MOL  # Convert to eV
                
                # Add to surface energy 
                if self.species[adsorbate].dissociative:
                    # For dissociative species, divide by 2
                    adsorbate_contribution += (coverage * eads_ev / 2) / area_per_atom
                else:
                    # For molecular species
                    adsorbate_contribution += (coverage * eads_ev) / area_per_atom
        
        return gamma_clean + adsorbate_contribution
    
    def calculate_equilibrium(self, temperature_range: Tuple[float, float, float],
                            pressures: Union[float, Dict[str, float]], 
                            facets: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Calculate equilibrium coverage and surface energies over temperature range.
        
        Args:
            temperature_range: (start, stop, step) for temperature in K
            pressures: Pressure(s) in Pa
            facets: List of facets to calculate (None = all available)
            
        Returns:
            DataFrame with results
        """
        if facets is None:
            facets = self.surface_props.get_available_facets(self.metal)
        
        temperatures = np.arange(*temperature_range)
        
        # Initialize data dictionary
        data = {'Temp': []}
        
        for facet in facets:
            for adsorbate in self.adsorbates:
                data[f'Cov_{adsorbate}_{facet}'] = []
                data[f'Eads_{adsorbate}_{facet}'] = []
            data[f'Gamma_{facet}'] = []
        
        # Calculate for each temperature
        for T in temperatures:
            # Add temperature only once
            if len(data['Temp']) == 0 or data['Temp'][-1] != T:
                data['Temp'].append(T)
            
            for facet in facets:
                # Calculate equilibrium coverage
                coverages = self.calculate_equilibrium_coverage(facet, T, pressures)
                
                # Store coverage and adsorption energies
                for adsorbate in self.adsorbates:
                    coverage = coverages.get(adsorbate, 0.0)
                    data[f'Cov_{adsorbate}_{facet}'].append(coverage)
                    
                    if coverage > 0:
                        eads = self._ads_energy(adsorbate, facet, coverages) / EV_TO_KJ_MOL
                    else:
                        eads = 0.0
                    data[f'Eads_{adsorbate}_{facet}'].append(eads)
                
                # Calculate surface energy
                gamma = self.calculate_surface_energy(facet, coverages, T, pressures)
                data[f'Gamma_{facet}'].append(gamma)
        
        return pd.DataFrame(data)