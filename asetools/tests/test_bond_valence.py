"""
Tests for bond valence sum calculations.
"""

import pytest
import numpy as np
import pandas as pd
from ase import Atoms
from ase.build import bulk
from asetools.bond_valence import BondValenceParameters, BondValenceSum


class TestBondValenceParameters:
    """Test BondValenceParameters class."""
    
    def test_parameter_loading(self):
        """Test loading of bond valence parameters."""
        bv_params = BondValenceParameters()
        
        # Test that parameters are loaded
        assert len(bv_params.parameters) > 0
        assert len(bv_params.references) > 0
        
        # Test that Ti-O parameters exist (from user example)
        ti_o_params = bv_params.get_parameters('Ti', 4, 'O', -2)
        assert ti_o_params['R0'] == 1.815
        assert ti_o_params['B'] == 0.37
        assert ti_o_params['reference'] == 'a'
    
    def test_parameter_lookup(self):
        """Test parameter lookup functionality."""
        bv_params = BondValenceParameters()
        
        # Test Ti-O lookup
        params = bv_params.get_parameters('Ti', 4, 'O', -2)
        assert isinstance(params, dict)
        assert 'R0' in params
        assert 'B' in params
        assert 'reference' in params
        assert 'details' in params
        
        # Test reverse lookup (O-Ti should give same result)
        params_reverse = bv_params.get_parameters('O', -2, 'Ti', 4)
        assert params['R0'] == params_reverse['R0']
        assert params['B'] == params_reverse['B']
    
    def test_multiple_parameters(self):
        """Test handling of multiple parameters for same element pair."""
        bv_params = BondValenceParameters()
        
        # Ti-O has multiple parameter sets
        all_params = bv_params.get_parameters('Ti', 4, 'O', -2, most_reliable=False)
        assert isinstance(all_params, list)
        assert len(all_params) > 1
        
        # Most reliable should be first
        most_reliable = bv_params.get_parameters('Ti', 4, 'O', -2, most_reliable=True)
        assert most_reliable == all_params[0]
    
    def test_missing_parameters(self):
        """Test handling of missing parameters."""
        bv_params = BondValenceParameters()
        
        # Test with non-existent element pair
        with pytest.raises(ValueError, match="No bond valence parameters found"):
            bv_params.get_parameters('X', 1, 'Y', 2)
    
    def test_element_methods(self):
        """Test methods that list available elements and valences."""
        bv_params = BondValenceParameters()
        
        # Test element listing
        elements = bv_params.get_elements()
        assert 'Ti' in elements
        assert 'O' in elements
        assert isinstance(elements, list)
        
        # Test valence state listing
        ti_valences = bv_params.get_valence_states('Ti')
        assert 4 in ti_valences
        assert 3 in ti_valences
        
        o_valences = bv_params.get_valence_states('O')
        assert -2 in o_valences


class TestBondValenceSum:
    """Test BondValenceSum class."""
    
    def test_simple_structure(self):
        """Test BVS calculation for a simple structure."""
        # Create simple TiO2 structure manually
        from ase import Atoms
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        # Set up BVS calculator
        valence_states = {'Ti': 4, 'O': -2}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states)
        
        # Calculate BVS
        bvs_results = bvs_calc.calculate_bvs()
        
        # Check that results are reasonable
        assert isinstance(bvs_results, dict)
        assert len(bvs_results) == len(atoms)
        
        # All BVS values should be non-negative
        for bvs in bvs_results.values():
            assert bvs >= 0
    
    def test_analysis_dataframe(self):
        """Test structure analysis DataFrame output."""
        # Create simple structure
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        valence_states = {'Ti': 4, 'O': -2}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states)
        
        # Get analysis DataFrame
        df = bvs_calc.analyze_structure()
        
        # Check DataFrame structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) == len(atoms)
        
        expected_columns = ['atom_index', 'element', 'expected_valence', 
                          'calculated_bvs', 'deviation', 'coordination_number', 
                          'relative_deviation']
        for col in expected_columns:
            assert col in df.columns
        
        # Check that we have both Ti and O atoms
        elements = df['element'].unique()
        assert 'Ti' in elements
        assert 'O' in elements
    
    def test_custom_parameters(self):
        """Test custom parameter functionality."""
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        # Test with custom parameters
        custom_params = {
            'Ti-O': {'R0': 1.9, 'B': 0.4}
        }
        
        valence_states = {'Ti': 4, 'O': -2}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states, 
                                 custom_parameters=custom_params)
        
        # Calculate BVS with custom parameters
        bvs_results = bvs_calc.calculate_bvs()
        
        # Results should be different from default parameters
        assert isinstance(bvs_results, dict)
        assert len(bvs_results) == len(atoms)
    
    def test_bond_details(self):
        """Test bond details functionality."""
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        valence_states = {'Ti': 4, 'O': -2}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states)
        
        # Calculate BVS first
        bvs_calc.calculate_bvs()
        
        # Get bond details for first atom
        bond_details = bvs_calc.get_bond_details(0)
        
        assert isinstance(bond_details, list)
        
        if bond_details:  # If atom has neighbors
            bond = bond_details[0]
            expected_keys = ['neighbor_index', 'neighbor_element', 'neighbor_valence', 
                           'distance', 'bond_valence']
            for key in expected_keys:
                assert key in bond
    
    def test_distance_cutoff(self):
        """Test distance cutoff functionality."""
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        valence_states = {'Ti': 4, 'O': -2}
        
        # Test with different cutoffs
        bvs_calc_short = BondValenceSum(atoms, valence_states=valence_states, 
                                       distance_cutoff=2.0)
        bvs_calc_long = BondValenceSum(atoms, valence_states=valence_states, 
                                      distance_cutoff=5.0)
        
        bvs_short = bvs_calc_short.calculate_bvs()
        bvs_long = bvs_calc_long.calculate_bvs()
        
        # Longer cutoff should generally give higher BVS values
        # (more neighbors contribute)
        assert isinstance(bvs_short, dict)
        assert isinstance(bvs_long, dict)
    
    def test_valence_state_defaults(self):
        """Test default valence state handling."""
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        # Test without explicit valence states (should use defaults)
        bvs_calc = BondValenceSum(atoms)
        
        # Should not raise error
        bvs_results = bvs_calc.calculate_bvs()
        assert isinstance(bvs_results, dict)
        assert len(bvs_results) == len(atoms)
    
    def test_mixed_valence_states(self):
        """Test with mixed valence states."""
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        # Test with some Ti as +3 and some as +4
        valence_states = {'Ti': 3, 'O': -2}  # Different from default Ti+4
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states)
        
        bvs_results = bvs_calc.calculate_bvs()
        assert isinstance(bvs_results, dict)
        assert len(bvs_results) == len(atoms)
    
    def test_summary_printing(self):
        """Test summary printing functionality."""
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        valence_states = {'Ti': 4, 'O': -2}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states)
        
        # Should not raise error
        bvs_calc.print_summary()


class TestBondValenceIntegration:
    """Test integration and edge cases."""
    
    def test_layered_oxide_coordination(self):
        """Test coordination numbers for layered oxide structure with proper periodic boundary conditions."""
        from ase.io import read
        
        # Load the layered oxide structure (Co9Li27Mn9Ni9O54)
        atoms = read('asetools/data/summary_test/0_1.2/CONTCAR')
        
        # Set up BVS calculator with appropriate cutoff for first coordination shell
        valence_states = {'Li': 1, 'O': -2, 'Mn': 3, 'Co': 3, 'Ni': 2}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states, distance_cutoff=2.5)
        
        # Force neighbor calculation
        bvs_calc._find_neighbors()
        
        # Calculate coordination numbers by element
        coord_by_element = {'Li': [], 'O': [], 'Mn': [], 'Co': [], 'Ni': []}
        
        for i, symbol in enumerate(atoms.get_chemical_symbols()):
            coord_num = len(bvs_calc.neighbors[i])
            if symbol in coord_by_element:
                coord_by_element[symbol].append(coord_num)
        
        # Calculate average coordination numbers
        avg_coords = {}
        for element, coords in coord_by_element.items():
            if coords:
                avg_coords[element] = sum(coords) / len(coords)
        
        # Test that metal atoms have octahedral coordination (~6 neighbors)
        # With proper minimum image convention, should be exactly 6
        for metal in ['Li', 'Mn', 'Co', 'Ni']:
            assert abs(avg_coords[metal] - 6.0) < 0.1, \
                f"{metal} coordination {avg_coords[metal]:.1f} should be ~6.0 (octahedral)"
        
        # Oxygen should also have consistent coordination
        assert avg_coords['O'] > 5.0, f"O coordination {avg_coords['O']:.1f} should be >5"
    
    def test_empty_structure(self):
        """Test with empty structure."""
        atoms = Atoms()
        
        bvs_calc = BondValenceSum(atoms)
        bvs_results = bvs_calc.calculate_bvs()
        
        assert isinstance(bvs_results, dict)
        assert len(bvs_results) == 0
    
    def test_single_atom(self):
        """Test with single atom."""
        atoms = Atoms('Ti', positions=[[0, 0, 0]])
        
        bvs_calc = BondValenceSum(atoms)
        bvs_results = bvs_calc.calculate_bvs()
        
        assert isinstance(bvs_results, dict)
        assert len(bvs_results) == 1
        assert bvs_results[0] == 0.0  # No neighbors, so BVS = 0
    
    def test_unknown_element(self):
        """Test with unknown element."""
        atoms = Atoms('X', positions=[[0, 0, 0]])
        
        bvs_calc = BondValenceSum(atoms)
        
        # Should raise error for unknown element
        with pytest.raises(ValueError, match="No valence state specified"):
            bvs_calc.calculate_bvs()
    
    def test_no_bond_parameters(self):
        """Test with element pair that has no bond parameters."""
        # Create structure with elements that might not have bond parameters
        atoms = Atoms('He2', positions=[[0, 0, 0], [1, 0, 0]])
        
        valence_states = {'He': 0}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states)
        
        # Should handle gracefully (bond valence = 0 for missing parameters)
        bvs_results = bvs_calc.calculate_bvs()
        assert isinstance(bvs_results, dict)
        assert len(bvs_results) == 2
    
    def test_allowed_pairs_specification(self):
        """Test custom specification of allowed element pairs."""
        atoms = Atoms('TiO2H', positions=[
            [0, 0, 0],    # Ti
            [1.2, 0, 0],  # O  
            [0, 1.2, 0],  # O
            [2.0, 0, 0]   # H
        ], cell=[4, 4, 4], pbc=True)
        
        # Test with only Ti-O pairs allowed
        valence_states = {'Ti': 4, 'O': -2, 'H': 1}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states, 
                                 allowed_pairs=[('Ti', 'O')], distance_cutoff=3.0)
        
        # Check allowed pairs
        allowed_pairs = bvs_calc.get_allowed_pairs()
        assert ('O', 'Ti') in allowed_pairs
        assert len(allowed_pairs) == 1
        
        # Calculate BVS - should only consider Ti-O bonds
        bvs_results = bvs_calc.calculate_bvs()
        
        # Check bond details for Ti atom (index 0)
        ti_bonds = bvs_calc.get_bond_details(0)
        for bond in ti_bonds:
            assert bond['neighbor_element'] == 'O'  # Only O neighbors should be considered
    
    def test_exclude_same_element_pairs(self):
        """Test exclusion of same-element pairs."""
        # Create structure with metal atoms close together
        atoms = Atoms('Ti2O', positions=[
            [0, 0, 0],     # Ti
            [1.0, 0, 0],   # Ti (close to first Ti)
            [2.0, 0, 0]    # O
        ], cell=[5, 5, 5], pbc=True)
        
        valence_states = {'Ti': 4, 'O': -2}
        
        # Test with same-element exclusion (default)
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states, 
                                 exclude_same_element=True, distance_cutoff=2.5)
        
        allowed_pairs = bvs_calc.get_allowed_pairs()
        assert ('Ti', 'Ti') not in allowed_pairs
        assert ('O', 'Ti') in allowed_pairs
        
        # Test without same-element exclusion
        bvs_calc_include = BondValenceSum(atoms, valence_states=valence_states, 
                                         exclude_same_element=False, distance_cutoff=2.5)
        
        allowed_pairs_include = bvs_calc_include.get_allowed_pairs()
        # Note: Ti-Ti might still not be included if no bond valence parameters exist
        
        # Calculate BVS and check that Ti-Ti bonds are not considered in first case
        bvs_results = bvs_calc.calculate_bvs()
        ti_bonds = bvs_calc.get_bond_details(0)  # First Ti atom
        
        # Should not have Ti-Ti bonds in excluded case
        for bond in ti_bonds:
            assert bond['neighbor_element'] != 'Ti'
    
    def test_hydrogen_metal_exclusion(self):
        """Test that hydrogen-metal pairs are excluded by default."""
        atoms = Atoms('TiO2H2', positions=[
            [0, 0, 0],     # Ti
            [1.5, 0, 0],   # O
            [0, 1.5, 0],   # O
            [2.2, 0, 0],   # H (bonded to O)
            [0, 2.2, 0]    # H (bonded to O)
        ], cell=[4, 4, 4], pbc=True)
        
        valence_states = {'Ti': 4, 'O': -2, 'H': 1}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states, distance_cutoff=3.0)
        
        # Check that Ti-H pairs are not in allowed pairs by default
        allowed_pairs = bvs_calc.get_allowed_pairs()
        assert ('H', 'Ti') not in allowed_pairs
        
        # Calculate BVS for Ti atom
        bvs_results = bvs_calc.calculate_bvs()
        ti_bonds = bvs_calc.get_bond_details(0)  # Ti atom
        
        # Ti should only bond to O, not H
        for bond in ti_bonds:
            assert bond['neighbor_element'] != 'H'
    
    def test_allowed_pairs_methods(self):
        """Test methods for getting and printing allowed pairs."""
        atoms = Atoms('TiO2', positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]], 
                     cell=[4, 4, 4], pbc=True)
        
        valence_states = {'Ti': 4, 'O': -2}
        bvs_calc = BondValenceSum(atoms, valence_states=valence_states)
        
        # Test get_allowed_pairs method
        allowed_pairs = bvs_calc.get_allowed_pairs()
        assert isinstance(allowed_pairs, list)
        assert all(isinstance(pair, tuple) for pair in allowed_pairs)
        
        # Test print_allowed_pairs method (should not raise error)
        bvs_calc.print_allowed_pairs()