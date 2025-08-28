import unittest
import numpy as np
import pandas as pd
import tempfile
import os
from asetools.ab_initio_thermodynamics import (
    AdsorbateSpecies, SurfaceProperties, InterpolationModel, 
    LatticeGasModel, ThermodynamicsCalculator
)


class TestAdsorbateSpecies(unittest.TestCase):
    
    def test_default_co_species(self):
        co = AdsorbateSpecies('CO')
        self.assertEqual(co.name, 'CO')
        self.assertFalse(co.dissociative)
        
        # Test entropy calculation
        entropy = co.entropy(298, 100000)  # 298 K, 1 bar
        self.assertIsInstance(entropy, float)
        self.assertGreater(entropy, 0)
    
    def test_default_o_species(self):
        o = AdsorbateSpecies('O')
        self.assertEqual(o.name, 'O')
        self.assertTrue(o.dissociative)
    
    def test_custom_species(self):
        custom = AdsorbateSpecies('H', entropy_params={'a': 100, 'b': 0.15}, dissociative=True)
        self.assertTrue(custom.dissociative)
        self.assertEqual(custom.entropy_params['a'], 100)


class TestSurfaceProperties(unittest.TestCase):
    
    def test_default_properties(self):
        props = SurfaceProperties()
        
        # Test available metals
        metals = props.get_available_metals()
        self.assertIn('Pd', metals)
        self.assertIn('Ru(fcc)', metals)
        
        # Test surface energy retrieval
        gamma_pd111 = props.get_surface_energy('Pd', '111')
        self.assertEqual(gamma_pd111, 0.12500)
        
        # Test surface area retrieval
        area_pd111 = props.get_surface_area('Pd', '111')
        self.assertEqual(area_pd111, 6.8589)
        
        # Test facet availability
        pd_facets = props.get_available_facets('Pd')
        self.assertIn('111', pd_facets)
        self.assertIn('100', pd_facets)


class TestInterpolationModel(unittest.TestCase):
    
    def setUp(self):
        # Create a temporary CSV file for testing
        self.test_data = pd.DataFrame({
            'Catalyst': ['Ru(fcc)', 'Ru(fcc)', 'Ru(fcc)'],
            'facet': ['111', '111', '111'],
            'Ads': ['CO', 'CO', 'CO'],
            'cov': [0.0, 0.5, 1.0],
            'Eads': [-1.0, -1.2, -1.5]
        })
        
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        self.test_data.to_csv(self.temp_file.name, index=False)
        self.temp_file.close()
    
    def tearDown(self):
        os.unlink(self.temp_file.name)
    
    def test_load_and_interpolate(self):
        model = InterpolationModel(self.temp_file.name)
        
        # Test exact values
        eads_0 = model.get_adsorption_energy('Ru(fcc)', '111', 'CO', 0.0)
        self.assertEqual(eads_0, -1.0)
        
        eads_1 = model.get_adsorption_energy('Ru(fcc)', '111', 'CO', 1.0)
        self.assertEqual(eads_1, -1.5)
        
        # Test interpolation
        eads_mid = model.get_adsorption_energy('Ru(fcc)', '111', 'CO', 0.25)
        self.assertAlmostEqual(eads_mid, -1.1, places=6)


class TestLatticeGasModel(unittest.TestCase):
    
    def setUp(self):
        # Example parameters from Pd notebook
        self.interaction_params = {
            'Pd': {
                '111': {
                    'coordination': 6,
                    'zero_coverage_energies': {'CO': -1.7341, 'O': -1.2138},
                    'self_interactions': {'CO': -0.16808325, 'O': -0.179494083},
                    'cross_interactions': {('CO', 'O'): -0.133321854}
                }
            }
        }
    
    def test_lattice_gas_calculation(self):
        model = LatticeGasModel(self.interaction_params)
        
        # Test CO adsorption energy with some coverage
        coverages = {'CO': 0.5, 'O': 0.2}
        eads_co = model.get_adsorption_energy('Pd', '111', 'CO', coverages)
        
        # Should be different from zero-coverage value due to interactions
        self.assertNotEqual(eads_co, -1.7341)
        self.assertIsInstance(eads_co, float)
        
        # Test O adsorption (dissociative, should be multiplied by 2)
        eads_o = model.get_adsorption_energy('Pd', '111', 'O', coverages)
        self.assertIsInstance(eads_o, float)


class TestThermodynamicsCalculator(unittest.TestCase):
    
    def setUp(self):
        # Create test CSV data
        self.test_data = pd.DataFrame({
            'Catalyst': ['Ru(fcc)'] * 3,
            'facet': ['111'] * 3,
            'Ads': ['CO'] * 3,
            'cov': [0.0, 0.5, 1.0],
            'Eads': [-1.0, -1.2, -1.5]
        })
        
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        self.test_data.to_csv(self.temp_file.name, index=False)
        self.temp_file.close()
        
        self.calc = ThermodynamicsCalculator('Ru(fcc)', ['CO'])
    
    def tearDown(self):
        os.unlink(self.temp_file.name)
    
    def test_initialization(self):
        self.assertEqual(self.calc.metal, 'Ru(fcc)')
        self.assertEqual(self.calc.adsorbates, ['CO'])
        self.assertIn('CO', self.calc.species)
    
    def test_load_interpolation_model(self):
        self.calc.load_interpolation_model(self.temp_file.name)
        self.assertEqual(self.calc.model_type, 'interpolation')
        self.assertIsNotNone(self.calc.model)
    
    def test_equilibrium_coverage_calculation(self):
        self.calc.load_interpolation_model(self.temp_file.name)
        
        # Test at high temperature (should give low coverage)
        coverages_high_T = self.calc.calculate_equilibrium_coverage('111', 1000, 14000)
        self.assertIn('CO', coverages_high_T)
        self.assertLessEqual(coverages_high_T['CO'], 1.0)
        self.assertGreaterEqual(coverages_high_T['CO'], 0.0)
        
        # Test at low temperature (should give higher coverage)
        coverages_low_T = self.calc.calculate_equilibrium_coverage('111', 300, 14000)
        self.assertGreaterEqual(coverages_low_T['CO'], coverages_high_T['CO'])
    
    def test_surface_energy_calculation(self):
        self.calc.load_interpolation_model(self.temp_file.name)
        
        coverages = {'CO': 0.5}
        gamma = self.calc.calculate_surface_energy('111', coverages, 500, 14000)
        
        # Should be different from clean surface energy
        clean_gamma = self.calc.surface_props.get_surface_energy('Ru(fcc)', '111')
        self.assertNotEqual(gamma, clean_gamma)
        self.assertIsInstance(gamma, float)
    
    def test_full_equilibrium_calculation(self):
        self.calc.load_interpolation_model(self.temp_file.name)
        
        # Test temperature range calculation
        results = self.calc.calculate_equilibrium(
            temperature_range=(400, 500, 50),
            pressures=14000,
            facets=['111']
        )
        
        self.assertIsInstance(results, pd.DataFrame)
        self.assertIn('Temp', results.columns)
        self.assertIn('Cov_CO_111', results.columns)
        self.assertIn('Gamma_111', results.columns)
        self.assertEqual(len(results), 2)  # 400, 450


class TestMultiAdsorbateSystem(unittest.TestCase):
    
    def setUp(self):
        self.interaction_params = {
            'Pd': {
                '111': {
                    'coordination': 6,
                    'zero_coverage_energies': {'CO': -1.7341, 'O': -1.2138},
                    'self_interactions': {'CO': -0.16808325, 'O': -0.179494083},
                    'cross_interactions': {('CO', 'O'): -0.133321854}
                }
            }
        }
        
        self.calc = ThermodynamicsCalculator('Pd', ['CO', 'O'])
    
    def test_multi_adsorbate_equilibrium(self):
        self.calc.load_lattice_gas_model(self.interaction_params)
        
        pressures = {'CO': 10000, 'O': 14000}
        
        # Test equilibrium calculation
        coverages = self.calc.calculate_equilibrium_coverage('111', 500, pressures)
        
        self.assertIn('CO', coverages)
        self.assertIn('O', coverages)
        
        # Total coverage should be ≤ 1
        total_coverage = sum(coverages.values())
        self.assertLessEqual(total_coverage, 1.0)
        self.assertGreaterEqual(total_coverage, 0.0)
    
    def test_competitive_adsorption_dataframe(self):
        self.calc.load_lattice_gas_model(self.interaction_params)
        
        pressures = {'CO': 10000, 'O': 14000}
        
        results = self.calc.calculate_equilibrium(
            temperature_range=(400, 600, 100),
            pressures=pressures,
            facets=['111']
        )
        
        self.assertIsInstance(results, pd.DataFrame)
        self.assertIn('Cov_CO_111', results.columns)
        self.assertIn('Cov_O_111', results.columns)
        self.assertIn('Eads_CO_111', results.columns)
        self.assertIn('Eads_O_111', results.columns)


if __name__ == '__main__':
    unittest.main()