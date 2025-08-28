#!/usr/bin/env python
"""
Example usage of the ab initio thermodynamics module.

This script demonstrates both interpolation and lattice gas models
for calculating equilibrium adsorbate coverage and surface energies.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from asetools.ab_initio_thermodynamics import (
    ThermodynamicsCalculator, AdsorbateSpecies, SurfaceProperties
)

def example_single_adsorbate_interpolation():
    """
    Example 1: Single adsorbate (CO) using interpolation model.
    Similar to the Ru notebook approach.
    """
    print("=== Example 1: Single Adsorbate with Interpolation Model ===")
    
    # Create example data (normally you'd have this from DFT calculations)
    ru_co_data = pd.DataFrame({
        'Catalyst': ['Ru(fcc)', 'Ru(fcc)', 'Ru(fcc)', 'Ru(fcc)'] * 4,
        'facet': ['111'] * 4 + ['100'] * 4 + ['111'] * 4 + ['100'] * 4,
        'Ads': ['CO'] * 16,
        'cov': [0.0, 0.33, 0.67, 1.0] * 4,
        'Eads': [
            # Ru(fcc) 111 facet
            -1.03, -1.15, -1.35, -1.56,
            # Ru(fcc) 100 facet  
            -1.50, -1.60, -1.65, -1.56,
            # Repeat for other combinations
            -1.03, -1.15, -1.35, -1.56,
            -1.50, -1.60, -1.65, -1.56
        ]
    })
    
    # Initialize calculator and load data directly (no CSV file needed!)
    calc = ThermodynamicsCalculator('Ru(fcc)', ['CO'])
    calc.load_interpolation_model(ru_co_data)  # Pass DataFrame directly
    
    # Calculate equilibrium over temperature range
    results = calc.calculate_equilibrium(
        temperature_range=(273, 1000, 10),
        pressures=14000,  # 14 kPa CO pressure
        facets=['111', '100']
    )
    
    print(f"Calculated {len(results)} temperature points")
    print("Sample results:")
    print(results[['Temp', 'Cov_CO_111', 'Cov_CO_100', 'Gamma_111', 'Gamma_100']].head())
    
    return results


def example_multi_adsorbate_lattice_gas():
    """
    Example 2: Multi-adsorbate system using lattice gas model.
    Similar to the Pd notebook approach with CO + O competitive adsorption.
    """
    print("\n=== Example 2: Multi-Adsorbate with Lattice Gas Model ===")
    
    # Define interaction parameters from Pd notebook
    interaction_params = {
        'Pd': {
            '111': {
                'coordination': 6,
                'zero_coverage_energies': {'CO': -1.7341, 'O': -1.2138},
                'self_interactions': {'CO': -0.16808325, 'O': -0.179494083},
                'cross_interactions': {('CO', 'O'): -0.133321854}
            },
            '100': {
                'coordination': 4,
                'zero_coverage_energies': {'CO': -1.6660, 'O': -1.0635},
                'self_interactions': {'CO': -0.14940725, 'O': -0.19239175},
                'cross_interactions': {('CO', 'O'): -0.149512563}
            }
        }
    }
    
    # Initialize calculator for CO + O system
    calc = ThermodynamicsCalculator('Pd', ['CO', 'O'])
    calc.load_lattice_gas_model(interaction_params)
    
    # Define pressures (CO: 10 kPa, O2: 14 kPa)
    pressures = {'CO': 10000, 'O': 14000}
    
    # Calculate equilibrium
    results = calc.calculate_equilibrium(
        temperature_range=(373, 998, 25),
        pressures=pressures,
        facets=['111', '100']
    )
    
    print(f"Calculated {len(results)} temperature points")
    print("Sample results:")
    print(results[['Temp', 'Cov_CO_111', 'Cov_O_111', 'Cov_CO_100', 'Cov_O_100']].head())
    
    return results


def example_custom_adsorbate_species():
    """
    Example 3: Using custom adsorbate species with specific entropy parameters.
    """
    print("\n=== Example 3: Custom Adsorbate Species ===")
    
    # Create custom H species with specific entropy parameters
    h_species = AdsorbateSpecies(
        'H', 
        entropy_params={'a': 95.0, 'b': 0.16},  # Custom parameters
        dissociative=True  # H2 dissociation
    )
    
    print(f"H species dissociative: {h_species.dissociative}")
    
    # Test entropy calculation
    entropy_298 = h_species.entropy(298, 100000)
    entropy_500 = h_species.entropy(500, 100000)
    
    print(f"H entropy at 298 K: {entropy_298:.4f} kJ/(mol·K)")
    print(f"H entropy at 500 K: {entropy_500:.4f} kJ/(mol·K)")


def example_surface_properties():
    """
    Example 4: Working with surface properties database.
    """
    print("\n=== Example 4: Surface Properties ===")
    
    props = SurfaceProperties()
    
    print("Available metals:", props.get_available_metals())
    print("Pd facets:", props.get_available_facets('Pd'))
    print("Ru(fcc) facets:", props.get_available_facets('Ru(fcc)'))
    
    # Compare surface energies
    print("\nSurface energies (eV/Å²):")
    for metal in ['Pd', 'Ru(fcc)', 'Pt']:
        if metal in props.get_available_metals():
            facets = props.get_available_facets(metal)
            for facet in facets:
                gamma = props.get_surface_energy(metal, facet)
                area = props.get_surface_area(metal, facet)
                print(f"  {metal}({facet}): γ = {gamma:.4f} eV/Å², area = {area:.3f} Å²/atom")


def plot_results(results, title="Thermodynamics Results"):
    """
    Simple plotting function for results visualization.
    """
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot coverage vs temperature
        for col in results.columns:
            if 'Cov_' in col:
                ax1.plot(results['Temp'], results[col], '-o', markersize=3, label=col)
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('Coverage')
        ax1.set_title('Coverage vs Temperature')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot surface energy vs temperature
        for col in results.columns:
            if 'Gamma_' in col:
                ax2.plot(results['Temp'], results[col], '-o', markersize=3, label=col)
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('Surface Energy (eV/Å²)')
        ax2.set_title('Surface Energy vs Temperature')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle(title)
        plt.tight_layout()
        plt.show()
        
    except ImportError:
        print("Matplotlib not available, skipping plots")


def main():
    """Run all examples."""
    print("Ab Initio Thermodynamics Module Examples")
    print("=" * 50)
    
    # Example 1: Single adsorbate interpolation
    try:
        results1 = example_single_adsorbate_interpolation()
        plot_results(results1, "Single Adsorbate (CO) - Interpolation Model")
    except Exception as e:
        print(f"Example 1 failed: {e}")
    
    # Example 2: Multi-adsorbate lattice gas
    try:
        results2 = example_multi_adsorbate_lattice_gas()
        plot_results(results2, "Multi-Adsorbate (CO + O) - Lattice Gas Model")
    except Exception as e:
        print(f"Example 2 failed: {e}")
    
    # Example 3: Custom species
    try:
        example_custom_adsorbate_species()
    except Exception as e:
        print(f"Example 3 failed: {e}")
    
    # Example 4: Surface properties
    try:
        example_surface_properties()
    except Exception as e:
        print(f"Example 4 failed: {e}")


if __name__ == '__main__':
    main()