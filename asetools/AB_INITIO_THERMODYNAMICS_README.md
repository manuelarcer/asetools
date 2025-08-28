# Ab Initio Thermodynamics Module

This module provides functionality for calculating equilibrium adsorbate coverage and surface energies under varying temperature and pressure conditions using ab initio thermodynamics principles.

## Features

- **Single and Multi-Adsorbate Systems**: Support for both simple single-adsorbate and competitive multi-adsorbate adsorption
- **Two Calculation Models**:
  - **Interpolation Model**: Simple coverage-dependent adsorption energies from CSV data
  - **Lattice Gas Model**: Interaction-based coverage dependence with adsorbate-adsorbate interactions
- **Flexible Species Definition**: Customizable adsorbate properties (entropy parameters, dissociation behavior)
- **Surface Property Database**: Built-in database for common metals and facets
- **Robust Optimization**: Multiple initial guesses and convergence checking for equilibrium calculations

## Based on Literature

The module implements approaches from:
- Zhu, Meng, and Gao, *J. Phys. Chem. C* **121** (2017) 5629-5634
- Standard Langmuir isotherm and competitive adsorption theory

## Quick Start

### Single Adsorbate with Interpolation Model

```python
from asetools import ThermodynamicsCalculator
import pandas as pd

# Create adsorption energy data
data = pd.DataFrame({
    'Catalyst': ['Ru(fcc)'] * 4,
    'facet': ['111'] * 4,
    'Ads': ['CO'] * 4,
    'cov': [0.0, 0.33, 0.67, 1.0],
    'Eads': [-1.03, -1.15, -1.35, -1.56]
})
data.to_csv('co_ads_data.csv', index=False)

# Initialize calculator
calc = ThermodynamicsCalculator('Ru(fcc)', ['CO'])
calc.load_interpolation_model('co_ads_data.csv')

# Calculate equilibrium
results = calc.calculate_equilibrium(
    temperature_range=(273, 1000, 50),
    pressures=14000,  # 14 kPa CO
    facets=['111']
)

print(results[['Temp', 'Cov_CO_111', 'Gamma_111']].head())
```

### Multi-Adsorbate with Lattice Gas Model

```python
from asetools import ThermodynamicsCalculator

# Define interaction parameters
interaction_params = {
    'Pd': {
        '111': {
            'coordination': 6,
            'zero_coverage_energies': {'CO': -1.7341, 'O': -1.2138},
            'self_interactions': {'CO': -0.16808325, 'O': -0.179494083},
            'cross_interactions': {('CO', 'O'): -0.133321854}
        }
    }
}

# Initialize multi-adsorbate calculator
calc = ThermodynamicsCalculator('Pd', ['CO', 'O'])
calc.load_lattice_gas_model(interaction_params)

# Calculate competitive adsorption
results = calc.calculate_equilibrium(
    temperature_range=(373, 1000, 50),
    pressures={'CO': 10000, 'O': 14000},  # kPa
    facets=['111']
)

print(results[['Temp', 'Cov_CO_111', 'Cov_O_111', 'Gamma_111']].head())
```

## Data Formats

### CSV Input for Interpolation Model

Required columns: `['Catalyst', 'facet', 'Ads', 'cov', 'Eads']`

```csv
Catalyst,facet,Ads,cov,Eads
Ru(fcc),111,CO,0.0,-1.03
Ru(fcc),111,CO,0.33,-1.15
Ru(fcc),111,CO,0.67,-1.35
Ru(fcc),111,CO,1.0,-1.56
```

### Lattice Gas Interaction Parameters

```python
interaction_params = {
    'Metal': {
        'facet': {
            'coordination': z_value,                    # Coordination number
            'zero_coverage_energies': {                 # E_0 values in eV
                'species': energy_eV
            },
            'self_interactions': {                      # w_ii parameters in eV
                'species': w_value_eV
            },
            'cross_interactions': {                     # w_ij parameters in eV
                ('species1', 'species2'): w_value_eV
            }
        }
    }
}
```

## Built-in Surface Database

The module includes surface energies and areas for:
- **Ru(fcc)**: facets 100, 111
- **Ru(hcp)**: facets 100, 101  
- **Pd**: facets 111, 100, 110, 120, 311
- **Pt**: facets 111, 100, 101, 331, 311

Custom surface properties can be provided via JSON file:

```json
{
  "surface_energies": {
    "Metal": {"facet": energy_eV_per_A2}
  },
  "surface_areas": {
    "Metal": {"facet": area_A2_per_atom}
  }
}
```

## Adsorbate Species

### Default Species
- **CO**: Molecular adsorption, entropy parameters from literature
- **O**: Dissociative from O₂, entropy parameters from literature
- **H**: Dissociative from H₂ (placeholder parameters)
- **N**: Dissociative from N₂ (placeholder parameters)

### Custom Species

```python
from asetools import AdsorbateSpecies

# Custom entropy parameters
custom_species = AdsorbateSpecies(
    'H',
    entropy_params={'a': 95.0, 'b': 0.16},
    dissociative=True
)
```

## Output Format

Results are returned as pandas DataFrame with columns:
- `Temp`: Temperature in K
- `Cov_{species}_{facet}`: Coverage for each species/facet combination
- `Eads_{species}_{facet}`: Adsorption energy in eV for each species/facet
- `Gamma_{facet}`: Surface energy in eV/Å² for each facet

## API Reference

### ThermodynamicsCalculator

Main calculation class:

- `__init__(metal, adsorbates, surface_properties=None)`
- `load_interpolation_model(csv_file)`
- `load_lattice_gas_model(interaction_params)`
- `calculate_equilibrium_coverage(facet, temperature, pressures)`
- `calculate_surface_energy(facet, coverages, temperature, pressures)`
- `calculate_equilibrium(temperature_range, pressures, facets=None)`

### AdsorbateSpecies

Handles adsorbate properties:
- `__init__(name, entropy_params=None, dissociative=None)`
- `entropy(temperature, pressure)`

### SurfaceProperties

Manages surface energy database:
- `__init__(surface_data_file=None)`
- `get_surface_energy(metal, facet)`
- `get_surface_area(metal, facet)`
- `get_available_metals()`
- `get_available_facets(metal)`

## Examples

See `examples/ab_initio_thermodynamics_example.py` for complete working examples.

## Testing

Run tests with:
```bash
python -m pytest asetools/tests/test_ab_initio_thermodynamics.py -v
```

## Notes

- Temperatures should be provided in Kelvin
- Pressures should be provided in Pascals
- Adsorption energies should be in eV (negative values for favorable adsorption)
- Surface energies are in eV/Å²
- The module handles both molecular (CO) and dissociative (O₂→2O) adsorption automatically
- Convergence is achieved through multiple random initial guesses with SLSQP optimization