# API Reference

Module-level Python API for ASEtools.

All subpackages are lazily loaded — `import asetools` does not eagerly import everything. Access subpackages via `asetools.<subpackage>` or import directly.

---

## `asetools.analysis`

VASP output analysis: convergence checking, energy/force extraction, parameter parsing.

**Location:** `asetools/analysis/vasp.py`

### Key Functions

| Function | Description |
|----------|-------------|
| `check_outcar_convergence(outcar, verbose=False)` | Check if VASP run converged. Returns `(bool, str)` — convergence status and VASP version. |
| `check_energy_and_maxforce(outcar, magmom=False, verbose=False)` | Extract final energy (eV) and max force (eV/Å). Optionally returns magnetic moments. |
| `extract_magnetic_moments(outcar, atoms_list)` | Get magnetic moments for specified atom indices. |
| `get_parameter_from_run(outcar, parameter='ISIF', check_converg=True)` | Extract a specific INCAR parameter value from OUTCAR. |
| `extract_comprehensive_metadata(outcar)` | Extract full metadata dict (structure, parameters, POTCAR, timing). |
| `find_initial_structure(path, pattern='POSCAR*')` | Locate initial structure file in a calculation directory. |

### Symmetry Submodule

**Location:** `asetools/analysis/symmetry.py`

Symmetry analysis using spglib.

---

## `asetools.electronic`

Density of States analysis with orbital projections and band center calculations.

**Location:** `asetools/electronic/doscar.py`

### `DOS` Class (recommended)

Object-oriented interface for DOSCAR analysis.

```python
from asetools.electronic.doscar import DOS

dos = DOS('DOSCAR')
```

| Method | Description |
|--------|-------------|
| `DOS(doscar_file)` | Parse DOSCAR file |
| `calculate_band_center(atoms, orbitals='all-d', spin_treatment='sum')` | d-band/p-band center. `spin_treatment`: `'sum'`, `'separate'`, `'up'`, `'down'` |
| `get_fermi_energy()` | Fermi energy in eV |
| `plot_dos(...)` | Plot total or projected DOS |
| `get_pdos(atom, orbital)` | Get projected DOS for specific atom/orbital |

### Legacy Functions

| Function | Description |
|----------|-------------|
| `extract_dos(doscarfile)` | Parse complete DOSCAR into dict |
| `extract_pdos_perstate(data, atoms, states)` | Project DOS onto s/p/d state groups |
| `extract_pdos_perorbital(data, atoms, orbitals)` | Detailed orbital projections (t2g, eg, etc.) |
| `extract_fermi_e(doscarfile)` | Extract Fermi energy |
| `calculate_band_center(doscarfile, atoms, orbitals)` | Band center calculation (functional API) |

### Orbital Specifications

| Orbital | Components |
|---------|------------|
| `'s'` | s orbital |
| `'p_states'` | px + py + pz |
| `'d_states'` | all 5 d orbitals |
| `'all-d'` | all d (same as `d_states`) |
| `'t2g'` | dxy + dxz + dyz |
| `'eg'` | dz² + dx²-y² |

---

## `asetools.structure`

Surface science and structural validation tools.

### `SurfaceAnalyzer`

**Location:** `asetools/structure/adsorbate.py`

Surface analysis and automated adsorbate placement.

```python
from asetools.structure.adsorbate import SurfaceAnalyzer

analyzer = SurfaceAnalyzer(atoms)
sites = analyzer.find_surface_neighbors()
analyzer.add_adsorbate_to_mid(sites[0], adsorbate='H', z_off=1.5)
```

| Method | Description |
|--------|-------------|
| `SurfaceAnalyzer(atoms)` | Initialize with ASE Atoms object |
| `find_surface_neighbors()` | Identify surface sites and neighbors |
| `add_adsorbate_to_mid(site, adsorbate, z_off)` | Place adsorbate at a surface site |
| `adsneighdistances` | Distances to neighboring surface atoms |

### `BondValenceSum`

**Location:** `asetools/structure/bond_valence.py`

Bond valence sum analysis using Brown's equation and the 2020 parameter database.

```python
from asetools.structure.bond_valence import BondValenceSum

bvs = BondValenceSum(atoms, valence_states={'Ti': 4, 'O': -2})
df = bvs.analyze_structure()
```

| Method / Parameter | Description |
|-------------------|-------------|
| `BondValenceSum(atoms, valence_states, ...)` | Initialize calculator |
| `auto_determine_valence=True` | Auto-optimize metal valence states |
| `per_atom_valence=True` | Allow mixed valences per element |
| `allowed_pairs` | Restrict to specific element pairs |
| `exclude_same_element=True` | Skip same-element bonds |
| `analyze_structure()` | Full analysis → DataFrame |
| `print_valence_optimization_summary()` | Print optimization details |

---

## `asetools.electrochemistry`

Applied potential calculations with Fermi shift corrections.

**Location:** `asetools/electrochemistry/appliedpotential.py`

```python
from asetools.electrochemistry.appliedpotential import (
    extract_corrected_energy_fermie,
    fit_data,
    get_energy_at_givenpotential,
)
```

| Function | Description |
|----------|-------------|
| `extract_corrected_energy_fermie(folders, calc_zero)` | Extract potential-dependent energies with Fermi corrections |
| `fit_data(X, Y, fit_type='polynomial', order=2)` | Fit energy vs. potential data |
| `get_energy_at_givenpotential(results, desiredU=0.)` | Interpolate energy at a specific potential |

---

## `asetools.pathways`

Reaction pathway methods for transition state searches.

### NEB

**Location:** `asetools/pathways/neb.py`

| Function | Description |
|----------|-------------|
| `extract_neb_data(folder_path, final=False)` | Extract NEB energies from image directories |

### Dimer

**Location:** `asetools/pathways/dimer.py`

Dimer method implementation for saddle point searches. Requires POSCAR in working directory.

---

## `asetools.thermodynamics`

Ab initio thermodynamics for equilibrium coverage and surface energies.

**Location:** `asetools/thermodynamics/ab_initio.py`

```python
from asetools.thermodynamics.ab_initio import (
    ThermodynamicsCalculator,
    AdsorbateSpecies,
    SurfaceProperties,
    InterpolationModel,
    LatticeGasModel,
)
```

These classes are also accessible at top level: `from asetools import ThermodynamicsCalculator`.

| Class | Description |
|-------|-------------|
| `ThermodynamicsCalculator(metal, adsorbates)` | Main calculator for coverage and surface energy |
| `AdsorbateSpecies(name, entropy_params, dissociative)` | Define adsorbate properties |
| `SurfaceProperties(surface_data_file)` | Surface energy/area database |
| `InterpolationModel` | Coverage-dependent energies from CSV data |
| `LatticeGasModel` | Interaction-based lattice gas model |

See [ab_initio_thermodynamics.md](ab_initio_thermodynamics.md) for detailed documentation.

---

## `asetools.workflow`

YAML-based multi-stage workflow runner for VASP calculations.

**Location:** `asetools/workflow/calculatorsetuptools.py`

```python
from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML
```

| Class/Function | Description |
|---------------|-------------|
| `VASPConfigurationFromYAML(config_file, system=None)` | Load VASP config from YAML |
| `run_workflow(atoms, config, workflow_name)` | Execute multi-stage workflow |
| `make_calculator(config, overrides=None)` | Create VASP calculator from config |

### YAML Structure

```yaml
basic_config:
  xc: PBE
  encut: 400
  kpts: [4, 4, 4]

workflows:
  optimization:
    stages:
      - name: "rough_opt"
        ediff: 1e-4
        nsw: 100
        constraints:          # optional Hookean constraints
          type: hookean
          config_file: pairs.json
          spring_constant: 20.0
        steps:
          - name: optimization
            optimizer: FIRE
            optimizer_kwargs: { fmax: 0.02 }
```

See [constraints_quickstart.md](constraints_quickstart.md) for constraint configuration.

---

## `asetools.database`

ASE database integration with duplicate detection.

**Location:** `asetools/database/databases.py`

```python
from asetools.database.databases import add_config_to_db, check_if_exists_in_db, db_to_pandas
```

| Function | Description |
|----------|-------------|
| `add_config_to_db(db, outcar, idname=None)` | Add structure to ASE database |
| `check_if_exists_in_db(db, atoms)` | Check for duplicate structures |
| `db_to_pandas(db)` | Convert database to pandas DataFrame |

---

## `asetools.plotting`

Energy profile and potential energy surface plotting.

**Location:** `asetools/plotting/plots.py`

| Function | Description |
|----------|-------------|
| `plot_nebs(list_dfs, ...)` | Plot multiple NEB pathways |
| `add_line_to_pes(ax, data)` | Add a line to a PES diagram |

---

## `asetools.parsers`

Low-level VASP output parsing.

**Location:** `asetools/parsers/vasp_outcar.py`

Direct OUTCAR parsing utilities used by higher-level analysis functions.

---

## Top-Level Imports

For convenience, key thermodynamics classes are re-exported at the package level:

```python
from asetools import ThermodynamicsCalculator, AdsorbateSpecies, SurfaceProperties
```

All other imports use the full subpackage path:

```python
from asetools.analysis import check_outcar_convergence
from asetools.electronic.doscar import DOS
from asetools.structure.bond_valence import BondValenceSum
```
