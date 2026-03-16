# API Reference

Python API for all ASEtools subpackages, documented from source code.

All subpackages must be imported explicitly. `import asetools` does not eagerly import everything.

---

## `asetools.analysis`

VASP output analysis: convergence, energy/force extraction, parameter parsing, metadata.

**Location:** `asetools/analysis/vasp.py`

**Import:**
```python
from asetools.analysis import (
    check_outcar_convergence,
    check_energy_and_maxforce,
    extract_magnetic_moments,
    get_parameter_from_run,
    classify_calculation_type,
    find_initial_structure,
    extract_comprehensive_metadata,
)
```

### Functions

#### `check_outcar_convergence(outcar, verbose=False)`

Check if a VASP run converged.

```python
def check_outcar_convergence(outcar: str, verbose: bool = False) -> Tuple[bool, str]
```

- **`outcar`** — path to OUTCAR file
- **`verbose`** — print convergence details to stdout
- **Returns** — `(converged: bool, vasp_version: str)` where `vasp_version` is `"vasp5"`, `"vasp6"`, or `""`

Handles optimization jobs (IBRION=1/2/3), single-point jobs (NSW≤1), and ASE optimizer jobs (IBRION=-1, NSW>0).

---

#### `check_energy_and_maxforce(outcar, magmom=False, verbose=False)`

Extract final energy and maximum force from OUTCAR.

```python
def check_energy_and_maxforce(
    outcar: str, magmom: bool = False, verbose: bool = False
) -> Union[Tuple[float, float], Tuple[float, float, float]]
```

- **`outcar`** — path to OUTCAR file
- **`magmom`** — if `True`, also return total magnetic moment
- **Returns** — `(energy, maxforce)` or `(energy, maxforce, magmom)` in eV / eV·Å⁻¹

---

#### `extract_magnetic_moments(outcar, listatoms, verbose=False)`

Get per-atom magnetic moments for specified atom indices.

```python
def extract_magnetic_moments(
    outcar: str, listatoms: List[int], verbose: bool = False
) -> List[float]
```

- **`listatoms`** — list of 0-based atom indices
- **Returns** — list of magnetic moments (rounded to 2 decimal places)

---

#### `get_parameter_from_run(outcar, check_converg=True, parameter="ISIF")`

Extract a single VASP parameter value from OUTCAR.

```python
def get_parameter_from_run(
    outcar: str, check_converg: bool = True, parameter: str = "ISIF"
) -> Tuple[Union[int, float, str, None], Union[bool, str]]
```

- **`parameter`** — INCAR parameter name to search for (e.g., `"ENCUT"`, `"EDIFFG"`, `"GGA"`)
- **`check_converg`** — if `True`, also check convergence before extracting
- **Returns** — `(value, convergence)` where value is int, float, str, or None

Raises `ValueError` if parameter not found and `check_converg=True`.

---

#### `classify_calculation_type(outcar_path, calc_dir)`

Classify a VASP calculation based on OUTCAR parameters and auxiliary files.

```python
def classify_calculation_type(outcar_path: str, calc_dir: str) -> str
```

- **Returns** — one of: `"dimer"`, `"neb"`, `"ase-optimizer"`, `"finite-diff"`, `"md"`, `"single-point"`, `"cell-relax"`, `"optimization"`, `"unknown"`

Classification priority: DIMER.log → ASE optimizer logs → IMAGES tag → NFREE > 0 → IBRION=0 → NSW≤1 → ISIF.

---

#### `find_initial_structure(calc_dir, pattern="*.vasp")`

Locate an initial structure file in a calculation directory.

```python
def find_initial_structure(calc_dir: str, pattern: str = "*.vasp") -> str
```

- **Returns** — absolute path to the initial structure file
- **Raises** `FileNotFoundError` — if no match for pattern and no POSCAR fallback found

Falls back to POSCAR if the glob pattern finds nothing.

---

#### `extract_comprehensive_metadata(outcar_path, incar_path=None, potcar_path=None)`

Extract comprehensive metadata from a VASP calculation.

```python
def extract_comprehensive_metadata(
    outcar_path: str,
    incar_path: Optional[str] = None,
    potcar_path: Optional[str] = None,
) -> Dict[str, Any]
```

- **Returns** — dictionary with keys:
  - `CalcType` (str), `Formula` (str)
  - `ENCUT`, `KSPACING`, `EDIFF`, `EDIFFG` (float)
  - `GGA` (str)
  - `IBRION`, `ISPIN`, `NSW`, `ISIF`, `NFREE` (int)
  - `Energy` (float, eV), `MaxForce` (float, eV/Å), `TotMagMom` (float)
  - `Converged` (bool), `VASPVersion` (str)
  - `INCAR_full` (str or None), `POTCAR_info` (list of TITEL strings)

---

### `asetools.analysis.symmetry`

Symmetry analysis using spglib.

**Location:** `asetools/analysis/symmetry.py`

**Import:** `from asetools.analysis.symmetry import SymmetryAnalyzer`

Requires optional dependency: `pip install -e ".[symmetry]"` (spglib).

---

## `asetools.electronic`

Density of States analysis with orbital projections and band center calculations.

**Location:** `asetools/electronic/doscar.py`

**Import:** `from asetools.electronic.doscar import DOS`

### `DOS` class

Parse a VASP DOSCAR file and expose analysis and plotting methods.

```python
class DOS:
    def __init__(self, doscarfile: str)
```

**Attributes after initialization:**
- `natoms` (int) — number of atoms
- `has_partial_dos` (bool) — True if DOSCAR contains PDOS (LORBIT >= 10)
- `fermi_energy` (float) — Fermi energy in eV
- `energy` (np.ndarray) — energy grid, referenced to Fermi energy
- `dos_up` (np.ndarray) — spin-up total DOS
- `dos_down` (np.ndarray) — spin-down total DOS (negative, mirror convention)
- `total_dos` (np.ndarray) — total DOS (sum of spin channels)
- `data` (dict) — raw parsed data; per-atom data under `data["at-N"]`

**Per-atom PDOS channels** (in `data["at-N"]`): `s+`, `s-`, `py+`, `py-`, `pz+`, `pz-`, `px+`, `px-`, `dxy+`, `dxy-`, `dyz+`, `dyz-`, `dz2+`, `dz2-`, `dxz+`, `dxz-`, `dx2+`, `dx2-`.

#### `get_pdos_by_orbitals(atoms, orbital)`

Get summed PDOS over selected atoms for a given orbital group.

```python
def get_pdos_by_orbitals(
    atoms: List[int], orbital: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]
```

- **`atoms`** — list of 0-based atom indices
- **`orbital`** — one of: `"all-s"`, `"all-p"`, `"all-d"`, `"t2g"`, `"eg"`, `"all"`
- **Returns** — `(energy, pdos_up, pdos_down)`

---

#### `plot_total_dos(ax=None)`

Plot total spin-up and spin-down DOS on axes.

```python
def plot_total_dos(self, ax=None)
```

---

#### `plot_multi_atom_pdos(atoms, orbitals, ax, same_color_spins=False, colors=None, linewidth=1.5)`

Plot summed PDOS for a list of atoms.

```python
def plot_multi_atom_pdos(
    self,
    atoms: List[int],
    orbitals: str,
    ax,
    same_color_spins: bool = False,
    colors=None,
    linewidth: float = 1.5,
)
```

---

#### `calculate_band_center(atoms, orbitals, energy_range=None, spin_treatment="combined")`

Calculate the band center (first moment of the DOS).

```python
def calculate_band_center(
    self,
    atoms: List[int],
    orbitals: str,
    energy_range: Optional[Tuple[float, float]] = None,
    spin_treatment: str = "combined",
) -> Union[float, Dict[str, float]]
```

- **`energy_range`** — `(min_eV, max_eV)` for integration; `None` uses full range
- **`spin_treatment`** — `"combined"` (default), `"separate"`, `"up"`, `"down"`
- **Returns** — float when `spin_treatment != "separate"`, dict `{"up": float, "down": float}` when `"separate"`

---

## `asetools.structure`

Surface analysis and structural validation tools.

### `SurfaceAnalyzer`

**Location:** `asetools/structure/adsorbate.py`

**Import:** `from asetools.structure.adsorbate import SurfaceAnalyzer`

Surface analysis and automated adsorbate placement.

```python
sa = SurfaceAnalyzer(atoms)
```

---

### `BondValenceSum`

**Location:** `asetools/structure/bond_valence.py`

**Import:** `from asetools.structure.bond_valence import BondValenceSum, BondValenceParameters`

Bond valence sum analysis using Brown's equation and the bvparm2020.cif parameter database (included in `asetools/data/`).

```python
from asetools.structure.bond_valence import BondValenceSum

bvs = BondValenceSum(atoms, valence_states={"Ni": 3, "O": -2})
df = bvs.analyze_structure()
```

---

## `asetools.electrochemistry`

Applied potential calculations with Fermi energy corrections.

**Location:** `asetools/electrochemistry/appliedpotential.py`

**Import:**
```python
from asetools.electrochemistry.appliedpotential import (
    extract_corrected_energy_fermie,
    fit_data,
    get_energy_at_givenpotential,
)
```

---

## `asetools.pathways`

Reaction pathway methods for transition state searches.

### NEB

**Location:** `asetools/pathways/neb.py`

NEB calculation support and energy extraction from image directories.

### Dimer

**Location:** `asetools/pathways/dimer.py`

Dimer method for saddle point searches. Requires `vasp-interactive`: `pip install -e ".[interactive]"`.

---

## `asetools.thermodynamics`

Ab initio thermodynamics for equilibrium coverage and surface energies.

**Location:** `asetools/thermodynamics/ab_initio.py`

**Import:**
```python
from asetools.thermodynamics.ab_initio import (
    ThermodynamicsCalculator,
    AdsorbateSpecies,
    SurfaceProperties,
    InterpolationModel,
    LatticeGasModel,
)
```

Key classes are also re-exported at the top level:
```python
from asetools import ThermodynamicsCalculator, AdsorbateSpecies, SurfaceProperties
```

| Class | Description |
|-------|-------------|
| `ThermodynamicsCalculator` | Main calculator for equilibrium coverage and surface energy |
| `AdsorbateSpecies` | Define adsorbate properties (name, entropy params, dissociative flag) |
| `SurfaceProperties` | Surface energy/area database backed by `surface_properties.json` |
| `InterpolationModel` | Coverage-dependent energies from CSV data |
| `LatticeGasModel` | Interaction-based lattice gas model |

See [ab_initio_thermodynamics.md](ab_initio_thermodynamics.md) for detailed documentation.

---

## `asetools.workflow`

YAML-based VASP workflow configuration and execution.

### `VASPConfigurationFromYAML`

**Location:** `asetools/workflow/calculatorsetuptools.py`

Load VASP calculator configuration from a YAML file.

```python
from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML

config = VASPConfigurationFromYAML("config.yaml", system="NCA")
```

```python
class VASPConfigurationFromYAML:
    def __init__(self, config_file: str, system: str = "default")
```

**Attributes:**
- `config` (dict) — full parsed YAML
- `basic_config` (dict) — contents of `basic:` key
- `workflows` (dict) — contents of `workflows:` key
- `globals` (dict) — contents of `globals:` key
- `initial_magmom_data` (dict) — magnetic moment data for this system

**Properties:**
- `system_config` → dict — configuration for the selected system (from `systems:` key); returns `{}` if not found

**Methods:**
- `initial_magmom()` → dict — returns `magmom` or `initial_magmom` from system config, or `{}` if absent

---

**Required YAML top-level keys:** `basic`, `systems`, `workflows`, `globals`

```yaml
basic:
  xc: PBE
  encut: 520
  ...

systems:
  NCA:
    magmom:
      Ni: 2.0
      Co: 3.0
      ...

workflows:
  bulk_opt:
    stages:
      - name: relax
        ...

globals:
  VASP_PP_PATH: /path/to/potentials
```

---

### Module-level functions (`calculatorsetuptools.py`)

#### `load_yaml_config(config_file)`

```python
def load_yaml_config(config_file: str) -> dict
```

Load and validate a YAML configuration file. Raises `KeyError` if required keys are missing.

---

#### `verify_configuration_keys(cfg)`

```python
def verify_configuration_keys(cfg: dict) -> None
```

Check that `cfg` contains `basic`, `systems`, `workflows`, `globals`. Raises `KeyError` if any are missing.

---

#### `deep_update(base, override)`

```python
def deep_update(base: dict, override: dict) -> dict
```

Recursively merge `override` into `base`. For dict-valued keys, recurses; for all others, replaces. Modifies `base` in place and returns it.

---

#### `setup_initial_magmom(atoms, magmom_data)`

```python
def setup_initial_magmom(atoms, magmom_data) -> atoms
```

Set initial magnetic moments on an ASE Atoms object.

- **`magmom_data`** accepts:
  - `dict` — element-based, e.g., `{"Ni": 2.0, "O": 0.0}`; atoms not in dict get 0.0
  - `list/tuple/np.ndarray` — per-atom values; padded with 0.0 if shorter than `len(atoms)`; raises `ValueError` if longer
  - `None` or empty — sets all to 0.0
- Raises `TypeError` for other types.
- **Returns** — modified `atoms` object

---

### `ConstraintManager`

**Location:** `asetools/workflow/constraints.py`

**Import:** `from asetools.workflow.constraints import ConstraintManager`

Manage ASE constraints (Hookean springs and FixAtoms) with JSON configuration support.

```python
class ConstraintManager:
    def __init__(self, distance_factor: float = 1.134)
```

`distance_factor` converts sum of covalent radii to equilibrium constraint distance (default 1.134 ≈ 13.4% buffer, calibrated for O-H bonds).

#### `load_constraint_config(json_file)`

```python
def load_constraint_config(self, json_file: str) -> Dict[str, Any]
```

Load JSON config. Raises `FileNotFoundError` or `json.JSONDecodeError`.

---

#### `calculate_bond_distance(atoms, idx1, idx2)`

```python
def calculate_bond_distance(self, atoms: Atoms, idx1: int, idx2: int) -> float
```

Compute equilibrium bond distance as `(r1 + r2) * distance_factor` using covalent radii.

---

#### `apply_hookean_from_pairs(atoms, pairs, k=20.0, distance_factor=None)`

```python
def apply_hookean_from_pairs(
    self,
    atoms: Atoms,
    pairs: List[Tuple[int, int]],
    k: float = 20.0,
    distance_factor: Optional[float] = None,
) -> List[Hookean]
```

Create Hookean constraints for a list of atom index pairs. Raises `ValueError` for out-of-range indices. Returns list of `Hookean` objects.

---

#### `get_existing_fix_indices(atoms)`

```python
def get_existing_fix_indices(self, atoms: Atoms) -> List[int]
```

Extract atom indices from existing `FixAtoms` constraints on the atoms object.

---

#### `merge_constraints(atoms, new_constraints)`

```python
def merge_constraints(self, atoms: Atoms, new_constraints: List) -> None
```

Merge new constraints with existing `FixAtoms` constraints. Rebuilds the constraint list preserving fixed atoms, then adds `new_constraints`. Modifies `atoms` in place. Raises `RuntimeError` if no constraints exist and none are provided.

---

#### `apply_from_json(atoms, json_file, k=20.0, distance_factor=None)`

```python
def apply_from_json(
    self, atoms: Atoms, json_file: str, k: float = 20.0,
    distance_factor: Optional[float] = None
) -> None
```

Main convenience method. Loads JSON, creates Hookean constraints, and merges with existing constraints. JSON must have a `"pairs"` key (list of `[idx1, idx2]`). Optional `"metadata"` key can contain `"spring_constant"` and `"distance_factor"`.

Raises `RuntimeError` if `"pairs"` key is missing or empty.

---

#### `apply_stage_constraints(atoms, constraint_config)`

```python
def apply_stage_constraints(self, atoms: Atoms, constraint_config: Dict[str, Any]) -> None
```

Apply constraints from a workflow stage configuration dict. Expected keys: `type` (must be `"hookean"`), `config_file` (required), `spring_constant` (optional), `distance_factor` (optional).

---

## `asetools.database`

ASE database integration with duplicate detection.

**Location:** `asetools/database/databases.py`

**Import:**
```python
from asetools.database.databases import check_if_exists_in_db, add_config_to_db, db_to_pandas
```

### Functions

#### `check_if_exists_in_db(db, atoms)`

```python
def check_if_exists_in_db(db: Any, atoms: Atoms) -> Tuple[bool, Optional[int]]
```

Check if a structure already exists in an ASE database by comparing forces.

- **Returns** — `(in_db: bool, row_id: Optional[int])`

---

#### `add_config_to_db(db, outcar, idname=None, update=False)`

```python
def add_config_to_db(
    db: Any, outcar: str, idname: Optional[str] = None, update: bool = False
) -> Optional[int]
```

Add a converged VASP configuration to an ASE database.

- **`idname`** — name to assign; defaults to the parent directory name of `outcar`
- **`update`** — if `True` and already in DB, update the existing entry
- Skips non-converged calculations (calls `check_outcar_convergence`)
- **Returns** — database row ID or None

---

#### `db_to_pandas(db, columns=None)`

```python
def db_to_pandas(db: Any, columns: Optional[List[str]] = None) -> pd.DataFrame
```

Convert an ASE database to a pandas DataFrame.

- **`columns`** — list of column names to extract; defaults to `["name", "id", "energy", "free_energy", "magmom"]`
- **Returns** — `pd.DataFrame`

---

## `asetools.plotting`

Energy profile and potential energy surface (PES) plotting.

**Location:** `asetools/plotting/plots.py`

**Import:**
```python
from asetools.plotting.plots import add_line_to_pes, beautify_pes_plot, validate_columns
```

Default column names expected in input DataFrames: `Label`, `Type-Conf`, `E`, `nPCET`.

### Functions

#### `add_line_to_pes(ax, data, energy_col='E', type_col='Type-Conf', c='k', label=None, indexes=None, col=None, style='-', lw=2, lw_connector=0.5)`

Add a line (energy profile) to a PES axes. Returns `ax`.

---

#### `beautify_pes_plot(ax, xlim=None, ylim=None, zero=True, leg=True, fs=12, data=None, label_col='Label', type_col='Type-Conf', npcet_col='nPCET', show_labels=True, show_npcet=False, indexes=None, frame=True, y_decimals=2)`

Apply formatting and annotations to a PES plot. Returns `ax`.

---

#### `validate_columns(data, required_cols, optional_cols=[])`

```python
def validate_columns(data, required_cols, optional_cols=[]) -> bool
```

Check that a DataFrame contains the required and optional columns.

---

## `asetools.parsers`

Low-level VASP OUTCAR parsing with an extensible class hierarchy.

**Location:** `asetools/parsers/vasp_outcar.py`

**Import:**
```python
from asetools.parsers.vasp_outcar import (
    outcarchunks,
    OutcarChunkParser,
    OutcarHeaderParser,
    DefaultParsersContainer,
)
```

### Class Hierarchy

```
VaspPropertyParser (ABC)
├── SimpleProperty
├── VaspChunkPropertyParser
│   └── (concrete: Stress, Cell, PositionsAndForces, Magmom, Magmoms, EFermi, Energy, Kpoints)
└── VaspHeaderPropertyParser
    └── (concrete: Spinpol, SpeciesTypes, IonsPerSpecies, KpointHeader)
```

### Top-level Parsers

- **`OutcarChunkParser`** — parses individual ionic step chunks
- **`OutcarHeaderParser`** — parses the OUTCAR header section
- **`DefaultParsersContainer`** — pre-configured container with all standard parsers

### Iterator

#### `outcarchunks(fd, chunk_parser, header_parser)`

```python
def outcarchunks(fd, chunk_parser, header_parser) -> Iterator[OUTCARChunk]
```

Iterate over ionic step chunks in an OUTCAR file.

### Concrete Parsers

| Class | Extracts |
|-------|---------|
| `Spinpol` | ISPIN value |
| `SpeciesTypes` | Element species names |
| `IonsPerSpecies` | Number of ions per species |
| `KpointHeader` | k-point header information |
| `Stress` | Stress tensor |
| `Cell` | Unit cell vectors |
| `PositionsAndForces` | Atomic positions and forces |
| `Magmom` | Total magnetic moment |
| `Magmoms` | Per-atom magnetic moments |
| `EFermi` | Fermi energy |
| `Energy` | Total energy |
| `Kpoints` | k-point list and weights |
