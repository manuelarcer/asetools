# CLI Reference

ASEtools ships 17 command-line tools, all registered as console scripts via `pyproject.toml`. After `pip install -e .`, they are available on `PATH`.

---

## Analysis Tools

### `getenergy`

Quick convergence check and energy extraction from a single OUTCAR.

```bash
getenergy [OUTCAR]
```

**Output:** Convergence status, max force (eV/Å), and total energy (eV).

**Example:**
```
Converged   MaxForce      Energy
True        0.015     -234.5678
```

---

### `summaryfolders`

Batch analysis across multiple calculation directories. Scans subdirectories for OUTCAR files, extracts convergence, energy, forces, and optionally magnetic moments.

```bash
summaryfolders [-m] [-f] [-p] [--folders dir1 dir2 ...]
```

| Flag | Description |
|------|-------------|
| `-m`, `--magmom` | Extract and display magnetic moments |
| `-f`, `--fast` | Fast mode (reduced parsing) |
| `-p`, `--pyatoms` | Pyatoms multi-step mode (reads step_* folders) |
| `--folders` | Specific directories to process |

**Output:** Formatted table (pandas DataFrame) with results per calculation.

---

### `comparevaspparam`

Side-by-side comparison of VASP parameters between two OUTCAR files.

```bash
comparevaspparam OUTCAR1 OUTCAR2
```

Compares INCAR parameters, k-points, basis set sizes, and other runtime settings. Highlights differences.

---

### `inputparam-from-outcar`

Extract INCAR parameters from an OUTCAR file, organized by category (electronic, ionic, DOS, etc.).

```bash
inputparam-from-outcar OUTCAR
```

Categories: StartParam, Electronic, IonicRelax, System, DOS, ElectronicRelax, Write, DipoleCorr, ExCorr, LinearResp.

---

### `outcar-extract`

Extract a specific ionic step frame from an OUTCAR and write it as a structure file.

```bash
outcar-extract OUTCAR [index] [output_file] [format]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `index` | `-1` (last) | Frame index to extract |
| `output_file` | `final.vasp` | Output file path |
| `format` | auto-detect | ASE format string (e.g., `extxyz`) |

**Examples:**
```bash
outcar-extract OUTCAR                        # last frame → final.vasp
outcar-extract OUTCAR 5                      # 5th frame → final.vasp
outcar-extract OUTCAR 5 freq/CONTCAR.vasp    # 5th frame → specific file
outcar-extract OUTCAR -1 final.extxyz extxyz # with explicit format
```

---

## Electronic Structure

### `plotdos`

Interactive DOS/PDOS plotting from VASP DOSCAR files. Supports total DOS, partial DOS with orbital projections, band center markers, and spin-resolved plots.

```bash
plotdos [DOSCAR] [options]
```

| Flag | Description |
|------|-------------|
| `--atoms SPEC` | Atom indices (`0,2,5` or `0-5` or `all`) |
| `--orbitals TYPE` | Orbital type: `s`, `p`, `d`, `all-d`, `t2g`, `eg`, etc. |
| `--states TYPE` | State groups: `s_states`, `p_states`, `d_states` |
| `--range MIN,MAX` | Energy range for plot (eV) |
| `--band-center` | Show band center marker |
| `--spin separate` | Spin-resolved plotting |
| `--save FILE` | Save plot to file instead of displaying |
| `--fermi-shift` | Shift energies relative to Fermi level |

**Examples:**
```bash
plotdos                                       # total DOS
plotdos DOSCAR --atoms 0-5 --orbitals all-d   # d-band PDOS for atoms 0-5
plotdos --atoms 0 --orbitals t2g,eg --band-center  # crystal field splitting
plotdos --save dos.png --range -5,5           # save with energy window
```

---

## Thermochemistry

### `thermochem`

Calculate Gibbs free energy corrections from vibrational frequency analysis (VASP frequency calculation OUTCAR).

```bash
thermochem OUTCAR [options]
```

| Flag | Description |
|------|-------------|
| `--temp T` | Temperature in K (default: 298.15) |
| `--press P` | Pressure in Pa (default: 101325) |
| `--atoms INDICES` | Specific atoms for partial Hessian |
| `--gas` | Treat as ideal gas molecule |
| `--geometry TYPE` | For gas: `linear` or `nonlinear` |
| `--spin S` | Spin multiplicity for gas-phase |
| `--symmetry N` | Symmetry number for gas-phase |

Uses ASE's `HarmonicThermo` (surfaces) or `IdealGasThermo` (gas-phase molecules).

---

### `gibbsFE`

Alias/variant of `thermochem` with identical functionality. Calculate corrections to electronic energy for Gibbs free energy.

```bash
gibbsFE OUTCAR [options]
```

Same flags as `thermochem`.

---

## Structure Manipulation

### `reorder_atoms`

Reorder atoms in a structure file by chemical symbol and optionally by z-coordinate. Preserves `FixAtoms` constraints.

```bash
reorder_atoms INPUT [options]
```

| Flag | Description |
|------|-------------|
| `--order EL1 EL2 ...` | Element ordering (e.g., `Cu O H`) |
| `--z-order MODE` | `top-bottom` or `bottom-top` within each element |
| `--output FILE` | Output file (default: overwrite input) |
| `--format FMT` | ASE output format |

**Example:**
```bash
reorder_atoms POSCAR --order Cu O H --z-order top-bottom
```

---

### `asegui`

Launch ASE GUI for structure visualization. Includes preprocessing to handle corrupted files with stray forward slashes.

```bash
asegui FILE [FILE2 ...]
```

---

### `remove-slashes`

Remove all forward slashes from a file. Utility for repairing corrupted VASP files.

```bash
remove-slashes INPUT [OUTPUT]
```

If `OUTPUT` is omitted, the input file is overwritten in place.

---

## Electrochemistry

### `genpotentialdepcalc`

Generate a set of calculation directories for potential-dependent DFT calculations. Creates folders with modified electron counts (NELECT) for applied potential studies.

```bash
genpotentialdepcalc [options]
```

Sets up folders with incrementally modified electron counts based on a reference POSCAR and valence electron configuration.

---

## Database

### `vasp2db`

Extract VASP calculation data from multiple directories into a pandas DataFrame, saved as a pickle file. Extracts structures, parameters, INCAR settings, POTCAR info.

```bash
vasp2db --paths DIR1 DIR2 ... --output DATABASE.pkl [options]
```

| Flag | Description |
|------|-------------|
| `--paths` | Calculation directories to process |
| `--output` | Output pickle file |
| `-r`, `--relative-energy` | Calculate relative energies |
| `-v`, `--verbose` | Verbose output |
| `--skip-errors` | Continue on errors |

**Examples:**
```bash
vasp2db --paths CuCu_2 NiCu_1 --output database.pkl
vasp2db --paths */ --output all_calcs.pkl --relative-energy -v
```

---

## File Management

### `vaspbackup`

Create backup of VASP calculation files with optional compression.

```bash
vaspbackup BACKUPNAME [-c] [-p PATTERN]
```

| Flag | Description |
|------|-------------|
| `-c` | Compress backed-up files with gzip |
| `-p PATTERN` | Custom glob pattern for files to backup |

Copies CONTCAR, POSCAR, KPOINTS, INCAR, OUTCAR, vasprun.xml, and CIF files to a named backup directory.

---

### `vasp2arc`

Convert VASP output files to ARC (DMol3) format.

```bash
vasp2arc INPUT [OUTPUT]
```

If `OUTPUT` is omitted, defaults to `INPUT.arc`.

---

## Visualization

### `view-outcars`

Batch visualization of structures from OUTCAR files across multiple directories.

```bash
view-outcars [options]
```

| Flag | Description |
|------|-------------|
| `--converged` | Only show converged calculations |
| `--folders DIR1 DIR2` | Specific directories |
| `--pyatoms` | Pyatoms directory structure |

Opens ASE GUI for each structure.

---

## Charge Analysis

### `calculate-bader-charges`

Run Bader charge analysis on VASP output. Extracts ZVAL from OUTCAR, reads Bader output (ACF.dat), and calculates atomic charges.

```bash
calculate-bader-charges [options]
```

| Flag | Description |
|------|-------------|
| `-s FILE` | Structure file (default: CONTCAR) |
| `-o FILE` | OUTCAR file (default: OUTCAR) |
| `--format FMT` | Output format: `table` (default) or `csv` |

Requires Bader analysis tools (`chgsum.pl`, `bader`) for initial ACF.dat generation.

---

## Summary Table

| Command | Purpose |
|---------|---------|
| `getenergy` | Quick energy/convergence from OUTCAR |
| `summaryfolders` | Batch analysis across directories |
| `comparevaspparam` | Compare VASP parameters between runs |
| `inputparam-from-outcar` | Extract INCAR parameters from OUTCAR |
| `outcar-extract` | Extract frames from OUTCAR |
| `plotdos` | DOS/PDOS plotting |
| `thermochem` | Gibbs free energy corrections |
| `gibbsFE` | Gibbs free energy corrections (alias) |
| `reorder_atoms` | Reorder atoms by element/z-coordinate |
| `asegui` | ASE GUI launcher |
| `remove-slashes` | File repair (remove `/` characters) |
| `genpotentialdepcalc` | Generate potential-dependent calc folders |
| `vasp2db` | Extract VASP data to pickle DataFrame |
| `vaspbackup` | Backup VASP files |
| `vasp2arc` | Convert VASP → ARC format |
| `view-outcars` | Batch structure visualization |
| `calculate-bader-charges` | Bader charge analysis |
