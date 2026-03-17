# CLI Reference

ASEtools ships 17 command-line tools, all registered as console scripts in `pyproject.toml`. After `pip install -e .` they are available on `PATH`.

---

## Analysis Tools

### `getenergy`

Quick convergence check and energy/force extraction from a single OUTCAR.

```bash
getenergy [OUTCAR]
```

- `OUTCAR` — path to OUTCAR file. Defaults to `OUTCAR` in the current directory if omitted.

**Output:** Three-column table: Converged, MaxForce (eV/Å), Energy (eV).

```
Converged   MaxForce      Energy
True           0.015   -234.5678
```

---

### `summaryfolders`

Batch analysis across subdirectories. Scans for OUTCAR files, extracts convergence status, ENCUT, EDIFFG (target fmax), max force, and energy. Caches results to `summary.log`.

```bash
summaryfolders [-m] [-f] [-p]
```

| Flag | Long flag | Description |
|------|-----------|-------------|
| `-m` | `--magmom` | Extract and display total magnetic moments |
| `-f` | `--fast` | Fast mode: reads `vasp.out` or `out.txt` instead of OUTCAR |
| `-p` | `--pyatoms` | Pyatoms mode: analyze multi-step pyatoms calculations (reads `step_*/OUTCAR`) |

**Output columns (standard mode):** Config, ISIF, Converged, ENCUT, Target-fmax, MaxForce, Energy, Rel.E (and MagMom with `-m`).

**Output columns (pyatoms mode):** Config, E_step0, E_step1, ..., MagMom, Converged, Rel.E. Also writes `summary_pyatoms.log`.

Detects ASE optimizer log files (BFGS.log, FIRE.log, LBFGS.log, etc.) and uses their convergence status when present.

---

### `comparevaspparam`

Side-by-side comparison of VASP parameters between two OUTCAR files. Interactive: prompts for file paths via stdin.

```bash
comparevaspparam
```

No command-line arguments. Prompts:
```
Enter path to first OUTCAR file:
Enter path to second OUTCAR file:
```

Compares all parameters found in both files and prints differences.

---

### `inputparam-from-outcar`

Extract VASP parameters from an OUTCAR file, grouped by category.

```bash
inputparam-from-outcar OUTCAR
```

| Argument | Description |
|----------|-------------|
| `outcar_file` | Path to OUTCAR file (required positional) |

**Categories reported:** StartParam (PREC, ISTART, ICHARG, ISPIN), Electronic (ENCUT, NELM, EDIFF, LREAL, LMAXMIX, VOSKOWN), IonicRelax (EDIFFG, NSW, IBRION, NFREE, ISIF, ISYM, LCORR, POTIM), System (POMASS, ZVAL, RWIGS, VCA, NELECT, NUPDOWN), DOS (EMIN, EMAX, EFERMI, ISMEAR, IALGO, SIGMA), ElectronicRelax (IALGO), Write (LWAVE, LCHARG, LVTOT, LVHAR, LELF, LORBIT), DipoleCorr (LMONO, LDIPOL, IDIPOL, EPSILON), ExCorr (GGA, VOSKOWN, LHFCALC, LHFONE, AEXX), LinearResp (LEPSILON).

---

### `outcar-extract`

Extract a structure frame from an OUTCAR and write it to a file.

```bash
outcar-extract OUTCAR [index] [output_file] [format]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `OUTCAR` | — | Path to OUTCAR file (required) |
| `index` | `-1` (last frame) | Frame index to extract |
| `output_file` | `final.vasp` | Output file path |
| `format` | auto-detect from extension | ASE format string (e.g., `extxyz`) |

**Examples:**
```bash
outcar-extract OUTCAR                              # last frame → final.vasp
outcar-extract OUTCAR 5                            # 5th frame → final.vasp
outcar-extract OUTCAR -1 final.extxyz extxyz       # with explicit format
```

---

## Electronic Structure

### `plotdos`

Plot total DOS or partial DOS (PDOS) from a VASP DOSCAR file.

```bash
plotdos [options]
```

**Input/output:**

| Flag | Default | Description |
|------|---------|-------------|
| `--doscar FILE` | `DOSCAR` | Path to DOSCAR file |
| `-o FILE`, `--output FILE` | interactive | Output file path; format determined by extension (png, pdf, svg, etc.) |
| `--save-data FILE` | — | Save plot data to CSV (`.csv`) or tab-separated (`.txt`) |

**What to plot:**

| Flag | Description |
|------|-------------|
| `--atoms SPEC` | Atom indices: `0`, `0,2,5`, `0-5`, or `all` |
| `--orbitals TYPE` | Orbital(s): `s`, `p`, `d`, `t2g`, `eg`, `all` (comma-separated for multiple) |
| `--total` | Plot total DOS (overrides `--atoms`/`--orbitals`) |
| `--overlay` | Overlay multiple orbitals on the same axes (use with multiple `--orbitals`) |

When neither `--atoms` nor `--orbitals` is given, total DOS is plotted.

**Appearance:**

| Flag | Default | Description |
|------|---------|-------------|
| `--xlim MIN,MAX` | full range | Energy range in eV, e.g., `-5,5` |
| `--ylim MIN,MAX` | auto | DOS range |
| `--title TEXT` | auto | Plot title |
| `--dpi N` | `300` | Resolution for saved figure |
| `--figsize W,H` | `8,6` | Figure size in inches |
| `--linewidth N` | `1.5` | Line width |
| `--colors C1,C2,...` | auto | Custom colors (comma-separated) |
| `--same-color-spins` | — | Use same color for spin-up and spin-down |
| `--no-legend` | — | Disable legend |

**Analysis:**

| Flag | Default | Description |
|------|---------|-------------|
| `--band-center` | — | Calculate and print band center (requires `--atoms` and `--orbitals`) |
| `--energy-range MIN,MAX` | full range | Energy range for band center integration |
| `--spin-treatment MODE` | `combined` | Spin handling: `combined`, `separate`, `up`, `down` |

**Orbital names:** `s`, `p`, `d` are the primary interface (internally mapped to `all-s`, `all-p`, `all-d`). Also accepted: `t2g` (dxy+dyz+dxz), `eg` (dz2+dx2-y2), `all`.

**Examples:**
```bash
plotdos                                         # total DOS
plotdos --atoms 0 --orbitals d                  # d PDOS for atom 0
plotdos --atoms 0-5 --orbitals d                # d PDOS summed over atoms 0-5
plotdos --atoms all --orbitals s,p,d --overlay  # overlay s/p/d for all atoms
plotdos --atoms 0 --orbitals t2g,eg --overlay   # crystal field splitting
plotdos --atoms 0 --orbitals d --band-center    # d-band center
plotdos --atoms all --orbitals d --band-center --spin-treatment separate  # spin-resolved band center
plotdos --atoms 0 --orbitals d --output dos.png --dpi 600
plotdos --atoms all --orbitals p,d --overlay --save-data dos_data.csv
```

---

## Thermochemistry

### `thermochem`

Calculate thermochemistry corrections (ZPE, CpT, entropy) from a VASP vibrational frequency OUTCAR. Uses ASE's `HarmonicThermo` for surface adsorbates.

```bash
thermochem OUTCAR [options]
```

| Argument/Flag | Default | Description |
|---------------|---------|-------------|
| `outcar` | — | Path to OUTCAR file (required positional) |
| `--writevib {y,n}` | `n` | Write vibration trajectory files |
| `--temp T` | `298` | Temperature in K |
| `--gas BOOL` | `False` | Use `IdealGasThermo` for gas-phase molecules |
| `--symnum N` | `2` | Symmetry number (for `--gas` mode) |
| `--geom {monoatomic,linear,nonlinear}` | `linear` | Molecular geometry (for `--gas` mode) |
| `--pressure P` | `1.0` | Gas pressure in bar (for `--gas` mode) |
| `--select-atoms SPEC` | — | Atom indices for mode selection, e.g., `0-4` or `0,1,2` |
| `--threshold N` | `0.5` | Minimum character fraction (0–1) for mode selection |

**Output:** Frequency table, ZPE, thermal correction (CpT), entropy (S), harmonic limit via ASE HarmonicThermo. Results also logged to `FREQ_ANALYSIS_YYYYMMDD_HHMM.log`.

**Examples:**
```bash
thermochem OUTCAR                              # all modes, T=298 K
thermochem OUTCAR --temp 400                   # at 400 K
thermochem OUTCAR --select-atoms 0-4           # only modes on atoms 0-4
thermochem OUTCAR --select-atoms 0,1,2 --threshold 0.7
thermochem OUTCAR --gas True --geom linear --symnum 1   # gas-phase CO
```

---

### `gibbsFE`

Same functionality as `thermochem`. Calculate corrections to electronic energy for Gibbs free energy.

```bash
gibbsFE OUTCAR [options]
```

Accepts the same arguments as `thermochem`.

---

## Structure Manipulation

### `reorder-atoms`

Reorder atoms in a structure file by element symbol and/or z-coordinate.

```bash
reorder-atoms input_file [options]
```

| Argument/Flag | Default | Description |
|---------------|---------|-------------|
| `input_file` | — | Input structure file (required positional) |
| `-o/--order EL1 EL2 ...` | — | Element ordering, e.g., `Cu O H` |
| `-z/--z-order {top-bottom,bottom-top}` | — | Sort within each element by z-coordinate |
| `-out/--output FILE` | `reordered_POSCAR` | Output file path |
| `-f/--format FMT` | `vasp` | ASE output format |

**Examples:**
```bash
reorder-atoms POSCAR --order Cu O H
reorder-atoms POSCAR --order Cu O H --z-order top-bottom
reorder-atoms POSCAR --z-order bottom-top --output sorted.vasp
```

---

### `asegui`

Open a VASP structure file in ASE GUI. Strips stray `/` characters from element symbols before loading.

```bash
asegui FILENAME
```

| Argument | Description |
|----------|-------------|
| `FILENAME` | Path to structure file (POSCAR, CONTCAR, OUTCAR, etc.) |

---

### `remove-slashes`

Remove all `/` characters from a file. Used to repair VASP files with corrupted element symbols.

```bash
remove-slashes input_file [output_file]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `input_file` | — | Input file (required) |
| `output_file` | overwrites input | Output file path |

---

## Electrochemistry

### `genpotentialdepcalc`

Generate calculation directories for potential-dependent DFT calculations by creating folders with modified electron counts (NELECT).

```bash
genpotentialdepcalc poscar pythonfile [--lower N] [--higher N] [--step N]
```

| Argument/Flag | Default | Description |
|---------------|---------|-------------|
| `poscar` | — | Path to POSCAR file (required positional) |
| `pythonfile` | — | Path to Python submission script (required positional) |
| `--lower N` | `-1` | Lower bound of potential range |
| `--higher N` | `1` | Upper bound of potential range |
| `--step N` | `0.2` | Step size for potential range |

Hardcoded valence electrons: Cu=11, Zn=12, C=4, O=6, H=1.

---

## Database

### `vasp2db`

Extract VASP calculation data from multiple directories into a pandas DataFrame saved as a pickle file. Extracts structures, VASP parameters, energies, convergence status, INCAR content, and POTCAR info.

```bash
vasp2db --paths DIR1 DIR2 ... --output DATABASE.pkl [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--paths DIR ...` | — | Calculation directories (required, one or more) |
| `-o/--output FILE` | — | Output pickle file path (required) |
| `--initial-pattern GLOB` | `*.vasp` | Glob pattern for initial structure files |
| `-r/--relative-energy` | — | Add `Rel.E` column (relative to minimum energy) |
| `-v/--verbose` | — | Print progress details |
| `--skip-errors` | — | Skip folders with errors instead of stopping |

**DataFrame columns:** Path, AbsPath, InitialStructure (ASE Atoms), FinalStructure (ASE Atoms), CalcType, Formula, ENCUT, KSPACING, EDIFF, EDIFFG, GGA, IBRION, ISPIN, NSW, ISIF, NFREE, Energy, TotMagMom, MaxForce, Converged, VASPVersion, INCAR_full, POTCAR_info, Rel.E (if `-r`).

**Loading the database:**
```python
import pickle
with open("database.pkl", "rb") as f:
    df = pickle.load(f)
print(df[["Path", "CalcType", "Energy", "Converged"]])
```

**Examples:**
```bash
vasp2db --paths CuCu_2 NiCu_1 NiCu_2 --output database.pkl
vasp2db --paths */ --output all_calcs.pkl --relative-energy -v
vasp2db --paths calc1 calc2 --output db.pkl --skip-errors
```

---

## File Management

### `vaspbackup`

Backup VASP calculation files into a named subdirectory, compressing large files with gzip.

```bash
vaspbackup backupname
```

| Argument | Description |
|----------|-------------|
| `backupname` | Name for the backup directory (required positional) |

Copies CONTCAR, POSCAR, KPOINTS, INCAR, and other VASP files; compresses OUTCAR and vasprun.xml with gzip.

---

### `vasp2arc`

Convert a VASP structure or trajectory file to DMol3 ARC format.

```bash
vasp2arc input_file [output_file]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `input_file` | — | Input VASP file (required positional) |
| `output_file` | auto (`.arc` extension) | Output ARC file (optional positional) |

---

## Visualization

### `view-outcars`

Visualize the last ionic configuration from OUTCAR files across subdirectories.

```bash
view-outcars [--converged] [--folders DIR ...] [--pyatoms]
```

| Flag | Description |
|------|-------------|
| `--converged` | Only show/open converged calculations |
| `-f/--folders DIR ...` | Specific directories to check (default: all subdirectories) |
| `--pyatoms` | Pyatoms directory layout |

Opens all found structures together in a single ASE GUI viewer.

---

## Charge Analysis

### `calculate-bader-charges`

Perform Bader charge analysis from VASP output files. Reads ZVAL from OUTCAR, reads Bader output (ACF.dat), and computes final atomic charges.

```bash
calculate-bader-charges [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-s/--structure FILE` | `CONTCAR` | Structure file |
| `-o/--outcar FILE` | `OUTCAR` | OUTCAR file |
| `-a/--acf FILE` | `ACF.dat` | Bader ACF.dat output file |
| `--output FILE` | `bader_charges.txt` | Results output file |
| `--format {table,csv,json}` | `table` | Output format |
| `--skip-bader` | — | Skip running bader, use existing ACF.dat |
| `-q/--quiet` | — | Suppress non-essential output |
| `-v/--verbose` | — | Verbose output |

Requires external Bader analysis tools (`chgsum.pl`, `bader`) to be installed and on PATH. The tool runs them automatically unless `--skip-bader` is passed (which expects an existing ACF.dat).

---

## Summary Table

| Command | Purpose |
|---------|---------|
| `getenergy` | Quick energy/convergence/force from OUTCAR |
| `summaryfolders` | Batch analysis across calculation directories |
| `comparevaspparam` | Interactively compare VASP parameters between two OUTCARs |
| `inputparam-from-outcar` | Extract categorized VASP parameters from OUTCAR |
| `outcar-extract` | Extract a frame from OUTCAR to a structure file |
| `plotdos` | Total or partial DOS/PDOS plotting with band center |
| `thermochem` | Thermochemistry corrections (ZPE, CpT, entropy) |
| `gibbsFE` | Gibbs free energy corrections (same as thermochem) |
| `reorder-atoms` | Reorder atoms by element and/or z-coordinate |
| `asegui` | Open VASP structure file in ASE GUI |
| `remove-slashes` | Remove `/` characters from a file |
| `genpotentialdepcalc` | Generate potential-dependent calculation folders |
| `vasp2db` | Extract VASP calculations to pandas DataFrame pickle |
| `vaspbackup` | Backup and compress VASP output files |
| `vasp2arc` | Convert VASP structure to DMol3 ARC format |
| `view-outcars` | Batch structure visualization from OUTCARs |
| `calculate-bader-charges` | Bader charge analysis |
