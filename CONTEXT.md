# ASEtools — Project Context

## Overview
Python toolkit for computational materials science with VASP and ASE. 10.6k lines, 17 CLI tools, domain-organized subpackages.

## Repo Location
`/home/manuelarcer/.openclaw/workspace/repos/asetools`

## Package Structure
```
asetools/
├── analysis/          # VASP output analysis (convergence, energy, forces, symmetry)
├── electronic/        # DOS, PDOS, orbital projections, band center
├── structure/         # SurfaceAnalyzer, bond valence
├── electrochemistry/  # Applied potential calculations
├── pathways/          # NEB, dimer method
├── thermodynamics/    # Ab initio thermodynamics
├── workflow/          # YAML-based multi-stage workflow runner
├── database/          # ASE database integration
├── plotting/          # PES and energy profile plotting
├── cli/               # 17 CLI scripts
├── parsers/           # OUTCAR parser
└── data/              # Runtime data (bvparm2020.cif, surface_properties.json)
```

## Development
```bash
pip install -e ".[dev]"   # Editable install with test deps
pytest                     # Run tests (tests/ directory)
```

## Improvement Roadmap

### A — Console script entry points ✅ DONE (2026-03-09)
Registered all 17 CLI tools as `[project.scripts]` in pyproject.toml. After `pip install`, all commands are on PATH. Branch: `feature/console-entry-points`, commit `3e74bb1`.

### B — Fix top-level exports ✅ DONE (2026-03-10)
Added lazy subpackage loading via `__getattr__` to top-level `__init__.py`. All 11 subpackages accessible via `asetools.<subpackage>` without eager imports. Backward-compatible thermodynamics re-exports preserved. Added `__version__` and module docstring. Branch: `feature/fix-top-level-exports`, commit `3e32836`.

### C — Test infrastructure ✅ DONE (2026-03-11)
- Created .venv, installed dev deps (pytest, spglib)
- Fixed np.trapz → np.trapezoid for NumPy 2.0+ compatibility (doscar.py)
- Added skip decorators for dimer tests requiring POSCAR in working directory
- Added conftest.py with shared fixtures
- Added test_imports.py (lazy loading, version, subpackages)
- Added test_electronic.py (DOS class, band center, Fermi energy)
- Added test_structure.py (SurfaceAnalyzer, BondValence imports)
- Added test_electrochemistry.py (constants, polynomial, electron counting)
- Result: 150 passed, 2 skipped, 0 failures
- Branch: feature/test-infrastructure, commit 2a73728

### D — Type hints ✅ DONE (2026-03-12)
Added type annotations to all public APIs across every module. Modules annotated:
analysis/vasp.py, electrochemistry/appliedpotential.py, structure/adsorbate.py,
database/databases.py, pathways/neb.py, pathways/dimer.py, plotting/plots.py.
Added py.typed marker (PEP 561). Fixed mutable default arg in plot_nebs().
Modules already typed (unchanged): electronic/doscar.py, thermodynamics/ab_initio.py,
analysis/symmetry.py, parsers/vasp_outcar.py, structure/bond_valence.py, workflow/*.py.
Branch: feature/type-hints, commit 40f8da5.

### E — Code quality tooling ✅ DONE (2026-03-13)
Added ruff config to pyproject.toml with rule sets E/W/F/I/UP/B/SIM/RUF. Auto-fixed 672 issues (whitespace, import sorting, f-strings, pyupgrade). Manual fixes for unused vars and collapsible ifs. Legacy patterns (bare except, unicode math, mutable class defaults) documented as ignored with rationale. ruff format deferred (44 files would change — could be a follow-up). Added ruff to dev deps. Branch: feature/code-quality-tooling, commit d59093b.

### F — Documentation refresh ✅ DONE (2026-03-16)
Added docs/index.md (package overview, structure, installation, quick start), docs/cli_reference.md (all 17 CLI tools with usage, flags, examples), and docs/api_reference.md (module-level Python API for all 11 subpackages). Existing docs (thermodynamics, constraints, VASP guide) linked from index. Branch: feature/documentation-refresh, commit e257de8.

### G — CI setup ✅ DONE (2026-03-17)
GitHub Actions workflow for pytest + ruff on push/PR. Lint job (ruff check + format --check, Python 3.12) and test matrix (pytest, Python 3.9–3.12). Fixed remaining 105 lint issues as part of this task. Branch: feature/ci-setup, commit 2154a6e.

### H — Error handling cleanup ✅ DONE (2026-03-19)
Replaced all bare `except:` with specific exception types across 8 files. Added `raise ... from exc` chaining to all re-raises in except clauses (11 instances). Removed E722 and B904 from ruff ignore list. Branch: feature/error-handling-cleanup, commit d0f49bb.

### I — Merge feature branches to main ✅ DONE (2026-03-19)
All feature branches merged to main. `feature/missing-cli-entries` changes (4 CLI entry points, version 0.2.0, requires-python >=3.9, modernized type hints) already in main.

### J — Test coverage expansion ✅ DONE (2026-03-20)
Added 48 new tests across 3 test files:
- test_parsers_utils.py: 18 tests for OUTCAR parser utility functions
- test_plots.py: 21 tests for PES plotting (validate_columns, add_line_to_pes, beautify_pes_plot)
- test_analysis_convergence.py: 9 tests for check_outcar_convergence with synthetic OUTCARs
Coverage: 29% → 32% (198 tests total). Branch: feature/test-coverage-expansion, commit ec56d76.

### K — Further test coverage (CLI tools, database, workflow) ✅ DONE (2026-03-23)
Added 65 new tests across 2 test files:
- test_cli_tools.py: 36 tests for remove_slashes, reorder_atoms, comparevaspparam,
  vaspbackup, gibbsFE helpers (vib extraction, corrections, atom index parsing,
  mode character, shift), summaryfolders
- test_workflow_utils.py: 29 tests for YAML config loading/validation, deep_update,
  setup_initial_magmom (dict/list/numpy/None/error cases), VASPConfigurationFromYAML,
  configure_logging, database functions (check_if_exists, db_to_pandas)
Coverage: 32% → 35% (215 tests total). Branch: feature/test-coverage-cli-workflow, commit 76c60bc.

### L — ruff format (auto-format all files) ✅ DONE (2026-03-24)
Verified all 67 files already pass `ruff format --check` (no formatting changes needed — resolved incrementally during prior tasks). Fixed one unused import in test_workflow_utils.py. Added .coverage to .gitignore. Branch: feature/ruff-format, commit d04ef6c.

### M — Core module test coverage (electronic, parsers, database) ✅ DONE (2026-03-24)
Added 96 new tests across 3 files:
- test_parsers_outcar.py: 33 tests for OUTCAR parser utilities, header parsers,
  exception classes, read_constraints_from_file, integration with real OUTCAR
- test_doscar_extended.py: 55 tests for DOS properties, PDOS extraction,
  band center edge cases, legacy functions, plotting smoke tests
- test_database_extended.py: 8 tests for check_if_exists_in_db and db_to_pandas
Coverage: 35% → 38% overall (doscar 57%→85%, parser 33%→51%).
Branch: feature/core-test-coverage, commit 9508d91.

### N — Package publishing prep (README, LICENSE, metadata) ✅ DONE (2026-03-25)
Added MIT LICENSE file, MANIFEST.in for sdist, complete pyproject.toml metadata (classifiers, keywords, urls, readme, license fields). Fixed README badge (3.9+), added MIT badge. All 311 tests pass. Branch: feature/package-publishing-prep, commit 5b4b0fa.

### O — Deprecation warning cleanup ✅ DONE (2026-03-26)
Addressed both deprecation warnings from test runs:
- **scipy.odr** (deprecated 1.17, removal 1.19): Created `_odr_compat.py` compatibility
  layer using `odrpack` package (recommended replacement). Falls back to scipy.odr
  with warning suppressed. Updated both appliedpotential.py files.
- **spglib OLD_ERROR_HANDLING**: Set `spglib.error.OLD_ERROR_HANDLING = False` at all
  import points. Updated `_get_symmetry_dataset` to catch `SpglibError`.
- Added `odrpack` to dependencies. All 311 tests pass, 0 deprecation warnings.
Branch: feature/deprecation-warning-cleanup, commit 69d64c7.

### S — Push local work to origin
Pending — 12+ commits ahead of remote. Needs owner approval.

### T — CHANGELOG.md ✅ DONE (2026-04-02)
Added CHANGELOG.md following Keep a Changelog format. Covers full 0.2.0 release
(all changes from tasks A–R) and placeholder 0.1.0 entry for pre-modernization code.
Branch: feature/changelog, commit 9b42181.

## Priority Order
A → B → C → D → E → F → G → H → I → J → K → L → M → N → O → P → Q → R → S → T (one per morning dev session)

## Potential Future Tasks
- **U** — Publish to PyPI (test first, then production)
- **V** — Incrementally remove mypy per-module ignore_errors (tighten type checking)

## Notes
- Owner: Juan M Arce-Ramos
- Venv at `.venv` (Python 3.14, all dev deps installed)
- Tests use real VASP output files in `tests/data/`
- CLAUDE.md exists with guidance for Claude Code sessions
