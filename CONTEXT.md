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

### I — Merge feature branches to main
Pending: `feature/missing-cli-entries` (4 missing CLI entry points, version bump to 0.2.0, requires-python >=3.9, modernized type hints). Needs Juan's review.

## Priority Order
A → B → C → D → E → F → G → H (one per morning dev session)

## Notes
- Owner: Juan M Arce-Ramos
- No venv currently set up in repo — create one at `.venv`
- Tests use real VASP output files in `tests/data/`
- CLAUDE.md exists with guidance for Claude Code sessions
