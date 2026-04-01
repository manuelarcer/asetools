# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-03-31

Major modernization release: package restructured for proper installation, fully
linted, typed, documented, tested (311+ tests), and CI-enabled.

### Added

- **Console script entry points** for all 17 CLI tools — installable via `pip install`
  and available on PATH without manual `scripts/` management (A)
- **Lazy subpackage loading** via `__getattr__` in top-level `__init__.py` — all 11
  subpackages accessible as `asetools.<subpackage>` without eager imports (B)
- **`__version__`** attribute and module docstring at package level (B)
- **Test infrastructure** with pytest, conftest.py fixtures, and skip decorators for
  environment-dependent tests (C)
- **Type annotations** on all public APIs with `py.typed` marker (PEP 561) (D)
- **ruff** configuration (E/W/F/I/UP/B/SIM/RUF rule sets) and auto-formatted
  codebase (E, L)
- **Documentation**: `docs/index.md` (overview, installation, quick start),
  `docs/cli_reference.md` (all 17 CLI tools), `docs/api_reference.md` (full
  Python API for all 11 subpackages) (F)
- **GitHub Actions CI** with pytest (Python 3.9–3.12) and ruff lint/format checks (G)
- **mypy** type checking in CI with pragmatic per-module configuration (R)
- **MIT LICENSE** file, `MANIFEST.in`, and complete `pyproject.toml` metadata
  (classifiers, keywords, urls, readme, license) (N)
- **ODR compatibility layer** (`_odr_compat.py`) using `odrpack` package as
  replacement for deprecated `scipy.odr` (O)
- **311+ tests** across all major modules:
  - Imports, electronic (DOS, band center), structure, electrochemistry (C)
  - Parsers, plotting, convergence analysis (J)
  - CLI tools, workflow utils, database (K)
  - OUTCAR parser, DOSCAR extended, database extended (M)
  - Electrochemistry ODR compat and applied potential (P)
  - NEB and dimer pathways (Q)

### Fixed

- `np.trapz` → `np.trapezoid` for NumPy 2.0+ compatibility (C)
- Mutable default argument in `plot_nebs()` (D)
- 672 lint issues (whitespace, import sorting, f-strings, pyupgrade patterns) (E)
- 105 additional lint issues for CI compliance (G)
- All bare `except:` replaced with specific exception types (8 files) (H)
- `raise ... from exc` chaining added to all re-raises (11 instances) (H)
- `spglib` `OLD_ERROR_HANDLING` set to `False` at all import points;
  `_get_symmetry_dataset` catches `SpglibError` (O)
- `fit_data` spline path crash (accessing `.beta` on `UnivariateSpline`) (P)
- `interpolate_new_x` wrong coefficient order with `np.polyval` (P)
- Missing `raise` in `calculate_band_center` found by mypy (R)
- Missing `return` in `get_cluster_around_adsorbate` found by mypy (R)
- Implicit `Optional` parameters in `vasp_outcar.py` (R)
- Dict type inference in `calculate_bader_charges.py` (R)

### Changed

- Minimum Python version: `>=3.9` (was unspecified)
- Version bumped from 0.1.0 to 0.2.0
- `odrpack` added to dependencies (replaces deprecated `scipy.odr`)
- `mypy`, `types-PyYAML` added to dev dependencies
- `ruff` added to dev dependencies

## [0.1.0] — Pre-2026

Initial package with VASP analysis tools, electrochemistry module, NEB/dimer
pathways, thermodynamics, workflow runner, database integration, and CLI scripts.
Developed as a research utility by Juan M Arce-Ramos.
