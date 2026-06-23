# Modular Stage Composition + MLIP Pre-optimization

**Date:** 2026-06-23
**Status:** Draft — awaiting implementation on `feature/modular-stage-composition`
**Scope:** New YAML section, new public Python API, one new subprocess runner. Additive — existing `workflows:` path untouched.

## Problem

The workflow manager defines a multi-stage VASP procedure entirely inside the YAML `workflows:` section. Each named workflow is a fixed, fully-specified list of stages. Composing a new procedure means authoring a new named workflow, and the existing file shows the cost: `relax`, `relax_low_potim`, `relax_last_step`, and `relax_ase_opt` are near-duplicates differing in one parameter (a `potim`, a constraint, an ENCUT), and `relax_ase_opt`'s 340 → 420 → 500 ENCUT ramp is the *same* stage written out three times. Every objective spawns another copy.

Two things are wanted:

1. **Compose the procedure in the Python submission script**, from reusable, parameterized building blocks, while the canonical *final* optimization parameters stay authoritative in YAML.
2. **Allow MLIP optimizations as pre-optimization stages**, driven by the `mlip_platform` package, so a cheap warm-start can precede the DFT stages.

## Glossary (authoritative: `CONTEXT.md`)

- **Stage** — unit of restart, backup, and constraint application. Writes `STAGE_{name}_DONE`. Owns its constraints and overrides. The unit a submission script composes.
- **Step** — sub-phase *within* a stage (shared geometry/constraints/VaspInteractive context). Not an independent restart unit.
- **Stage template** — reusable, parameterized stage body in the new `stage_templates:` YAML section.
- **Production stage** — the canonical final optimization, from the template named `production` (default).
- **Engine** — `vasp` (default) or `mlip`; what performs a stage's optimization.
- **MLIP stage** — `engine: mlip` stage; relaxes via `mlip_platform` in a subprocess, in a separate **MLIP env**.
- **MLIP env registry** — `globals.mlip_envs:` map from env key to interpreter path.

## Decisions

### D1 — `stage_templates:` section (additive)

New top-level YAML key, parallel to `workflows:`. Optional; **not** added to the required-keys check in `verify_configuration_keys`. Maps a template name to a stage body. The existing `workflows:` section and every named workflow (`relax*`, `dos`, `bader`, `freq`, `dimer`, the spec'd refactor) are untouched and keep working. Old calculation directories keep resuming against their named workflows. Coexistence now; the `relax*` family is documented as superseded and deprecated once legacy projects finish.

```yaml
stage_templates:
  encut_ramp:                    # a VASP stage template
    constraints:
      type: hookean
      config_file: hookean_c_pairs.json
      spring_constant: 20.0
    steps:
      - name: single_point
        overrides: { nsw: 0, nelm: 800, ispin: 1, lcharg: true, ediff: 0.0001 }
      - name: optimization_bfgs
        overrides: { nsw: 500, lcharg: true, ediff: 0.0001 }
        optimizer: BFGS
        optimizer_kwargs: { maxstep: 0.2 }

  production:                    # the canonical final optimization (no constraints block)
    steps:
      - name: single_point
        overrides: { nsw: 0, nelm: 800, ispin: 1, lcharg: true, encut: 500 }
      - name: optimization
        overrides: { encut: 500 }
```

### D2 — `cfg.stage()` factory

Method on `VASPConfigurationFromYAML`. Instantiates a template into a concrete stage dict of **the exact shape `_run_stage` already consumes** (`{name, constraints?, steps:[...]}`), plus a private `_template` key (provenance; ignored at execution).

```python
cfg.stage(template, *, name, overrides=None, optimizer_kwargs=None, constraints=None)
```

- `name` — **required, explicit**. Becomes the `STAGE_{name}_DONE` sentinel, backup suffix, restart key. No auto-naming (auto-generated names silently change when a parameter is tweaked on resubmit, re-running completed stages).
- `overrides` / `optimizer_kwargs` / `constraints` — **explicit buckets**, each **deep-merged** onto the template's baked-in values (patch only what differs). Applied to **every step** in the instance.
- Unknown template name → clear `KeyError`.
- Overrides are **unrestricted** for now (no declared safe-set; can touch any VASP tag).

### D3 — Precedence

VASP kwargs reaching the calculator, lowest → highest:

```
basic → system → template step overrides → per-instance stage overrides → run_overrides
```

`run_overrides` is the most specific and always wins (see [ADR 0003](../../adr/0003-run-overrides-is-authoritative.md)); template/per-instance step values are defaults it can override run-wide.

> Note: this supersedes the original design, which placed `run_overrides` *below* step-level values for backward compatibility (ADR 0001). That position made a run-wide override silently lose to a per-step value and was reversed by ADR 0003.

### D4 — Two-door execution API

```python
run_stages(atoms, cfg, *, stages, production='production', run_overrides=None, dry_run=False, magmoms=None)
run_workflow(atoms, cfg, *, workflow_name, run_overrides=None, dry_run=False, magmoms=None)
```

- `run_stages` is the **engine**: runs an explicit ordered list of resolved stage dicts. The magmom / atom-reorder / reference-file machinery and the per-stage restart/backup logic currently at the top of `run_workflow` move here so both paths share them.
- `run_workflow` becomes a **thin wrapper**: looks up `cfg.workflows[workflow_name]['stages']`, calls `run_stages(..., production=None)`. No logic duplicated.
- The two doors are kept distinct (not one function with mutually-exclusive args) so a submission script self-documents which mode it is in.

Submission script:

```python
PREOPTIMIZATION = [
    cfg.stage('encut_ramp', name='OPT_340', overrides={'encut': 340}, optimizer_kwargs={'fmax': 0.05}),
    cfg.stage('encut_ramp', name='OPT_420', overrides={'encut': 420}, optimizer_kwargs={'fmax': 0.02}),
]
run_stages(atoms, cfg, stages=PREOPTIMIZATION + [cfg.stage('production', name='OPT_500')])
```

### D5 — Production-stage warning

`run_stages` emits a **non-blocking warning** if the `production`-template stage (name from `production=`, default `'production'`) is **not the last element** of `stages`. Covers both "production absent" and "production followed by another stage". `production=None` disables the check (e.g. a pure-DOS modular pipeline). The named-workflow path passes `production=None`.

```
⚠ Production stage 'production' was not the final stage — workflow ends on 'OPT_420'.
  The last optimization may not use canonical YAML parameters.
```

Out of scope: warning on `run_overrides` contamination of a present production stage (R1).

### D6 — Constraints

- A template body may carry a `constraints:` block (default).
- `cfg.stage(constraints=...)` **deep-merges** onto it (e.g. `constraints={'spring_constant': 30.0}` keeps `type`/`config_file`).
- `config_file` resolves **relative to the calculation directory** (cwd), as today.
- To run a stage unconstrained, use a template with **no** `constraints:` block (the production template). No `constraints=False` escape hatch for now.

### D7 — MLIP stage (engine: mlip)

MLIP packages pin mutually-incompatible torch stacks and live in per-MLIP envs (see `mlip_platform` ADR 0001), so an MLIP stage **cannot** run in the VASP job env. It executes as a **subprocess** in a separate MLIP env. `asetools` is a declared dependency of `mlip_platform`, so the same `ConstraintManager` code is available inside every MLIP env.

Template shape:

```yaml
stage_templates:
  uma_preopt:
    engine: mlip
    mlip: uma-s-1p2          # MLIP tag passed to setup_calculator
    env: uma                 # key into globals.mlip_envs
    optimizer: bfgs
    fmax: 0.10
    max_steps: 300
    relax_cell: false
    device: auto
    constraints:             # re-applied in-subprocess from the same JSON
      type: hookean
      config_file: hookean_c_pairs.json
      spring_constant: 20.0

globals:
  mlip_envs:
    uma:  /path/to/uma-env/bin/python
    mace: /path/to/mace-env/bin/python
```

The env registry centralizes interpreter paths so a template says `env: uma` and only `globals.mlip_envs` changes across machines (laptop / ASPIRE2A / cos-cluster). Resolves to an explicit interpreter path — no `conda activate` subshell.

#### MLIP runner (new)

A small entry point added to `asetools` (importable in the MLIP env), invoked as:

```
<interpreter> -m asetools.workflow.mlip_runner \
    --structure CONTCAR \
    --mlip uma-s-1p2 --optimizer bfgs --fmax 0.10 --max-steps 300 \
    --device auto [--relax-cell] \
    --constraints-json hookean_c_pairs.json --spring-constant 20.0 [--distance-factor ...] \
    --name OPT_MLIP
```

Runner contract:

1. `read()` the current structure (carries `FixAtoms` via selective dynamics).
2. If constraints given: `ConstraintManager().apply_stage_constraints(atoms, cfg)` — rebuilds Hookean springs from the same JSON the VASP stages use. `run_optimization` honors `atoms.constraints` automatically (constraint-adjusted forces).
3. `setup_calculator(atoms, mlip, uma_task/mace_head, device)`.
4. `run_optimization(atoms, optimizer, fmax, max_steps, relax_cell, output_dir=cwd)`.
5. Write **`CONTCAR`** (relaxed structure) and back up as `CONTCAR_{name}` plus `opt.log` / `opt.traj` / `opt_convergence.csv`.
6. Exit non-zero if `run_optimization` returned not-converged, so the manager refuses to mark the stage `DONE` (parity with VASP stages raising on non-convergence).

The manager builds the command from the resolved stage dict, runs it via `subprocess`, checks the exit code, then writes `STAGE_{name}_DONE`.

### D8 — Structure handoff and the ordering invariant

**Invariant: MLIP stages are pre-optimization only — every MLIP stage precedes every VASP stage.** Rationale: MLIP is a cheap warm-start; nothing downgrades to an MLIP after paying for DFT.

Given the invariant, the MLIP runner writing a `CONTCAR` is sufficient and `load_structure` is **left unchanged**: a stale `OUTCAR` can only come from a VASP stage, and no VASP stage runs before an MLIP stage, so when the post-MLIP stage loads its structure there is no `OUTCAR` on disk and `load_structure` falls through to the MLIP `CONTCAR`. Restart is consistent: the last completed stage at any resume point is either the MLIP stage (no `OUTCAR` → `CONTCAR` read) or a later VASP stage (its `OUTCAR` is freshest).

`run_stages` enforces the invariant with a **hard error at validation time** (before any compute): an `engine: mlip` stage that appears after an `engine: vasp` stage raises. Unlike the production warning, this is a correctness failure (silently loading the wrong structure), so it fails fast.

## Backwards compatibility

- `workflows:` section, all named workflows, and the existing `run_workflow(workflow_name=...)` signature unchanged. `workflow_name` becomes keyword-only is **not** required — keep the existing positional/keyword call working.
- `stage_templates:` and `globals.mlip_envs:` are optional; configs without them load exactly as before.
- Existing calculation directories keep resuming against named workflows.
- No migration tooling.

## Non-goals

- Removing/collapsing the `relax*` family (deprecate later).
- `run_overrides`-contamination warning on a present production stage (R1).
- A declared safe-set restricting which keys per-instance overrides may touch.
- In-process MLIP (rejected — see ADR 0002).
- MLIP stages anywhere except before all VASP stages.

## Acceptance criteria

1. `cfg.stage(template, name=..., overrides=..., optimizer_kwargs=..., constraints=...)` returns a resolved stage dict matching `_run_stage`'s expected shape, with deep-merge precedence per D3 and a `_template` provenance key.
2. `run_stages(atoms, cfg, stages=[...])` runs an explicit list; `run_workflow(workflow_name=...)` delegates to it with `production=None`; existing named workflows behave identically (verified against `dos`, `relax_ase_opt`).
3. Production warning fires iff the `production` stage is not last; suppressed by `production=None`.
4. MLIP stage: `run_stages` builds and runs the subprocess against the resolved env interpreter, the runner re-applies Hookean from JSON, writes `CONTCAR` + `CONTCAR_{name}`, and a non-converged MLIP run leaves no `STAGE_*_DONE`.
5. Ordering guard: composing `[..., vasp_stage, mlip_stage, ...]` raises at validation time before any stage runs.
6. `dry_run=True` enumerates the composed stage names (VASP and MLIP) in order with no compute and no subprocess launch.
7. A config with no `stage_templates:` / `globals.mlip_envs:` still loads and runs named workflows.
8. Unit tests cover `cfg.stage` merge/precedence, production-warning logic, ordering-guard error, and MLIP command construction (mock subprocess). Live MLIP smoke test is user-executed (out of scope for the implementer).

## Out-of-scope live test (user-executed)

In a scratch directory: compose `[uma_preopt, OPT_340, OPT_420, production]`, run, kill after `OPT_340`, resubmit, confirm only `OPT_420` + `production` run and the structure carried forward is the MLIP-then-VASP-refined one. Separately confirm a constrained `uma_preopt` honors the Hookean springs (compare against an unconstrained run).
