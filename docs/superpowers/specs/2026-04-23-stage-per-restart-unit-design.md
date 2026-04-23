# Stage-per-Restart-Unit YAML Refactor

**Date:** 2026-04-23
**Status:** Draft — awaiting user review
**Scope:** YAML-only. No Python code changes.

## Problem

The `workflow.manager` tracks progress at stage granularity via `STAGE_{name}_DONE` sentinel files written by `_mark_done` at the end of `_run_stage`. On resubmission, `stages_to_run` filters out completed stages, and `_run_stage` runs the stage's steps from the first one.

Workflows that pack multiple independent SCF calculations into a single stage therefore cannot be resumed mid-pipeline. The `dos` workflow is the canonical case: one stage (`DOS`) with three steps (`non-spinpolarized` → `spinpolarized` → `non-self-consistent`). If the job is killed during step 3, no sentinel is written, and on resubmit all three steps re-run from scratch — wasting the completed steps' wall-time.

The same pattern exists in `bader` (two steps) and `freq` (two steps: non-spinpolarized warm-start + ibrion=5 Hessian).

## Decision

Promote each independent SCF into its own stage. The existing stage-level restart machinery then produces correct resume behaviour with **zero code changes**.

### Why this works

- `stages_to_run` skips stages with a `STAGE_{name}_DONE` sentinel. Splitting one stage into N gives N independent restart points.
- `backup_output_files` runs per stage, so each completed step's `OUTCAR`/`INCAR`/`CONTCAR`/`OSZICAR`/`POSCAR` get suffixed with the stage name — strictly more inspection artifacts than today's single stage-level backup, which only captures the last step.
- `CHGCAR` is **not** in the backup list. It stays unsuffixed on disk between stages, so downstream stages (`icharg=2`, `icharg=11`, or `icharg=1` implicit) read the most recent CHGCAR as before.
- `load_structure` reloads atoms from `OUTCAR` (index=-1) between stages. All three workflows in scope use `nsw=0` on every step, so geometry is identical across stages — reloading is a no-op in practice.
- Magmom reference handling (`.asetools_magmom_reference.json`) runs once at the top of `run_workflow`, independent of how stages are cut.

### Correctness check for `freq` mid-Hessian kill

For `ibrion=5`, each finite-difference displacement is written to `OUTCAR` as an ionic step, so `read("OUTCAR", index=-1)` would return a displaced structure rather than the equilibrium. This is only a concern if `load_structure` is called after a Hessian stage is killed — but in `run_workflow`, `load_structure` is only called **between stages** (after a stage completes). On a kill-and-resume, stage 2 is not yet `DONE`, so `run_workflow` enters the loop with the caller-supplied initial `atoms` (equilibrium POSCAR) and stage 1's completion is skipped via the sentinel. VASP writes a fresh POSCAR from `atoms` and restarts ibrion=5 from scratch. Stage 1 is preserved. No regression.

## Scope

### YAML files to refactor

| File | Workflows to refactor |
|---|---|
| `/Users/juar/computers/AOR_NiOH2/asetools.manager/aor.yaml` | `dos`, `bader`, `freq` |
| `asetools/workflow/sample_yaml/aor.yaml` | `dos`, `bader`, `freq` |
| `asetools/workflow/sample_yaml/aor_with_ase_optimizers.yaml` | `freq` |

Untouched in these files: `dos-w-chgcar` (single step, already one stage), `relax*`, `dimer`, `gas`, `sp_lreal_false`, `bulk-relax-keep-shape`.

### Workflows explicitly not refactored

- **`relax*`** — the `single_point + optimization` step pair is a meaningful grouping: shared ENCUT, shared constraints, atomic restart unit. Step 1 is a short warm-start that can be tolerated on resubmit.
- **`dimer`** — `single_point + dimer` shares the VaspInteractive context required by the ASE dimer optimizer. Cannot be split without code changes.

## Proposed stage names

| Workflow | New stage names | Rationale |
|---|---|---|
| `dos` | `DOS_NONSPIN`, `DOS_SPIN_SCF`, `DOS_NONSC` | workflow-prefixed; role-descriptive |
| `bader` | `BADER_NONSPIN`, `BADER_SPIN` | same |
| `freq` | `FREQ_NONSPIN`, `FREQ_HESSIAN` | `HESSIAN` reflects what ibrion=5 actually computes |

Each new stage contains one step with the existing overrides moved verbatim. Example for `dos`:

```yaml
dos:
  stages:
    - name: DOS_NONSPIN
      steps:
        - name: non-spinpolarized
          overrides: { nsw: 0, nelm: 800, ispin: 1, lcharg: true, encut: 500, prec: Accurate, lreal: false, kspacing: 0.25, kpar: 4, npar: 2 }
    - name: DOS_SPIN_SCF
      steps:
        - name: spinpolarized
          overrides: { nsw: 0, nelm: 800, ispin: 2, lcharg: true, encut: 500, laechg: true, prec: Accurate, lreal: false, kspacing: 0.25, icharg: 2, kpar: 4, npar: 2 }
    - name: DOS_NONSC
      steps:
        - name: non-self-consistent
          overrides: { nsw: 0, nelm: 800, ispin: 2, lcharg: false, encut: 500, prec: Accurate, lreal: false, kspacing: 0.25, icharg: 11, nedos: 3001, emin: -20.0, emax: 10.0, kpar: 4, npar: 2 }
```

## Documentation addition

Append a short section to `docs/vasp_calculation_guide.md` (new subsection: *Stages vs. steps: when to split*). Guidance:

> A **stage** is the unit of restart, backup, and constraint application. A **step** is a sub-phase within a stage that shares geometry, constraints, and (for ASE-optimizer workflows) a single VaspInteractive context.
>
> Use separate **stages** when each sub-calculation is independent enough that completing one is worth preserving across job resubmissions — e.g., the non-spin SCF, spin SCF, and non-self-consistent DOS steps.
>
> Use **steps within one stage** when sub-phases genuinely belong together: a warm-start single point followed by an optimization that reads its WAVECAR, or any workflow using an ASE optimizer that needs a single VaspInteractive context across all steps.

## Backwards compatibility

- Existing completed runs have `STAGE_DOS_DONE` / `STAGE_BADER_DONE` / `STAGE_FREQ_DONE` sentinels. Under the refactored YAML these names are no longer referenced, so `stages_to_run` will see the new stages as unfinished. For already-completed directories this is harmless (the user has their output) but would trigger re-runs if the workflow is re-invoked there. Acceptable: these directories don't need re-runs.
- In-flight runs (partially completed old-style stages) will restart from scratch under the new YAML — same as they would today.
- No migration tooling is required.

## Non-goals

- Intra-step restart (WAVECAR/CHGCAR-aware resume of a killed SCF). Addressable separately if step kills become a pain point.
- Changes to `relax*` or `dimer` grouping.
- Any change to `asetools/workflow/manager.py`.

## Acceptance criteria

1. `aor.yaml` (user's working copy, outside the repo) has `dos`, `bader`, `freq` expressed as multi-stage workflows with one step per stage, using the names in the stage-names table.
2. `asetools/workflow/sample_yaml/aor.yaml` and `asetools/workflow/sample_yaml/aor_with_ase_optimizers.yaml` updated to match.
3. `docs/vasp_calculation_guide.md` has a new *Stages vs. steps* subsection.
4. No Python files changed.
5. Each refactored YAML is verified to parse and load cleanly via `VASPConfigurationFromYAML` (importable construction, no KeyError / schema error).
6. A `dry_run=True` invocation of `run_workflow` on each refactored workflow enumerates the expected new stage names in order and produces no errors.

Live cluster smoke test (user-executed, out of scope for the implementer): run the refactored `dos` workflow in a scratch directory, stop after stage 2 completes, resubmit, confirm only stage 3 runs and produces `DOSCAR` with the expected `NEDOS=3001` grid.
