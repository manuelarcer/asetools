# 0003 — run_overrides is authoritative (beats step overrides)

**Status**: accepted (2026-06-24)
**Context spec**: [docs/superpowers/specs/2026-06-23-modular-stage-composition-design.md](../superpowers/specs/2026-06-23-modular-stage-composition-design.md)
**Supersedes**: ADR 0001, precedence decision (final bullet of its Consequences)

## Context

ADR 0001 fixed the parameter precedence as `basic → system → run_overrides → template step → per-instance`, deliberately keeping `run_overrides` *below* step-level values "for backward compatibility", which means it was explicitly *not* a global override hammer.

In practice this surprised users. A run-wide override is the obvious tool for forcing a single tag across an entire procedure (e.g. `run_workflow(atoms, cfg, "cellopt", run_overrides={"kpar": 1})` on a node where the YAML step pins `kpar: 16`), but under the old precedence the per-step override silently won and the run-wide value was ignored. The log even showed both — `run_overrides` applied at calculator build, then the step override applied after — making the loss invisible unless the resulting `INCAR` was inspected. There is no natural use for "run-wide override that loses to a step default": if a caller passes `run_overrides`, they mean it.

## Decision

`run_overrides` is the **final word**. It is re-applied *after* the per-step `overrides` in both execution paths, so it beats any per-step value for any VASP tag.

New precedence (low to high):

```
basic → system → template step overrides → per-instance stage overrides → run_overrides
```

Implementation in `asetools/workflow/manager.py`:

- **Regular VASP path** (`_run_step`): after `atoms.calc.set(**overrides)`, re-apply `atoms.calc.set(**run_overrides)` when `run_overrides` is non-empty. `run_overrides` is threaded into `_run_step` from `_run_stage`.
- **VaspInteractive path** (`_run_stage_with_vaspinteractive`): the per-step parameter set is built by a new `_layer_step_params(cfg, run_overrides, overrides)` helper that ends with `params.update(run_overrides or {})`.

Both paths are covered by precedence unit tests in `tests/test_manager.py` (`TestRunOverridesPrecedence`).

## Considered options

- **Keep ADR 0001 precedence (run_overrides below step)** — rejected. It is the documented cause of the bug; users cannot force a tag run-wide without editing every step.
- **Config-level workaround only** (move the contested tag from the step into `basic` so `run_overrides` already wins) — already used live for the `kpar` case in a project YAML, but it is per-tag and per-project. It does not fix the general semantics, so it is not a substitute for the engine change.
- **A declared "global" override channel separate from `run_overrides`** — rejected as over-engineering. `run_overrides` is already the run-wide channel; giving it authoritative precedence is the least surprising behavior.

## Consequences

- **Behavior change for a shared package.** `asetools` is used by other projects (e.g. Battery, AOR). Any workflow that relied — knowingly or not — on a per-step override beating a passed `run_overrides` changes outcome. Audit: such reliance is implausible (it requires passing `run_overrides` *and* wanting it ignored), but it is a real semantic change and is recorded here for that reason.
- The R1 "`run_overrides` contamination of a present production stage" concern from ADR 0001 / the spec is now sharper: `run_overrides` overrides the production stage's own values too. That remains the caller's responsibility; the production-stage warning is unchanged and still does not inspect contents.
- Logs already announce `Overriding run parameters with: {...}`; the applied value now matches what lands in `INCAR`.

## Reversibility

Easy in code (drop the trailing re-apply / `_layer_step_params` final update), but reverting would re-break the user expectation this fixes and would itself be a behavior change for any config written after this ADR. Treat the new precedence as the contract.
