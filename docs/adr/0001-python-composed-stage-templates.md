# 0001 — Python-composed stage templates with a two-door run API

**Status**: accepted (2026-06-23)
**Context spec**: [docs/superpowers/specs/2026-06-23-modular-stage-composition-design.md](../superpowers/specs/2026-06-23-modular-stage-composition-design.md)

## Context

A multi-stage VASP procedure was defined entirely in the YAML `workflows:` section: each named workflow is a fixed, fully-specified list of stages, selected from the submission script by name. Authoring a new procedure means copying a whole named workflow. The result was duplication — `relax`, `relax_low_potim`, `relax_last_step` differ by a single parameter, and `relax_ase_opt`'s 340 → 420 → 500 ENCUT ramp is one stage written three times. The goal: compose the procedure in the Python submission script from reusable, parameterized building blocks, while the canonical final-optimization parameters stay authoritative in YAML.

## Decision

Introduce **stage templates**: a new top-level `stage_templates:` YAML section (parallel to `workflows:`, optional, not in the required-keys check). A `cfg.stage(template, *, name, overrides=, optimizer_kwargs=, constraints=)` factory instantiates a template into a concrete stage dict of the shape the manager already consumes, with each bucket **deep-merged** onto the template defaults (patch only what differs), applied to every step. `name` is required and explicit — it is the restart/sentinel/backup key.

Execution uses **two doors**: a `run_stages(atoms, cfg, stages=[...], production='production')` engine that runs an explicit ordered list, and `run_workflow(workflow_name=...)` reduced to a thin wrapper that looks up the named workflow's stages and delegates to `run_stages(..., production=None)`. The engine warns (non-blocking) if the production stage is not run last.

## Considered options

- **YAML library, Python selects+orders only (no overrides)** — rejected. Every parameter variant still needs its own YAML entry, so the ENCUT ramp and the `relax*` family stay duplicated. No real win over the status quo.
- **Python builds stages inline (factory with all params in the script)** — rejected. Parameter values escape the YAML, contradicting "final/canonical parameters live in YAML" and scattering authority across submission scripts.
- **YAML templates + Python per-instance overrides (chosen)** — one template, instantiated N times with the differing parameter at the call site; production stays canonical in YAML. Kills the duplication while keeping YAML authoritative.
- **One `run_workflow` with mutually-exclusive `workflow_name=` / `stages=` args** — rejected in favor of two named functions. A single function with "pass exactly one of these" is a footgun and obscures which mode a call is in; two verbs make each submission-script call self-documenting. Internally there is still one engine, so no logic is duplicated.
- **Auto-generated instance names (hash of overrides)** — rejected. The name is the restart key; an auto name silently changes when a parameter is tweaked on resubmit, re-running completed stages. Explicit required `name=` keeps the restart unit stable and visible.

## Consequences

- A second way to run a workflow now exists. The named-workflow path is unchanged and fully backward-compatible; old calculation directories keep resuming. The `relax*` family is documented as superseded and will be deprecated once legacy projects finish.
- The magmom / atom-reorder / reference-file machinery moves from `run_workflow` into `run_stages` so both paths share it. This touches the resume path and is covered by tests.
- Per-instance overrides are **unrestricted** for now — they can touch any VASP tag, including ones (`ediffg`, constraints) that a careful user would not vary on a pre-opt stage. A declared safe-set was deferred. The production-stage warning is the only guardrail, and it only catches a missing/misplaced production stage, not a contaminated-but-present one (`run_overrides` leak, R1).
- Precedence is fixed as `basic → system → run_overrides → template step → per-instance`. `run_overrides` keeps its existing (below-step) position for compatibility, which means it is *not* a global override hammer.

## Reversibility

Moderately hard. `cfg.stage`, `run_stages`, and the `stage_templates:` schema become public surface that submission scripts and saved YAML depend on. Removing them later would break those scripts. This is why it is recorded.
