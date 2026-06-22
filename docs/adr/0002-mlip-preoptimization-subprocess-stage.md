# 0002 — MLIP pre-optimization as a subprocess stage

**Status**: accepted (2026-06-23)
**Context spec**: [docs/superpowers/specs/2026-06-23-modular-stage-composition-design.md](../superpowers/specs/2026-06-23-modular-stage-composition-design.md)
**Related**: `mlip_platform` ADR 0001 (per-MLIP envs)

## Context

We want a cheap MLIP relaxation as a pre-optimization stage ahead of the DFT stages, driven by the `mlip_platform` package. The blocker: `mlip_platform` ADR 0001 establishes that MLIP packages (`mace-torch`, `fairchem-core`, `sevenn`, `chgnet`) pin mutually-incompatible `torch`/`e3nn` stacks and must live one-per-env, and that `setup_calculator()` does in-process ASE imports assuming the correct MLIP env is already active. The workflow manager runs inside a VASP/asetools job env that cannot also host an MLIP torch stack without reintroducing exactly the corruption ADR 0001 prevents.

## Decision

An MLIP stage (`engine: mlip` in a stage template) executes as a **subprocess in a separate MLIP env**, not as an in-process import. The manager hands the subprocess the current structure file plus the stage's MLIP and constraint parameters; a new `asetools.workflow.mlip_runner` entry point — importable in the MLIP env because `asetools` is a declared dependency of `mlip_platform` — reads the structure, re-applies the stage's constraints from the same JSON via `ConstraintManager`, attaches the calculator via `setup_calculator`, and runs `mlip_platform.core.optimize.run_optimization`. The env is resolved through a `globals.mlip_envs:` registry (env key → interpreter path).

The runner writes a `CONTCAR`. **MLIP stages are pre-optimization only** — every MLIP stage must precede every VASP stage — and `run_stages` enforces this with a hard error at validation time.

## Considered options

- **In-process import of `run_optimization`** — rejected. Requires the MLIP torch stack inside the VASP job env, violating `mlip_platform` ADR 0001 and corrupting the VASP env. Only viable on a purpose-built combined env, which we do not have on the target machines.
- **Subprocess via the stock `optimize run` CLI** — rejected for constrained stages. The CLI only `read()`s a structure and optimizes; Hookean springs have no POSCAR representation and would be silently dropped. A constraint-preserving serialization (`.traj`/`.json`) was considered as a fix but is more fragile than re-applying from the JSON.
- **Subprocess via a new `asetools` MLIP runner that re-applies constraints from JSON (chosen)** — robust. `run_optimization` honors `atoms.constraints` automatically, and the same `ConstraintManager` code runs in the MLIP env, so Hookean pre-optimization works without serializing the constraint object.
- **Per-template conda name / interpreter path instead of a `globals` registry** — rejected. Interpreter paths differ per machine (laptop / ASPIRE2A / cos-cluster); a registry keeps templates portable so only `globals.mlip_envs` changes when relocating.
- **Allow MLIP stages anywhere; make the manager track each stage's output structure explicitly** — rejected as unnecessary. MLIP is a warm-start; there is no real scenario where an MLIP stage follows a DFT stage. Constraining MLIP to pre-opt-only means a stale `OUTCAR` can never precede an MLIP stage, so `load_structure` needs no change and the runner writing a `CONTCAR` is sufficient. The ordering guard makes the invariant safe.

## Consequences

- Adding MLIP support touches three places: the `engine: mlip` branch in `run_stages` (build + launch the subprocess, check exit code), the new `asetools.workflow.mlip_runner`, and the `globals.mlip_envs:` schema.
- The structure crosses the process boundary as a file; the constraint *config* crosses as arguments. `FixAtoms` survives via selective dynamics in the structure file; Hookean is rebuilt from JSON in the runner.
- `load_structure` is unchanged. Correctness depends on the pre-opt-only invariant, which is enforced by a hard error (an `engine: mlip` stage after an `engine: vasp` stage). If that invariant is ever relaxed, the stale-`OUTCAR` handoff problem returns and explicit per-stage output tracking becomes necessary.
- Non-convergence is signaled by a non-zero subprocess exit code; the manager refuses to write `STAGE_{name}_DONE`, matching VASP-stage behavior.
- The user maintains the per-MLIP envs and their interpreter paths; the platform does not abstract env management (consistent with `mlip_platform` ADR 0001).

## Reversibility

Hard. The `engine: mlip` template schema, the `globals.mlip_envs:` registry, the runner CLI contract, and the pre-opt-only invariant all become things saved configs and submission scripts depend on. Recorded for that reason.
