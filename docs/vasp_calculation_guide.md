# VASP/DFT Calculation Guide

Personal reference for computational materials science calculations and lessons learned.

## Purpose

This document serves as a personal knowledge base for VASP calculations, collecting tips, best practices, and solutions discovered through experience. It will be expanded over time as new techniques and insights are gained.

The focus is on practical, actionable advice that has proven useful in real calculations, particularly for challenging systems and convergence issues.

## Spin-Polarized Calculations with Magnetic Materials

Description: Guidelines and tips for handling magnetic systems in VASP, particularly for challenging convergence cases and DFT+U calculations.

### Key Challenges
- Achieving convergence to desired magnetic moments in atoms
- Handling strongly correlated systems with DFT+U
- Dealing with magnetic frustration and multiple magnetic states
- Convergence issues in complex magnetic structures

### Tips and Best Practices

#### Two-Step Convergence Approach
For difficult magnetic systems, especially with DFT+U, use a staged approach:

1. **First step**: Non-spin-polarized single-point energy calculation
   ```bash
   # INCAR settings for step 1
   NSW = 0          # Single point calculation
   NELM = 800       # High number of electronic steps
   ISPIN = 1        # Non-spin-polarized
   LCHARG = .TRUE.  # Save charge density to CHGCAR
   ```

2. **Second step**: Spin-polarized calculation with initial magnetic moment
   ```bash
   # INCAR settings for step 2
   ICHARG = 1       # Read charge density from CHGCAR
   ISPIN = 2        # Spin-polarized
   MAGMOM = ...     # Set initial magnetic moments
   # Include DFT+U parameters if needed
   LDAU = .TRUE.
   LDAUTYPE = 2
   LDAUL = ...
   LDAUU = ...
   ```

**Why this works**: The first step provides a good starting charge density without the complexity of magnetic degrees of freedom, making the second step more likely to converge to the desired magnetic state.

#### Convergence Parameter Tuning
When standard SCF convergence fails, adjust magnetic mixing parameters as recommended in the [VASP wiki](https://www.vasp.at/wiki/index.php/AMIX_MAG):

```bash
# Conservative mixing for difficult magnetic systems
AMIX = 0.2        # Linear mixing parameter
BMIX = 0.0001     # Cutoff wave vector for Kerker mixing
AMIX_MAG = 0.8    # Linear mixing parameter for magnetization
BMIX_MAG = 0.0001 # Cutoff wave vector for magnetic mixing
```

**When to use**: 
- Oscillating magnetic moments during SCF
- Convergence failures in strongly correlated systems
- Systems with competing magnetic states

#### Initial Magnetic Moment Specification
Proper MAGMOM initialization is crucial for convergence:

```bash
# Example for a system with Fe (4 μB) and O (0 μB) atoms
MAGMOM = 4*4.0 8*0.0  # 4 Fe atoms, 8 O atoms

# For mixed valence systems, specify per atom
MAGMOM = 4.5 3.5 0.0 0.0  # Different Fe oxidation states
```

**Tips**:
- Start with reasonable magnetic moments based on expected oxidation states
- For transition metals, use typical values (Fe: 4-5 μB, Co: 3 μB, Ni: 2 μB)
- Set non-magnetic atoms (O, C, etc.) to 0.0
- For unknown systems, try different initial values if convergence fails

### Common Issues and Solutions

| Problem | Solution |
|---------|----------|
| Oscillating magnetic moments | Reduce AMIX and AMIX_MAG |
| Wrong magnetic ground state | Try different MAGMOM initialization |
| SCF not converging | Increase NELM, use two-step approach |
| Unexpected magnetic moments | Check for charge transfer, verify DFT+U parameters |

### Additional Notes

- Always check final magnetic moments in OUTCAR to verify they match expectations
- For antiferromagnetic systems, initialize with alternating spin directions
- Consider using NUPDOWN for systems with known net magnetic moment
- Monitor convergence in OSZICAR for both energy and magnetic moments

---

## Modular stage composition

Two ways to run a multi-stage procedure:

1. **Named workflow** (original) — the whole stage list is defined in the YAML
   `workflows:` section and selected by name:
   ```python
   run_workflow(atoms, cfg, workflow_name="relax_ase_opt")
   ```
2. **Modular composition** (new) — reusable `stage_templates:` are instantiated
   and ordered in the Python submission script, while the canonical final
   ("production") parameters stay authoritative in the YAML:
   ```python
   from asetools.workflow.manager import run_stages

   PREOPTIMIZATION = [
       cfg.stage("encut_ramp", name="OPT_340",
                 overrides={"encut": 340}, optimizer_kwargs={"fmax": 0.05}),
       cfg.stage("encut_ramp", name="OPT_420",
                 overrides={"encut": 420}, optimizer_kwargs={"fmax": 0.02}),
   ]
   PRODUCTION = cfg.stage("production", name="OPT_500")
   run_stages(atoms, cfg, stages=PREOPTIMIZATION + [PRODUCTION])
   ```
   One template instantiated N times replaces N near-duplicate named workflows.

Both paths share the same engine, restart sentinels (`STAGE_{name}_DONE`),
backups, and constraint handling. The named-workflow path is unchanged; existing
configs keep working.

See `asetools/workflow/sample_yaml/modular_stage_composition.yaml` and
`submit_modular_example.py` for a complete example, and the design spec at
`docs/superpowers/specs/2026-06-23-modular-stage-composition-design.md`.

### `cfg.stage()` — instantiating a template

`cfg.stage(template, *, name, overrides=None, optimizer_kwargs=None, constraints=None)`
returns a resolved stage. `name` is **required** and becomes the restart/backup
key. `overrides`, `optimizer_kwargs`, and `constraints` are **deep-merged** onto
the template (patch only what differs) and applied to every step. Precedence,
lowest to highest:

```
basic -> system -> run_overrides -> template step values -> per-instance values
```

### Production-stage warning

`run_stages` logs a non-blocking warning if the stage from the template named
`production` (override via `production=`) is not the last stage run — catching a
forgotten or misplaced final optimization. Pass `production=None` to disable it
(e.g. a pure-DOS modular pipeline). The named-workflow path disables it
automatically.

### Stages vs. steps

A **stage** is the unit of restart, backup, and constraint application (one
`STAGE_{name}_DONE` sentinel). A **step** is a sub-phase within a stage that
shares geometry, constraints, and (for ASE-optimizer workflows) one
VaspInteractive context — e.g. a warm-start single point followed by an
optimization that reads its WAVECAR. Use separate stages when completing one is
worth preserving across resubmissions; use steps when sub-phases genuinely
belong together.

## MLIP pre-optimization stages

A stage template with `engine: mlip` runs a machine-learning interatomic
potential (via the `mlip_platform` package) as a cheap warm-start before the DFT
stages:

```yaml
stage_templates:
  uma_preopt:
    engine: mlip
    mlip: uma-s-1p2          # MLIP tag for setup_calculator
    env: uma                 # key into globals.mlip_envs
    optimizer: bfgs
    fmax: 0.10
    max_steps: 300
    device: auto
    constraints:             # re-applied inside the MLIP env from the same JSON
      type: hookean
      config_file: hookean_c_pairs.json
      spring_constant: 20.0

globals:
  mlip_envs:
    uma:  /path/to/uma-env/bin/python
    mace: /path/to/mace-env/bin/python
```

Because MLIP packages pin mutually incompatible torch stacks (see `mlip_platform`
ADR 0001), an MLIP stage **cannot** run in the VASP job env. It executes as a
**subprocess** in the MLIP env named by `env:` (resolved through
`globals.mlip_envs`), invoking `asetools.workflow.mlip_runner`. That runner reads
the current structure, re-applies the stage's constraints from the same JSON the
VASP stages use, attaches the MLIP calculator, optimizes, and writes a `CONTCAR`
the next VASP stage picks up.

**Invariant — MLIP stages are pre-optimization only:** every MLIP stage must
precede every VASP stage. `run_stages` raises a hard error at validation time if
an MLIP stage is composed after a VASP stage (an MLIP stage produces no OUTCAR,
so a later MLIP relaxation would be silently discarded on the next structure
load). A non-converged MLIP run exits non-zero and no `STAGE_*_DONE` is written.

---

*This section will be expanded with additional calculation types and techniques as they are encountered and mastered.*