# Hookean Constraints - Quick Start Guide

## What is This?

The new constraint system eliminates repetitive boilerplate code for applying Hookean spring constraints in VASP calculations. Perfect for proton-coupled electron transfer (PCET) reactions and maintaining O-H bond distances.

## Quick Start (2 minutes)

### 1. Your JSON file (unchanged)

`proton_mappings.json`:
```json
{
  "pairs": [
    [114, 179],
    [116, 180],
    [112, 175]
  ]
}
```

### 2. Add to your YAML config

```yaml
workflows:
  relax_with_constraints:
    stages:
      - name: OPT_WITH_HOOKEAN
        constraints:
          type: hookean
          config_file: proton_mappings.json
          spring_constant: 20.0
        steps:
          - name: optimization
            overrides: { nsw: 0 }
            optimizer: FIRE
            optimizer_kwargs: { fmax: 0.02 }
```

### 3. Simplified Python script

```python
#!/usr/bin/env python

from asetools.manager.manager import load_structure, run_workflow
from asetools.manager.calculatorsetuptools import VASPConfigurationFromYAML
from asetools.manager.logger import configure_logging

configure_logging(file_prefix="run")
cfg = VASPConfigurationFromYAML(config_file='config.yaml', system='CuNiOOH')
atoms = load_structure('POSCAR')

run_workflow(atoms, cfg=cfg, workflow_name='relax_with_constraints')
```

**That's it!** No more 65 lines of constraint boilerplate.

---

## What Happened to Your Old Code?

**Before** (65+ lines per script):
```python
# Load mapping
with open('proton_map.json') as fh:
    mapping = json.load(fh)

# Build Hookean constraints
DISTANCE_FACTOR = 1.134
hooks = []
for idx1, idx2 in mapping.get('pairs', []):
    i = int(idx1)
    j = int(idx2)
    symbol1 = atoms[i].symbol
    symbol2 = atoms[j].symbol
    radius1 = covalent_radii[atoms[i].number]
    radius2 = covalent_radii[atoms[j].number]
    r0 = (radius1 + radius2) * DISTANCE_FACTOR
    k = 20.0
    hooks.append((i, j, k, r0))

# Get existing FixAtoms
existing_constraints = atoms.constraints
fix_indices = []
if existing_constraints:
    if isinstance(existing_constraints, FixAtoms):
        fix_indices = list(existing_constraints.index)
    elif isinstance(existing_constraints, (list, tuple)):
        for c in existing_constraints:
            if isinstance(c, FixAtoms):
                fix_indices.extend(c.index)

# Merge constraints
try:
    new_constraints = []
    if fix_indices:
       new_constraints.append(FixAtoms(indices=sorted(set(fix_indices))))
    if not hooks:
       raise RuntimeError("No atom pairs found")
    for (i, j, k, r0) in hooks:
       new_constraints.append(Hookean(a1=i, a2=j, k=k, rt=r0))
    atoms.set_constraint(new_constraints)
except Exception as e:
    raise RuntimeError(f"Failed: {e}")
```

**After** (in YAML config):
```yaml
constraints:
  type: hookean
  config_file: proton_map.json
  spring_constant: 20.0
```

---

## Key Features

✅ **Automatic r0 calculation** from covalent radii
✅ **Preserves FixAtoms** constraints (bottom layer stays fixed)
✅ **YAML-driven** configuration (no Python code changes)
✅ **Stage-specific** constraint parameters
✅ **Comprehensive logging** of all constraint applications
✅ **Error handling** with clear messages

---

## Common Use Cases

### Different constraint strengths per stage

```yaml
workflows:
  staged_optimization:
    stages:
      - name: COARSE
        constraints: { type: hookean, config_file: pairs.json, spring_constant: 10.0 }
        steps: [...]

      - name: FINE
        constraints: { type: hookean, config_file: pairs.json, spring_constant: 30.0 }
        steps: [...]
```

### Direct Python API (if needed)

```python
from asetools.constraints import ConstraintManager

cm = ConstraintManager()
cm.apply_from_json(atoms, 'proton_mappings.json', k=20.0)
```

---

## Parameters

- **spring_constant** (k): Force constant in eV/Å² (typical: 10-50)
- **distance_factor**: Multiplier for covalent radii sum (default: 1.134 for O-H)
- **config_file**: Path to JSON file with atom pairs

---

## Documentation

- **Full migration guide**: `CONSTRAINTS_MIGRATION_GUIDE.md`
- **Complete YAML example**: `asetools/manager/sample_yaml/aor_with_hookean_constraints.yaml`
- **Module documentation**: `asetools/constraints.py`
- **Test examples**: `asetools/tests/test_constraints.py`

---

## Verification

Run integration tests:
```bash
python test_constraint_integration.py
```

All 81 tests should pass.

---

## Need Help?

1. Check logs for constraint application messages
2. Verify JSON file format: `{"pairs": [[idx1, idx2], ...]}`
3. Ensure `constraints` is at stage level (not step level)
4. Check atom indices are valid for your structure

---

**Status**: ✅ Production ready, all tests passing
