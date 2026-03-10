# Constraint Management Integration - Migration Guide

This guide shows how to migrate from manual constraint setup to the new integrated constraint management system.

## Overview

The new `asetools.constraints` module provides:
- **Simplified workflow**: Constraint configuration in YAML instead of Python code
- **JSON-based constraint definitions**: Easy to version control and modify
- **Automatic constraint merging**: Preserves existing FixAtoms constraints
- **Consistent API**: Same approach across all calculations

---

## Before & After Comparison

### OLD APPROACH (Manual in each script)

```python
#!/usr/bin/env python

from ase.constraints import FixAtoms, Hookean
from ase.data import covalent_radii
from asetools.manager.manager import load_structure, run_workflow
from asetools.manager.calculatorsetuptools import VASPConfigurationFromYAML
import json

# Load structure
cfg = VASPConfigurationFromYAML(config_file='config.yaml', system='CuNiOOH')
atoms = load_structure('POSCAR')

# ---- BEGIN 65 LINES OF BOILERPLATE CODE ----
# Load mapping
with open('proton_map_O_vac.json') as fh:
    mapping = json.load(fh)

# Build Hookean constraints
DISTANCE_FACTOR = 1.134
hooks = []

for idx1, idx2 in mapping.get('pairs', []):
    i = int(idx1)
    j = int(idx2)

    # Get atomic symbols and radii
    symbol1 = atoms[i].symbol
    symbol2 = atoms[j].symbol
    radius1 = covalent_radii[atoms[i].number]
    radius2 = covalent_radii[atoms[j].number]

    # Calculate r0
    r0 = (radius1 + radius2) * DISTANCE_FACTOR
    k = 20.0

    print(f"Constraint: {symbol1}({i})—{symbol2}({j}), r0={r0:.3f} Å")
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
    raise RuntimeError(f"Failed to apply constraints: {e}")
# ---- END BOILERPLATE CODE ----

# Run workflow
run_workflow(atoms, cfg=cfg, workflow_name='relax_ase_opt')
```

**Problems:**
- 65+ lines of repetitive code in every script
- Easy to make mistakes when copying
- Hard to maintain consistency across calculations
- Difficult to modify constraint parameters

---

### NEW APPROACH (Integrated with YAML)

#### Option 1: YAML Configuration (Recommended)

**Python script** (simplified to ~10 lines):
```python
#!/usr/bin/env python

from asetools.manager.manager import load_structure, run_workflow
from asetools.manager.calculatorsetuptools import VASPConfigurationFromYAML
from asetools.manager.logger import configure_logging

configure_logging(file_prefix="run")
cfg = VASPConfigurationFromYAML(config_file='config.yaml', system='CuNiOOH')
atoms = load_structure('POSCAR')

# Constraints are automatically applied from YAML!
run_workflow(atoms, cfg=cfg, workflow_name='relax_with_constraints')
```

**YAML configuration** (`config.yaml`):
```yaml
workflows:
  relax_with_constraints:
    stages:
      - name: OPT_WITH_HOOKEAN
        constraints:
          type: hookean
          config_file: proton_map_O_vac.json
          spring_constant: 20.0
          distance_factor: 1.134
        steps:
          - name: optimization
            overrides: { nsw: 0 }
            optimizer: FIRE
            optimizer_kwargs:
              fmax: 0.02
```

**JSON file** (`proton_map_O_vac.json` - unchanged):
```json
{
  "pairs": [
    [114, 179],
    [116, 180],
    [112, 175]
  ]
}
```

---

#### Option 2: Direct Python API

For cases where you need more control:

```python
from asetools.constraints import ConstraintManager

# Simple one-liner
cm = ConstraintManager()
cm.apply_from_json(atoms, 'proton_map_O_vac.json', k=20.0)

# Or with custom parameters
cm = ConstraintManager(distance_factor=1.2)
cm.apply_from_json(atoms, 'proton_map_O_vac.json', k=25.0)
```

---

## Migration Steps

### Step 1: No Changes Needed to JSON Files

Your existing JSON files work as-is:
```json
{
  "pairs": [[114, 179], [116, 180], [112, 175], ...]
}
```

Optional: Add metadata for documentation:
```json
{
  "pairs": [[114, 179], [116, 180], ...],
  "metadata": {
    "description": "O-H proton mappings for NiOOH surface",
    "spring_constant": 20.0,
    "distance_factor": 1.134
  }
}
```

### Step 2: Add Constraints Section to YAML

In your workflow YAML, add a `constraints` section to the stage:

```yaml
workflows:
  your_workflow_name:
    stages:
      - name: YOUR_STAGE
        constraints:                        # NEW
          type: hookean                     # NEW
          config_file: proton_map.json      # NEW
          spring_constant: 20.0             # NEW
          distance_factor: 1.134            # NEW (optional)
        steps:
          - name: optimization
            # ... existing step configuration ...
```

### Step 3: Simplify Python Script

Remove all the constraint boilerplate code. Your script becomes:

```python
#!/usr/bin/env python

from asetools.manager.manager import load_structure, run_workflow
from asetools.manager.calculatorsetuptools import VASPConfigurationFromYAML
from asetools.manager.logger import configure_logging

configure_logging(file_prefix="run")
cfg = VASPConfigurationFromYAML(config_file='config.yaml', system='CuNiOOH')
atoms = load_structure('POSCAR')

run_workflow(atoms, cfg=cfg, workflow_name='your_workflow_name')
```

That's it!

---

## Advanced Features

### Different Constraints Per Stage

```yaml
workflows:
  staged_optimization:
    stages:
      # Loose constraints for initial relaxation
      - name: COARSE
        constraints:
          type: hookean
          config_file: constraints.json
          spring_constant: 10.0  # Weak spring
        steps:
          - name: coarse_opt
            optimizer: FIRE
            optimizer_kwargs: { fmax: 0.05 }

      # Tight constraints for final optimization
      - name: FINE
        constraints:
          type: hookean
          config_file: constraints.json
          spring_constant: 30.0  # Strong spring
        steps:
          - name: fine_opt
            optimizer: BFGS
            optimizer_kwargs: { fmax: 0.01 }
```

### Scanning Spring Constants

```yaml
workflows:
  scan_k_values:
    stages:
      - name: OPT_K10
        constraints: { type: hookean, config_file: pairs.json, spring_constant: 10.0 }
        steps: [...]

      - name: OPT_K20
        constraints: { type: hookean, config_file: pairs.json, spring_constant: 20.0 }
        steps: [...]

      - name: OPT_K30
        constraints: { type: hookean, config_file: pairs.json, spring_constant: 30.0 }
        steps: [...]
```

### Direct Python API for Complex Cases

If you need constraint logic not covered by JSON:

```python
from asetools.constraints import ConstraintManager

cm = ConstraintManager(distance_factor=1.2)

# Apply from JSON
cm.apply_from_json(atoms, 'constraints.json', k=25.0)

# Or build pairs programmatically
pairs = identify_proton_pairs(atoms)  # Your custom function
hookean_constraints = cm.apply_hookean_from_pairs(atoms, pairs, k=20.0)
cm.merge_constraints(atoms, hookean_constraints)
```

---

## Benefits Summary

| Aspect | OLD Approach | NEW Approach |
|--------|-------------|-------------|
| **Lines of code** | 65+ per script | ~10 per script |
| **Configuration** | Hard-coded | YAML-based |
| **Consistency** | Manual | Automatic |
| **Maintainability** | Low | High |
| **Error-prone** | Yes | No |
| **Version control** | Difficult | Easy |
| **Testing** | None | Comprehensive |

---

## Troubleshooting

### "Missing 'config_file' in constraint configuration"
**Solution:** Add `config_file` to your constraints section in YAML

### "File not found: constraints.json"
**Solution:** Ensure the JSON file path is relative to your working directory, or use absolute path

### "Unsupported constraint type"
**Solution:** Currently only `type: hookean` is supported. Other types coming in Phase 2.

### Constraints not being applied
**Solution:** Check logs for constraint application messages. Verify `constraints` is at the stage level, not step level.

---

## Example: Complete Migration

See `asetools/manager/sample_yaml/aor_with_hookean_constraints.yaml` for a complete working example with:
- Multiple workflow patterns
- Different constraint strengths
- Spring constant scanning
- Detailed comments and documentation

---

## Questions?

- Check the comprehensive YAML example: `asetools/manager/sample_yaml/aor_with_hookean_constraints.yaml`
- Run integration tests: `python test_constraint_integration.py`
- Read module docstrings: `asetools/constraints.py`
- Check test cases: `asetools/tests/test_constraints.py`
