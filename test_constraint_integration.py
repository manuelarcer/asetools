#!/usr/bin/env python
"""
Integration test for constraint support in asetools workflow manager.

This script demonstrates the complete workflow with Hookean constraints.
"""

import json
import tempfile
from pathlib import Path

from ase.build import molecule
from ase.constraints import FixAtoms
from asetools.constraints import ConstraintManager

def test_basic_constraint_application():
    """Test 1: Basic constraint application from JSON"""
    print("=" * 60)
    print("Test 1: Basic Constraint Application")
    print("=" * 60)

    # Create test molecule
    atoms = molecule('H2O2')
    print(f"\nCreated H2O2 molecule with {len(atoms)} atoms")
    print(f"Atoms: {atoms.get_chemical_symbols()}")

    # Create temporary JSON file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        config = {
            "pairs": [[0, 1], [2, 3]],  # O-H pairs
            "metadata": {
                "description": "O-H bonds in H2O2",
                "spring_constant": 20.0
            }
        }
        json.dump(config, f)
        json_file = f.name

    # Apply constraints
    cm = ConstraintManager()
    cm.apply_from_json(atoms, json_file, k=20.0)

    print(f"\nConstraints applied: {len(atoms.constraints)}")
    for i, constraint in enumerate(atoms.constraints):
        print(f"  {i+1}. {type(constraint).__name__}")

    # Cleanup
    Path(json_file).unlink()
    print("\n✅ Test 1 PASSED\n")

def test_constraint_with_fixatoms():
    """Test 2: Constraints with existing FixAtoms"""
    print("=" * 60)
    print("Test 2: Constraints with Existing FixAtoms")
    print("=" * 60)

    # Create test molecule
    atoms = molecule('H2O2')

    # Add FixAtoms constraint (fix first atom)
    atoms.set_constraint(FixAtoms(indices=[0]))
    print(f"\nAdded FixAtoms constraint for atom 0")

    # Create temporary JSON file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        config = {"pairs": [[1, 2], [2, 3]]}
        json.dump(config, f)
        json_file = f.name

    # Apply Hookean constraints
    cm = ConstraintManager()
    cm.apply_from_json(atoms, json_file, k=15.0)

    print(f"\nTotal constraints: {len(atoms.constraints)}")
    for i, constraint in enumerate(atoms.constraints):
        constraint_type = type(constraint).__name__
        print(f"  {i+1}. {constraint_type}")
        if constraint_type == 'FixAtoms':
            print(f"     Fixed atoms: {list(constraint.index)}")

    # Verify FixAtoms is preserved
    fixatoms_count = sum(1 for c in atoms.constraints if isinstance(c, FixAtoms))
    assert fixatoms_count == 1, "FixAtoms constraint should be preserved"

    # Cleanup
    Path(json_file).unlink()
    print("\n✅ Test 2 PASSED\n")

def test_stage_constraint_config():
    """Test 3: Stage constraint configuration format"""
    print("=" * 60)
    print("Test 3: Stage Constraint Configuration")
    print("=" * 60)

    atoms = molecule('H2O')

    # Create temporary JSON file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        config = {"pairs": [[0, 1], [0, 2]]}  # O bonded to both H
        json.dump(config, f)
        json_file = f.name

    # Simulate stage configuration from YAML
    stage_config = {
        'type': 'hookean',
        'config_file': json_file,
        'spring_constant': 25.0,
        'distance_factor': 1.15
    }

    print(f"\nStage configuration:")
    for key, value in stage_config.items():
        print(f"  {key}: {value}")

    # Apply constraints as would happen in workflow
    cm = ConstraintManager()
    cm.apply_stage_constraints(atoms, stage_config)

    print(f"\nConstraints applied: {len(atoms.constraints)}")
    assert len(atoms.constraints) == 2, "Should have 2 Hookean constraints"

    # Cleanup
    Path(json_file).unlink()
    print("\n✅ Test 3 PASSED\n")

def test_covalent_radii_calculation():
    """Test 4: Bond distance calculation from covalent radii"""
    print("=" * 60)
    print("Test 4: Covalent Radii Bond Distance Calculation")
    print("=" * 60)

    atoms = molecule('H2O')
    cm = ConstraintManager(distance_factor=1.134)

    # Calculate O-H bond distance
    r0 = cm.calculate_bond_distance(atoms, 0, 1)  # O-H

    print(f"\nH2O molecule:")
    print(f"  Atom 0: {atoms[0].symbol}")
    print(f"  Atom 1: {atoms[1].symbol}")
    print(f"  Calculated r0: {r0:.3f} Å")
    print(f"  Distance factor: {cm.distance_factor}")

    # Typical O-H bond length is ~0.96 Å
    # With factor 1.134, we expect ~1.10 Å
    assert 1.0 < r0 < 1.2, f"O-H bond distance {r0:.3f} Å seems incorrect"

    print("\n✅ Test 4 PASSED\n")

def test_real_world_scenario():
    """Test 5: Real-world scenario similar to user's AOR case"""
    print("=" * 60)
    print("Test 5: Real-World AOR-like Scenario")
    print("=" * 60)

    # Simulate a larger system (like NiOOH surface)
    from ase.build import bulk
    from ase import Atom

    # Create simple slab
    slab = bulk('Ni', 'fcc', a=3.52).repeat((3, 3, 4))

    # Fix bottom 2 layers
    layer_height = 3.52 / 2
    fix_indices = [i for i in range(len(slab)) if slab[i].position[2] < 2 * layer_height]
    slab.set_constraint(FixAtoms(indices=fix_indices))

    print(f"\nCreated Ni slab:")
    print(f"  Total atoms: {len(slab)}")
    print(f"  Fixed atoms: {len(fix_indices)} (bottom layers)")

    # Add some O and H atoms on top (simulate adsorbates)
    slab.append(Atom('O', position=[1.76, 1.76, 14.0]))
    slab.append(Atom('H', position=[1.76, 1.76, 15.0]))
    slab.append(Atom('O', position=[5.28, 1.76, 14.0]))
    slab.append(Atom('H', position=[5.28, 1.76, 15.0]))

    o1_idx = len(slab) - 4
    h1_idx = len(slab) - 3
    o2_idx = len(slab) - 2
    h2_idx = len(slab) - 1

    print(f"  Added adsorbates:")
    print(f"    O at index {o1_idx}")
    print(f"    H at index {h1_idx}")
    print(f"    O at index {o2_idx}")
    print(f"    H at index {h2_idx}")

    # Create constraint file for O-H pairs
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        config = {
            "pairs": [[o1_idx, h1_idx], [o2_idx, h2_idx]],
            "metadata": {
                "description": "O-H pairs on NiOOH surface",
                "spring_constant": 20.0,
                "distance_factor": 1.134
            }
        }
        json.dump(config, f)
        json_file = f.name

    # Apply constraints (preserving FixAtoms)
    cm = ConstraintManager()
    cm.apply_from_json(slab, json_file, k=20.0)

    print(f"\nFinal constraint setup:")
    print(f"  Total constraints: {len(slab.constraints)}")

    fixatoms_found = False
    hookean_count = 0
    for constraint in slab.constraints:
        if isinstance(constraint, FixAtoms):
            fixatoms_found = True
            print(f"  - FixAtoms: {len(constraint.index)} atoms")
        else:
            hookean_count += 1

    print(f"  - Hookean: {hookean_count} constraints")

    assert fixatoms_found, "FixAtoms should be preserved"
    assert hookean_count == 2, "Should have 2 Hookean constraints"

    # Cleanup
    Path(json_file).unlink()
    print("\n✅ Test 5 PASSED - Real-world scenario works!\n")

def main():
    """Run all integration tests"""
    print("\n" + "=" * 60)
    print("ASETOOLS CONSTRAINT INTEGRATION TESTS")
    print("=" * 60 + "\n")

    try:
        test_basic_constraint_application()
        test_constraint_with_fixatoms()
        test_stage_constraint_config()
        test_covalent_radii_calculation()
        test_real_world_scenario()

        print("=" * 60)
        print("✅ ALL INTEGRATION TESTS PASSED!")
        print("=" * 60)
        print("\nConstraint support is ready for use in workflows.")
        print("See: asetools/manager/sample_yaml/aor_with_hookean_constraints.yaml")
        print("\n")

    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0

if __name__ == '__main__':
    exit(main())
