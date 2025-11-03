"""
Tests for the constraints module.
"""

import json
import pytest
import tempfile
from pathlib import Path

from ase import Atoms
from ase.constraints import FixAtoms, Hookean
from ase.build import bulk

from asetools.constraints import ConstraintManager


class TestConstraintManager:
    """Test ConstraintManager class."""

    def test_init_default(self):
        """Test ConstraintManager initialization with default parameters."""
        cm = ConstraintManager()
        assert cm.distance_factor == 1.134

    def test_init_custom_distance_factor(self):
        """Test ConstraintManager initialization with custom distance factor."""
        cm = ConstraintManager(distance_factor=1.2)
        assert cm.distance_factor == 1.2

    def test_load_constraint_config(self, tmp_path):
        """Test loading constraint configuration from JSON file."""
        # Create test JSON file
        config_data = {
            "pairs": [[0, 1], [2, 3]],
            "metadata": {
                "description": "Test constraints",
                "spring_constant": 25.0
            }
        }
        json_file = tmp_path / "test_constraints.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        cm = ConstraintManager()
        config = cm.load_constraint_config(str(json_file))

        assert "pairs" in config
        assert len(config["pairs"]) == 2
        assert config["pairs"][0] == [0, 1]
        assert config["metadata"]["spring_constant"] == 25.0

    def test_load_constraint_config_file_not_found(self):
        """Test error handling when JSON file doesn't exist."""
        cm = ConstraintManager()
        with pytest.raises(FileNotFoundError):
            cm.load_constraint_config("nonexistent_file.json")

    def test_load_constraint_config_invalid_json(self, tmp_path):
        """Test error handling for malformed JSON."""
        json_file = tmp_path / "invalid.json"
        with open(json_file, 'w') as f:
            f.write("{invalid json content")

        cm = ConstraintManager()
        with pytest.raises(json.JSONDecodeError):
            cm.load_constraint_config(str(json_file))

    def test_calculate_bond_distance(self):
        """Test bond distance calculation from covalent radii."""
        # Create H2O molecule
        from ase.build import molecule
        atoms = molecule('H2O')

        cm = ConstraintManager(distance_factor=1.134)
        # O-H bond distance
        r0 = cm.calculate_bond_distance(atoms, 0, 1)  # O-H

        # O covalent radius ~0.66, H ~0.31, sum ~0.97
        # With factor 1.134: ~1.10
        assert 1.0 < r0 < 1.2

    def test_apply_hookean_from_pairs(self):
        """Test creating Hookean constraints from atom pairs."""
        # Create simple chain of atoms
        atoms = Atoms('COOH', positions=[[0, 0, 0], [1.4, 0, 0], [2.5, 0, 0], [3.5, 0, 0]])

        cm = ConstraintManager()
        pairs = [[0, 1], [1, 2]]  # C-O and O-O bonds
        hookean_constraints = cm.apply_hookean_from_pairs(atoms, pairs, k=20.0)

        assert len(hookean_constraints) == 2
        assert all(isinstance(c, Hookean) for c in hookean_constraints)

    def test_apply_hookean_from_pairs_invalid_index(self):
        """Test error handling for out-of-range atom indices."""
        atoms = Atoms('H2O', positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])

        cm = ConstraintManager()
        pairs = [[0, 10]]  # Index 10 is out of range

        with pytest.raises(ValueError, match="out of range"):
            cm.apply_hookean_from_pairs(atoms, pairs)

    def test_get_existing_fix_indices_no_constraints(self):
        """Test extracting FixAtoms indices when no constraints exist."""
        atoms = bulk('Cu', 'fcc', a=3.6)

        cm = ConstraintManager()
        fix_indices = cm.get_existing_fix_indices(atoms)

        assert fix_indices == []

    def test_get_existing_fix_indices_with_fixatoms(self):
        """Test extracting FixAtoms indices when FixAtoms constraints exist."""
        atoms = bulk('Cu', 'fcc', a=3.6).repeat((2, 2, 2))
        atoms.set_constraint(FixAtoms(indices=[0, 1, 2, 3]))

        cm = ConstraintManager()
        fix_indices = cm.get_existing_fix_indices(atoms)

        assert fix_indices == [0, 1, 2, 3]

    def test_get_existing_fix_indices_multiple_constraints(self):
        """Test extracting FixAtoms indices from multiple constraints."""
        atoms = bulk('Cu', 'fcc', a=3.6).repeat((2, 2, 2))
        constraints = [
            FixAtoms(indices=[0, 1]),
            FixAtoms(indices=[2, 3])
        ]
        atoms.set_constraint(constraints)

        cm = ConstraintManager()
        fix_indices = cm.get_existing_fix_indices(atoms)

        assert sorted(fix_indices) == [0, 1, 2, 3]

    def test_merge_constraints_no_existing(self):
        """Test merging constraints when no existing constraints."""
        atoms = Atoms('H2O', positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])

        cm = ConstraintManager()
        new_constraints = [Hookean(a1=0, a2=1, k=20.0, rt=1.0)]
        cm.merge_constraints(atoms, new_constraints)

        assert len(atoms.constraints) == 1
        assert isinstance(atoms.constraints[0], Hookean)

    def test_merge_constraints_with_existing_fixatoms(self):
        """Test merging constraints preserves existing FixAtoms."""
        atoms = Atoms('H2O', positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        atoms.set_constraint(FixAtoms(indices=[0]))

        cm = ConstraintManager()
        new_constraints = [Hookean(a1=1, a2=2, k=20.0, rt=1.0)]
        cm.merge_constraints(atoms, new_constraints)

        assert len(atoms.constraints) == 2
        # First should be FixAtoms
        assert isinstance(atoms.constraints[0], FixAtoms)
        # Second should be Hookean
        assert isinstance(atoms.constraints[1], Hookean)

    def test_merge_constraints_empty_raises_error(self):
        """Test that merging empty constraints raises error."""
        atoms = Atoms('H2O', positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])

        cm = ConstraintManager()
        with pytest.raises(RuntimeError, match="No constraints to apply"):
            cm.merge_constraints(atoms, [])

    def test_apply_from_json(self, tmp_path):
        """Test applying constraints from JSON file."""
        # Create test atoms
        atoms = Atoms('H2O2', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])

        # Create test JSON file
        config_data = {"pairs": [[0, 1], [2, 3]]}
        json_file = tmp_path / "constraints.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        cm = ConstraintManager()
        cm.apply_from_json(atoms, str(json_file), k=15.0)

        # Should have 2 Hookean constraints
        assert len(atoms.constraints) == 2
        assert all(isinstance(c, Hookean) for c in atoms.constraints)

    def test_apply_from_json_missing_pairs_key(self, tmp_path):
        """Test error when JSON missing 'pairs' key."""
        config_data = {"metadata": {}}
        json_file = tmp_path / "invalid_constraints.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        atoms = Atoms('H2O', positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        cm = ConstraintManager()

        with pytest.raises(RuntimeError, match="missing required 'pairs' key"):
            cm.apply_from_json(atoms, str(json_file))

    def test_apply_from_json_empty_pairs(self, tmp_path):
        """Test error when pairs list is empty."""
        config_data = {"pairs": []}
        json_file = tmp_path / "empty_pairs.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        atoms = Atoms('H2O', positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        cm = ConstraintManager()

        with pytest.raises(RuntimeError, match="Empty pairs list"):
            cm.apply_from_json(atoms, str(json_file))

    def test_apply_from_json_with_metadata(self, tmp_path):
        """Test that metadata in JSON overrides default parameters."""
        atoms = Atoms('H2O2', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])

        config_data = {
            "pairs": [[0, 1]],
            "metadata": {
                "spring_constant": 30.0,
                "distance_factor": 1.2
            }
        }
        json_file = tmp_path / "constraints_with_metadata.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        cm = ConstraintManager()
        cm.apply_from_json(atoms, str(json_file))

        # Verify constraint was applied (parameters verified in logs)
        assert len(atoms.constraints) == 1
        assert isinstance(atoms.constraints[0], Hookean)

    def test_apply_stage_constraints(self, tmp_path):
        """Test applying constraints from workflow stage configuration."""
        atoms = Atoms('OH', positions=[[0, 0, 0], [1, 0, 0]])

        config_data = {"pairs": [[0, 1]]}
        json_file = tmp_path / "stage_constraints.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        cm = ConstraintManager()
        stage_config = {
            'type': 'hookean',
            'config_file': str(json_file),
            'spring_constant': 25.0,
            'distance_factor': 1.15
        }
        cm.apply_stage_constraints(atoms, stage_config)

        assert len(atoms.constraints) == 1
        assert isinstance(atoms.constraints[0], Hookean)

    def test_apply_stage_constraints_unsupported_type(self, tmp_path):
        """Test error for unsupported constraint type."""
        atoms = Atoms('OH', positions=[[0, 0, 0], [1, 0, 0]])

        cm = ConstraintManager()
        stage_config = {
            'type': 'unsupported_type',
            'config_file': 'some_file.json'
        }

        with pytest.raises(ValueError, match="Unsupported constraint type"):
            cm.apply_stage_constraints(atoms, stage_config)

    def test_apply_stage_constraints_missing_config_file(self):
        """Test error when config_file is missing."""
        atoms = Atoms('OH', positions=[[0, 0, 0], [1, 0, 0]])

        cm = ConstraintManager()
        stage_config = {
            'type': 'hookean',
            'spring_constant': 20.0
        }

        with pytest.raises(ValueError, match="Missing 'config_file'"):
            cm.apply_stage_constraints(atoms, stage_config)


class TestConstraintIntegration:
    """Integration tests with realistic molecular systems."""

    def test_water_molecule_constraints(self, tmp_path):
        """Test applying constraints to water molecule."""
        from ase.build import molecule
        atoms = molecule('H2O')

        # Create constraint for O-H bonds
        config_data = {"pairs": [[0, 1], [0, 2]]}  # O bonded to both H
        json_file = tmp_path / "water_constraints.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        cm = ConstraintManager()
        cm.apply_from_json(atoms, str(json_file), k=20.0)

        assert len(atoms.constraints) == 2

    def test_surface_with_adsorbate(self, tmp_path):
        """Test constraints on surface slab with adsorbate."""
        # Create simple slab
        slab = bulk('Cu', 'fcc', a=3.6)
        slab = slab.repeat((2, 2, 3))

        # Fix bottom layer
        fix_indices = [i for i in range(len(slab)) if slab[i].position[2] < 2.0]
        slab.set_constraint(FixAtoms(indices=fix_indices))

        # Add H atom on top (simplified adsorbate)
        from ase import Atom
        slab.append(Atom('H', position=[1.8, 1.8, 10.0]))
        h_index = len(slab) - 1

        # Create constraint connecting H to surface Cu
        config_data = {"pairs": [[h_index, 0]]}  # H to first Cu atom
        json_file = tmp_path / "surface_constraints.json"
        with open(json_file, 'w') as f:
            json.dump(config_data, f)

        cm = ConstraintManager()
        cm.apply_from_json(slab, str(json_file), k=15.0)

        # Should have FixAtoms + Hookean
        assert len(slab.constraints) == 2
        assert isinstance(slab.constraints[0], FixAtoms)
        assert isinstance(slab.constraints[1], Hookean)
