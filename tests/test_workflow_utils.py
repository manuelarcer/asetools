"""Tests for workflow utility modules (calculatorsetuptools, logger)."""

import logging
import os
import shutil
import tempfile

import numpy as np
import pytest
import yaml
from ase import Atoms
from ase.build import molecule


class TestLoadYamlConfig:
    """Test YAML configuration loading."""

    def setup_method(self):
        self.test_dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.test_dir)

    def _write_config(self, data, name="config.yaml"):
        path = os.path.join(self.test_dir, name)
        with open(path, "w") as f:
            yaml.dump(data, f)
        return path

    def _minimal_config(self, **overrides):
        cfg = {
            "basic": {"encut": 500, "ediff": 1e-5},
            "systems": {"default": None},
            "workflows": {"relax": {"stages": []}},
            "globals": {"vasp_pp_path": "/path"},
        }
        cfg.update(overrides)
        return cfg

    def test_load_valid_config(self):
        from asetools.workflow.calculatorsetuptools import load_yaml_config

        path = self._write_config(self._minimal_config())
        cfg = load_yaml_config(path)
        assert "basic" in cfg
        assert "systems" in cfg
        assert cfg["basic"]["encut"] == 500

    def test_load_config_missing_key(self):
        from asetools.workflow.calculatorsetuptools import load_yaml_config

        path = self._write_config({"basic": {}, "systems": {}})
        with pytest.raises(KeyError, match="workflows"):
            load_yaml_config(path)

    def test_load_config_file_not_found(self):
        from asetools.workflow.calculatorsetuptools import load_yaml_config

        with pytest.raises(FileNotFoundError):
            load_yaml_config("/nonexistent/config.yaml")


class TestVerifyConfigurationKeys:
    def test_all_keys_present(self):
        from asetools.workflow.calculatorsetuptools import verify_configuration_keys

        cfg = {"basic": {}, "systems": {}, "workflows": {}, "globals": {}}
        verify_configuration_keys(cfg)  # Should not raise

    def test_missing_basic(self):
        from asetools.workflow.calculatorsetuptools import verify_configuration_keys

        with pytest.raises(KeyError, match="basic"):
            verify_configuration_keys({"systems": {}, "workflows": {}, "globals": {}})

    def test_missing_globals(self):
        from asetools.workflow.calculatorsetuptools import verify_configuration_keys

        with pytest.raises(KeyError, match="globals"):
            verify_configuration_keys({"basic": {}, "systems": {}, "workflows": {}})


class TestDeepUpdate:
    def test_flat_override(self):
        from asetools.workflow.calculatorsetuptools import deep_update

        base = {"a": 1, "b": 2}
        result = deep_update(base, {"b": 3, "c": 4})
        assert result == {"a": 1, "b": 3, "c": 4}

    def test_nested_override(self):
        from asetools.workflow.calculatorsetuptools import deep_update

        base = {"a": {"x": 1, "y": 2}, "b": 3}
        result = deep_update(base, {"a": {"y": 99}})
        assert result["a"]["x"] == 1
        assert result["a"]["y"] == 99
        assert result["b"] == 3

    def test_nested_replace_with_non_dict(self):
        from asetools.workflow.calculatorsetuptools import deep_update

        base = {"a": {"x": 1}}
        result = deep_update(base, {"a": 42})
        assert result["a"] == 42

    def test_empty_override(self):
        from asetools.workflow.calculatorsetuptools import deep_update

        base = {"a": 1}
        result = deep_update(base, {})
        assert result == {"a": 1}

    def test_empty_base(self):
        from asetools.workflow.calculatorsetuptools import deep_update

        result = deep_update({}, {"a": 1})
        assert result == {"a": 1}


class TestSetupInitialMagmom:
    def test_dict_based_magmom(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = molecule("H2O")
        setup_initial_magmom(atoms, {"O": 2.0, "H": 0.5})
        magmoms = atoms.get_initial_magnetic_moments()
        # O is index 0 in H2O molecule
        for i, atom in enumerate(atoms):
            if atom.symbol == "O":
                assert magmoms[i] == pytest.approx(2.0)
            else:
                assert magmoms[i] == pytest.approx(0.5)

    def test_dict_missing_element(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = molecule("H2O")
        setup_initial_magmom(atoms, {"O": 2.0})
        magmoms = atoms.get_initial_magnetic_moments()
        for i, atom in enumerate(atoms):
            if atom.symbol == "H":
                assert magmoms[i] == pytest.approx(0.0)

    def test_list_based_magmom(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = Atoms("CuO", positions=[[0, 0, 0], [1, 0, 0]])
        setup_initial_magmom(atoms, [3.0, 1.0])
        magmoms = atoms.get_initial_magnetic_moments()
        assert magmoms[0] == pytest.approx(3.0)
        assert magmoms[1] == pytest.approx(1.0)

    def test_list_shorter_pads_with_zero(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = Atoms("CuOH", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        setup_initial_magmom(atoms, [3.0])
        magmoms = atoms.get_initial_magnetic_moments()
        assert magmoms[0] == pytest.approx(3.0)
        assert magmoms[1] == pytest.approx(0.0)
        assert magmoms[2] == pytest.approx(0.0)

    def test_list_too_long_raises(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = Atoms("Cu", positions=[[0, 0, 0]])
        with pytest.raises(ValueError, match="longer than"):
            setup_initial_magmom(atoms, [1.0, 2.0, 3.0])

    def test_none_sets_zero(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = Atoms("CuO", positions=[[0, 0, 0], [1, 0, 0]])
        setup_initial_magmom(atoms, None)
        magmoms = atoms.get_initial_magnetic_moments()
        assert all(m == pytest.approx(0.0) for m in magmoms)

    def test_numpy_array_magmom(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = Atoms("CuO", positions=[[0, 0, 0], [1, 0, 0]])
        setup_initial_magmom(atoms, np.array([2.0, 1.0]))
        magmoms = atoms.get_initial_magnetic_moments()
        assert magmoms[0] == pytest.approx(2.0)

    def test_invalid_type_raises(self):
        from asetools.workflow.calculatorsetuptools import setup_initial_magmom

        atoms = Atoms("Cu", positions=[[0, 0, 0]])
        with pytest.raises(TypeError, match="must be dict, list, or None"):
            setup_initial_magmom(atoms, "invalid")


class TestVASPConfigurationFromYAML:
    """Test VASPConfigurationFromYAML class."""

    def setup_method(self):
        self.test_dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.test_dir)

    def _write_config(self, data):
        path = os.path.join(self.test_dir, "config.yaml")
        with open(path, "w") as f:
            yaml.dump(data, f)
        return path

    def test_basic_init(self):
        from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML

        cfg = {
            "basic": {"encut": 500},
            "systems": {"default": None},
            "workflows": {"relax": {}},
            "globals": {"vasp_pp_path": "/path"},
        }
        path = self._write_config(cfg)
        vc = VASPConfigurationFromYAML(path)
        assert vc.basic_config["encut"] == 500

    def test_system_config_none(self):
        from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML

        cfg = {
            "basic": {},
            "systems": {"default": None},
            "workflows": {},
            "globals": {},
        }
        path = self._write_config(cfg)
        vc = VASPConfigurationFromYAML(path)
        assert vc.system_config == {}

    def test_system_config_missing_system(self):
        from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML

        cfg = {
            "basic": {},
            "systems": {"NCA": {"ispin": 2}},
            "workflows": {},
            "globals": {},
        }
        path = self._write_config(cfg)
        vc = VASPConfigurationFromYAML(path, system="nonexistent")
        assert vc.system_config == {}

    def test_initial_magmom_from_system(self):
        from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML

        cfg = {
            "basic": {},
            "systems": {"NCA": {"magmom": {"Ni": 2.0, "Co": 3.0}}},
            "workflows": {},
            "globals": {},
        }
        path = self._write_config(cfg)
        vc = VASPConfigurationFromYAML(path, system="NCA")
        assert vc.initial_magmom_data == {"Ni": 2.0, "Co": 3.0}

    def test_initial_magmom_empty(self):
        from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML

        cfg = {
            "basic": {},
            "systems": {"default": {"ispin": 2}},
            "workflows": {},
            "globals": {},
        }
        path = self._write_config(cfg)
        vc = VASPConfigurationFromYAML(path)
        assert vc.initial_magmom_data == {}


class TestConfigureLogging:
    """Test the workflow logger."""

    def setup_method(self):
        self.test_dir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.test_dir)

    def teardown_method(self):
        os.chdir(self.orig_dir)
        # Clean up handlers to avoid affecting other tests
        root = logging.getLogger()
        for handler in root.handlers[:]:
            root.removeHandler(handler)
            handler.close()
        shutil.rmtree(self.test_dir)

    def test_creates_logfile(self):
        from asetools.workflow.logger import configure_logging

        configure_logging(file_prefix="test")
        log_files = [f for f in os.listdir(".") if f.startswith("test_") and f.endswith(".log")]
        assert len(log_files) >= 1

    def test_custom_prefix(self):
        from asetools.workflow.logger import configure_logging

        configure_logging(file_prefix="myrun")
        log_files = [f for f in os.listdir(".") if f.startswith("myrun_") and f.endswith(".log")]
        assert len(log_files) >= 1


class TestDatabaseFunctions:
    """Test database module functions."""

    def test_check_if_exists_empty_db(self):
        from ase.db import connect

        from asetools.database.databases import check_if_exists_in_db

        db_path = os.path.join(tempfile.mkdtemp(), "test.db")
        db = connect(db_path)
        atoms = molecule("H2O")
        # Need forces set
        from ase.calculators.singlepoint import SinglePointCalculator

        calc = SinglePointCalculator(atoms, energy=-10.0, forces=np.zeros((3, 3)))
        atoms.calc = calc

        exists, idx = check_if_exists_in_db(db, atoms)
        assert exists is False
        assert idx is None

    def test_db_to_pandas_empty(self):
        from ase.db import connect

        from asetools.database.databases import db_to_pandas

        db_path = os.path.join(tempfile.mkdtemp(), "test.db")
        db = connect(db_path)
        df = db_to_pandas(db, columns=["id"])
        assert len(df) == 0
        assert "id" in df.columns

    def test_db_to_pandas_with_data(self):
        from ase.calculators.singlepoint import SinglePointCalculator
        from ase.db import connect

        from asetools.database.databases import db_to_pandas

        db_path = os.path.join(tempfile.mkdtemp(), "test.db")
        db = connect(db_path)

        atoms = molecule("H2O")
        calc = SinglePointCalculator(atoms, energy=-10.5, forces=np.zeros((3, 3)))
        atoms.calc = calc
        db.write(atoms, name="water")

        df = db_to_pandas(db, columns=["name", "id", "energy"])
        assert len(df) == 1
        assert df.iloc[0]["name"] == "water"
        assert df.iloc[0]["energy"] == pytest.approx(-10.5)
