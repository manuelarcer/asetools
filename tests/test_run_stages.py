"""Tests for the run_stages engine: ordering guard, production warning, and
the run_workflow delegation wrapper."""

import logging
import os
import shutil
import tempfile

import pytest
import yaml
from ase import Atoms

from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML


def _vasp_stage(name):
    return {"name": name, "engine": "vasp", "_template": "t",
            "steps": [{"name": "sp", "overrides": {"nsw": 0}}]}


def _mlip_stage(name):
    return {"name": name, "engine": "mlip", "_template": "uma_preopt",
            "mlip": "uma-s-1p2", "env": "uma", "optimizer": "bfgs",
            "fmax": 0.1, "max_steps": 50}


class TestOrderingGuard:
    def test_mlip_after_vasp_raises(self):
        from asetools.workflow.manager import _validate_stage_ordering

        stages = [_vasp_stage("V1"), _mlip_stage("M1")]
        with pytest.raises(ValueError, match="MLIP"):
            _validate_stage_ordering(stages)

    def test_mlip_before_vasp_is_allowed(self):
        from asetools.workflow.manager import _validate_stage_ordering

        stages = [_mlip_stage("M1"), _vasp_stage("V1"), _vasp_stage("V2")]
        _validate_stage_ordering(stages)  # must not raise

    def test_all_vasp_is_allowed(self):
        from asetools.workflow.manager import _validate_stage_ordering

        _validate_stage_ordering([_vasp_stage("V1"), _vasp_stage("V2")])

    def test_multiple_mlip_then_vasp_allowed(self):
        from asetools.workflow.manager import _validate_stage_ordering

        _validate_stage_ordering(
            [_mlip_stage("M1"), _mlip_stage("M2"), _vasp_stage("V1")]
        )


class TestProductionWarning:
    def test_warns_when_production_not_last(self, caplog):
        from asetools.workflow.manager import _warn_if_production_not_last

        stages = [
            {"name": "OPT_340", "_template": "encut_ramp"},
            {"name": "OPT_500", "_template": "production"},
            {"name": "OPT_420", "_template": "encut_ramp"},
        ]
        with caplog.at_level(logging.WARNING):
            _warn_if_production_not_last(stages, "production")
        assert any("production" in r.message.lower() for r in caplog.records)

    def test_warns_when_production_absent(self, caplog):
        from asetools.workflow.manager import _warn_if_production_not_last

        stages = [{"name": "OPT_340", "_template": "encut_ramp"}]
        with caplog.at_level(logging.WARNING):
            _warn_if_production_not_last(stages, "production")
        assert any("production" in r.message.lower() for r in caplog.records)

    def test_no_warning_when_production_last(self, caplog):
        from asetools.workflow.manager import _warn_if_production_not_last

        stages = [
            {"name": "OPT_340", "_template": "encut_ramp"},
            {"name": "OPT_500", "_template": "production"},
        ]
        with caplog.at_level(logging.WARNING):
            _warn_if_production_not_last(stages, "production")
        assert not any("production" in r.message.lower() for r in caplog.records)

    def test_disabled_when_production_none(self, caplog):
        from asetools.workflow.manager import _warn_if_production_not_last

        stages = [{"name": "OPT_340", "_template": "encut_ramp"}]
        with caplog.at_level(logging.WARNING):
            _warn_if_production_not_last(stages, None)
        assert not any("production" in r.message.lower() for r in caplog.records)


class TestRunStagesDryRun:
    """Engine orchestration verified with dry_run (no VASP execution)."""

    def setup_method(self):
        self.dir = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.dir)

    def teardown_method(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.dir)

    def _cfg(self):
        data = {
            "basic": {"encut": 500, "ediff": 1e-5},
            "systems": {"default": None},
            "workflows": {"relax": {"stages": [
                {"name": "WF_A", "steps": [{"name": "sp", "overrides": {"nsw": 0}}]},
            ]}},
            "globals": {"initial_conf_pattern": "POSCAR"},
        }
        path = os.path.join(self.dir, "config.yaml")
        with open(path, "w") as f:
            yaml.dump(data, f)
        return VASPConfigurationFromYAML(config_file=path, system="default")

    def test_dry_run_marks_each_stage_done_in_order(self):
        from asetools.workflow.manager import run_stages

        cfg = self._cfg()
        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]], cell=[10, 10, 10])
        stages = [_vasp_stage("S1"), _vasp_stage("S2")]
        run_stages(atoms, cfg, stages=stages, dry_run=True, production=None)
        assert os.path.exists("STAGE_S1_DONE")
        assert os.path.exists("STAGE_S2_DONE")

    def test_dry_run_skips_already_done_stage(self):
        from asetools.workflow.manager import run_stages

        cfg = self._cfg()
        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]], cell=[10, 10, 10])
        with open("STAGE_S1_DONE", "w") as f:
            f.write("done\n")
        # S1 already done; running should not error and should create S2
        run_stages(atoms, cfg, stages=[_vasp_stage("S1"), _vasp_stage("S2")],
                   dry_run=True, production=None)
        assert os.path.exists("STAGE_S2_DONE")

    def test_run_workflow_delegates_to_engine(self):
        from asetools.workflow.manager import run_workflow

        cfg = self._cfg()
        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]], cell=[10, 10, 10])
        run_workflow(atoms, cfg, workflow_name="relax", dry_run=True)
        assert os.path.exists("STAGE_WF_A_DONE")


class TestMlipStageBranch:
    def setup_method(self):
        self.dir = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.dir)

    def teardown_method(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.dir)

    def _cfg(self):
        data = {
            "basic": {"encut": 500},
            "systems": {"default": None},
            "workflows": {"relax": {"stages": []}},
            "globals": {
                "initial_conf_pattern": "POSCAR",
                "mlip_envs": {"uma": "/fake/uma-env/bin/python"},
            },
        }
        path = os.path.join(self.dir, "config.yaml")
        with open(path, "w") as f:
            yaml.dump(data, f)
        return VASPConfigurationFromYAML(config_file=path, system="default")

    def _mlip_stage_with_constraints(self):
        s = _mlip_stage("M1")
        s["device"] = "auto"
        s["relax_cell"] = False
        s["constraints"] = {
            "type": "hookean",
            "config_file": "pairs.json",
            "spring_constant": 20.0,
        }
        return s

    def test_build_mlip_command_contains_expected_flags(self):
        from asetools.workflow.manager import _build_mlip_command

        cmd = _build_mlip_command("/x/py", "POSCAR", self._mlip_stage_with_constraints())
        assert cmd[:4] == ["/x/py", "-m", "asetools.workflow.mlip_runner", "--structure"]
        # flag/value pairs present
        for flag, val in [
            ("--structure", "POSCAR"), ("--name", "M1"), ("--mlip", "uma-s-1p2"),
            ("--optimizer", "bfgs"), ("--fmax", "0.1"), ("--max-steps", "50"),
            ("--device", "auto"), ("--constraints-json", "pairs.json"),
            ("--spring-constant", "20.0"),
        ]:
            assert cmd[cmd.index(flag) + 1] == val
        assert "--relax-cell" not in cmd  # False → flag omitted

    def test_relax_cell_true_adds_flag(self):
        from asetools.workflow.manager import _build_mlip_command

        s = self._mlip_stage_with_constraints()
        s["relax_cell"] = True
        cmd = _build_mlip_command("/x/py", "POSCAR", s)
        assert "--relax-cell" in cmd

    def test_missing_env_key_raises(self):
        from asetools.workflow.manager import _run_mlip_stage

        cfg = self._cfg()
        s = _mlip_stage("M1")
        s["env"] = "not_registered"
        with open("POSCAR", "w") as f:
            f.write("dummy\n")
        with pytest.raises(KeyError, match="not_registered"):
            _run_mlip_stage(cfg, s, dry_run=False)

    def test_dry_run_marks_done_without_launching(self, monkeypatch):
        import asetools.workflow.manager as mgr

        def _boom(*a, **k):
            raise AssertionError("subprocess must not run in dry_run")

        monkeypatch.setattr(mgr.subprocess, "run", _boom)
        cfg = self._cfg()
        with open("POSCAR", "w") as f:
            f.write("dummy\n")
        mgr._run_mlip_stage(cfg, _mlip_stage("M1"), dry_run=True)
        assert os.path.exists("STAGE_M1_DONE")

    def test_nonzero_exit_raises_and_no_sentinel(self, monkeypatch):
        import asetools.workflow.manager as mgr

        class _Result:
            returncode = 1

        monkeypatch.setattr(mgr.subprocess, "run", lambda *a, **k: _Result())
        cfg = self._cfg()
        with open("POSCAR", "w") as f:
            f.write("dummy\n")
        with pytest.raises(RuntimeError, match="converge"):
            mgr._run_mlip_stage(cfg, _mlip_stage("M1"), dry_run=False)
        assert not os.path.exists("STAGE_M1_DONE")

    def test_zero_exit_marks_done(self, monkeypatch):
        import asetools.workflow.manager as mgr

        class _Result:
            returncode = 0

        monkeypatch.setattr(mgr.subprocess, "run", lambda *a, **k: _Result())
        cfg = self._cfg()
        with open("POSCAR", "w") as f:
            f.write("dummy\n")
        mgr._run_mlip_stage(cfg, _mlip_stage("M1"), dry_run=False)
        assert os.path.exists("STAGE_M1_DONE")

    def test_run_stages_dry_run_mixed_pipeline(self):
        from asetools.workflow.manager import run_stages

        cfg = self._cfg()
        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]], cell=[10, 10, 10])
        with open("POSCAR", "w") as f:
            f.write("dummy\n")
        stages = [_mlip_stage("M1"), _vasp_stage("V1")]
        run_stages(atoms, cfg, stages=stages, dry_run=True, production=None)
        assert os.path.exists("STAGE_M1_DONE")
        assert os.path.exists("STAGE_V1_DONE")
