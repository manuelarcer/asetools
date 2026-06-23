"""Tests for the modular stage-composition factory (cfg.stage())."""

import os
import shutil
import tempfile

import pytest
import yaml

from asetools.workflow.calculatorsetuptools import VASPConfigurationFromYAML


def _write(data, directory, name="config.yaml"):
    path = os.path.join(directory, name)
    with open(path, "w") as f:
        yaml.dump(data, f)
    return path


def _config_dict(**overrides):
    cfg = {
        "basic": {"encut": 500, "ediff": 1e-5},
        "systems": {"default": None},
        "workflows": {"relax": {"stages": []}},
        "globals": {"vasp_pp_path": "/path"},
        "stage_templates": {
            "encut_ramp": {
                "constraints": {
                    "type": "hookean",
                    "config_file": "hookean_c_pairs.json",
                    "spring_constant": 20.0,
                },
                "steps": [
                    {
                        "name": "single_point",
                        "overrides": {"nsw": 0, "ispin": 1, "encut": 340},
                    },
                    {
                        "name": "optimization_bfgs",
                        "overrides": {"nsw": 500, "encut": 340},
                        "optimizer": "BFGS",
                        "optimizer_kwargs": {"maxstep": 0.2},
                    },
                ],
            },
            "production": {
                "steps": [
                    {"name": "optimization", "overrides": {"encut": 500}},
                ],
            },
            "uma_preopt": {
                "engine": "mlip",
                "mlip": "uma-s-1p2",
                "env": "uma",
                "optimizer": "bfgs",
                "fmax": 0.10,
                "max_steps": 300,
                "relax_cell": False,
                "device": "auto",
                "constraints": {
                    "type": "hookean",
                    "config_file": "hookean_c_pairs.json",
                    "spring_constant": 20.0,
                },
            },
        },
    }
    cfg.update(overrides)
    return cfg


class TestStageFactory:
    def setup_method(self):
        self.dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.dir)

    def _cfg(self, **overrides):
        path = _write(_config_dict(**overrides), self.dir)
        return VASPConfigurationFromYAML(config_file=path, system="default")

    def test_resolves_vasp_template_with_name_and_provenance(self):
        cfg = self._cfg()
        stage = cfg.stage("encut_ramp", name="OPT_340")
        assert stage["name"] == "OPT_340"
        assert stage["_template"] == "encut_ramp"
        assert stage.get("engine", "vasp") == "vasp"
        assert [s["name"] for s in stage["steps"]] == ["single_point", "optimization_bfgs"]

    def test_override_merges_into_every_step(self):
        cfg = self._cfg()
        stage = cfg.stage("encut_ramp", name="OPT_420", overrides={"encut": 420})
        for step in stage["steps"]:
            assert step["overrides"]["encut"] == 420
        # untouched template values survive
        assert stage["steps"][0]["overrides"]["ispin"] == 1
        assert stage["steps"][1]["overrides"]["nsw"] == 500

    def test_optimizer_kwargs_merge_into_steps(self):
        cfg = self._cfg()
        stage = cfg.stage(
            "encut_ramp", name="OPT_340", optimizer_kwargs={"fmax": 0.05}
        )
        opt_step = stage["steps"][1]
        assert opt_step["optimizer_kwargs"]["fmax"] == 0.05
        # template optimizer_kwargs preserved
        assert opt_step["optimizer_kwargs"]["maxstep"] == 0.2

    def test_constraints_deep_merge_patches_only_changed_field(self):
        cfg = self._cfg()
        stage = cfg.stage(
            "encut_ramp", name="OPT_340", constraints={"spring_constant": 30.0}
        )
        c = stage["constraints"]
        assert c["spring_constant"] == 30.0
        assert c["type"] == "hookean"
        assert c["config_file"] == "hookean_c_pairs.json"

    def test_unknown_template_raises_keyerror(self):
        cfg = self._cfg()
        with pytest.raises(KeyError, match="nonexistent"):
            cfg.stage("nonexistent", name="X")

    def test_name_is_required(self):
        cfg = self._cfg()
        with pytest.raises(TypeError):
            cfg.stage("encut_ramp")  # missing required name=

    def test_does_not_mutate_template(self):
        cfg = self._cfg()
        cfg.stage("encut_ramp", name="OPT_340", overrides={"encut": 999})
        again = cfg.stage("encut_ramp", name="OPT_OTHER")
        assert again["steps"][0]["overrides"]["encut"] == 340

    def test_mlip_template_sets_engine_and_flat_fields(self):
        cfg = self._cfg()
        stage = cfg.stage("uma_preopt", name="MLIP_PRE")
        assert stage["engine"] == "mlip"
        assert stage["mlip"] == "uma-s-1p2"
        assert stage["env"] == "uma"
        assert stage["fmax"] == 0.10
        assert stage["max_steps"] == 300

    def test_mlip_optimizer_kwargs_override_flat_fields(self):
        cfg = self._cfg()
        stage = cfg.stage(
            "uma_preopt", name="MLIP_PRE", optimizer_kwargs={"fmax": 0.2, "max_steps": 100}
        )
        assert stage["fmax"] == 0.2
        assert stage["max_steps"] == 100

    def test_mlip_rejects_vasp_overrides(self):
        cfg = self._cfg()
        with pytest.raises(ValueError, match="overrides"):
            cfg.stage("uma_preopt", name="MLIP_PRE", overrides={"encut": 340})

    def test_missing_stage_templates_section_raises_keyerror(self):
        path = _write(
            {
                "basic": {},
                "systems": {"default": None},
                "workflows": {"relax": {"stages": []}},
                "globals": {},
            },
            self.dir,
            name="bare.yaml",
        )
        cfg = VASPConfigurationFromYAML(config_file=path, system="default")
        with pytest.raises(KeyError):
            cfg.stage("encut_ramp", name="X")
