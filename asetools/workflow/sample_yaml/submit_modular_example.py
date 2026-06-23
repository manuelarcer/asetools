#!/usr/bin/env python
"""Example submission script for the modular stage-composition path.

Composes the optimization procedure here in Python from reusable stage
templates defined in modular_stage_composition.yaml, instead of selecting a
single fully specified named workflow. The canonical final ("production")
parameters stay authoritative in the YAML.

Pairs with: modular_stage_composition.yaml
Docs:       docs/vasp_calculation_guide.md ("Modular stage composition")
"""

import logging

from asetools.workflow.calculatorsetuptools import (
    VASPConfigurationFromYAML,
    setup_initial_magmom,
)
from asetools.workflow.logger import configure_logging
from asetools.workflow.manager import load_structure, run_stages

yaml_config_path = "modular_stage_composition.yaml"
system = "CuNiOOH"

configure_logging(file_prefix="run")
logger = logging.getLogger(__name__)

cfg = VASPConfigurationFromYAML(config_file=yaml_config_path, system=system)
atoms = load_structure(cfg.globals["initial_conf_pattern"])
setup_initial_magmom(atoms, cfg.initial_magmom_data)

# --- Compose the procedure ---------------------------------------------------
# Optional cheap MLIP warm-start (pre-optimization only; must precede VASP).
PREOPTIMIZATION = [
    cfg.stage("uma_preopt", name="OPT_MLIP"),
    cfg.stage(
        "encut_ramp",
        name="OPT_340",
        overrides={"encut": 340},
        optimizer_kwargs={"fmax": 0.05},
    ),
    cfg.stage(
        "encut_ramp",
        name="OPT_420",
        overrides={"encut": 420},
        optimizer_kwargs={"fmax": 0.02},
    ),
]

# Canonical final optimization — no overrides, YAML values stay authoritative.
PRODUCTION = cfg.stage("production", name="OPT_500")

# run_stages warns if the production-template stage is not run last.
run_stages(atoms, cfg, stages=PREOPTIMIZATION + [PRODUCTION], dry_run=False)
