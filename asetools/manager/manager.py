# Description of the script:

import logging
import os
import shutil
import sys
import glob
from ase.io import read
from .calculatorsetuptools import VASPConfigurationFromYAML, deep_update
from ase.calculators.vasp import Vasp
from ase import Atoms

logger = logging.getLogger(__name__)    

def make_calculator(cfg: VASPConfigurationFromYAML, run_overrides: dict = None) -> Vasp:
    if run_overrides is None:
        run_overrides = {}
    else:
        logger.info(f" ** Overriding run parameters with: {run_overrides}")
    vasp_kwargs = deep_update(
        deep_update(cfg.basic_config.copy(), cfg.system_config),
        run_overrides
        )
    calc = Vasp(**vasp_kwargs)

    if 'magmom' in cfg.system_config:
        logger.warning(f' ** WARNING **  There is a "magmom" key in the system config, '
                       f' make sure you run setup_initial_magmom() from asetools/manager/calculatorsetuptools.py ')

    logger.info(" ** Calculator created")
    return calc

def run_workflow(atoms: Atoms, calc: Vasp, cfg: VASPConfigurationFromYAML, workflow_name: str):
    stages_todo = stages_to_run(cfg, workflow_name)

    if not stages_todo:
        logger.info(f"-->  Workflow '{workflow_name}' is already completed, nothing to do  <--")
        return

    for stage in cfg.workflows[workflow_name]['stages']:
        if stage['name'] not in stages_todo:
            logger.info(f"Skipping STAGE: {stage['name']}, already done")
            continue
        logger.info(f"Running STAGE: {stage['name']}")
        for step in stage['steps']:
            step_overrides = step.get('overrides', {})
            logger.info(f"  – Running STEP: * {step['name']} * with overrides: {step_overrides}")
            atoms.calc = calc
            atoms.calc.set(**step_overrides)
            atoms.get_potential_energy()  # This will trigger the VASP calculation
            logger.info(f"    STEP {step['name']} completed.")

            done_file = f"STAGE_{stage['name']}_DONE"
            with open(done_file, 'w') as f:
                f.write(f"Stage {stage['name']} completed successfully.\n")
    logger.info(f"-->  Workflow '{workflow_name}' completed successfully  <--")

def load_structure(pattern_initial_default: str = 'POSCAR'):
    """Load the structure from CONTCAR or initial configuration."""
    if os.path.exists('CONTCAR'):
        logger.info("Loading structure from CONTCAR")
        atoms = read('CONTCAR', format='vasp')
    else:
        matches = glob.glob(pattern_initial_default)
        if not matches:
            logger.warning(f"NO CONTCAR and NO match for pattern: {pattern_initial_default}")
            sys.exit(1)
        atoms = read(matches[0])
        logger.info(f"Loading structure from: {matches[0]}")
    return atoms

def stages_to_run(cfg: VASPConfigurationFromYAML, workflow_name: str = 'default') -> list:
    to_run = []
    for stage in cfg.workflows[workflow_name]['stages']:
        stage_name = stage['name']
        # Each stage is considered DONE when a file named 'STAGE_{name}_DONE' exists
        done_file = f"STAGE_{stage_name}_DONE"
        if os.path.exists(done_file):
            logger.info(f"Stage {stage_name} is DONE, skipping.")
            continue
        else:
            logger.info(f"Stage {stage_name} is NOT DONE, adding to run list.")
            to_run.append(stage_name)
    return to_run
        
def backup_output_files(name='backup'):
    for fname in ('POSCAR','INCAR','CONTCAR','OUTCAR','OSZICAR'):
        dst = f"{fname}_{name}"
        try:
            shutil.copy(fname, dst)
        except FileNotFoundError:
            logger.warning(f"File {fname} not found, skipping backup.")
    logger.info(f"Backup of output files completed with prefix '{name}'.")