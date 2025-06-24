# Description of the script:

import logging
import os
import shutil
import sys
import glob
from ase.io import read
from .calculatorsetuptools import VASPConfigurationFromYAML, deep_update, setup_initial_magmom
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

def run_workflow(atoms: Atoms, cfg: VASPConfigurationFromYAML, workflow_name: str, run_overrides: dict = None, dry_run: bool = False):
    
    # get the initial magmom from the config to be used at each stage
    initial_magmom = cfg.initial_magmom     # {} if not defined

    to_run = stages_to_run(cfg, workflow_name)
    stages = cfg.workflows[workflow_name]['stages']
    
    if not to_run:
        logger.info(f"-->  Workflow '{workflow_name}' is already completed, nothing to do  <--")
        return

    for stage in stages:
        if stage['name'] not in to_run:
            logger.info(f"Skipping STAGE: {stage['name']}, already done")
            continue
        _run_stage(atoms, cfg, stage, run_overrides, dry_run, initial_magmom=initial_magmom)
    logger.info(f"-->  Workflow '{workflow_name}' completed successfully  <--")
    

def _run_stage(atoms: Atoms, cfg: VASPConfigurationFromYAML, stage: dict, run_overrides: dict, dry_run: bool, initial_magmom: dict):
    # run_overrides is different from the overrides in each step
    name  = stage['name']
    steps = stage['steps']
    logger.info(f"Running STAGE: {stage['name']}")
    for step in steps:
        # Make calculator for each step
        logger.info('  * Setting up calculator from config')
        calc = make_calculator(cfg, run_overrides=run_overrides)
        atoms.calc = calc
        atoms = setup_initial_magmom(atoms, magmom_dict=initial_magmom)
        _run_step(atoms, step, dry_run)
    
    backup_output_files(name=name)
    logger.info(f" -- ✅ Stage '{name}' completed and backed up")
    _mark_done(name)

def _run_step(atoms: Atoms, step: dict, dry_run: bool):
    name      = step['name']
    overrides = step.get('overrides', {})
    logger.info(f"  • Step '{name}' overrides={overrides}")
    atoms.calc.set(**overrides)
    if dry_run:
        logger.info("    (dry-run, skipping calculation)")
        return
    
    # trigger VASP run
    atoms.get_potential_energy()
    logger.info(f"    ✔ Step '{name}' done")
    

def _mark_done(step_name: str):
    path = f"STAGE_{step_name}_DONE"
    with open(path, 'w') as f:
        f.write(f"{step_name} completed\n")


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