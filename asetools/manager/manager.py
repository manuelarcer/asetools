# Description of the script:

import logging
import os
import shutil
import sys
import glob
from ase.io import read
from .calculatorsetuptools import VASPConfigurationFromYAML, deep_update, setup_initial_magmom
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive 
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS, GPMin, MDMin, QuasiNewton

logger = logging.getLogger(__name__)    

def make_calculator(cfg: VASPConfigurationFromYAML, run_overrides: dict = None):
    if run_overrides is None:
        run_overrides = {}
    else:
        logger.info(f" ** Overriding run parameters with: {run_overrides}")
    
    # Determine calculator type from configuration
    calculator_type = cfg.globals.get('calculator_type', 'vasp_interactive')
    
    vasp_kwargs = deep_update(
        deep_update(cfg.basic_config.copy(), cfg.system_config),
        run_overrides
    )
    
    if calculator_type.lower() == 'vasp_interactive':
        calc = VaspInteractive(**vasp_kwargs)
        logger.info(" ** VaspInteractive calculator created")
    elif calculator_type.lower() == 'vasp':
        calc = Vasp(**vasp_kwargs)
        logger.info(" ** Regular Vasp calculator created")
    else:
        raise ValueError(f"Unknown calculator type: {calculator_type}. Use 'vasp_interactive' or 'vasp'")

    return calc

def run_workflow(atoms: Atoms, cfg: VASPConfigurationFromYAML, workflow_name: str, run_overrides: dict = None, dry_run: bool = False):
    
    # get the initial magmom from the config to be used at each stage
    initial_magmom = cfg.initial_magmom_data     # {} if not defined

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
    optimizer_name = step.get('optimizer', None)
    optimizer_kwargs = step.get('optimizer_kwargs', {})
    
    logger.info(f"  • Step '{name}' overrides={overrides}")
    atoms.calc.set(**overrides)
    if dry_run:
        logger.info("    (dry-run, skipping calculation)")
        return
    
    # Check if we should use ASE optimizer (VaspInteractive only)
    if optimizer_name and isinstance(atoms.calc, VaspInteractive):
        logger.info(f"    Using ASE optimizer: {optimizer_name}")
        _run_with_ase_optimizer(atoms, optimizer_name, optimizer_kwargs)
    else:
        # trigger VASP run (regular optimization or single point)
        atoms.get_potential_energy()
    
    logger.info(f"    ✔ Step '{name}' done")


def _run_with_ase_optimizer(atoms: Atoms, optimizer_name: str, optimizer_kwargs: dict):
    """Run optimization using ASE optimizers with VaspInteractive calculator."""
    
    # Map optimizer names to classes
    optimizer_map = {
        'bfgs': BFGS,
        'fire': FIRE,
        'lbfgs': LBFGS,
        'gpmin': GPMin,
        'mdmin': MDMin,
        'quasinewton': QuasiNewton,
    }
    
    optimizer_name_lower = optimizer_name.lower()
    if optimizer_name_lower not in optimizer_map:
        raise ValueError(f"Unknown optimizer: {optimizer_name}. Available: {list(optimizer_map.keys())}")
    
    OptClass = optimizer_map[optimizer_name_lower]
    
    # Default parameters
    default_kwargs = {'logfile': f'{optimizer_name.upper()}.log'}
    default_kwargs.update(optimizer_kwargs)
    
    logger.info(f"    Initializing {optimizer_name} optimizer with kwargs: {default_kwargs}")
    
    # Create and run optimizer
    opt = OptClass(atoms, **default_kwargs)
    
    # Set default convergence criteria if not specified
    if 'fmax' not in default_kwargs:
        fmax = 0.02  # Default force convergence
        logger.info(f"    Using default fmax={fmax} eV/Å")
        opt.run(fmax=fmax)
    else:
        opt.run()


def _mark_done(step_name: str):
    path = f"STAGE_{step_name}_DONE"
    with open(path, 'w') as f:
        f.write(f"{step_name} completed\n")


def load_structure(pattern_initial_default: str = 'POSCAR') -> Atoms:
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