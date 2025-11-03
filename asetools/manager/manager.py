# Description of the script:

import logging
import os
import shutil
import sys
import glob
import numpy as np
from ase.io import read
from .calculatorsetuptools import VASPConfigurationFromYAML, deep_update, setup_initial_magmom
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS, GPMin, MDMin, QuasiNewton
from ase.mep import MinModeTranslate
from ..analysis import check_outcar_convergence

logger = logging.getLogger(__name__)    

def make_calculator(cfg: VASPConfigurationFromYAML, run_overrides: dict = None):
    if run_overrides is None:
        run_overrides = {}
    else:
        logger.info(f" ** Overriding run parameters with: {run_overrides}")
    
    # Determine calculator type from configuration
    calculator_type = cfg.globals.get('calculator_type', 'vasp')
    
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


def _make_step_calculator(cfg: VASPConfigurationFromYAML, step: dict, run_overrides: dict = None):
    """
    Create calculator for a specific step, choosing between Vasp and VaspInteractive
    based on step requirements:
    
    - VaspInteractive: when step has 'optimizer' field (ASE optimizers)
    - Regular Vasp: for single point calculations and VASP internal optimization
    """
    if run_overrides is None:
        run_overrides = {}
    
    # Check if step requires VaspInteractive
    needs_interactive = step.get('optimizer') is not None
    
    # Build calculator parameters
    vasp_kwargs = deep_update(
        deep_update(cfg.basic_config.copy(), cfg.system_config),
        run_overrides
    )
    
    if needs_interactive:
        calc = VaspInteractive(**vasp_kwargs)
        logger.info(" ** VaspInteractive calculator created for ASE optimizer")
    else:
        calc = Vasp(**vasp_kwargs)
        logger.info(" ** Regular Vasp calculator created")
    
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
        # Make calculator for each step based on step requirements
        logger.info('  * Setting up calculator from config')
        
        # Check if step needs VaspInteractive with proper context management
        if step.get('optimizer') is not None:
            _run_step_with_vaspinteractive(atoms, cfg, step, run_overrides, initial_magmom, dry_run)
        else:
            # Regular Vasp calculator
            calc = _make_step_calculator(cfg, step, run_overrides)
            atoms.calc = calc
            atoms = setup_initial_magmom(atoms, magmom_dict=initial_magmom)
            _run_step(atoms, step, dry_run)

    # Check convergence before marking as done
    if not dry_run:
        convergence, vasp_version = check_outcar_convergence('OUTCAR', verbose=False)
        if not convergence:
            logger.error(f" ❌ Stage '{name}' did NOT converge - STAGE_*_DONE file will NOT be created")
            logger.error(f"    The calculation likely hit NSW limit or failed to meet convergence criteria")
            raise RuntimeError(f"Stage '{name}' did not converge. Fix the issue and re-run.")

    backup_output_files(name=name)
    logger.info(f" -- ✅ Stage '{name}' completed, converged, and backed up")
    _mark_done(name)


def _run_step_with_vaspinteractive(atoms: Atoms, cfg: VASPConfigurationFromYAML, step: dict, run_overrides: dict, initial_magmom: dict, dry_run: bool):
    """
    Run a step that requires VaspInteractive using the documented with-clause pattern.
    This ensures proper process management and prevents hanging.
    """
    # Build VaspInteractive parameters
    vasp_kwargs = deep_update(
        deep_update(cfg.basic_config.copy(), cfg.system_config),
        run_overrides or {}
    )
    
    # Apply step overrides and VaspInteractive parameter fixes
    overrides = step.get('overrides', {})
    overrides = _fix_vaspinteractive_params(overrides, using_ase_optimizer=True)
    
    # Merge overrides into vasp_kwargs
    vasp_kwargs.update(overrides)
    
    logger.info(f" ** Creating VaspInteractive with parameters for ASE optimizer")
    logger.info(f"    VaspInteractive parameters: {overrides}")
    
    if dry_run:
        logger.info("    (dry-run, skipping VaspInteractive calculation)")
        return
    
    # Use the documented with-clause pattern
    with VaspInteractive(**vasp_kwargs) as calc:
        atoms.calc = calc
        atoms = setup_initial_magmom(atoms, magmom_dict=initial_magmom)
        
        # Run ASE optimizer
        optimizer_name = step.get('optimizer')
        optimizer_kwargs = step.get('optimizer_kwargs', {})
        
        logger.info(f"    Running ASE optimizer {optimizer_name} with VaspInteractive")
        _run_with_ase_optimizer(atoms, optimizer_name, optimizer_kwargs)
        
    logger.info(f"    ✔ VaspInteractive step '{step['name']}' completed and finalized")


def _run_step(atoms: Atoms, step: dict, dry_run: bool):
    name      = step['name']
    overrides = step.get('overrides', {})
    optimizer_name = step.get('optimizer', None)
    optimizer_kwargs = step.get('optimizer_kwargs', {})
    
    logger.info(f"  • Step '{name}' overrides={overrides}")
    
    # For VaspInteractive, we need to ensure proper parameter handling
    if isinstance(atoms.calc, VaspInteractive):
        # Handle VaspInteractive-specific parameter conflicts
        overrides = _fix_vaspinteractive_params(overrides, optimizer_name is not None)
        logger.info(f"    VaspInteractive adjusted overrides: {overrides}")
    
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
        logger.info("    Starting VASP calculation...")
        try:
            energy = atoms.get_potential_energy()
            logger.info(f"    Calculation completed, energy: {energy:.6f} eV")
        except Exception as e:
            logger.error(f"    VASP calculation failed: {e}")
            raise
    
    logger.info(f"    ✔ Step '{name}' done")


def _fix_vaspinteractive_params(overrides: dict, using_ase_optimizer: bool = False) -> dict:
    """
    Fix parameter conflicts for VaspInteractive calculator based on documentation:
    
    1. IBRION doesn't enable VASP's internal optimizer - VaspInteractive handles this
    2. NSW should be >= max steps for ASE optimizers (default 2000) 
    3. For single point (NSW=0): disable interactive mode to avoid stdin hang
    """
    fixed_overrides = overrides.copy()
    
    # Critical fix for single point calculations hanging at stdin
    if overrides.get('nsw', 0) == 0:
        # Single point - disable interactive mode to avoid "reading from stdin" hang
        fixed_overrides['ibrion'] = -1  # No ionic optimization
        logger.debug("    Single point: set IBRION=-1 to avoid VaspInteractive stdin hang")
    
    # For ASE optimizers, ensure proper NSW and IBRION settings
    if using_ase_optimizer:
        # Ensure NSW is high enough for ASE optimizer steps (VaspInteractive requirement)
        if 'nsw' not in fixed_overrides or fixed_overrides.get('nsw', 0) < 100:
            fixed_overrides['nsw'] = 1000  # VaspInteractive default
        fixed_overrides['ibrion'] = -1   # VaspInteractive handles optimization externally
        logger.debug(f"    ASE optimizer: NSW={fixed_overrides['nsw']}, IBRION=-1")
    
    return fixed_overrides


def _run_with_ase_optimizer(atoms: Atoms, optimizer_name: str, optimizer_kwargs: dict):
    """
    Run optimization using ASE optimizers with VaspInteractive calculator.
    Uses the documented with-clause pattern for proper process management.
    """
    
    # Map optimizer names to classes
    optimizer_map = {
        'bfgs': BFGS,
        'fire': FIRE,
        'lbfgs': LBFGS,
        'gpmin': GPMin,
        'mdmin': MDMin,
        'quasinewton': QuasiNewton,
        'dimer': MinModeTranslate,
    }
    
    optimizer_name_lower = optimizer_name.lower()
    if optimizer_name_lower not in optimizer_map:
        raise ValueError(f"Unknown optimizer: {optimizer_name}. Available: {list(optimizer_map.keys())}")

    # Special handling for dimer optimizer
    if optimizer_name_lower == 'dimer':
        return _run_with_dimer_optimizer(atoms, optimizer_kwargs)

    OptClass = optimizer_map[optimizer_name_lower]
    
    # Separate initialization kwargs from run kwargs
    init_kwargs = {'logfile': f'{optimizer_name.upper()}.log'}
    run_kwargs = {}
    
    # Define which parameters go to run() vs __init__() for each optimizer
    run_only_params = {'fmax', 'steps'}  # These always go to run()
    
    # Optimizer-specific parameter routing
    optimizer_init_params = {
        'bfgs': {'maxstep', 'alpha', 'damping'},
        'fire': {'dt', 'maxmove', 'dtmax', 'Nmin', 'finc', 'fdec', 'astart', 'fa'},
        'lbfgs': {'maxstep', 'memory', 'damping', 'alpha'},
    }
    
    # Route parameters correctly
    current_optimizer_init_params = optimizer_init_params.get(optimizer_name_lower, set())
    
    for key, value in optimizer_kwargs.items():
        if key in run_only_params:
            run_kwargs[key] = value
        elif key in current_optimizer_init_params:
            init_kwargs[key] = value
        else:
            # Default to run kwargs for unknown parameters
            run_kwargs[key] = value
    
    logger.info(f"    Initializing {optimizer_name} optimizer with init_kwargs: {init_kwargs}")
    logger.info(f"    Running optimization with run_kwargs: {run_kwargs}")
    
    # Set default convergence criteria if not specified
    if 'fmax' not in run_kwargs:
        run_kwargs['fmax'] = 0.02  # Default force convergence
        logger.info(f"    Using default fmax={run_kwargs['fmax']} eV/Å")
    
    # Use the documented pattern: VaspInteractive should already be initialized
    # Just create and run the ASE optimizer with the existing calculator
    if not isinstance(atoms.calc, VaspInteractive):
        raise RuntimeError("Expected VaspInteractive calculator for ASE optimization")
    
    logger.info(f"    Running ASE {optimizer_name} optimization...")
    
    # Create and run optimizer with existing VaspInteractive calculator
    opt = OptClass(atoms, **init_kwargs)
    opt.run(**run_kwargs)
    
    logger.info(f"    ASE {optimizer_name} optimization completed")


def _run_with_dimer_optimizer(atoms: Atoms, optimizer_kwargs: dict):
    """
    Run dimer optimization using ASE dimer method with VaspInteractive calculator.
    Handles MODECAR file detection and MinModeAtoms setup.
    """
    from ..dimer import setup_dimer_atoms, validate_dimer_kwargs

    # Validate and separate dimer-specific kwargs
    dimer_control_kwargs, optimizer_init_kwargs, optimizer_run_kwargs = validate_dimer_kwargs(optimizer_kwargs)

    logger.info(f"    Setting up dimer calculation")
    logger.info(f"    Dimer control kwargs: {dimer_control_kwargs}")
    logger.info(f"    Optimizer init kwargs: {optimizer_init_kwargs}")
    logger.info(f"    Optimizer run kwargs: {optimizer_run_kwargs}")

    # Check for VaspInteractive calculator
    if not isinstance(atoms.calc, VaspInteractive):
        raise RuntimeError("Expected VaspInteractive calculator for dimer optimization")

    # Handle displacement vector
    displacement_vector = None
    if 'displacement_vector' in optimizer_kwargs:
        displacement_vector = np.array(optimizer_kwargs['displacement_vector'])
        if displacement_vector.shape != (len(atoms), 3):
            raise ValueError(f"Displacement vector shape {displacement_vector.shape} doesn't match atoms shape ({len(atoms)}, 3)")

    # Set up MinModeAtoms with dimer control
    d_atoms = setup_dimer_atoms(atoms, displacement_vector, dimer_control_kwargs)

    # Set default convergence criteria if not specified
    if 'fmax' not in optimizer_run_kwargs:
        optimizer_run_kwargs['fmax'] = 0.01  # Stricter default for dimer
        logger.info(f"    Using default dimer fmax={optimizer_run_kwargs['fmax']} eV/Å")

    # Set default logfile if not specified
    if 'logfile' not in optimizer_init_kwargs:
        optimizer_init_kwargs['logfile'] = 'DIMER.log'

    logger.info(f"    Running dimer optimization...")

    # Create and run dimer optimizer
    from ase.mep import MinModeTranslate
    dimer_opt = MinModeTranslate(d_atoms, **optimizer_init_kwargs)
    dimer_opt.run(**optimizer_run_kwargs)

    logger.info(f"    Dimer optimization completed")

    # Save final MODECAR file with converged eigenvector
    try:
        from ..dimer import write_modecar
        if hasattr(d_atoms, 'eigenvector') and d_atoms.eigenvector is not None:
            write_modecar(d_atoms.eigenvector, atoms, 'MODECAR_final')
            logger.info("    Saved final eigenvector to MODECAR_final")
    except Exception as e:
        logger.warning(f"    Could not save final MODECAR: {e}")

    # Log convergence information
    try:
        from ..dimer import check_dimer_convergence, extract_saddle_point_info
        conv_info = check_dimer_convergence(d_atoms)
        if conv_info.get('converged', False):
            saddle_info = extract_saddle_point_info(d_atoms)
            logger.info(f"    ✓ Converged to saddle point: E={saddle_info.get('energy', 'N/A'):.6f} eV")
            logger.info(f"    Final eigenvalue: {saddle_info.get('eigenvalue', 'N/A')}")
        else:
            logger.warning(f"    Dimer may not have converged properly")
            logger.warning(f"    Final max force: {conv_info.get('max_force', 'N/A'):.4f} eV/Å")
    except Exception as e:
        logger.warning(f"    Could not check convergence: {e}")


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