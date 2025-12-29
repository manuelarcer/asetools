# Description of the script:

import logging
import os
import shutil
import sys
import glob
import json
import numpy as np
from ase.io import read
from .calculatorsetuptools import VASPConfigurationFromYAML, deep_update, setup_initial_magmom
from ase.calculators.vasp import Vasp
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS, GPMin, MDMin, QuasiNewton
from ase.mep import MinModeTranslate
from ..analysis import check_outcar_convergence

logger = logging.getLogger(__name__)

# Reference file for atom order tracking across workflow restarts
MAGMOM_REFERENCE_FILE = '.asetools_magmom_reference.json'


def _save_magmom_reference(symbols, magmoms):
    """
    Save reference atom order and magmoms to file for workflow restart.

    Args:
        symbols: List of chemical symbols in original order
        magmoms: Original magmom list
    """
    data = {
        'reference_symbols': symbols,
        'original_magmoms': magmoms,
        'version': '1.0'
    }
    with open(MAGMOM_REFERENCE_FILE, 'w') as f:
        json.dump(data, f, indent=2)
    logger.info(f"  * Saved atom order reference to {MAGMOM_REFERENCE_FILE}")


def _load_magmom_reference():
    """
    Load reference atom order and magmoms from file.

    Returns:
        tuple: (reference_symbols, original_magmoms) or (None, None) if file doesn't exist
    """
    if not os.path.exists(MAGMOM_REFERENCE_FILE):
        return None, None

    try:
        with open(MAGMOM_REFERENCE_FILE, 'r') as f:
            data = json.load(f)
        logger.info(f"  * Loaded atom order reference from {MAGMOM_REFERENCE_FILE}")
        return data['reference_symbols'], data['original_magmoms']
    except Exception as e:
        logger.warning(f"  ⚠ Could not load {MAGMOM_REFERENCE_FILE}: {e}")
        return None, None


def _remove_magmom_reference():
    """Remove reference file after successful workflow completion."""
    if os.path.exists(MAGMOM_REFERENCE_FILE):
        os.remove(MAGMOM_REFERENCE_FILE)
        logger.debug(f"  * Removed {MAGMOM_REFERENCE_FILE}")


def _detect_atom_reordering(atoms, reference_symbols=None):
    """
    Detect if atoms have been reordered compared to a reference configuration.

    Args:
        atoms: Current ASE Atoms object
        reference_symbols: List of chemical symbols in original order

    Returns:
        tuple: (reordered, sort_indices, resort_indices)
            - reordered: True if atoms were reordered
            - sort_indices: Mapping from original to current order
            - resort_indices: Mapping from current back to original order
    """
    current_symbols = list(atoms.get_chemical_symbols())

    if reference_symbols is None or len(reference_symbols) != len(current_symbols):
        # Cannot detect reordering without reference
        return False, None, None

    # Check if order changed
    if current_symbols == reference_symbols:
        return False, None, None

    # Atoms were reordered - find the mapping
    # ASE VASP calculator sorts by element, then preserves relative order within each element
    from collections import defaultdict

    # Build mapping: for each element, track original indices
    element_indices_original = defaultdict(list)
    for i, symbol in enumerate(reference_symbols):
        element_indices_original[symbol].append(i)

    # Build mapping: for each element, track current indices
    element_indices_current = defaultdict(list)
    for i, symbol in enumerate(current_symbols):
        element_indices_current[symbol].append(i)

    # Create sort_indices: original_index -> current_index
    sort_indices = [None] * len(current_symbols)
    for symbol in element_indices_original.keys():
        orig_list = element_indices_original[symbol]
        curr_list = element_indices_current[symbol]
        for orig_pos, orig_idx in enumerate(orig_list):
            curr_idx = curr_list[orig_pos]
            sort_indices[orig_idx] = curr_idx

    # Create resort_indices: current_index -> original_index
    resort_indices = [None] * len(current_symbols)
    for orig_idx, curr_idx in enumerate(sort_indices):
        resort_indices[curr_idx] = orig_idx

    return True, sort_indices, resort_indices


def _reorder_magmom_list(magmom_list, sort_indices):
    """
    Reorder a magmom list to match reordered atoms.

    Args:
        magmom_list: Original magmom list
        sort_indices: Mapping from original to current order

    Returns:
        Reordered magmom list
    """
    if sort_indices is None:
        return magmom_list

    reordered = [0.0] * len(magmom_list)
    for orig_idx, curr_idx in enumerate(sort_indices):
        reordered[curr_idx] = magmom_list[orig_idx]

    return reordered


def _log_magmom_summary(atoms, prefix="  *"):
    """Log summary of magnetic moments on atoms object."""
    try:
        magmoms = atoms.get_initial_magnetic_moments()

        if len(magmoms) == 0:
            logger.info(f"{prefix} No magnetic moments set")
            return

        min_mag = magmoms.min()
        max_mag = magmoms.max()
        sum_mag = magmoms.sum()
        nonzero_count = np.sum(np.abs(magmoms) > 0.01)

        logger.info(f"{prefix} Magnetic moments: min={min_mag:.3f}, max={max_mag:.3f}, "
                   f"sum={sum_mag:.3f} μB")
        if nonzero_count > 0:
            logger.info(f"{prefix} Non-zero moments on {nonzero_count}/{len(magmoms)} atoms")
            # Show distribution of non-zero values
            nonzero_vals = magmoms[np.abs(magmoms) > 0.01]
            unique_vals = np.unique(np.round(nonzero_vals, 2))
            if len(unique_vals) <= 5:
                logger.info(f"{prefix} Unique moment values: {list(unique_vals)}")
        else:
            logger.info(f"{prefix} All magnetic moments are zero")

    except Exception as e:
        logger.debug(f"{prefix} Could not extract magnetic moment summary: {e}")


def _log_calculator_params(calc, prefix="  *"):
    """Log critical VASP calculator parameters."""
    try:
        params = calc.parameters if hasattr(calc, 'parameters') else {}

        # Critical parameters to log
        critical_params = {
            'ibrion': params.get('ibrion', 'default'),
            'nsw': params.get('nsw', 'default'),
            'ispin': params.get('ispin', 'default'),
            'ediff': params.get('ediff', 'default'),
            'ediffg': params.get('ediffg', 'default'),
        }

        logger.info(f"{prefix} VASP parameters: IBRION={critical_params['ibrion']}, "
                   f"NSW={critical_params['nsw']}, ISPIN={critical_params['ispin']}")
        logger.info(f"{prefix} Convergence: EDIFF={critical_params['ediff']}, "
                   f"EDIFFG={critical_params['ediffg']}")

        # Check for MAGMOM in calculator
        if 'magmom' in params:
            magmom_in_calc = params['magmom']
            if isinstance(magmom_in_calc, (list, tuple)):
                logger.info(f"{prefix} MAGMOM in calculator: {len(magmom_in_calc)} values")
            else:
                logger.info(f"{prefix} MAGMOM in calculator: {magmom_in_calc}")

    except Exception as e:
        logger.debug(f"{prefix} Could not extract calculator parameters: {e}")



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
        try:
            from vasp_interactive import VaspInteractive
        except ImportError:
            raise ImportError(
                "VaspInteractive calculator requested but 'vasp_interactive' package is not installed. "
                "Install it with: pip install vasp-interactive"
            )
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
        try:
            from vasp_interactive import VaspInteractive
        except ImportError:
            raise ImportError(
                "ASE optimizer requested but 'vasp_interactive' package is not installed. "
                "Install it with: pip install vasp-interactive"
            )
        calc = VaspInteractive(**vasp_kwargs)
        logger.info(" ** VaspInteractive calculator created for ASE optimizer")
    else:
        calc = Vasp(**vasp_kwargs)
        logger.info(" ** Regular Vasp calculator created")

    return calc

def run_workflow(atoms: Atoms, cfg: VASPConfigurationFromYAML, workflow_name: str, run_overrides: dict = None, dry_run: bool = False, magmoms=None):
    """
    Run a multi-stage VASP workflow.

    Args:
        atoms: ASE Atoms object to optimize
        cfg: Configuration from YAML file
        workflow_name: Name of workflow to run (from YAML)
        run_overrides: Dict of VASP parameters to override
        dry_run: If True, skip actual calculations
        magmoms: Optional list of per-atom magnetic moments.
            - If provided, overrides YAML-based element magmoms
            - If shorter than len(atoms), remaining atoms set to 0.0
            - If longer, raises ValueError
            - Example: magmoms = 32*[0] + 64*[2.0] + 16*[0]
    """

    # Determine magnetic moment source: runtime parameter overrides YAML config
    if magmoms is not None:
        initial_magmom = magmoms  # Use runtime list
        logger.info(f"Using runtime magmoms parameter (list with {len(magmoms)} values)")
    else:
        initial_magmom = cfg.initial_magmom_data  # Use YAML dict
        if initial_magmom:
            logger.info(f"Using YAML-based element magmoms: {initial_magmom}")
        else:
            logger.info("No magnetic moments specified (YAML or runtime)")

    # Handle reference atom order for reordering detection (only needed for list-based magmoms)
    reference_symbols = None
    if isinstance(initial_magmom, (list, tuple, np.ndarray)):
        # Check if reference file exists (workflow restart)
        saved_ref, saved_magmoms = _load_magmom_reference()

        if saved_ref is not None:
            # Workflow restart: use saved reference
            reference_symbols = saved_ref
            # Verify saved magmoms match what user provided
            if list(saved_magmoms) != list(initial_magmom):
                logger.warning(
                    f"  ⚠ Magmom list changed from previous run!\n"
                    f"    Saved: {saved_magmoms[:5]}...\n"
                    f"    Current: {list(initial_magmom)[:5]}...\n"
                    f"    Using current magmoms but will remap based on saved atom order"
                )
        else:
            # First run: store current atom order as reference
            reference_symbols = list(atoms.get_chemical_symbols())
            _save_magmom_reference(reference_symbols, list(initial_magmom))
            logger.info(f"  * Stored reference atom order for magmom mapping ({len(reference_symbols)} atoms)")

    to_run = stages_to_run(cfg, workflow_name)
    stages = cfg.workflows[workflow_name]['stages']

    if not to_run:
        logger.info(f"-->  Workflow '{workflow_name}' is already completed, nothing to do  <--")
        return

    for stage in stages:
        if stage['name'] not in to_run:
            logger.info(f"Skipping STAGE: {stage['name']}, already done")
            continue
        _run_stage(atoms, cfg, stage, run_overrides, dry_run,
                   initial_magmom=initial_magmom, reference_symbols=reference_symbols)

    # Workflow completed successfully - clean up reference file
    if isinstance(initial_magmom, (list, tuple, np.ndarray)):
        _remove_magmom_reference()

    logger.info(f"-->  Workflow '{workflow_name}' completed successfully  <--")
    

def _run_stage(atoms: Atoms, cfg: VASPConfigurationFromYAML, stage: dict, run_overrides: dict, dry_run: bool, initial_magmom, reference_symbols=None):
    # run_overrides is different from the overrides in each step
    name  = stage['name']
    steps = stage['steps']
    logger.info(f"Running STAGE: {stage['name']}")

    # Detect and handle atom reordering for list-based magmoms
    magmom_to_apply = initial_magmom
    if isinstance(initial_magmom, (list, tuple, np.ndarray)) and reference_symbols is not None:
        reordered, sort_idx, _ = _detect_atom_reordering(atoms, reference_symbols)
        if reordered:
            logger.warning(f"  ⚠ Atoms were reordered by VASP - remapping magmom list to match current order")
            magmom_to_apply = _reorder_magmom_list(list(initial_magmom), sort_idx)
            logger.info(f"  * Magmom list remapped for reordered structure")

    # Apply constraints if specified in stage configuration
    if 'constraints' in stage:
        from ..constraints import ConstraintManager
        logger.info(f"  * Applying constraints from stage configuration")
        cm = ConstraintManager()
        cm.apply_stage_constraints(atoms, stage['constraints'])

    # Track if any step used ASE optimizer (for proper convergence checking)
    ase_optimizer_converged = None

    for step in steps:
        # Make calculator for each step based on step requirements
        logger.info('  * Setting up calculator from config')

        # Check if step needs VaspInteractive with proper context management
        if step.get('optimizer') is not None:
            opt_converged = _run_step_with_vaspinteractive(atoms, cfg, step, run_overrides, magmom_to_apply, dry_run)
            if opt_converged is not None:
                ase_optimizer_converged = opt_converged
        else:
            # Regular Vasp calculator
            calc = _make_step_calculator(cfg, step, run_overrides)
            atoms.calc = calc
            _log_calculator_params(calc, prefix="  *")
            atoms = setup_initial_magmom(atoms, magmom_to_apply)
            _log_magmom_summary(atoms, prefix="  *")
            _run_step(atoms, step, dry_run)

    # Check convergence before marking as done
    if not dry_run:
        # If an ASE optimizer was used, trust its convergence status
        # (ASE optimizers use constraint-adjusted forces, not raw VASP forces)
        if ase_optimizer_converged is not None:
            convergence = ase_optimizer_converged
            if convergence:
                logger.info(f"  ✓ ASE optimizer converged successfully")
            else:
                logger.error(f" ❌ Stage '{name}' did NOT converge - STAGE_*_DONE file will NOT be created")
                logger.error(f"    ASE optimizer did not reach convergence criteria")
                raise RuntimeError(f"Stage '{name}' did not converge. Fix the issue and re-run.")
        else:
            # No ASE optimizer used, check OUTCAR convergence (standard VASP optimization)
            convergence, vasp_version = check_outcar_convergence('OUTCAR', verbose=False)

            # Log convergence details if successful
            if convergence:
                try:
                    final_atoms = read('OUTCAR', format='vasp-out', index=-1)
                    forces = final_atoms.get_forces()
                    max_force = np.linalg.norm(forces, axis=1).max()
                    logger.info(f"  ✓ VASP converged: max_force={max_force:.4f} eV/Å")
                except Exception:
                    logger.info(f"  ✓ VASP converged")

            if not convergence:
                logger.error(f" ❌ Stage '{name}' did NOT converge - STAGE_*_DONE file will NOT be created")
                logger.error(f"    The calculation likely hit NSW limit or failed to meet convergence criteria")
                raise RuntimeError(f"Stage '{name}' did not converge. Fix the issue and re-run.")

    backup_output_files(name=name)

    # Stage summary including magnetic moments
    if not dry_run:
        try:
            final_atoms = read('OUTCAR', format='vasp-out', index=-1)
            final_magmoms = final_atoms.get_magnetic_moments()
            logger.info(f" -- Stage '{name}' magnetic summary:")
            logger.info(f"    Total magnetic moment: {final_magmoms.sum():.3f} μB")
            logger.info(f"    Max |moment|: {np.abs(final_magmoms).max():.3f} μB")
        except Exception:
            logger.debug(f" -- Could not extract final magnetic moments from OUTCAR")

    logger.info(f" -- ✅ Stage '{name}' completed, converged, and backed up")
    _mark_done(name)


def _run_step_with_vaspinteractive(atoms: Atoms, cfg: VASPConfigurationFromYAML, step: dict, run_overrides: dict, initial_magmom, dry_run: bool):
    """
    Run a step that requires VaspInteractive using the documented with-clause pattern.
    This ensures proper process management and prevents hanging.

    Returns:
        bool or None: Optimizer convergence status (True/False) or None if dry_run
    """
    try:
        from vasp_interactive import VaspInteractive
    except ImportError:
        raise ImportError(
            "ASE optimizer requested but 'vasp_interactive' package is not installed. "
            "Install it with: pip install vasp-interactive"
        )

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
        return None

    # Use the documented with-clause pattern
    with VaspInteractive(**vasp_kwargs) as calc:
        atoms.calc = calc
        _log_calculator_params(calc, prefix="    ")
        atoms = setup_initial_magmom(atoms, initial_magmom)
        _log_magmom_summary(atoms, prefix="    ")

        # Run ASE optimizer
        optimizer_name = step.get('optimizer')
        optimizer_kwargs = step.get('optimizer_kwargs', {})

        logger.info(f"    Running ASE optimizer {optimizer_name} with VaspInteractive")
        converged = _run_with_ase_optimizer(atoms, optimizer_name, optimizer_kwargs)

    logger.info(f"    ✔ VaspInteractive step '{step['name']}' completed and finalized")
    return converged


def _run_step(atoms: Atoms, step: dict, dry_run: bool):
    name      = step['name']
    overrides = step.get('overrides', {})
    optimizer_name = step.get('optimizer', None)
    optimizer_kwargs = step.get('optimizer_kwargs', {})
    
    logger.info(f"  • Step '{name}' overrides={overrides}")

    # For VaspInteractive, we need to ensure proper parameter handling
    try:
        from vasp_interactive import VaspInteractive
        is_vasp_interactive = isinstance(atoms.calc, VaspInteractive)
    except ImportError:
        is_vasp_interactive = False

    if is_vasp_interactive:
        # Handle VaspInteractive-specific parameter conflicts
        overrides = _fix_vaspinteractive_params(overrides, optimizer_name is not None)
        logger.info(f"    VaspInteractive adjusted overrides: {overrides}")

    atoms.calc.set(**overrides)
    
    if dry_run:
        logger.info("    (dry-run, skipping calculation)")
        return

    # Check if we should use ASE optimizer (VaspInteractive only)
    try:
        from vasp_interactive import VaspInteractive
        is_vasp_interactive = isinstance(atoms.calc, VaspInteractive)
    except ImportError:
        is_vasp_interactive = False

    if optimizer_name and is_vasp_interactive:
        logger.info(f"    Using ASE optimizer: {optimizer_name}")
        _run_with_ase_optimizer(atoms, optimizer_name, optimizer_kwargs)
    else:
        # trigger VASP run (regular optimization or single point)
        logger.info("    Starting VASP calculation...")
        try:
            energy = atoms.get_potential_energy()
            logger.info(f"    Calculation completed, energy: {energy:.6f} eV")

            # Report final magnetic moments from VASP calculation
            try:
                final_magmoms = atoms.get_magnetic_moments()
                total_mag = final_magmoms.sum()
                max_mag = np.abs(final_magmoms).max()
                logger.info(f"    Final magnetic moments: total={total_mag:.3f} μB, "
                           f"max={max_mag:.3f} μB")
            except Exception:
                logger.debug("    Could not extract final magnetic moments")

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

    Returns:
        bool: True if optimizer converged, False otherwise
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
        'fire': {'dt', 'maxstep', 'maxmove', 'dtmax', 'Nmin', 'finc', 'fdec', 'astart', 'fa'},
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
    logger.info(f"    Running ASE {optimizer_name} optimization...")

    # Verify VaspInteractive calculator is present
    try:
        from vasp_interactive import VaspInteractive
        if not isinstance(atoms.calc, VaspInteractive):
            raise RuntimeError("Expected VaspInteractive calculator for ASE optimization")
    except ImportError:
        raise ImportError(
            "ASE optimizer requested but 'vasp_interactive' package is not installed. "
            "Install it with: pip install vasp-interactive"
        )

    # Create and run optimizer with existing VaspInteractive calculator
    opt = OptClass(atoms, **init_kwargs)
    converged = opt.run(**run_kwargs)

    logger.info(f"    ASE {optimizer_name} optimization completed")

    # Check convergence status
    # ASE optimizers return True if converged, False if hit step limit
    if converged:
        logger.info(f"    ✓ Optimizer converged to fmax={run_kwargs['fmax']} eV/Å")
    else:
        logger.warning(f"    ⚠ Optimizer did not converge (hit step limit)")

    return converged


def _run_with_dimer_optimizer(atoms: Atoms, optimizer_kwargs: dict):
    """
    Run dimer optimization using ASE dimer method with VaspInteractive calculator.
    Handles MODECAR file detection and MinModeAtoms setup.

    Returns:
        bool: True if dimer converged, False otherwise
    """
    from ..dimer import setup_dimer_atoms, validate_dimer_kwargs

    # Validate and separate dimer-specific kwargs
    dimer_control_kwargs, optimizer_init_kwargs, optimizer_run_kwargs = validate_dimer_kwargs(optimizer_kwargs)

    logger.info(f"    Setting up dimer calculation")
    logger.info(f"    Dimer control kwargs: {dimer_control_kwargs}")
    logger.info(f"    Optimizer init kwargs: {optimizer_init_kwargs}")
    logger.info(f"    Optimizer run kwargs: {optimizer_run_kwargs}")

    # Check for VaspInteractive calculator
    try:
        from vasp_interactive import VaspInteractive
        if not isinstance(atoms.calc, VaspInteractive):
            raise RuntimeError("Expected VaspInteractive calculator for dimer optimization")
    except ImportError:
        raise ImportError(
            "Dimer optimizer requested but 'vasp_interactive' package is not installed. "
            "Install it with: pip install vasp-interactive"
        )

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
    converged = dimer_opt.run(**optimizer_run_kwargs)

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
    dimer_converged = False
    try:
        from ..dimer import check_dimer_convergence, extract_saddle_point_info
        conv_info = check_dimer_convergence(d_atoms)
        dimer_converged = conv_info.get('converged', False)
        if dimer_converged:
            saddle_info = extract_saddle_point_info(d_atoms)
            logger.info(f"    ✓ Converged to saddle point: E={saddle_info.get('energy', 'N/A'):.6f} eV")
            logger.info(f"    Final eigenvalue: {saddle_info.get('eigenvalue', 'N/A')}")
        else:
            logger.warning(f"    Dimer may not have converged properly")
            logger.warning(f"    Final max force: {conv_info.get('max_force', 'N/A'):.4f} eV/Å")
    except Exception as e:
        logger.warning(f"    Could not check convergence: {e}")

    # Use ASE optimizer convergence if dimer-specific check not available
    if converged is not None:
        return converged
    return dimer_converged


def _mark_done(step_name: str):
    path = f"STAGE_{step_name}_DONE"
    with open(path, 'w') as f:
        f.write(f"{step_name} completed\n")


def load_structure(pattern_initial_default: str = 'POSCAR') -> Atoms:
    """Load the structure from OUTCAR (preferred), CONTCAR, or initial configuration.

    The function attempts to load structures in the following priority order:
    1. OUTCAR (last ionic step, index=-1) - Most reliably updated during calculations
    2. CONTCAR - Fallback for compatibility with existing workflows
    3. Pattern-matched file (default: POSCAR) - Initial structure for first stage

    This approach ensures reliable workflow continuation even when CONTCAR is not
    properly updated before job crashes.

    Args:
        pattern_initial_default: Glob pattern for initial structure file (default: 'POSCAR')

    Returns:
        Atoms: ASE Atoms object containing the loaded structure
    """
    # Try OUTCAR first (most reliable)
    if os.path.exists('OUTCAR'):
        try:
            atoms = read('OUTCAR', format='vasp-out', index=-1)
            logger.info("Loading structure from OUTCAR (last configuration)")
            return atoms
        except Exception as e:
            logger.warning(f"Failed to read OUTCAR: {e}. Trying CONTCAR...")

    # Fallback to CONTCAR
    if os.path.exists('CONTCAR'):
        try:
            atoms = read('CONTCAR', format='vasp')
            logger.info("Loading structure from CONTCAR")
            return atoms
        except Exception as e:
            logger.warning(f"Failed to read CONTCAR: {e}. Trying initial pattern...")

    # Final fallback to initial configuration
    matches = glob.glob(pattern_initial_default)
    if not matches:
        logger.error(f"No OUTCAR, CONTCAR, or match for pattern: {pattern_initial_default}")
        sys.exit(1)

    atoms = read(matches[0])
    logger.info(f"Loading structure from initial file: {matches[0]}")
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