# Description of the script:

import logging
import os
import shutil
import sys
import glob
from ase.io import read
import yaml
from .calculatorsetuptools import VASPConfigurationFromYAML

logger = logging.getLogger(__name__)    

def backup_output_files(name='backup'):
    for fname in ('POSCAR','INCAR','CONTCAR','OUTCAR','OSZICAR'):
        dst = f"{fname}_{name}"
        try:
            shutil.copy(fname, dst)
        except FileNotFoundError:
            logger.warning(f"File {fname} not found, skipping backup.")

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

# TODO: DO WE NEED THIS?
def load_last_logger(logger_basename: str = 'logger'):
    """Load the last logger configuration from the log file."""
    log_files = glob.glob(f"{logger_basename}*")
    if not log_files:
        logger.warning(f"No log files found for logger: {logger_basename}")
        return None
    # Sort log files by modification time
    log_files.sort(key=os.path.getmtime)
    with open(log_files[-1], 'r') as file:
        return file.read()

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
            logger.info(f"STAGE {stage_name} is NOT DONE, adding to run list.")
            to_run.append(stage_name)
    return to_run
        