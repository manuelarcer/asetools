# ASE/calculator setup utilities

import logging
import yaml

logger = logging.getLogger(__name__)

class VASPConfigurationFromYAML:
    def __init__(self, config_file: str, system: str = 'default'):
        self.config = load_yaml_config(config_file)
        self.system = system
        self.basic_config = self.config['basic']
        self.workflows = self.config['workflows']
        self.globals = self.config['globals']
        
        self.initial_magmom_data = self.initial_magmom()
        
        # Remove the magmom shorthand so it won't sneak into system_config
        systems = self.config.get('systems', {})
        if system in systems and systems[system] is not None:
            systems[system].pop('magmom', None)
            

    @property
    def system_config(self) -> dict:
        try:
            system_dict = self.config['systems'][self.system]
            if system_dict is None:
                logger.warning(f"System configuration for '{self.system}' is None. Returning an empty dictionary.")
                return {}
        except KeyError:
            logger.warning(f"System configuration for '{self.system}' not found. Returning an empty dictionary.")
            return {}
        return system_dict

    def initial_magmom(self) -> dict:
        system_config = self.system_config
        if system_config is None:
            return {}
        if 'magmom' in system_config:
            return system_config['magmom']
        elif 'initial_magmom' in system_config:
            return system_config['initial_magmom']
        else:
            return {}

def load_yaml_config(config_file: str) -> dict:
    """
    Load a YAML configuration file and return its contents as a dictionary.
    The YAML file should contain the following 1st level keys:
           - basic: for basic VASP settings
           - systems: for system-specific settings, eg, NCA, LFP
           - workflows: to specify the type of job, eg. bulk_opt, slab_opt, etc
               - each workflow should have 'stages' containing a name and a dict with overrides
           - globals: to define settings like VASP_PP_PATH or initial_conf_pattern
    """
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)
    verify_configuration_keys(cfg)   # Verify keys
    return cfg

def verify_configuration_keys(cfg: dict) -> None:
    required_keys = ['basic', 'systems', 'workflows', 'globals']
    for key in required_keys:
        if key not in cfg:
            raise KeyError(f"Configuration file is missing required key: '{key}'")

def deep_update(base: dict, override: dict):
    """
    For each (k, v) in override:
      - if base[k] is also a dict, recurse into it
      - else, replace base[k] with v
    """
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base

def setup_initial_magmom(atoms, magmom_dict: dict):
    """
    Set the initial magnetic moments based on the provided dictionary.
    The dictionary should be provided in the vasp_parameters.yaml file
    under the 'systems' section
    """
    for at in atoms:
        if at.symbol in  magmom_dict:
            at.magmom = magmom_dict[at.symbol]
        else:
            at.magmom = 0.0
    return atoms

