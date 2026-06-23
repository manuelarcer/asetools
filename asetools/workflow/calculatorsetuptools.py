# ASE/calculator setup utilities

import copy
import logging

import numpy as np
import yaml

logger = logging.getLogger(__name__)


class VASPConfigurationFromYAML:
    def __init__(self, config_file: str, system: str = "default"):
        self.config = load_yaml_config(config_file)
        self.system = system
        self.basic_config = self.config["basic"]
        self.workflows = self.config["workflows"]
        self.globals = self.config["globals"]
        # stage_templates is optional (modular stage composition); empty if absent
        self.stage_templates = self.config.get("stage_templates", {})

        self.initial_magmom_data = self.initial_magmom()

        # Remove the magmom shorthand so it won't sneak into system_config
        systems = self.config.get("systems", {})
        if system in systems and systems[system] is not None:
            systems[system].pop("magmom", None)

    @property
    def system_config(self) -> dict:
        try:
            system_dict = self.config["systems"][self.system]
            if system_dict is None:
                logger.warning(
                    f"System configuration for '{self.system}' is None. Returning an empty dictionary."
                )
                return {}
        except KeyError:
            logger.warning(
                f"System configuration for '{self.system}' not found. Returning an empty dictionary."
            )
            return {}
        return system_dict

    def initial_magmom(self) -> dict:
        system_config = self.system_config
        if system_config is None:
            return {}
        if "magmom" in system_config:
            return system_config["magmom"]
        elif "initial_magmom" in system_config:
            return system_config["initial_magmom"]
        else:
            return {}

    # MLIP stage fields copied verbatim from a template into the resolved stage.
    _MLIP_FIELDS = (
        "mlip",
        "env",
        "optimizer",
        "fmax",
        "max_steps",
        "relax_cell",
        "device",
        "uma_task",
        "mace_head",
    )

    def stage(
        self,
        template: str,
        *,
        name: str,
        overrides: dict = None,
        optimizer_kwargs: dict = None,
        constraints: dict = None,
    ) -> dict:
        """Instantiate a stage template into a concrete, resolved stage dict.

        The returned dict matches the shape the workflow manager consumes
        (``{name, engine, constraints?, steps?}``) plus a private ``_template``
        provenance key (ignored at execution, used for the production-stage
        warning and the MLIP-ordering guard).

        Per-instance ``overrides`` / ``optimizer_kwargs`` / ``constraints`` are
        deep-merged onto the template's baked-in values (patch only what
        differs). For VASP stages they apply to *every* step. For MLIP stages
        (``engine: mlip``) ``overrides`` is rejected (VASP-only); tune the MLIP
        optimization via ``optimizer_kwargs`` (e.g. ``{'fmax': 0.2}``).

        Args:
            template: Key into the ``stage_templates`` YAML section.
            name: Required, explicit stage name (sentinel / backup / restart key).
            overrides: VASP parameter patch, merged into each step's overrides.
            optimizer_kwargs: Optimizer-parameter patch.
            constraints: Patch onto the template's constraint block.

        Raises:
            KeyError: If ``template`` is not in ``stage_templates``.
            ValueError: If ``overrides`` is passed to an MLIP stage.
        """
        if template not in self.stage_templates:
            raise KeyError(
                f"Stage template '{template}' not found in 'stage_templates' section "
                f"(available: {sorted(self.stage_templates)})"
            )

        body = copy.deepcopy(self.stage_templates[template])
        engine = body.get("engine", "vasp")
        resolved = {"name": name, "_template": template, "engine": engine}

        # Constraints: template default deep-merged with the per-instance patch.
        tmpl_constraints = body.get("constraints")
        if tmpl_constraints is not None or constraints is not None:
            merged = copy.deepcopy(tmpl_constraints) if tmpl_constraints else {}
            if constraints:
                deep_update(merged, copy.deepcopy(constraints))
            resolved["constraints"] = merged

        if engine == "mlip":
            if overrides:
                raise ValueError(
                    "overrides apply to VASP stages; tune an MLIP stage via "
                    "optimizer_kwargs (e.g. optimizer_kwargs={'fmax': 0.2})"
                )
            for key in self._MLIP_FIELDS:
                if key in body:
                    resolved[key] = body[key]
            if optimizer_kwargs:
                deep_update(resolved, copy.deepcopy(optimizer_kwargs))
            return resolved

        # VASP engine: deep-merge the per-instance patches into every step.
        steps = copy.deepcopy(body.get("steps", []))
        for step in steps:
            if overrides:
                deep_update(step.setdefault("overrides", {}), copy.deepcopy(overrides))
            if optimizer_kwargs:
                deep_update(
                    step.setdefault("optimizer_kwargs", {}),
                    copy.deepcopy(optimizer_kwargs),
                )
        resolved["steps"] = steps
        return resolved


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
    with open(config_file, "r") as file:
        cfg = yaml.safe_load(file)
    verify_configuration_keys(cfg)  # Verify keys
    return cfg


def verify_configuration_keys(cfg: dict) -> None:
    required_keys = ["basic", "systems", "workflows", "globals"]
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


def setup_initial_magmom(atoms, magmom_data):
    """
    Set initial magnetic moments based on provided data.

    Args:
        atoms: ASE Atoms object
        magmom_data: dict (element-based) OR list/array (per-atom) OR None
            - dict: {"Ni": 2.0, "Cu": 0.0} applies by element symbol
            - list: [1.0, 2.0, 3.0, ...] applies per atom index
            - None: sets all to 0.0

    Returns:
        atoms: ASE Atoms object with magnetic moments set

    Raises:
        ValueError: If list is longer than number of atoms
        TypeError: If magmom_data is not dict, list, or None
    """

    # Handle list-based (per-atom) magmoms
    if isinstance(magmom_data, (list, tuple, np.ndarray)):
        magmom_list = list(magmom_data)

        # Validate: error if too long
        if len(magmom_list) > len(atoms):
            raise ValueError(
                f"Magnetic moment list is longer than number of atoms. "
                f"Got {len(magmom_list)} magmoms for {len(atoms)} atoms."
            )

        # Pad with zeros if shorter (per user requirement)
        if len(magmom_list) < len(atoms):
            logger.info(
                f"Padding magnetic moments: {len(magmom_list)} values "
                f"provided for {len(atoms)} atoms, remaining set to 0.0"
            )
            magmom_list.extend([0.0] * (len(atoms) - len(magmom_list)))

        # Apply per-atom magmoms using ASE bulk setter
        atoms.set_initial_magnetic_moments(magmom_list)
        logger.info(f"Applied per-atom magnetic moments (list with {len(magmom_list)} values)")

    # Handle dict-based (element-based) magmoms - EXISTING BEHAVIOR
    elif isinstance(magmom_data, dict):
        for at in atoms:
            if at.symbol in magmom_data:
                at.magmom = magmom_data[at.symbol]
            else:
                at.magmom = 0.0
        logger.info(f"Applied element-based magnetic moments from dict: {magmom_data}")

    # Handle None or empty containers
    elif magmom_data is None or (hasattr(magmom_data, "__len__") and len(magmom_data) == 0):
        for at in atoms:
            at.magmom = 0.0
        logger.info("No magnetic moments specified, set all to 0.0")

    else:
        raise TypeError(f"magmom_data must be dict, list, or None. Got {type(magmom_data)}")

    return atoms
