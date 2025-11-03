"""
Constraint management module for ASEtools.

This module provides utilities for managing ASE constraints in VASP calculations,
including Hookean springs for maintaining bond distances and FixAtoms constraints.

Example usage:
    from asetools.constraints import ConstraintManager

    cm = ConstraintManager()
    cm.apply_from_json(atoms, 'proton_mappings.json', k=20.0)
"""

import json
import logging
from typing import List, Tuple, Optional, Dict, Any
from pathlib import Path

from ase import Atoms
from ase.constraints import FixAtoms, Hookean
from ase.data import covalent_radii

logger = logging.getLogger(__name__)


class ConstraintManager:
    """
    Manages ASE constraints with JSON configuration support.

    This class handles loading constraint configurations from JSON files,
    calculating appropriate constraint parameters, and applying constraints
    to ASE Atoms objects while preserving existing constraints.
    """

    def __init__(self, distance_factor: float = 1.134):
        """
        Initialize ConstraintManager.

        Parameters
        ----------
        distance_factor : float, optional
            Factor to convert sum of covalent radii to constraint distance.
            Based on O-H: r0=1.10 Å vs covalent_radii(O+H)=0.97 Å
            -> factor=1.134 (~13.4% buffer). Default is 1.134.
        """
        self.distance_factor = distance_factor

    def load_constraint_config(self, json_file: str) -> Dict[str, Any]:
        """
        Load constraint configuration from JSON file.

        Parameters
        ----------
        json_file : str
            Path to JSON configuration file

        Returns
        -------
        dict
            Constraint configuration dictionary

        Raises
        ------
        FileNotFoundError
            If JSON file doesn't exist
        json.JSONDecodeError
            If JSON file is malformed
        """
        json_path = Path(json_file)
        if not json_path.exists():
            raise FileNotFoundError(f"Constraint configuration file not found: {json_file}")

        try:
            with open(json_path, 'r') as f:
                config = json.load(f)
            logger.info(f"Loaded constraint configuration from {json_file}")
            return config
        except json.JSONDecodeError as e:
            raise json.JSONDecodeError(
                f"Failed to parse JSON file {json_file}: {e.msg}",
                e.doc, e.pos
            )

    def calculate_bond_distance(
        self,
        atoms: Atoms,
        idx1: int,
        idx2: int
    ) -> float:
        """
        Calculate equilibrium bond distance based on covalent radii.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object
        idx1 : int
            Index of first atom
        idx2 : int
            Index of second atom

        Returns
        -------
        float
            Equilibrium distance in Angstroms
        """
        radius1 = covalent_radii[atoms[idx1].number]
        radius2 = covalent_radii[atoms[idx2].number]
        r0 = (radius1 + radius2) * self.distance_factor
        return r0

    def apply_hookean_from_pairs(
        self,
        atoms: Atoms,
        pairs: List[Tuple[int, int]],
        k: float = 20.0,
        distance_factor: Optional[float] = None
    ) -> List[Hookean]:
        """
        Create Hookean spring constraints from atom index pairs.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object
        pairs : list of tuple
            List of (idx1, idx2) tuples specifying atom pairs to constrain
        k : float, optional
            Spring constant in eV/Ang^2. Default is 20.0.
        distance_factor : float, optional
            Override the instance distance_factor. If None, uses instance value.

        Returns
        -------
        list of Hookean
            List of Hookean constraint objects

        Raises
        ------
        ValueError
            If any atom index is out of range
        """
        if distance_factor is None:
            distance_factor = self.distance_factor

        hookean_constraints = []

        for idx1, idx2 in pairs:
            # Validate indices
            if idx1 < 0 or idx1 >= len(atoms):
                raise ValueError(f"Atom index {idx1} out of range [0, {len(atoms)-1}]")
            if idx2 < 0 or idx2 >= len(atoms):
                raise ValueError(f"Atom index {idx2} out of range [0, {len(atoms)-1}]")

            # Get atomic symbols and radii
            symbol1 = atoms[idx1].symbol
            symbol2 = atoms[idx2].symbol
            radius1 = covalent_radii[atoms[idx1].number]
            radius2 = covalent_radii[atoms[idx2].number]

            # Calculate equilibrium distance
            r0 = (radius1 + radius2) * distance_factor

            logger.info(
                f"Hookean constraint: {symbol1}({idx1})—{symbol2}({idx2}), "
                f"r0={r0:.3f} Å, k={k:.1f} eV/Å² "
                f"(radii: {radius1:.3f} + {radius2:.3f} = {radius1+radius2:.3f} Å)"
            )

            # Create Hookean constraint
            hookean_constraints.append(Hookean(a1=idx1, a2=idx2, k=k, rt=r0))

        return hookean_constraints

    def get_existing_fix_indices(self, atoms: Atoms) -> List[int]:
        """
        Extract existing FixAtoms constraint indices from atoms object.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object

        Returns
        -------
        list of int
            List of atom indices that are fixed
        """
        existing_constraints = atoms.constraints
        fix_indices = []

        if not existing_constraints:
            return fix_indices

        # Handle single constraint
        if isinstance(existing_constraints, FixAtoms):
            fix_indices = list(existing_constraints.index)
        # Handle list of constraints
        elif isinstance(existing_constraints, (list, tuple)):
            for constraint in existing_constraints:
                if isinstance(constraint, FixAtoms):
                    fix_indices.extend(constraint.index)

        return sorted(set(fix_indices))

    def merge_constraints(
        self,
        atoms: Atoms,
        new_constraints: List
    ) -> None:
        """
        Merge new constraints with existing constraints on atoms object.

        Preserves existing FixAtoms constraints and adds new constraints.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object (modified in place)
        new_constraints : list
            List of new constraint objects to add

        Raises
        ------
        RuntimeError
            If no constraints were applied
        """
        # Get existing fix indices
        fix_indices = self.get_existing_fix_indices(atoms)

        # Build combined constraint list
        combined_constraints = []

        # Add FixAtoms if there are any fixed atoms
        if fix_indices:
            combined_constraints.append(FixAtoms(indices=fix_indices))
            logger.info(f"Preserved {len(fix_indices)} FixAtoms constraints")

        # Add new constraints
        if not new_constraints:
            if not fix_indices:
                raise RuntimeError("No constraints to apply (neither existing nor new)")
            else:
                logger.warning("No new constraints provided, keeping only existing FixAtoms")
        else:
            combined_constraints.extend(new_constraints)
            logger.info(f"Added {len(new_constraints)} new constraints")

        # Apply combined constraints
        atoms.set_constraint(combined_constraints)
        logger.info(f"Total constraints applied: {len(combined_constraints)}")

    def apply_from_json(
        self,
        atoms: Atoms,
        json_file: str,
        k: float = 20.0,
        distance_factor: Optional[float] = None
    ) -> None:
        """
        Load and apply Hookean constraints from JSON configuration file.

        This is the main convenience method for applying constraints from
        a JSON file. It preserves existing FixAtoms constraints.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object (modified in place)
        json_file : str
            Path to JSON configuration file containing "pairs" key
        k : float, optional
            Spring constant in eV/Ang^2. Default is 20.0.
        distance_factor : float, optional
            Override the instance distance_factor. If None, uses instance value.

        Raises
        ------
        RuntimeError
            If JSON file has no "pairs" key or pairs list is empty
        """
        # Load configuration
        config = self.load_constraint_config(json_file)

        # Extract pairs
        if 'pairs' not in config:
            raise RuntimeError(
                f"JSON file {json_file} missing required 'pairs' key. "
                "Expected format: {{\"pairs\": [[idx1, idx2], ...]}}"
            )

        pairs = config['pairs']
        if not pairs:
            raise RuntimeError(f"Empty pairs list in {json_file}")

        logger.info(f"Applying {len(pairs)} Hookean constraints from {json_file}")

        # Override spring constant if specified in JSON metadata
        if 'metadata' in config and 'spring_constant' in config['metadata']:
            k = config['metadata']['spring_constant']
            logger.info(f"Using spring constant from JSON: k={k:.1f} eV/Å²")

        # Override distance factor if specified in JSON metadata
        if distance_factor is None:
            if 'metadata' in config and 'distance_factor' in config['metadata']:
                distance_factor = config['metadata']['distance_factor']
                logger.info(f"Using distance factor from JSON: {distance_factor:.3f}")
            else:
                distance_factor = self.distance_factor

        # Create Hookean constraints
        hookean_constraints = self.apply_hookean_from_pairs(
            atoms, pairs, k=k, distance_factor=distance_factor
        )

        # Merge with existing constraints
        self.merge_constraints(atoms, hookean_constraints)

        logger.info(f"Successfully applied constraints from {json_file}")

    def apply_stage_constraints(
        self,
        atoms: Atoms,
        constraint_config: Dict[str, Any]
    ) -> None:
        """
        Apply constraints from stage configuration (YAML workflow).

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object (modified in place)
        constraint_config : dict
            Constraint configuration from workflow stage

        Expected format:
            {
                'type': 'hookean',
                'config_file': 'path/to/pairs.json',
                'spring_constant': 20.0,
                'distance_factor': 1.134
            }
        """
        constraint_type = constraint_config.get('type', 'hookean')

        if constraint_type != 'hookean':
            raise ValueError(
                f"Unsupported constraint type: {constraint_type}. "
                "Currently only 'hookean' is supported."
            )

        config_file = constraint_config.get('config_file')
        if not config_file:
            raise ValueError("Missing 'config_file' in constraint configuration")

        k = constraint_config.get('spring_constant', 20.0)
        distance_factor = constraint_config.get('distance_factor', None)

        self.apply_from_json(atoms, config_file, k=k, distance_factor=distance_factor)
