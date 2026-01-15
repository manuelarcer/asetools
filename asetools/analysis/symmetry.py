"""
Symmetry analysis for atomic structures using spglib.

This module provides functionality for analyzing symmetry in bulk and surface structures,
with a focus on identifying equivalent atomic sites for adsorption studies.

The spglib package is an optional dependency and will be imported only when
symmetry analysis features are used.

Classes:
    SymmetryAnalyzer: Main class for symmetry analysis of ASE Atoms objects
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional

from ase import Atoms

logger = logging.getLogger(__name__)


def _get_spglib():
    """Import and return spglib, raising helpful error if not installed."""
    try:
        import spglib
        return spglib
    except ImportError:
        raise ImportError(
            "spglib is required for symmetry analysis but is not installed. "
            "Install it with: pip install spglib"
        )


class SymmetryAnalyzer:
    """
    Analyze symmetry of atomic structures for identifying equivalent sites.

    This class provides functionality to:
    - Determine space group and point group symmetry
    - Identify symmetry operations
    - Group atoms by symmetry equivalence
    - Check if two specific atoms are symmetrically equivalent
    - Handle surface slabs with broken symmetry

    Uses spglib for symmetry detection (optional dependency).

    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object to analyze
    symprec : float, optional
        Symmetry precision for spglib (default: 1e-5)
    angle_tolerance : float, optional
        Angle tolerance in degrees for symmetry detection (default: -1.0,
        meaning automatic)
    is_surface : bool, optional
        If True, analyze only the surface layer for symmetry (default: False)
    surface_threshold : float, optional
        Distance threshold in Angstrom for surface atom detection when
        is_surface=True (default: 2.0)

    Examples
    --------
    >>> from ase.build import bulk
    >>> atoms = bulk('Cu', 'fcc', a=3.6)
    >>> analyzer = SymmetryAnalyzer(atoms)
    >>> print(analyzer.get_spacegroup())
    'Fm-3m (225)'

    >>> groups = analyzer.get_equivalent_groups()
    >>> # All Cu atoms in bulk FCC are equivalent

    >>> # Surface slab example
    >>> from ase.build import fcc111
    >>> slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
    >>> analyzer = SymmetryAnalyzer(slab, is_surface=True)
    >>> groups = analyzer.get_equivalent_groups()
    >>> # Surface atoms grouped by equivalence
    """

    def __init__(
        self,
        atoms: Atoms,
        symprec: float = 1e-5,
        angle_tolerance: float = -1.0,
        is_surface: bool = False,
        surface_threshold: float = 2.0
    ):
        """Initialize SymmetryAnalyzer."""
        self.spglib = _get_spglib()

        self.atoms = atoms
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.is_surface = is_surface
        self.surface_threshold = surface_threshold

        # Cached results
        self._symmetry_dataset = None
        self._surface_indices = None

    def _get_spglib_cell(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Convert ASE Atoms to spglib cell format."""
        lattice = np.array(self.atoms.get_cell())
        positions = self.atoms.get_scaled_positions()
        numbers = self.atoms.get_atomic_numbers()
        return (lattice, positions, numbers)

    def _get_symmetry_dataset(self) -> dict:
        """Get or compute symmetry dataset from spglib.

        Returns a dict-like object with symmetry information. Uses attribute
        access for modern spglib versions to avoid deprecation warnings.
        """
        if self._symmetry_dataset is None:
            cell = self._get_spglib_cell()
            dataset = self.spglib.get_symmetry_dataset(
                cell,
                symprec=self.symprec,
                angle_tolerance=self.angle_tolerance
            )
            if dataset is None:
                logger.warning("spglib could not determine symmetry, using P1")
                self._symmetry_dataset = {
                    'number': 1,
                    'international': 'P1',
                    'equivalent_atoms': np.arange(len(self.atoms)),
                    'rotations': np.array([np.eye(3, dtype=int)]),
                    'translations': np.array([[0, 0, 0]])
                }
            else:
                # Use attribute access for modern spglib, fall back to dict for older versions
                self._symmetry_dataset = {
                    'number': getattr(dataset, 'number', None) or dataset['number'],
                    'international': getattr(dataset, 'international', None) or dataset['international'],
                    'equivalent_atoms': getattr(dataset, 'equivalent_atoms', None)
                                        if hasattr(dataset, 'equivalent_atoms') else dataset['equivalent_atoms'],
                    'rotations': getattr(dataset, 'rotations', None)
                                 if hasattr(dataset, 'rotations') else dataset['rotations'],
                    'translations': getattr(dataset, 'translations', None)
                                    if hasattr(dataset, 'translations') else dataset['translations'],
                }
        return self._symmetry_dataset

    def get_spacegroup(self) -> str:
        """
        Get space group symbol and number.

        Returns
        -------
        str
            Space group symbol with number, e.g., 'Fm-3m (225)'
        """
        dataset = self._get_symmetry_dataset()
        return f"{dataset['international']} ({dataset['number']})"

    def get_spacegroup_number(self) -> int:
        """
        Get space group number.

        Returns
        -------
        int
            Space group number (1-230)
        """
        dataset = self._get_symmetry_dataset()
        return dataset['number']

    def get_equivalent_atoms(self) -> np.ndarray:
        """
        Get array mapping each atom to its equivalence class.

        Returns
        -------
        np.ndarray
            Array where equivalent_atoms[i] gives the representative atom index
            for atom i. Atoms with the same value are symmetrically equivalent.
        """
        dataset = self._get_symmetry_dataset()
        return dataset['equivalent_atoms']

    def get_symmetry_operations(self) -> List[Tuple[np.ndarray, np.ndarray]]:
        """
        Get list of symmetry operations as (rotation, translation) pairs.

        Returns
        -------
        list of tuples
            Each tuple contains (rotation_matrix, translation_vector)
            where rotation is a 3x3 integer matrix and translation is
            a length-3 fractional coordinate vector.
        """
        dataset = self._get_symmetry_dataset()
        rotations = dataset['rotations']
        translations = dataset['translations']
        return [(rotations[i], translations[i]) for i in range(len(rotations))]

    def are_equivalent(self, index1: int, index2: int) -> bool:
        """
        Check if two atoms are symmetrically equivalent.

        Parameters
        ----------
        index1 : int
            Index of first atom
        index2 : int
            Index of second atom

        Returns
        -------
        bool
            True if atoms are symmetrically equivalent

        Raises
        ------
        IndexError
            If atom indices are out of range
        """
        n_atoms = len(self.atoms)
        if not (0 <= index1 < n_atoms and 0 <= index2 < n_atoms):
            raise IndexError(
                f"Atom indices must be in range [0, {n_atoms}), "
                f"got {index1} and {index2}"
            )

        equiv = self.get_equivalent_atoms()
        return equiv[index1] == equiv[index2]

    def get_equivalent_groups(
        self,
        element: Optional[str] = None,
        indices: Optional[List[int]] = None
    ) -> Dict[int, List[int]]:
        """
        Group atoms by symmetry equivalence.

        Parameters
        ----------
        element : str, optional
            Filter by element symbol (e.g., 'O' for oxygen atoms only)
        indices : list of int, optional
            Consider only specific atom indices (useful for surface atoms)

        Returns
        -------
        dict
            Dictionary mapping representative atom index to list of
            equivalent atom indices. Keys are the smallest index in each group.

        Examples
        --------
        >>> analyzer = SymmetryAnalyzer(atoms)
        >>> groups = analyzer.get_equivalent_groups(element='O')
        >>> # {0: [0, 3, 5], 1: [1, 2, 4]} means O atoms 0,3,5 are equivalent
        """
        equiv = self.get_equivalent_atoms()
        symbols = self.atoms.get_chemical_symbols()

        # Build groups
        groups: Dict[int, List[int]] = {}

        for i in range(len(self.atoms)):
            # Filter by element if specified
            if element is not None and symbols[i] != element:
                continue

            # Filter by indices if specified
            if indices is not None and i not in indices:
                continue

            rep = equiv[i]
            if rep not in groups:
                groups[rep] = []
            groups[rep].append(i)

        # Renumber groups to use smallest index as key
        renumbered = {}
        for rep, members in groups.items():
            min_idx = min(members)
            renumbered[min_idx] = sorted(members)

        return dict(sorted(renumbered.items()))

    def get_unique_sites(
        self,
        element: Optional[str] = None
    ) -> List[int]:
        """
        Get one representative atom from each equivalence class.

        Parameters
        ----------
        element : str, optional
            Filter by element symbol

        Returns
        -------
        list of int
            Atom indices, one per unique site
        """
        groups = self.get_equivalent_groups(element=element)
        return sorted(groups.keys())

    def get_multiplicity(self, index: int) -> int:
        """
        Get the multiplicity (number of equivalent sites) for an atom.

        Parameters
        ----------
        index : int
            Atom index

        Returns
        -------
        int
            Number of atoms equivalent to this one (including itself)
        """
        equiv = self.get_equivalent_atoms()
        return int(np.sum(equiv == equiv[index]))

    # === Surface-specific methods ===

    def get_surface_indices(self) -> List[int]:
        """
        Identify surface atom indices based on z-coordinate threshold.

        Returns
        -------
        list of int
            Indices of atoms in the surface layer
        """
        if self._surface_indices is None:
            positions = self.atoms.get_positions()
            z_max = positions[:, 2].max()
            self._surface_indices = [
                i for i, pos in enumerate(positions)
                if pos[2] > z_max - self.surface_threshold
            ]
        return self._surface_indices

    def get_surface_equivalent_groups(
        self,
        element: Optional[str] = None
    ) -> Dict[int, List[int]]:
        """
        Get equivalence groups for surface atoms only.

        This is particularly useful for adsorption studies where only
        surface sites are relevant.

        Parameters
        ----------
        element : str, optional
            Filter by element symbol

        Returns
        -------
        dict
            Dictionary mapping representative surface atom index to list of
            equivalent surface atom indices
        """
        surface_indices = self.get_surface_indices()
        return self.get_equivalent_groups(element=element, indices=surface_indices)

    def get_unique_adsorption_sites(
        self,
        element: Optional[str] = None
    ) -> List[int]:
        """
        Get unique surface sites for adsorption studies.

        Returns one representative atom from each equivalence class
        among surface atoms.

        Parameters
        ----------
        element : str, optional
            Filter by element symbol (e.g., 'O' for oxygen sites)

        Returns
        -------
        list of int
            Indices of unique adsorption sites
        """
        groups = self.get_surface_equivalent_groups(element=element)
        return sorted(groups.keys())

    # === Analysis methods ===

    def analyze(self) -> Dict:
        """
        Perform full symmetry analysis and return summary.

        Returns
        -------
        dict
            Analysis results including:
            - spacegroup: Space group symbol and number
            - spacegroup_number: Space group number
            - n_operations: Number of symmetry operations
            - n_equivalent_groups: Number of equivalence classes
            - equivalent_groups: Dict of equivalent atom groups
            - unique_sites: List of unique site indices
        """
        results = {
            'spacegroup': self.get_spacegroup(),
            'spacegroup_number': self.get_spacegroup_number(),
            'n_operations': len(self.get_symmetry_operations()),
            'n_equivalent_groups': len(self.get_equivalent_groups()),
            'equivalent_groups': self.get_equivalent_groups(),
            'unique_sites': self.get_unique_sites(),
        }

        if self.is_surface:
            results.update({
                'surface_indices': self.get_surface_indices(),
                'surface_equivalent_groups': self.get_surface_equivalent_groups(),
                'unique_adsorption_sites': self.get_unique_adsorption_sites(),
            })

        return results

    def print_summary(self):
        """Print a human-readable summary of symmetry analysis."""
        analysis = self.analyze()

        print("Symmetry Analysis Summary")
        print("=" * 40)
        print(f"Space group: {analysis['spacegroup']}")
        print(f"Number of symmetry operations: {analysis['n_operations']}")
        print(f"Number of equivalence classes: {analysis['n_equivalent_groups']}")
        print()

        # Group by element
        symbols = self.atoms.get_chemical_symbols()
        elements = sorted(set(symbols))

        for element in elements:
            groups = self.get_equivalent_groups(element=element)
            print(f"{element} atoms:")
            for rep, members in groups.items():
                multiplicity = len(members)
                print(f"  Site {rep}: {multiplicity} equivalent atoms {members}")
            print()

        if self.is_surface:
            print("Surface Analysis:")
            print(f"  Surface atoms: {len(analysis['surface_indices'])} atoms")
            print(f"  Unique adsorption sites: {analysis['unique_adsorption_sites']}")
