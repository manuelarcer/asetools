"""
Tests for symmetry analysis functionality.

Uses ASE-generated structures for integration testing.
"""

import pytest
import numpy as np

# Check if spglib is available for test skipping
try:
    import spglib
    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False

pytestmark = pytest.mark.skipif(
    not HAS_SPGLIB,
    reason="spglib not installed"
)


class TestSymmetryAnalyzer:
    """Test SymmetryAnalyzer class."""

    def test_bulk_fcc_symmetry(self):
        """Test symmetry detection for FCC bulk structure."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        # Create FCC Cu bulk
        atoms = bulk('Cu', 'fcc', a=3.6)
        analyzer = SymmetryAnalyzer(atoms)

        # FCC should have high symmetry (Fm-3m, #225)
        assert analyzer.get_spacegroup_number() == 225
        assert 'Fm-3m' in analyzer.get_spacegroup()

        # All atoms in primitive cell are equivalent
        groups = analyzer.get_equivalent_groups()
        assert len(groups) == 1

    def test_bulk_supercell(self):
        """Test that supercell atoms are all equivalent."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('Cu', 'fcc', a=3.6) * (2, 2, 2)
        analyzer = SymmetryAnalyzer(atoms)

        # All Cu atoms should be equivalent
        groups = analyzer.get_equivalent_groups()
        assert len(groups) == 1
        assert len(groups[0]) == len(atoms)

    def test_rocksalt_structure(self):
        """Test NaCl rocksalt with two unique sites."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('NaCl', 'rocksalt', a=5.64)
        analyzer = SymmetryAnalyzer(atoms)

        # Rocksalt should be Fm-3m (#225)
        assert analyzer.get_spacegroup_number() == 225

        # Two unique sites (Na and Cl)
        unique = analyzer.get_unique_sites()
        assert len(unique) == 2

        # Filter by element
        na_sites = analyzer.get_unique_sites(element='Na')
        cl_sites = analyzer.get_unique_sites(element='Cl')
        assert len(na_sites) == 1
        assert len(cl_sites) == 1

    def test_are_equivalent_method(self):
        """Test are_equivalent() method."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('Cu', 'fcc', a=3.6) * (2, 2, 2)
        analyzer = SymmetryAnalyzer(atoms)

        # All Cu atoms should be equivalent to each other
        for i in range(len(atoms)):
            for j in range(len(atoms)):
                assert analyzer.are_equivalent(i, j) == True

        # Test index bounds checking
        with pytest.raises(IndexError):
            analyzer.are_equivalent(0, len(atoms) + 1)

    def test_are_equivalent_different_elements(self):
        """Test that different elements are not equivalent."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('NaCl', 'rocksalt', a=5.64)
        analyzer = SymmetryAnalyzer(atoms)

        # Find Na and Cl indices
        symbols = atoms.get_chemical_symbols()
        na_idx = symbols.index('Na')
        cl_idx = symbols.index('Cl')

        # Na and Cl should not be equivalent
        assert analyzer.are_equivalent(na_idx, cl_idx) == False

    def test_get_multiplicity(self):
        """Test get_multiplicity() method."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('Cu', 'fcc', a=3.6)
        analyzer = SymmetryAnalyzer(atoms)

        # In primitive cell, only one atom
        mult = analyzer.get_multiplicity(0)
        assert mult == 1

        # In supercell, all equivalent
        atoms_super = atoms * (2, 2, 2)
        analyzer_super = SymmetryAnalyzer(atoms_super)
        mult_super = analyzer_super.get_multiplicity(0)
        assert mult_super == len(atoms_super)

    def test_get_equivalent_groups_element_filter(self):
        """Test element filtering in get_equivalent_groups."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('NaCl', 'rocksalt', a=5.64) * (2, 2, 2)
        analyzer = SymmetryAnalyzer(atoms)

        # Get only Na groups
        na_groups = analyzer.get_equivalent_groups(element='Na')
        for members in na_groups.values():
            for idx in members:
                assert atoms[idx].symbol == 'Na'

        # Get only Cl groups
        cl_groups = analyzer.get_equivalent_groups(element='Cl')
        for members in cl_groups.values():
            for idx in members:
                assert atoms[idx].symbol == 'Cl'


class TestSurfaceSymmetry:
    """Test surface-specific symmetry analysis."""

    def test_surface_indices(self):
        """Test surface atom detection."""
        from ase.build import fcc111
        from asetools.analysis import SymmetryAnalyzer

        # Create Pt(111) slab
        slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
        analyzer = SymmetryAnalyzer(slab, is_surface=True, surface_threshold=2.5)

        surface_indices = analyzer.get_surface_indices()

        # Should find top layer atoms
        assert len(surface_indices) > 0

        # All surface atoms should be near the top
        positions = slab.get_positions()
        z_max = positions[:, 2].max()
        for idx in surface_indices:
            assert positions[idx, 2] > z_max - 2.5

    def test_surface_equivalent_groups(self):
        """Test equivalence grouping for surface atoms."""
        from ase.build import fcc111
        from asetools.analysis import SymmetryAnalyzer

        slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
        analyzer = SymmetryAnalyzer(slab, is_surface=True)

        groups = analyzer.get_surface_equivalent_groups()

        # Should have grouped surface atoms
        assert len(groups) >= 1

        # All indices should be surface atoms
        surface_indices = set(analyzer.get_surface_indices())
        for members in groups.values():
            for idx in members:
                assert idx in surface_indices

    def test_unique_adsorption_sites(self):
        """Test identification of unique adsorption sites."""
        from ase.build import fcc111
        from asetools.analysis import SymmetryAnalyzer

        slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
        analyzer = SymmetryAnalyzer(slab, is_surface=True)

        unique = analyzer.get_unique_adsorption_sites()

        # Should return representative sites
        assert len(unique) >= 1

        # Each unique site should be a surface atom
        surface_indices = analyzer.get_surface_indices()
        for idx in unique:
            assert idx in surface_indices


class TestToleranceSettings:
    """Test symmetry precision and tolerance settings."""

    def test_symprec_effect(self):
        """Test effect of symmetry precision on detection."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('Cu', 'fcc', a=3.6)

        # Add small perturbation
        atoms.positions[0] += [0.001, 0.001, 0.001]

        # With tight precision, may lose symmetry
        analyzer_tight = SymmetryAnalyzer(atoms, symprec=1e-6)

        # With loose precision, should still find symmetry
        analyzer_loose = SymmetryAnalyzer(atoms, symprec=0.01)

        # Loose should have >= symmetry operations as tight
        ops_tight = len(analyzer_tight.get_symmetry_operations())
        ops_loose = len(analyzer_loose.get_symmetry_operations())
        assert ops_loose >= ops_tight


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_single_atom(self):
        """Test with single atom (trivial symmetry)."""
        from ase import Atoms
        from asetools.analysis import SymmetryAnalyzer

        atoms = Atoms('H', positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        analyzer = SymmetryAnalyzer(atoms)

        # Single atom has trivial equivalence
        groups = analyzer.get_equivalent_groups()
        assert len(groups) == 1
        assert groups[0] == [0]

    def test_molecule_with_pbc(self):
        """Test with molecule in periodic box."""
        from ase.build import molecule
        from asetools.analysis import SymmetryAnalyzer

        # Water molecule
        h2o = molecule('H2O')
        h2o.set_cell([10, 10, 10])
        h2o.pbc = True

        analyzer = SymmetryAnalyzer(h2o)

        # Should find some symmetry
        groups = analyzer.get_equivalent_groups()
        assert len(groups) >= 1

        # Two H atoms should be equivalent
        h_groups = analyzer.get_equivalent_groups(element='H')
        assert len(h_groups) == 1  # Both H atoms equivalent
        assert len(list(h_groups.values())[0]) == 2

    def test_analyze_method(self):
        """Test full analysis method."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('Cu', 'fcc', a=3.6)
        analyzer = SymmetryAnalyzer(atoms)

        results = analyzer.analyze()

        assert 'spacegroup' in results
        assert 'spacegroup_number' in results
        assert 'n_operations' in results
        assert 'equivalent_groups' in results
        assert 'unique_sites' in results

    def test_analyze_surface_mode(self):
        """Test analysis in surface mode."""
        from ase.build import fcc111
        from asetools.analysis import SymmetryAnalyzer

        slab = fcc111('Pt', size=(2, 2, 3), vacuum=10.0)
        analyzer = SymmetryAnalyzer(slab, is_surface=True)

        results = analyzer.analyze()

        # Should have surface-specific results
        assert 'surface_indices' in results
        assert 'surface_equivalent_groups' in results
        assert 'unique_adsorption_sites' in results

    def test_print_summary_no_error(self):
        """Test that print_summary runs without error."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('Cu', 'fcc', a=3.6)
        analyzer = SymmetryAnalyzer(atoms)

        # Should not raise
        analyzer.print_summary()

    def test_symmetry_operations(self):
        """Test get_symmetry_operations returns valid data."""
        from ase.build import bulk
        from asetools.analysis import SymmetryAnalyzer

        atoms = bulk('Cu', 'fcc', a=3.6)
        analyzer = SymmetryAnalyzer(atoms)

        operations = analyzer.get_symmetry_operations()

        # FCC has 48 symmetry operations
        assert len(operations) == 48

        # Each operation is a (rotation, translation) tuple
        for rot, trans in operations:
            assert rot.shape == (3, 3)
            assert trans.shape == (3,)


class TestSymmetryAvailable:
    """Test symmetry_available function."""

    def test_symmetry_available(self):
        """Test that symmetry_available returns True when spglib is installed."""
        from asetools.analysis import symmetry_available

        # Since we're running these tests, spglib must be available
        assert symmetry_available() is True
