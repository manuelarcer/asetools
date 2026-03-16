"""Tests for asetools.structure module."""

import pytest
from ase.build import fcc111


class TestSurfaceAnalyzer:
    """Tests for SurfaceAnalyzer class."""

    @pytest.fixture
    def slab_with_analyzer(self):
        """Create a slab and return it with a SurfaceAnalyzer."""
        from asetools.structure.adsorbate import SurfaceAnalyzer

        slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
        return slab, SurfaceAnalyzer(slab)

    def test_find_surface_atoms(self, slab_with_analyzer):
        slab, analyzer = slab_with_analyzer
        surface = analyzer.surface_indices
        assert len(surface) > 0
        # Surface atoms should be near the top z
        z_max = slab.positions[:, 2].max()
        for idx in surface:
            assert slab.positions[idx, 2] > z_max - 0.5

    def test_find_surface_neighbors(self, slab_with_analyzer):
        _, analyzer = slab_with_analyzer
        neighbors = analyzer.surface_neighbors
        # Should find triangular sites on fcc(111)
        assert len(neighbors) > 0
        for triplet in neighbors:
            assert len(triplet) == 3

    def test_midpoint_three_atoms(self, slab_with_analyzer):
        slab, analyzer = slab_with_analyzer
        if analyzer.surface_neighbors:
            triplet = analyzer.surface_neighbors[0]
            mid = analyzer.midpoint_three_atoms(triplet)
            assert mid.shape == (3,)
            # Midpoint z should be near surface level
            z_max = slab.positions[:, 2].max()
            assert mid[2] > z_max - 2.0

    def test_add_adsorbate_to_mid(self, slab_with_analyzer):
        slab, analyzer = slab_with_analyzer
        if analyzer.surface_neighbors:
            original_len = len(slab)
            triplet = analyzer.surface_neighbors[0]
            analyzer.add_adsorbate_to_mid(triplet, adsorbate="H", z_off=1.0)
            assert len(slab) == original_len + 1
            assert slab[-1].symbol == "H"


class TestBondValence:
    """Tests for bond valence module."""

    def test_import(self):
        from asetools.structure.bond_valence import BondValenceParameters, BondValenceSum

        assert BondValenceSum is not None
        assert BondValenceParameters is not None

    def test_bvparm_data_exists(self):
        """Ensure the BV parameter data file is accessible."""

        # The data file should be bundled
        from pathlib import Path

        data_file = Path(__file__).parent.parent / "asetools" / "data" / "bvparm2020.cif"
        assert data_file.exists(), "bvparm2020.cif not found in asetools/data/"
