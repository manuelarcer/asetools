"""Tests for asetools.electronic module."""

import pytest
import os
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


class TestDOSImport:
    """Test that DOS class can be imported and instantiated."""

    def test_import_dos_class(self):
        from asetools.electronic.doscar import DOS

        assert DOS is not None

    def test_import_legacy_functions(self):
        from asetools.electronic.doscar import (
            calculate_band_center,
            extract_fermi_e,
        )

        assert callable(calculate_band_center)
        assert callable(extract_fermi_e)


class TestDOSWithData:
    """Tests that require DOSCAR test data files."""

    @pytest.fixture
    def doscar_path(self):
        p = DATA_DIR / "DOSCAR_fe3o4_Feoct_2x2.doscar"
        if not p.exists():
            pytest.skip(f"Test data not found: {p}")
        return str(p)

    @pytest.fixture
    def dos(self, doscar_path):
        from asetools.electronic.doscar import DOS

        return DOS(doscar_path)

    def test_dos_init(self, dos):
        assert dos is not None
        assert hasattr(dos, "energy")
        assert hasattr(dos, "total_dos")

    def test_fermi_energy(self, dos):
        fe = dos.fermi_energy
        # Fermi energy should be a finite number
        assert isinstance(fe, float)
        import math

        assert math.isfinite(fe)

    def test_extract_fermi_e(self, doscar_path):
        from asetools.electronic.doscar import extract_fermi_e

        fe = extract_fermi_e(doscar_path)
        assert isinstance(fe, float)

    def test_band_center_d_orbitals(self, dos):
        result = dos.calculate_band_center([0], orbitals="all-d")
        assert isinstance(result, float)

    def test_band_center_p_orbitals(self, dos):
        result = dos.calculate_band_center([0], orbitals="all-p")
        assert isinstance(result, float)

    def test_band_center_spin_treatments(self, dos):
        combined = dos.calculate_band_center([0], orbitals="all-d", spin_treatment="combined")
        assert isinstance(combined, float)

    def test_band_center_multiple_atoms(self, dos):
        single = dos.calculate_band_center([0], orbitals="all-d")
        # Should not error with multiple atoms
        assert isinstance(single, float)
