"""Tests for asetools.electrochemistry module."""

import pytest
import numpy as np


class TestAppliedPotentialImports:
    """Test that electrochemistry module imports work."""

    def test_import_module(self):
        from asetools.electrochemistry import appliedpotential
        assert appliedpotential is not None

    def test_ushe_constant(self):
        from asetools.electrochemistry.appliedpotential import U_SHE
        assert U_SHE == pytest.approx(4.43)

    def test_valences_dict(self):
        from asetools.electrochemistry.appliedpotential import VALANCES
        assert "Cu" in VALANCES
        assert "H" in VALANCES
        assert VALANCES["H"] == 1


class TestPureFunctions:
    """Test pure functions that don't need VASP files."""

    def test_custom_polynomial(self):
        from asetools.electrochemistry.appliedpotential import custom_polynomial
        # y = fixed_constant + beta[0]*x + beta[1]*x^2
        beta = [2.0, 3.0]
        x = np.array([0.0, 1.0, 2.0])
        result = custom_polynomial(beta, x, fixed_constant=1.0)
        expected = np.array([1.0, 6.0, 17.0])  # 1 + 2x + 3x^2
        np.testing.assert_array_almost_equal(result, expected)

    def test_custom_polynomial_zero(self):
        from asetools.electrochemistry.appliedpotential import custom_polynomial
        beta = [1.0]
        result = custom_polynomial(beta, 0.0, fixed_constant=5.0)
        assert result == pytest.approx(5.0)

    def test_get_sum_electrons(self, tmp_path):
        """Test electron counting with a simple structure."""
        from ase import Atoms
        from ase.io import write
        from asetools.electrochemistry.appliedpotential import get_sum_electrons

        # Create a simple Cu2 structure
        atoms = Atoms("Cu2", positions=[[0, 0, 0], [2.55, 0, 0]], cell=[10, 10, 10])
        poscar = str(tmp_path / "POSCAR")
        write(poscar, atoms, format="vasp")
        result = get_sum_electrons(poscar)
        assert result == 22  # 2 Cu atoms × 11 valence electrons
