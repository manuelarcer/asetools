"""Extended tests for asetools.electronic.doscar module.

Covers: DOS properties, PDOS extraction, legacy functions, to_dataframe,
error paths, and individual atom PDOS.
"""

from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).parent / "data"
DOSCAR_FILE = DATA_DIR / "DOSCAR_fe3o4_Feoct_2x2.doscar"


@pytest.fixture
def dos():
    """Load the test DOSCAR file."""
    if not DOSCAR_FILE.exists():
        pytest.skip("DOSCAR test data not found")
    from asetools.electronic.doscar import DOS
    return DOS(str(DOSCAR_FILE))


# ── Basic properties ──────────────────────────────────────────────────


class TestDOSProperties:
    def test_energy_is_array(self, dos):
        assert isinstance(dos.energy, np.ndarray)
        assert len(dos.energy) > 0

    def test_dos_up_shape(self, dos):
        assert dos.dos_up.shape == dos.energy.shape

    def test_dos_down_shape(self, dos):
        assert dos.dos_down.shape == dos.energy.shape

    def test_dos_down_negative(self, dos):
        # dos_down is stored as negative values
        assert np.all(dos.dos_down <= 0)

    def test_total_dos_nonnegative(self, dos):
        assert np.all(dos.total_dos >= 0)

    def test_total_dos_equals_sum(self, dos):
        expected = dos.dos_up + np.abs(dos.dos_down)
        np.testing.assert_array_equal(dos.total_dos, expected)

    def test_fermi_energy_is_float(self, dos):
        assert isinstance(dos.fermi_energy, float)

    def test_has_partial_dos(self, dos):
        assert isinstance(dos.has_partial_dos, bool)

    def test_natoms_positive(self, dos):
        if dos.has_partial_dos:
            assert dos.natoms > 0


# ── to_dataframe ──────────────────────────────────────────────────────


class TestToDataframe:
    def test_returns_dataframe(self, dos):
        import pandas as pd
        df = dos.to_dataframe()
        assert isinstance(df, pd.DataFrame)

    def test_dataframe_columns(self, dos):
        df = dos.to_dataframe()
        assert "energy" in df.columns
        assert "dos_up" in df.columns
        assert "dos_down" in df.columns
        assert "total_dos" in df.columns

    def test_dataframe_length(self, dos):
        df = dos.to_dataframe()
        assert len(df) == len(dos.energy)


# ── PDOS extraction ──────────────────────────────────────────────────


class TestPDOSByStates:
    def test_s_states(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_states([0], ["s_states"])
        assert len(energy) == len(dos.energy)
        assert len(up) == len(dos.energy)
        assert len(down) == len(dos.energy)

    def test_p_states(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_states([0], ["p_states"])
        np.testing.assert_array_equal(energy, dos.energy)

    def test_d_states(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_states([0], ["d_states"])
        assert up.shape == dos.energy.shape

    def test_multiple_states(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_states([0], ["s_states", "p_states"])
        # Should be sum of both
        e_s, up_s, _ = dos.get_pdos_by_states([0], ["s_states"])
        e_p, up_p, _ = dos.get_pdos_by_states([0], ["p_states"])
        np.testing.assert_allclose(up, up_s + up_p)

    def test_no_partial_dos_raises(self):
        """Test error when no partial DOS data."""
        from asetools.electronic.doscar import DOS
        # Create a mock DOS without partial DOS
        dos_obj = DOS.__new__(DOS)
        dos_obj.has_partial_dos = False
        dos_obj.data = {"energy": np.array([1, 2, 3])}
        with pytest.raises(ValueError, match="No partial DOS"):
            dos_obj.get_pdos_by_states([0], ["s_states"])


class TestPDOSByOrbitals:
    def test_all_d_shorthand(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_orbitals([0], "all-d")
        assert len(up) == len(dos.energy)

    def test_all_p_shorthand(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_orbitals([0], "all-p")
        assert len(up) == len(dos.energy)

    def test_all_s_shorthand(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_orbitals([0], "all-s")
        assert len(up) == len(dos.energy)

    def test_t2g_shorthand(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_orbitals([0], "t2g")
        assert len(up) == len(dos.energy)

    def test_eg_shorthand(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_orbitals([0], "eg")
        assert len(up) == len(dos.energy)

    def test_all_shorthand(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        energy, up, down = dos.get_pdos_by_orbitals([0], "all")
        assert len(up) == len(dos.energy)

    def test_t2g_plus_eg_equals_all_d(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        _, up_d, down_d = dos.get_pdos_by_orbitals([0], "all-d")
        _, up_t2g, down_t2g = dos.get_pdos_by_orbitals([0], "t2g")
        _, up_eg, down_eg = dos.get_pdos_by_orbitals([0], "eg")
        np.testing.assert_allclose(up_d, up_t2g + up_eg)
        np.testing.assert_allclose(down_d, down_t2g + down_eg)


# ── Individual atom PDOS ─────────────────────────────────────────────


class TestIndividualAtomPDOS:
    def test_with_orbitals(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        result = dos.get_individual_atom_pdos([0, 1], orbitals="all-d")
        assert 0 in result
        assert 1 in result
        assert "energy" in result[0]
        assert "dos_up" in result[0]
        assert "dos_down" in result[0]

    def test_with_states(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        result = dos.get_individual_atom_pdos([0], states=["d_states"])
        assert 0 in result

    def test_both_raises(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        with pytest.raises(ValueError, match="Cannot specify both"):
            dos.get_individual_atom_pdos([0], states=["d_states"], orbitals="all-d")

    def test_neither_raises(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        with pytest.raises(ValueError, match="Must specify either"):
            dos.get_individual_atom_pdos([0])


# ── Band center edge cases ──────────────────────────────────────────


class TestBandCenterEdgeCases:
    def test_separate_spin(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        result = dos.calculate_band_center([0], orbitals="all-d", spin_treatment="separate")
        assert isinstance(result, dict)
        assert "up" in result
        assert "down" in result
        assert isinstance(result["up"], float)
        assert isinstance(result["down"], float)

    def test_spin_up_only(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        result = dos.calculate_band_center([0], orbitals="all-d", spin_treatment="up")
        assert isinstance(result, float)

    def test_spin_down_only(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        result = dos.calculate_band_center([0], orbitals="all-d", spin_treatment="down")
        assert isinstance(result, float)

    def test_energy_range(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        full = dos.calculate_band_center([0], orbitals="all-d")
        restricted = dos.calculate_band_center([0], orbitals="all-d", energy_range=(-5, 0))
        # Restricted range should give different result (unless all DOS is in that range)
        assert isinstance(restricted, float)

    def test_invalid_spin_treatment(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        with pytest.raises(ValueError, match="spin_treatment"):
            dos.calculate_band_center([0], orbitals="all-d", spin_treatment="invalid")

    def test_both_states_and_orbitals_raises(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        with pytest.raises(ValueError, match="Cannot specify both"):
            dos.calculate_band_center([0], orbitals="all-d", states=["d_states"])

    def test_neither_states_nor_orbitals_raises(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        with pytest.raises(ValueError, match="Must specify either"):
            dos.calculate_band_center([0])

    def test_with_states(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        result = dos.calculate_band_center([0], states=["d_states"])
        assert isinstance(result, float)


# ── Legacy functions ─────────────────────────────────────────────────


class TestLegacyFunctions:
    def test_extract_dos(self, dos):
        from asetools.electronic.doscar import extract_dos
        data = extract_dos(str(DOSCAR_FILE))
        assert isinstance(data, dict)
        assert "energy" in data
        assert "DOSup" in data
        assert "DOSdown" in data

    def test_extract_fermi_e(self):
        from asetools.electronic.doscar import extract_fermi_e
        fe = extract_fermi_e(str(DOSCAR_FILE))
        assert isinstance(fe, float)

    def test_extract_fermi_e_missing_file(self):
        from asetools.electronic.doscar import extract_fermi_e
        result = extract_fermi_e("/nonexistent/DOSCAR")
        assert result is None

    def test_extract_pdos_perstate(self, dos):
        from asetools.electronic.doscar import extract_dos, extract_pdos_perstate
        data = extract_dos(str(DOSCAR_FILE))
        if "at-0" not in data:
            pytest.skip("No partial DOS in test data")
        energy, up, down = extract_pdos_perstate(data, [0], ["d_states"])
        assert len(energy) > 0
        assert len(up) == len(energy)

    def test_extract_pdos_perorbital(self, dos):
        from asetools.electronic.doscar import extract_dos, extract_pdos_perorbital
        data = extract_dos(str(DOSCAR_FILE))
        if "at-0" not in data:
            pytest.skip("No partial DOS in test data")
        energy, up, down = extract_pdos_perorbital(data, [0], "all-d")
        assert len(energy) > 0

    def test_extract_pdos_perorbital_shorthands(self, dos):
        from asetools.electronic.doscar import extract_dos, extract_pdos_perorbital
        data = extract_dos(str(DOSCAR_FILE))
        if "at-0" not in data:
            pytest.skip("No partial DOS in test data")
        for shorthand in ["all-s", "all-p", "all-d", "all", "t2g", "eg"]:
            energy, up, down = extract_pdos_perorbital(data, [0], shorthand)
            assert len(energy) > 0, f"Failed for {shorthand}"

    def test_calculate_band_center_legacy(self, dos):
        from asetools.electronic.doscar import calculate_band_center
        result = calculate_band_center(str(DOSCAR_FILE), [0], orbitals="all-d")
        assert isinstance(result, float)


# ── _calculate_moment edge cases ─────────────────────────────────────


class TestCalculateMoment:
    def test_empty_dos_raises(self):
        from asetools.electronic.doscar import DOS
        dos_obj = DOS.__new__(DOS)
        with pytest.raises(ValueError, match="empty or all zeros"):
            dos_obj._calculate_moment(np.array([]), np.array([]))

    def test_all_zeros_raises(self):
        from asetools.electronic.doscar import DOS
        dos_obj = DOS.__new__(DOS)
        with pytest.raises(ValueError, match="empty or all zeros"):
            dos_obj._calculate_moment(np.array([1, 2, 3]), np.array([0, 0, 0]))

    def test_delta_function(self):
        """Band center of a delta-like peak should be at that energy."""
        from asetools.electronic.doscar import DOS
        dos_obj = DOS.__new__(DOS)
        energy = np.linspace(-10, 10, 1000)
        dos_vals = np.exp(-((energy - 2.0) ** 2) / 0.01)  # Sharp peak at 2.0
        result = dos_obj._calculate_moment(energy, dos_vals)
        assert abs(result - 2.0) < 0.1


# ── Plotting (just that they don't crash, not visual checks) ─────────


class TestPlotting:
    def test_plot_total_dos(self, dos):
        import matplotlib
        matplotlib.use("Agg")
        ax = dos.plot_total_dos()
        assert ax is not None

    def test_plot_pdos_by_states(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        import matplotlib
        matplotlib.use("Agg")
        ax = dos.plot_pdos_by_states([0], ["d_states"])
        assert ax is not None

    def test_plot_pdos_by_orbitals(self, dos):
        if not dos.has_partial_dos:
            pytest.skip("No partial DOS")
        import matplotlib
        matplotlib.use("Agg")
        ax = dos.plot_pdos_by_orbitals([0], "all-d")
        assert ax is not None
