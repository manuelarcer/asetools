"""Tests for asetools.parsers.vasp_outcar module.

Tests utility functions, parser base classes, and concrete header/chunk parsers
using synthetic data (no real OUTCAR files needed for most tests).
"""

import numpy as np
import pytest

from asetools.parsers.vasp_outcar import (
    NoNonEmptyLines,
    UnableToLocateDelimiter,
    _check_line,
    convert_vasp_outcar_stress,
    find_next_non_empty_line,
    search_lines,
)


# ── Utility functions ──────────────────────────────────────────────────


class TestCheckLine:
    """Test _check_line for VASP numeric formatting quirks."""

    def test_normal_line_unchanged(self):
        line = "  1.234  5.678  9.012"
        assert _check_line(line) == line

    def test_missing_space_between_numbers(self):
        # VASP sometimes omits space: "1.234-5.678" should become "1.234 -5.678"
        line = "1.234-5.678"
        result = _check_line(line)
        assert " -" in result
        assert result == "1.234 -5.678"

    def test_multiple_missing_spaces(self):
        line = "1.234-5.678-9.012"
        result = _check_line(line)
        assert result == "1.234 -5.678 -9.012"

    def test_no_digit_before_minus(self):
        # Negative number at start shouldn't be affected
        line = " -1.234  5.678"
        assert _check_line(line) == line

    def test_actual_negative_number(self):
        # Pure negative with space before: no digit-minus-digit pattern
        line = "  -1.234  -5.678"
        assert _check_line(line) == line


class TestFindNextNonEmptyLine:
    """Test find_next_non_empty_line."""

    def test_cursor_on_non_empty(self):
        lines = ["hello", "world"]
        assert find_next_non_empty_line(0, lines) == 0

    def test_skip_empty_lines(self):
        lines = ["", "  ", "hello"]
        assert find_next_non_empty_line(0, lines) == 2

    def test_skip_whitespace_only(self):
        lines = ["   ", "\t", "data"]
        assert find_next_non_empty_line(0, lines) == 2

    def test_start_from_middle(self):
        lines = ["a", "", "", "b"]
        assert find_next_non_empty_line(1, lines) == 3

    def test_raises_when_all_empty(self):
        lines = ["", "  ", ""]
        with pytest.raises(NoNonEmptyLines):
            find_next_non_empty_line(0, lines)

    def test_raises_when_rest_empty(self):
        lines = ["data", "", ""]
        with pytest.raises(NoNonEmptyLines):
            find_next_non_empty_line(1, lines)


class TestSearchLines:
    """Test search_lines."""

    def test_find_delimiter(self):
        lines = ["line1", "line2", "DELIMITER here", "line4"]
        assert search_lines("DELIMITER", 0, lines) == 2

    def test_find_from_cursor(self):
        lines = ["DELIM early", "skip", "DELIM later"]
        assert search_lines("DELIM", 1, lines) == 2

    def test_raises_when_not_found(self):
        lines = ["line1", "line2", "line3"]
        with pytest.raises(UnableToLocateDelimiter) as exc_info:
            search_lines("MISSING", 0, lines)
        assert exc_info.value.delimiter == "MISSING"

    def test_partial_match(self):
        lines = ["this has ISPIN in it"]
        assert search_lines("ISPIN", 0, lines) == 0


class TestConvertVaspOutcarStress:
    """Test stress tensor conversion."""

    def test_correct_shape(self):
        stress = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        result = convert_vasp_outcar_stress(stress)
        assert result.shape == (6,)

    def test_wrong_shape_raises(self):
        with pytest.raises(ValueError, match="wrong shape"):
            convert_vasp_outcar_stress([1.0, 2.0, 3.0])

    def test_reordering(self):
        # Input: xx yy zz xy yz zx -> Output: xx yy zz yz zx xy
        stress = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        result = convert_vasp_outcar_stress(stress)
        # After reorder [0,1,2,4,5,3]: xx yy zz yz zx xy = 10 20 30 50 60 40
        # Then negated and scaled by 0.1 GPa
        import ase.units
        expected_order = np.array([-10.0, -20.0, -30.0, -50.0, -60.0, -40.0]) * 1e-1 * ase.units.GPa
        np.testing.assert_allclose(result, expected_order)


# ── Parser classes ─────────────────────────────────────────────────────


class TestSpinpol:
    """Test Spinpol header parser."""

    def test_spin_polarized(self):
        from asetools.parsers.vasp_outcar import Spinpol
        parser = Spinpol()
        lines = ["   ISPIN  =      2    spin polarized calculation?"]
        assert parser.has_property(0, lines) is True
        result = parser.parse(0, lines)
        assert result["spinpol"] is True

    def test_non_spin_polarized(self):
        from asetools.parsers.vasp_outcar import Spinpol
        parser = Spinpol()
        lines = ["   ISPIN  =      1    spin polarized calculation?"]
        result = parser.parse(0, lines)
        assert result["spinpol"] is False

    def test_no_match(self):
        from asetools.parsers.vasp_outcar import Spinpol
        parser = Spinpol()
        lines = ["   ENCUT  =   400.0"]
        assert parser.has_property(0, lines) is False


class TestIonsPerSpecies:
    """Test IonsPerSpecies header parser."""

    def test_parse_ion_types(self):
        from asetools.parsers.vasp_outcar import IonsPerSpecies
        parser = IonsPerSpecies()
        lines = ["   ions per type =              32  31   2"]
        assert parser.has_property(0, lines) is True
        result = parser.parse(0, lines)
        assert result["ion_types"] == [32, 31, 2]

    def test_single_species(self):
        from asetools.parsers.vasp_outcar import IonsPerSpecies
        parser = IonsPerSpecies()
        lines = ["   ions per type =               8"]
        result = parser.parse(0, lines)
        assert result["ion_types"] == [8]


class TestSpeciesTypes:
    """Test SpeciesTypes header parser."""

    def test_parse_standard_potcar(self):
        from asetools.parsers.vasp_outcar import SpeciesTypes
        parser = SpeciesTypes()
        lines = ["   POTCAR:    PAW_PBE Ni 02Aug2007"]
        result = parser.parse(0, lines)
        assert "Ni" in result["species"]

    def test_parse_tagged_potcar(self):
        from asetools.parsers.vasp_outcar import SpeciesTypes
        parser = SpeciesTypes()
        lines = ["   POTCAR:    PAW_PBE Fe_pv 02Aug2007"]
        result = parser.parse(0, lines)
        assert "Fe" in result["species"]

    def test_parse_hydrogen_variant(self):
        from asetools.parsers.vasp_outcar import SpeciesTypes
        parser = SpeciesTypes()
        lines = ["   POTCAR:    PAW_PBE H1.25 07Sep2000"]
        result = parser.parse(0, lines)
        assert "H" in result["species"]

    def test_double_counting_removal(self):
        from asetools.parsers.vasp_outcar import SpeciesTypes
        parser = SpeciesTypes()
        # Simulate OUTCAR double-listing
        lines = [
            "   POTCAR:    PAW_PBE Fe_pv 02Aug2007",
            "   POTCAR:    PAW_PBE O 08Apr2002",
            "   POTCAR:    PAW_PBE Fe_pv 02Aug2007",
            "   POTCAR:    PAW_PBE O 08Apr2002",
        ]
        for i in range(4):
            parser.parse(i, lines)
        species = parser.get_species()
        assert species == ["Fe", "O"]


class TestKpointHeader:
    """Test KpointHeader parser."""

    def test_parse_nkpts_nbands(self):
        from asetools.parsers.vasp_outcar import KpointHeader
        parser = KpointHeader()
        # Realistic line format from actual OUTCAR
        lines = [
            "k-points           NKPTS =      4   k-points in BZ     NKDIM =      4   number of bands    NBANDS=     48",
            "",
            " k-points in reciprocal lattice and weights: automatic mesh",
            "   0.25  0.25  0.25    0.250",
            "   0.25  0.25  0.75    0.250",
            "   0.25  0.75  0.75    0.250",
            "   0.75  0.75  0.75    0.250",
            "",
        ]
        assert parser.has_property(0, lines) is True
        result = parser.parse(0, lines)
        assert result["nkpts"] == 4
        assert result["nbands"] == 48


# ── Exception classes ──────────────────────────────────────────────────


class TestExceptions:
    """Test custom exception classes."""

    def test_no_non_empty_lines(self):
        with pytest.raises(NoNonEmptyLines):
            raise NoNonEmptyLines("test message")

    def test_unable_to_locate_delimiter(self):
        exc = UnableToLocateDelimiter("DELIM", "not found")
        assert exc.delimiter == "DELIM"
        assert "not found" in str(exc)


# ── read_constraints_from_file ─────────────────────────────────────────


class TestReadConstraints:
    """Test read_constraints_from_file."""

    def test_returns_none_for_empty_dir(self, tmp_path):
        from asetools.parsers.vasp_outcar import read_constraints_from_file
        result = read_constraints_from_file(tmp_path)
        assert result is None

    def test_reads_from_contcar(self, tmp_path):
        from asetools.parsers.vasp_outcar import read_constraints_from_file
        from ase.io import write
        from ase import Atoms

        # Write a minimal POSCAR/CONTCAR
        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[5, 5, 5], pbc=True)
        write(tmp_path / "CONTCAR", atoms, format="vasp")

        result = read_constraints_from_file(tmp_path)
        # No constraints set, so should be empty list
        assert result == []

    def test_reads_from_poscar_if_no_contcar(self, tmp_path):
        from asetools.parsers.vasp_outcar import read_constraints_from_file
        from ase.io import write
        from ase import Atoms

        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[5, 5, 5], pbc=True)
        write(tmp_path / "POSCAR", atoms, format="vasp")

        result = read_constraints_from_file(tmp_path)
        assert result == []


# ── Integration with real test data ────────────────────────────────────


class TestWithRealOutcar:
    """Integration tests using real OUTCAR files from tests/data/."""

    @pytest.fixture
    def outcar_path(self):
        from pathlib import Path
        p = Path(__file__).parent / "data" / "OUTCAR"
        if not p.exists():
            pytest.skip("Test OUTCAR not found")
        return p

    @pytest.fixture
    def outcar_v6_path(self):
        from pathlib import Path
        p = Path(__file__).parent / "data" / "OUTCAR_vasp6"
        if not p.exists():
            pytest.skip("Test OUTCAR_vasp6 not found")
        return p

    def test_outcar_readable(self, outcar_path):
        lines = outcar_path.read_text().splitlines()
        assert len(lines) > 100

    def test_find_scf_delimiter(self, outcar_path):
        from asetools.parsers.vasp_outcar import _OUTCAR_SCF_DELIM
        lines = outcar_path.read_text().splitlines()
        found = False
        for line in lines:
            if _OUTCAR_SCF_DELIM in line:
                found = True
                break
        assert found, "SCF delimiter not found in OUTCAR"

    def test_spinpol_detection(self, outcar_path):
        from asetools.parsers.vasp_outcar import Spinpol
        parser = Spinpol()
        lines = outcar_path.read_text().splitlines()
        found = False
        for i, line in enumerate(lines):
            if parser.has_property(i, lines):
                result = parser.parse(i, lines)
                assert "spinpol" in result
                found = True
                break
        assert found, "ISPIN line not found in OUTCAR"

    def test_ions_per_species_detection(self, outcar_path):
        from asetools.parsers.vasp_outcar import IonsPerSpecies
        parser = IonsPerSpecies()
        lines = outcar_path.read_text().splitlines()
        for i, line in enumerate(lines):
            if parser.has_property(i, lines):
                result = parser.parse(i, lines)
                assert len(result["ion_types"]) > 0
                assert all(isinstance(n, int) for n in result["ion_types"])
                return
        pytest.fail("ions per type line not found in OUTCAR")

    def test_species_detection(self, outcar_path):
        from asetools.parsers.vasp_outcar import SpeciesTypes
        parser = SpeciesTypes()
        lines = outcar_path.read_text().splitlines()
        found_count = 0
        for i, line in enumerate(lines):
            if parser.has_property(i, lines):
                parser.parse(i, lines)
                found_count += 1
        assert found_count > 0
        species = parser.get_species()
        assert len(species) > 0

    def test_vasp6_outcar_readable(self, outcar_v6_path):
        lines = outcar_v6_path.read_text().splitlines()
        assert len(lines) > 100
