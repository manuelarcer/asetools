"""Tests for parsers/vasp_outcar.py utility functions."""

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


class TestCheckLine:
    """Tests for _check_line: fixes VASP numeric formatting where
    digits touch minus signs (e.g., '1.23-4.56' -> '1.23 -4.56')."""

    def test_normal_line_unchanged(self):
        line = "  1.234  -5.678  9.012"
        assert _check_line(line) == line

    def test_joined_numbers_separated(self):
        line = "1.234-5.678"
        result = _check_line(line)
        assert result == "1.234 -5.678"

    def test_multiple_joined_numbers(self):
        line = "1.23-4.56-7.89"
        result = _check_line(line)
        # Each digit-minus-digit gets a space
        assert " -" in result
        parts = result.split()
        assert len(parts) == 3

    def test_empty_line(self):
        assert _check_line("") == ""

    def test_no_digits_with_minus(self):
        line = "some text - more text"
        assert _check_line(line) == line


class TestFindNextNonEmptyLine:
    def test_current_line_non_empty(self):
        lines = ["hello", "world"]
        assert find_next_non_empty_line(0, lines) == 0

    def test_skip_empty_lines(self):
        lines = ["", "  ", "hello"]
        assert find_next_non_empty_line(0, lines) == 2

    def test_start_from_offset(self):
        lines = ["first", "", "third"]
        assert find_next_non_empty_line(1, lines) == 2

    def test_all_empty_raises(self):
        lines = ["", "  ", "\n"]
        with pytest.raises(NoNonEmptyLines):
            find_next_non_empty_line(0, lines)

    def test_cursor_at_end_raises(self):
        lines = ["hello"]
        with pytest.raises(NoNonEmptyLines):
            find_next_non_empty_line(1, lines)


class TestSearchLines:
    def test_find_delimiter(self):
        lines = ["line 1", "line 2", "DELIMITER here", "line 4"]
        assert search_lines("DELIMITER", 0, lines) == 2

    def test_find_from_offset(self):
        lines = ["DELIM", "other", "DELIM again"]
        assert search_lines("DELIM", 1, lines) == 2

    def test_not_found_raises(self):
        lines = ["line 1", "line 2"]
        with pytest.raises(UnableToLocateDelimiter) as exc_info:
            search_lines("MISSING", 0, lines)
        assert exc_info.value.delimiter == "MISSING"

    def test_partial_match(self):
        lines = ["FREE ENERGIE OF THE ION-ELECTRON SYSTEM"]
        assert search_lines("FREE ENERGIE", 0, lines) == 0


class TestConvertVaspOutcarStress:
    def test_valid_stress(self):
        stress = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        result = convert_vasp_outcar_stress(stress)
        assert result.shape == (6,)
        # Check reordering: [0,1,2,4,5,3] -> values at positions 0,1,2,4,5,3
        # Then negated and multiplied by 1e-1 * GPa
        import ase

        factor = -1e-1 * ase.units.GPa
        expected_order = [1.0, 2.0, 3.0, 5.0, 6.0, 4.0]
        np.testing.assert_allclose(result, np.array(expected_order) * factor)

    def test_wrong_shape_raises(self):
        with pytest.raises(ValueError, match="wrong shape"):
            convert_vasp_outcar_stress([1.0, 2.0, 3.0])

    def test_wrong_shape_2d_raises(self):
        with pytest.raises(ValueError, match="wrong shape"):
            convert_vasp_outcar_stress([[1, 2, 3], [4, 5, 6]])


class TestUnableToLocateDelimiter:
    def test_stores_delimiter(self):
        exc = UnableToLocateDelimiter("MY_DELIM", "not found")
        assert exc.delimiter == "MY_DELIM"
        assert "not found" in str(exc)
