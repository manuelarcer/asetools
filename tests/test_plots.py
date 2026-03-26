"""Tests for plotting/plots.py — PES plotting functions."""

import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from asetools.plotting.plots import (
    add_line_to_pes,
    beautify_pes_plot,
    validate_columns,
)

matplotlib.use("Agg")  # Non-interactive backend for testing


class TestValidateColumns:
    def test_all_present(self):
        df = pd.DataFrame({"A": [1], "B": [2], "C": [3]})
        assert validate_columns(df, ["A", "B"]) is True

    def test_missing_required_raises(self):
        df = pd.DataFrame({"A": [1]})
        with pytest.raises(ValueError, match="Missing required columns"):
            validate_columns(df, ["A", "MISSING"])

    def test_optional_warns(self):
        df = pd.DataFrame({"A": [1]})
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_columns(df, ["A"], optional_cols=["OPT"])
            assert len(w) == 1
            assert "Optional" in str(w[0].message)

    def test_empty_required_ok(self):
        df = pd.DataFrame({"A": [1]})
        assert validate_columns(df, []) is True


class TestAddLineToPes:
    @pytest.fixture
    def pes_data(self):
        return pd.DataFrame(
            {
                "E": [0.0, 0.5, -0.3, 0.2],
                "Type-Conf": ["M", "M", "M", "M"],
                "Label": ["IS", "TS1", "IM", "FS"],
            }
        )

    @pytest.fixture
    def ax(self):
        fig, ax = plt.subplots()
        yield ax
        plt.close(fig)

    def test_default_style(self, ax, pes_data):
        result = add_line_to_pes(ax, pes_data)
        assert result is ax
        # Should have plotted lines
        assert len(ax.lines) > 0

    def test_step_style(self, ax, pes_data):
        result = add_line_to_pes(ax, pes_data, style="step")
        assert result is ax
        assert len(ax.lines) > 0

    def test_invalid_style_raises(self, ax, pes_data):
        with pytest.raises(ValueError, match="style must be"):
            add_line_to_pes(ax, pes_data, style="invalid")

    def test_custom_color_and_label(self, ax, pes_data):
        add_line_to_pes(ax, pes_data, c="red", label="Test")
        # First line should have the label
        labels = [line.get_label() for line in ax.lines]
        assert "Test" in labels

    def test_with_transition_state(self, ax):
        data = pd.DataFrame(
            {
                "E": [0.0, 0.8, -0.3],
                "Type-Conf": ["M", "T", "M"],
            }
        )
        result = add_line_to_pes(ax, data)
        assert result is ax
        assert len(ax.lines) > 0

    def test_with_gap(self, ax):
        data = pd.DataFrame(
            {
                "E": [0.0, np.nan, -0.3],
                "Type-Conf": ["M", np.nan, "M"],
            }
        )
        result = add_line_to_pes(ax, data)
        assert result is ax

    def test_deprecated_col_param(self, ax, pes_data):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            add_line_to_pes(ax, pes_data, col="E")
            assert any("deprecated" in str(warning.message).lower() for warning in w)

    def test_custom_energy_col(self, ax):
        data = pd.DataFrame(
            {
                "energy": [0.0, 0.5],
                "Type-Conf": ["M", "M"],
            }
        )
        result = add_line_to_pes(ax, data, energy_col="energy")
        assert result is ax

    def test_missing_column_raises(self, ax):
        data = pd.DataFrame({"X": [1]})
        with pytest.raises(ValueError, match="Missing required"):
            add_line_to_pes(ax, data)


class TestBeautifyPesPlot:
    @pytest.fixture
    def ax(self):
        fig, ax = plt.subplots()
        ax.plot([0, 10], [0, 1])  # Dummy line
        yield ax
        plt.close(fig)

    def test_basic_beautify(self, ax):
        result = beautify_pes_plot(ax)
        assert result is ax

    def test_zero_line(self, ax):
        beautify_pes_plot(ax, zero=True)
        # Should have the original line + zero line
        assert len(ax.lines) >= 2

    def test_no_zero_line(self, ax):
        initial_lines = len(ax.lines)
        beautify_pes_plot(ax, zero=False)
        assert len(ax.lines) == initial_lines

    def test_limits(self, ax):
        beautify_pes_plot(ax, xlim=(0, 10), ylim=(-1, 1))
        assert ax.get_xlim() == (0, 10)
        assert ax.get_ylim() == (-1, 1)

    def test_no_frame(self, ax):
        beautify_pes_plot(ax, frame=False)
        assert not ax.spines["top"].get_visible()
        assert not ax.spines["right"].get_visible()
        assert not ax.spines["bottom"].get_visible()

    def test_with_labels(self, ax):
        data = pd.DataFrame(
            {
                "E": [0.0, 0.5],
                "Type-Conf": ["M", "M"],
                "Label": ["IS", "FS"],
            }
        )
        beautify_pes_plot(ax, data=data, show_labels=True)
        # Should have tick labels set
        ticks = ax.get_xticklabels()
        assert len(ticks) == 2

    def test_y_decimals(self, ax):
        beautify_pes_plot(ax, y_decimals=2)
        # Formatter should be set

    def test_legend(self, ax):
        ax.plot([0, 1], [0, 1], label="test")
        beautify_pes_plot(ax, leg=True)
        legend = ax.get_legend()
        assert legend is not None
