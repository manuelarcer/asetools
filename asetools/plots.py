#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import numpy as np

# Default column names
DEFAULT_LABEL_COL = 'Label'
DEFAULT_TYPE_COL = 'Type-Conf'
DEFAULT_ENERGY_COL = 'E'
DEFAULT_NPCET_COL = 'nPCET'


def validate_columns(data, required_cols, optional_cols=None):
    """
    Validate that required columns exist in the DataFrame.

    Parameters
    ----------
    data : pandas.DataFrame
        The DataFrame to validate
    required_cols : list
        List of column names that must exist
    optional_cols : list, optional
        List of column names that are optional (will warn if missing)

    Returns
    -------
    bool
        True if all required columns exist

    Raises
    ------
    ValueError
        If any required column is missing
    """
    missing = [col for col in required_cols if col not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. "
                        f"Available columns: {list(data.columns)}")

    if optional_cols:
        missing_optional = [col for col in optional_cols if col not in data.columns]
        if missing_optional:
            import warnings
            warnings.warn(f"Optional columns not found: {missing_optional}")

    return True


def add_line_to_pes(ax, data, energy_col=None, type_col=None, c='k', label=None,
                   indexes=None, col=None, style='default', lw=None, lw_connector=None):
    """
    Add a potential energy surface line to a matplotlib axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axis object to plot on
    data : pandas.DataFrame
        DataFrame containing the PES data
    energy_col : str, optional
        Column name for energy values (default: 'E')
    type_col : str, optional
        Column name for configuration type M/T (default: 'Type-Conf')
    c : str, optional
        Color of the line (default: 'k')
    label : str, optional
        Label for the legend
    indexes : list, optional
        Custom x-positions for intermediates
    col : str, optional
        Deprecated: use energy_col instead. Kept for backward compatibility.
    style : str, optional
        Plot style: 'default' (thick horizontal lines with diagonal connectors)
        or 'step' (horizontal and vertical lines only). Default: 'default'
    lw : float, optional
        Line width for energy level lines (default: 3.5 for 'default', 1.5 for 'step')
    lw_connector : float, optional
        Line width for connector lines (default: 0.75 for 'default', same as lw for 'step')

    Returns
    -------
    matplotlib.axes.Axes
        The axis object with the added line
    """
    # Handle backward compatibility with 'col' parameter
    if col is not None:
        import warnings
        warnings.warn("Parameter 'col' is deprecated, use 'energy_col' instead",
                     DeprecationWarning)
        if energy_col is None:
            energy_col = col

    # Set defaults
    if energy_col is None:
        energy_col = DEFAULT_ENERGY_COL
    if type_col is None:
        type_col = DEFAULT_TYPE_COL

    # Validate style parameter
    if style not in ('default', 'step'):
        raise ValueError(f"style must be 'default' or 'step', got '{style}'")

    # Validate required columns
    validate_columns(data, [energy_col, type_col])

    if style == 'step':
        return _add_line_step_style(ax, data, energy_col, type_col, c, label, indexes, lw)
    else:
        return _add_line_default_style(ax, data, energy_col, type_col, c, label, indexes, lw, lw_connector)


def _add_line_default_style(ax, data, energy_col, type_col, c, label, indexes, lw=None, lw_connector=None):
    """Default style: thick horizontal lines with diagonal connectors."""
    # Set default line widths
    if lw is None:
        lw = 3.5
    if lw_connector is None:
        lw_connector = 0.75

    count = 0   # Count intermediates
    gap = False
    for i, t in enumerate(data[type_col]):
        if pd.isnull(t):    # If row is empty, increase count and continue
            count += 1
            gap = True
            continue
        if t == 'M':
            if indexes is None:
                x = count*2 + 1
            else:
                x = indexes[count]*2 + 1
            if count == 0:
                ax.plot([x,x+1], [data[energy_col][i], data[energy_col][i]], color=c, linewidth=lw, label=label)
            else:
                ax.plot([x,x+1], [data[energy_col][i], data[energy_col][i]], color=c, linewidth=lw)
            count += 1      # Count increases with each plotted M point
            if i < (len(data) - 1) and data[type_col][i+1] == 'M':
                if gap:
                    ax.plot([x+1, x+2], [data[energy_col][i], data[energy_col][i+1]], '-', color=c, linewidth=lw_connector)
                    gap = False
                else:
                    ax.plot([x+1, x+2], [data[energy_col][i], data[energy_col][i+1]], '-', color=c, linewidth=lw_connector)
        if t == 'T':
            if indexes is None:
                x = count*2 + 0.5       # Count should be +1 from previous M point
                X_Y_Spline = make_interp_spline([x-0.5, x, x+0.5], [data[energy_col][i-1], data[energy_col][i], data[energy_col][i+1]], k=2)
                X_ = np.linspace(x-0.5, x+0.5, 50)
            else:       # Both X_Y_Spline and X_ variables should change if using "indexes" as parameter
                diff_x = indexes[count]*2 - indexes[count-1]*2 - 1      # Diff between point before and after
                x = indexes[count]*2 + 1 - diff_x / 2.
                print(diff_x, x)
                X_Y_Spline = make_interp_spline([x - diff_x/2., x, x + diff_x/2.], [data[energy_col][i-1], data[energy_col][i], data[energy_col][i+1]], k=2)
                X_ = np.linspace(x - diff_x/2., x + diff_x/2., 50)

            Y_ = X_Y_Spline(X_)
            ax.plot(X_, Y_, '-', color=c, linewidth=lw_connector)
    return ax


def _add_line_step_style(ax, data, energy_col, type_col, c, label, indexes, lw=None):
    """Step style: horizontal and vertical lines only (no diagonal connectors)."""
    # Set default line width
    if lw is None:
        lw = 1.5

    # Filter to only M points (minima), ignoring T points and gaps
    m_points = []
    count = 0
    for i, t in enumerate(data[type_col]):
        if pd.isnull(t):
            count += 1
            continue
        if t == 'M':
            if indexes is None:
                x = count * 2 + 1.5  # Center of the horizontal segment
            else:
                x = indexes[count] * 2 + 1.5
            m_points.append((x, data[energy_col][i]))
            count += 1

    if not m_points:
        return ax

    # Build step plot coordinates with proper vertical connectors
    # Each intermediate gets a horizontal segment of width 2, centered at x
    # The vertical transition happens at the midpoint between consecutive centers
    x_coords = []
    y_coords = []

    for idx, (x, y) in enumerate(m_points):
        if idx == 0:
            # First point: horizontal line from x-1 to x+1 (width 2, centered at x)
            x_coords.extend([x - 1, x + 1])
            y_coords.extend([y, y])
        else:
            prev_x = m_points[idx - 1][0]
            prev_y = y_coords[-1]
            # Midpoint between previous center and current center
            x_mid = (prev_x + x) / 2
            # Horizontal from previous end to midpoint at prev_y
            x_coords.append(x_mid)
            y_coords.append(prev_y)
            # Vertical at midpoint from prev_y to current y
            x_coords.append(x_mid)
            y_coords.append(y)
            # Horizontal from midpoint to current end at current y
            if idx == len(m_points) - 1:
                # Last point: extend to x+1
                x_coords.append(x + 1)
            else:
                # Not last: extend to midpoint with next
                x_coords.append(x + 1)
            y_coords.append(y)

    # Plot the step line
    ax.plot(x_coords, y_coords, color=c, linewidth=lw, label=label)

    return ax

def beautify_pes_plot(ax, xlim=None, ylim=None, zero=True, leg=False, fs=12,
                      data=None, label_col=None, type_col=None, npcet_col=None,
                      show_labels=False, show_npcet=False, indexes=None, frame=True):
    """
    Beautify a PES plot by removing unnecessary spines and adding optional labels.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axis object
    xlim : tuple, optional
        x-axis limits (min, max)
    ylim : tuple, optional
        y-axis limits (min, max)
    zero : bool, optional
        Draw horizontal line at y=0 (default: True)
    leg : bool, optional
        Show legend (default: False)
    fs : int, optional
        Font size (default: 12)
    data : pandas.DataFrame, optional
        DataFrame with PES data (required if show_labels or show_npcet is True)
    label_col : str, optional
        Column name for species labels (default: 'Label')
    type_col : str, optional
        Column name for M/T type (default: 'Type-Conf')
    npcet_col : str, optional
        Column name for nPCET values (default: 'nPCET')
    show_labels : bool, optional
        Show species labels on x-axis (default: False)
    show_npcet : bool, optional
        Show nPCET values below labels (default: False)
    indexes : list, optional
        Custom x-positions for intermediates (must match add_line_to_pes)
    frame : bool, optional
        Keep plot frame/spines visible (default: True). If False, removes
        top, right, and bottom spines.

    Returns
    -------
    matplotlib.axes.Axes
        The beautified axis object
    """
    # Set column name defaults
    if label_col is None:
        label_col = DEFAULT_LABEL_COL
    if type_col is None:
        type_col = DEFAULT_TYPE_COL
    if npcet_col is None:
        npcet_col = DEFAULT_NPCET_COL

    if not frame:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(axis='both', labelsize=fs)

    if zero:
        ax.plot(ax.get_xlim(), [0, 0], '--', color='lightgray', linewidth=1)
    if leg:
        ax.legend(fontsize=fs)

    # Add x-axis labels if requested
    if (show_labels or show_npcet) and data is not None:
        # Validate columns
        required = [type_col]
        if show_labels:
            required.append(label_col)
        if show_npcet:
            required.append(npcet_col)
        validate_columns(data, required)

        # Calculate tick positions and labels
        tick_positions = []
        tick_labels = []
        count = 0

        for i, t in enumerate(data[type_col]):
            if pd.isnull(t):
                count += 1
                continue
            if t == 'M':
                if indexes is None:
                    x = count * 2 + 1.5  # Center of the horizontal line
                else:
                    x = indexes[count] * 2 + 1.5
                tick_positions.append(x)

                # Build label
                label_parts = []
                if show_labels and label_col in data.columns:
                    label_parts.append(str(data[label_col][i]))
                if show_npcet and npcet_col in data.columns:
                    npcet_val = data[npcet_col][i]
                    if not pd.isnull(npcet_val):
                        label_parts.append(f"n={int(npcet_val)}")

                tick_labels.append('\n'.join(label_parts))
                count += 1

        # Set ticks and labels
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=fs-2)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', length=0)

    return ax