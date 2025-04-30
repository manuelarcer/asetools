#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import numpy as np

def add_line_to_pes(ax, data, c='k', label=None, indexes=None):
    # ax is the plt axis object
    # data is the pandas df in the format established
    # c is the color of the line
    count = 0   # Count intermediates
    gap = False
    for i, t in enumerate(data['Type-Conf']):
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
                plt.plot([x,x+1], [data.E[i], data.E[i]], color=c, linewidth=3.5, label=label)
            else:
                plt.plot([x,x+1], [data.E[i], data.E[i]], color=c, linewidth=3.5)
            count += 1      # Count increases with each plotted M point
            if i < ( len(data) - 1) and data['Type-Conf'][i+1] == 'M':
                if gap:
                    plt.plot([x+1, x+2], [data.E[i], data.E[i+1]], '-', color=c, linewidth=0.75)
                    gap = False
                else:
                    #print(x+1, x+2)
                    plt.plot([x+1, x+2], [data.E[i], data.E[i+1]], '-', color=c, linewidth=0.75)
        if t == 'T':
            if indexes is None:
                x = count*2 + 0.5       # Count should be +1 from previous M point
                X_Y_Spline = make_interp_spline([x-0.5, x, x+0.5], [data.E[i-1], data.E[i], data.E[i+1]], k=2)
                X_ = np.linspace(x-0.5, x+0.5, 50)
            else:       # Both X_Y_Spline and X_ variables should change if using "indexes" as parameter
                diff_x = indexes[count]*2 - indexes[count-1]*2 - 1      # Diff between point before and after
                x = indexes[count]*2 + 1 - diff_x / 2.
                print(diff_x, x)
                X_Y_Spline = make_interp_spline([x - diff_x/2., x, x + diff_x/2.], [data.E[i-1], data.E[i], data.E[i+1]], k=2)
                X_ = np.linspace(x - diff_x/2., x + diff_x/2., 50)
            
            Y_ = X_Y_Spline(X_)
            plt.plot(X_, Y_, '-', color=c, linewidth=0.75)
            #plt.plot(x, data.E[i], 'o', color=c)
    return ax

def beautify_pes_plot(ax, xlim=None, ylim=None, zero=True, leg=False, fs=12):
    # ax
    # xlim, ylim: the limits of the plot
    # zero: draw the horizontal line at Zero
    # leg: print legend

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(axis='both', labelsize=fs)
    #ax.set_yticks(fontsize=fs)
    if zero:
        ax.plot( ax.get_xlim(), [0,0], '--', color='lightgray', linewidth=1 )
    if leg:
        ax.legend(frameon=False, fontsize=fs)

    return ax