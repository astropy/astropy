# Licensed under a 3-clause BSD style license - see LICENSE.rst

""" This module contains changes to matplotlib rc parameters in order
to adopt an astropy plotting style.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib as mpl

def set_astropy_plot_style(version=1):
    """
    matplotlib configuration parameters for astropy plotting style.
    """

    if version == 1:

        params = {

            # Lines
            'lines.linewidth': 1.0,
            'lines.antialiased': True,

            # Patches
            'patch.linewidth': 0.5,
            'patch.facecolor': '#348ABD',
            'patch.edgecolor': '#EEEEEE',
            'patch.antialiased': True,

            # Font
            'font.family': 'monospace',
            'font.size': 12.0,
            'font.monospace': ['Andale Mono', 'Nimbus Mono L', 'Courier New', 'Courier', 'Fixed', 'Terminal', 'monospace'],

            # Axes
            'axes.facecolor': '#EEEEEE',
            'axes.edgecolor': '#BCBCBC',
            'axes.linewidth': 1,
            'axes.grid': True,
            'axes.titlesize': 'x-large',
            'axes.labelsize': 'large',
            'axes.labelcolor': '#555555',
            'axes.axisbelow': True,
            'axes.color_cycle': ['#348ABD', '#7A68A6', '#A60628', '#467821', '#CF4457', '#188487', '#E24A33'],
                # 348ABD : blue
                # 7A68A6 : purple
                # A60628 : red
                # 467821 : green
                # CF4457 : pink
                # 188487 : turquoise
                # E24A33 : orange

            # Ticks
            'xtick.major.size': 0,
            'xtick.minor.size': 0,
            'xtick.major.pad': 6,
            'xtick.minor.pad': 6,
            'xtick.color': '#555555',
            'xtick.direction': 'in',
            'ytick.major.size': 0,
            'ytick.minor.size': 0,
            'ytick.major.pad': 6,
            'ytick.minor.pad': 6,
            'ytick.color': '#555555',
            'ytick.direction': 'in',

            # Legend
            'legend.fancybox': True,

            # Figure
            'figure.figsize': [8, 6],
            'figure.facecolor': '0.85',
            'figure.edgecolor': '0.50',
            'figure.subplot.hspace': 0.5,

            # Keymap
            'keymap.fullscreen': 'f',
            'keymap.home': ['h', 'r', 'home'],
            'keymap.back': ['left', 'c', 'backspace'],
            'keymap.forward': ['right', 'v'],
            'keymap.pan': 'p',
            'keymap.zoom': 'o',
            'keymap.save': 's',
            'keymap.grid': 'g',
            'keymap.yscale': 'l',
            'keymap.xscale': ['L', 'k'],
            'keymap.all_axes': 'a',

            # Other
            'savefig.dpi': 72,
            'interactive': True,
            'toolbar': 'toolbar2',
            'timezone': 'UTC'
        }

    else:
        raise ValueError('Unrecognized astropy_plot_style verion')

    mpl.rcParams.update(params)
