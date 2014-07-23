# Licensed under a 3-clause BSD style license - see LICENSE.rst

""" This module contains changes to matplotlib rc parameters in order
to adopt an astropy plotting style.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib as mpl

def set_astropy_plot_style():
    """
    matplotlib configuration parameters for astropy plotting style.
    """

    # Lines
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['lines.antialiased'] = True

    # Patches
    mpl.rcParams['patch.linewidth'] = 0.5
    mpl.rcParams['patch.facecolor'] = '348ABD'
    mpl.rcParams['patch.edgecolor'] = 'eeeeee'
    mpl.rcParams['patch.antialiased'] = True

    # Font
    mpl.rcParams['font.family'] = 'monospace'
    mpl.rcParams['font.size'] = 12.0
    mpl.rcParams['font.monospace'] = ['Andale Mono', 'Nimbus Mono L', 'Courier New', 'Courier', 'Fixed', 'Terminal', 'monospace']

    # Axes
    mpl.rcParams['axes.facecolor'] = 'eeeeee'
    mpl.rcParams['axes.edgecolor'] = 'bcbcbc'
    mpl.rcParams['axes.linewidth'] = 1
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['axes.titlesize'] = 'x-large'
    mpl.rcParams['axes.labelsize'] = 'large'
    mpl.rcParams['axes.labelcolor'] = '555555'
    mpl.rcParams['axes.axisbelow'] = True
    mpl.rcParams['axes.color_cycle'] = ['348ABD', '7A68A6', 'A60628', '467821', 'CF4457', '188487', 'E24A33']
        # 348ABD : blue
        # 7A68A6 : purple
        # A60628 : red
        # 467821 : green
        # CF4457 : pink
        # 188487 : turquoise
        # E24A33 : orange

    # Ticks
    mpl.rcParams['xtick.major.size'] = 0
    mpl.rcParams['xtick.minor.size'] = 0
    mpl.rcParams['xtick.major.pad'] = 6
    mpl.rcParams['xtick.minor.pad'] = 6
    mpl.rcParams['xtick.color'] = '555555'
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.major.size'] = 0
    mpl.rcParams['ytick.minor.size'] = 0
    mpl.rcParams['ytick.major.pad'] = 6
    mpl.rcParams['ytick.minor.pad'] = 6
    mpl.rcParams['ytick.color'] = '555555'
    mpl.rcParams['ytick.direction'] = 'in'

    # Legend
    mpl.rcParams['legend.fancybox'] = True

    # Figure
    mpl.rcParams['figure.figsize'] = [8, 6]
    mpl.rcParams['figure.facecolor'] = '0.85'
    mpl.rcParams['figure.edgecolor'] = '0.50'
    mpl.rcParams['figure.subplot.hspace'] = 0.5

    # Keymap
    mpl.rcParams['keymap.fullscreen'] = 'f'
    mpl.rcParams['keymap.home'] = ['h', 'r', 'home']
    mpl.rcParams['keymap.back'] = ['left', 'c', 'backspace']
    mpl.rcParams['keymap.forward'] = ['right', 'v']
    mpl.rcParams['keymap.pan'] = 'p'
    mpl.rcParams['keymap.zoom'] = 'o'
    mpl.rcParams['keymap.save'] = 's'
    mpl.rcParams['keymap.grid'] = 'g'
    mpl.rcParams['keymap.yscale'] = 'l'
    mpl.rcParams['keymap.xscale'] = ['L', 'k']
    mpl.rcParams['keymap.all_axes'] = 'a'

    # Other
    mpl.rcParams['savefig.dpi'] = 72
    mpl.rcParams['interactive'] = True
    mpl.rcParams['toolbar'] = 'toolbar2'
    mpl.rcParams['timezone'] = 'UTC'
