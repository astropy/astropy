# Licensed under a 3-clause BSD style license - see LICENSE.rst

""" This module contains changes to matplotlib rc parameters in order
to adopt an astropy plotting style.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


def set_astropy_plot_style(version=1):
    """
    Set matplotlib parameters for astropy plotting style.

    This methods sets some :mod:`matplotlib` defaults. It is mostly here to
    allow a consistent plotting style in tutorials, but can be used to
    prepare any :mod:`matplotlib` figure.

    This method is only useful, if the module :mod:`matplotlib` can be imported.

    Parameters
    ----------
    version : integer
        Different default plotting styles for astropy are labeled with 
        version numbers. Currently, only one style (``version=1``) is
        implemented.
    """
    import matplotlib as mpl
    if version == 1:

        params = {

            # Lines
            'lines.linewidth': 1.7,
            'lines.antialiased': True,

            # Patches
            'patch.linewidth': 1.0,
            'patch.facecolor': '#348ABD',
            'patch.edgecolor': '#CCCCCC',
            'patch.antialiased': True,

            # images
            'image.cmap': 'gist_heat',
            'image.origin': 'upper',

            # Font
            'font.size': 12.0,

            # Axes
            'axes.facecolor': '#FFFFFF',
            'axes.edgecolor': '#AAAAAA',
            'axes.linewidth': 1.0,
            'axes.grid': True,
            'axes.titlesize': 'x-large',
            'axes.labelsize': 'large',
            'axes.labelcolor': '#FFFFFF',
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
            'xtick.color': '#565656',
            'xtick.direction': 'in',
            'ytick.major.size': 0,
            'ytick.minor.size': 0,
            'ytick.major.pad': 6,
            'ytick.minor.pad': 6,
            'ytick.color': '#565656',
            'ytick.direction': 'in',

            # Legend
            'legend.fancybox': True,
            'legend.loc': 'best',

            # Figure
            'figure.figsize': [8, 6],
            'figure.facecolor': '1.0',
            'figure.edgecolor': '0.50',
            'figure.subplot.hspace': 0.5,

            # Other
            'savefig.dpi': 72,
            'interactive': True,
            'toolbar': 'toolbar2',
            'timezone': 'UTC'
        }

    else:
        raise ValueError('Unrecognized astropy_plot_style version.')

    mpl.rcdefaults()
    mpl.rcParams.update(params)
