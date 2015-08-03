# Licensed under a 3-clause BSD style license - see LICENSE.rst

__doctest_requires__ = {'*': ['matplotlib']}

""" This module contains dictionaries that can be used to set a
matplotlib plotting style.
It is mostly here to allow a consistent plotting style in tutorials,
but can be used to prepare any matplotlib figure.

Using a matplotlib version > 1.4 you can do::

    >>> import matplotlib.pyplot as plt
    >>> from astropy.visualization import astropy_mpl_style
    >>> plt.style.use(astropy_mpl_style)

for older versions of matplotlib the following works::

    >>> import matplotlib as mpl
    >>> from astropy.visualization import astropy_mpl_style
    >>> mpl.rcParams.update(astropy_mpl_style)

This applies the astropy style on top of your existing matplotlib
default parameters. If you want an exactly reproducible plot (again, this
is useful if you are writing teaching material and you want the plot
to come out exactly the same, independent of the users configuration for
example), you should reset the matplotlib settings to the library defaults
*before* applying the astropy style, e.g.::

    >>> import matplotlib as mpl
    >>> from astropy.visualization import astropy_mpl_style
    >>> mpl.rcdefaults()
    >>> mpl.rcParams.update(astropy_mpl_style)
"""


astropy_mpl_style_1 = {

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
    'axes.color_cycle': [
        '#E24A33',   # orange
        '#348ABD',   # blue
        '#467821',   # green
        '#A60628',   # red
        '#7A68A6',   # purple
        '#CF4457',   # pink
        '#188487'    # turquoise
    ],

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
}
'''
Version 1 astropy plotting style for matplotlib.

This style improves some settings over the matplotlib default.
'''

astropy_mpl_style = astropy_mpl_style_1
'''
Most recent version of the astropy plotting style for matplotlib.

This style improves some settings over the matplotlib default.
'''
