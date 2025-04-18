# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This module contains dictionaries that can be used to set a matplotlib
# plotting style. It is no longer documented/recommended as of Astropy v3.0
# but is kept here for backward-compatibility.

__all__ = ["astropy_mpl_style", "astropy_mpl_style_1"]

# Version 1 astropy plotting style for matplotlib
astropy_mpl_style_1 = {
    # Lines
    "lines.linewidth": 1.7,
    "lines.antialiased": True,
    # Patches
    "patch.linewidth": 1.0,
    "patch.facecolor": "#348ABD",
    "patch.edgecolor": "#CCCCCC",
    "patch.antialiased": True,
    # Images
    "image.cmap": "gist_heat",
    "image.origin": "upper",
    # Font
    "font.size": 12.0,
    # Axes
    "axes.facecolor": "#FFFFFF",
    "axes.edgecolor": "#AAAAAA",
    "axes.linewidth": 1.0,
    "axes.grid": True,
    "axes.titlesize": "x-large",
    "axes.labelsize": "large",
    "axes.labelcolor": "k",
    "axes.axisbelow": True,
    # Ticks
    "xtick.major.size": 0,
    "xtick.minor.size": 0,
    "xtick.major.pad": 6,
    "xtick.minor.pad": 6,
    "xtick.color": "#565656",
    "xtick.direction": "in",
    "ytick.major.size": 0,
    "ytick.minor.size": 0,
    "ytick.major.pad": 6,
    "ytick.minor.pad": 6,
    "ytick.color": "#565656",
    "ytick.direction": "in",
    # Legend
    "legend.fancybox": True,
    "legend.loc": "best",
    # Figure
    "figure.figsize": [8, 6],
    "figure.facecolor": "1.0",
    "figure.edgecolor": "0.50",
    "figure.subplot.hspace": 0.5,
    # Other
    "savefig.dpi": 72,
}
color_cycle = [
    "#348ABD",  # blue
    "#7A68A6",  # purple
    "#A60628",  # red
    "#467821",  # green
    "#CF4457",  # pink
    "#188487",  # turquoise
    "#E24A33",
]  # orange

try:
    # This is a dependency of matplotlib, so should be present if matplotlib
    # is installed.
    from cycler import cycler

    astropy_mpl_style_1["axes.prop_cycle"] = cycler("color", color_cycle)
except ImportError:
    astropy_mpl_style_1["axes.color_cycle"] = color_cycle


astropy_mpl_style = astropy_mpl_style_1
"""The most recent version of the astropy plotting style."""
