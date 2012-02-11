# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The `nddata` package manages n-dimensional array based data (e.g. CCD
images, IFU Data, grid-based simulation data, ...)
"""

from astropy import setup_helpers

from .nddata import *

if setup_helpers.is_in_build_mode():
    pass
else:
    from .convolution.convolve import convolve
