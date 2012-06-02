# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The `nddata` subpackage provides the `~astropy.nddata.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU Data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.
"""

from astropy import setup_helpers

from .nddata import *

if not setup_helpers.is_in_build_mode():
    from .convolution.convolve import convolve,convolve_fft
    from .convolution.make_kernel import make_kernel
