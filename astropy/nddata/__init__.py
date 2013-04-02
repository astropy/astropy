# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The `nddata` subpackage provides the `~astropy.nddata.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU Data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.
"""

from .nddata import *
from .nduncertainty import *
from .flag_collection import *

try:
    # Not guaranteed available at setup time
    from .convolution.convolve import convolve,convolve_fft
    from .convolution.make_kernel import make_kernel
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
