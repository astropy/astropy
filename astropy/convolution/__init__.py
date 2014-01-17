# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .core import *
from .kernels import *
from .utils import discretize_model

try:
    # Not guaranteed available at setup time
    from .convolve import convolve, convolve_fft
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
