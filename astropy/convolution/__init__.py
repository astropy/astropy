# Licensed under a 3-clause BSD style license - see LICENSE.rst


from .core import *
from .kernels import *
from .utils import discretize_model

try:
    # Not guaranteed available at setup time
    from .convolve import convolve, convolve_fft, interpolate_replace_nans, convolve_models
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
