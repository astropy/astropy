# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for work with time series
data sets.

"""

from ._astropy_init import *

if not _ASTROPY_SETUP_:
    from .core import *  # noqa
    from .sampled import *  # noqa
    from .binned import *  # noqa
    from . import io  # noqa
    from .downsample import *  # noqa
