# Licensed under a 3-clause BSD style license - see LICENSE.rst
try:
    # The ERFA wrappers are not guaranteed available at setup time
    from .constants import *
    from .core import *
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
