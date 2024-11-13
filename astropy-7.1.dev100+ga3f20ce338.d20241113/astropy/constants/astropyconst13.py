# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants for Astropy v1.3 and earlier.
See :mod:`astropy.constants` for a complete listing of constants
defined in Astropy.
"""

from astropy.utils import find_current_module

from . import codata2010, iau2012
from . import utils as _utils

codata = codata2010
iaudata = iau2012

_utils._set_c(codata, iaudata, find_current_module())

# Clean up namespace
del find_current_module
del _utils
