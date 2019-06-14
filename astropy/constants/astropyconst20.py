# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants for Astropy v2.0.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
from astropy.utils import find_current_module
from . import utils as _utils
from . import codata2014, iau2015

_utils._set_c(codata2014, iau2015, find_current_module())

# Clean up namespace
del find_current_module
del _utils
