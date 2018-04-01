# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants for Astropy v1.3 and earlier.
See :mod:`astropy.constants` for a complete listing of constants
defined in Astropy.
"""
import inspect
from . import utils as _utils
from . import codata2010, iau2012

_utils._set_c(codata2010, iau2012, inspect.getmodule(inspect.currentframe()))

# Clean up namespace
del inspect
del _utils
