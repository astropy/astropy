# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants for Astropy v2.0.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
import inspect
from . import utils as _utils
from . import codata2014, iau2015

_utils._set_c(codata2014, iau2015, inspect.getmodule(inspect.currentframe()))

# Clean up namespace
del inspect
del _utils
