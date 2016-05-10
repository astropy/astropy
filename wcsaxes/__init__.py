# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
WCSAxes implementation
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    from .core import *
    from .coordinate_helpers import CoordinateHelper
    from .coordinates_map import CoordinatesMap
    from .wcs_wrapper import WCS
    from .patches import *
