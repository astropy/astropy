# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for celestial coordinates
of astronomical objects. It also contains a framework for conversions
between coordinate systems.
"""

from .errors import *
from .angles import *
from .baseframe import *
from .attributes import *
from .distances import *
from .earth import *
from .transformations import *
from .builtin_frames import *
from .name_resolve import *
from .matching import *
from .representation import *
from .sky_coordinate import *
from .funcs import *
from .calculation import *
from .solar_system import *

# This is for backwards-compatibility -- can be removed in v3.0 when the
# deprecation warnings are removed
from .attributes import (TimeFrameAttribute, QuantityFrameAttribute,
                         CartesianRepresentationFrameAttribute)

__doc__ += builtin_frames._transform_graph_docs + """

.. note::

    The ecliptic coordinate systems (added in Astropy v1.1) have not been
    extensively tested for accuracy or consistency with other implementations of
    ecliptic coordinates.  We welcome contributions to add such testing, but in
    the meantime, users who depend on consistency with other implementations may
    wish to check test inputs against good datasets before using Astropy's
    ecliptic coordinates.

"""
