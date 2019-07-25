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
