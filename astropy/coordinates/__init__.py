# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for celestial coordinates
of astronomical objects. It also contains a framework for conversions
between coordinate systems.
"""

from .errors import *
from .angles import *
from .coordsystems import *
from .distances import *
from .transformations import *
from .builtin_systems import *

__doc__ += builtin_systems._transform_graph_docs