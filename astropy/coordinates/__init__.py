# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for celestial coordinates
of astronomical objects. It also contains a framework for conversions
between coordinate systems.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .errors import *
from .angles import *
from .coordsystems import *
from .distances import *
from .transformations import *
from .builtin_systems import *
from .name_resolve import *
from .matching import *
from .old_builtin_systems_names import *  # TODO: remove this in next version, along with module file

__doc__ += builtin_systems._transform_graph_docs
