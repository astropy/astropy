# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for defining and converting
between different physical units.

This code is adapted from the `pynbody
<https://github.com/pynbody/pynbody>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
"""
# Lots of things to import - go from more basic to advanced, so that
# whatever advanced ones need generally has been imported already;
# this helps prevent circular imports and makes it easier to understand
# where most time is spent (e.g., using python -X importtime).
from .core import *
from .quantity import *

from . import si
from . import cgs
from . import astrophys
from . import photometric
from .function import units as function_units

from .si import *
from .astrophys import *
from .photometric import *
from .cgs import *
from .physical import *
from .function.units import *

from .equivalencies import *

from .function.core import *
from .function.logarithmic import *
from .function import magnitude_zero_points

from .decorators import *

del bases

# Enable the set of default units.  This notably does *not* include
# Imperial units.

set_enabled_units([si, cgs, astrophys, function_units, photometric])
