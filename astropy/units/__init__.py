# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for defining and converting
between different physical units.

This code is adapted from the `pynbody
<https://github.com/pynbody/pynbody>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from .core import *
from .quantity import *
from .decorators import *

from . import si
from . import cgs
from . import astrophys
from .function import units as function_units

from .si import *
from .astrophys import *
from .cgs import *
from .physical import *
from .function.units import *

from .equivalencies import *

from .function.core import *
from .function.logarithmic import *
from .function import magnitude_zero_points

del bases

# Enable the set of default units.  This notably does *not* include
# Imperial units.

set_enabled_units([si, cgs, astrophys, function_units])
