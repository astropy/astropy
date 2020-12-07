# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines miscellaneous units.  They are also
available in the `astropy.units` namespace.
"""


from . import si
from astropy.constants import si as _si
from .core import (UnitBase, def_unit, si_prefixes, binary_prefixes,
                   set_enabled_units)

# To ensure si units of the constants can be interpreted.
set_enabled_units([si])

import numpy as _numpy

_ns = globals()
