# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for defining and converting
between different physical units.

This code is adapted from the `pynbody
<https://github.com/pynbody/pynbody>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
"""

# Import order matters here -- circular dependencies abound!
# Lots of things to import - go from more basic to advanced, so that
# whatever advanced ones need generally has been imported already;
# this also makes it easier to understand where most time is spent
# (e.g., using python -X importtime).

# isort: off
from .core import *
from .quantity import *

from . import astrophys, cgs, misc, photometric, si
from .function import units as function_units

from .si import *
from .astrophys import *
from .photometric import *
from .cgs import *
from .physical import *
from .function.units import *
from .misc import *

from .equivalencies import *

from .function.core import *
from .function.logarithmic import *

from .decorators import *
from .structured import *
# isort: on

del bases

# Enable the set of default units.  This notably does *not* include
# Imperial units.

set_enabled_units([si, cgs, astrophys, function_units, misc, photometric])


# -------------------------------------------------------------------------

def __getattr__(attr):
    if attr == "littleh":
        from astropy.units.astrophys import littleh
        return littleh
    elif attr == "with_H0":
        from astropy.units.equivalencies import with_H0
        return with_H0

    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}.")
