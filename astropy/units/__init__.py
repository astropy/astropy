# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for defining and converting
between different physical units.

This code is adapted from the `pynbody
<https://github.com/pynbody/pynbody>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
"""

from . import (
    astrophys,
    cgs,
    core,
    decorators,
    misc,
    photometric,
    physical,
    quantity,
    si,
    structured,
)
from .astrophys import *
from .cgs import *
from .core import *
from .core import set_enabled_units
from .decorators import *
from .misc import *
from .photometric import *
from .physical import *
from .quantity import *
from .si import *
from .structured import *

# Import order matters here -- circular dependencies abound!
# isort: split
from . import equivalencies, function
from .equivalencies import *
from .function import *

__all__ = []
__all__ += core.__all__
__all__ += quantity.__all__
__all__ += decorators.__all__
__all__ += structured.__all__
__all__ += si.__all__
__all__ += cgs.__all__
__all__ += astrophys.__all__
__all__ += physical.__all__
__all__ += misc.__all__
__all__ += photometric.__all__
__all__ += structured.__all__
__all__ += equivalencies.__all__
__all__ += function.__all__
# We sort `__all__` for the docs. We use `globals` to avoid confusing any static
# analysis tools.
globals()["__all__"].sort()

# Enable the set of default units.  This notably does *not* include
# Imperial units.
set_enabled_units([si, cgs, astrophys, function.units, misc, photometric])
