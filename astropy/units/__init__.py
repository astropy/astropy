"""
Light-weight physical units module.

This code is adapted from the `pynbody
<http://code.google.com/p/pynbody/>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
"""

from .core import *
from .standard_units import *

# After importing standard_units, set the unit namespace to this
# top-level module so that new units are added here.
UnitBase._set_namespace(globals())
