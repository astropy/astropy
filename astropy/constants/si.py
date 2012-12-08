# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units. The constants
available (with approximate values) are:
"""
# This docstring is extended by __init__.py

# The values of the constants in this module are defined in
# `astropy/_constants/si.py` in order to avoid circular dependencies, since
# Constant depends on Quantity, Quantity depends on Unit, and Unit depends on
# the value of some of the constants, so they need to be kept separately.

from .constant import Constant
from .._constants import si as _si
from .._constants.definition import ConstantDefinition

# Define actual Quantity-based Constants. Need to use _c in the loop instead
# of c to avoid overwriting the speed of light constant.
for nm, val in sorted(_si.__dict__.items()):
    if isinstance(val, ConstantDefinition):
        _c = Constant(val.value, val.units, val.uncertainty, val.name, val.reference)
        locals()[nm] = _c

del _si, nm, val, _c
