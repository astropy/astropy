# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines the `Quantity` object, which represents a number with some
associated units. `Quantity` objects support operations like ordinary numbers,
but will deal with unit conversions internally.
"""

from . import quantity, utils
from .quantity import *
from .utils import *

# Extra imports for help for those who imported private atttributes
# from the old quantity.py
from .quantity import _UNIT_NOT_INITIALISED, conf
from .utils import _unquantify_allclose_arguments


__all__ = quantity.__all__ + utils.__all__
