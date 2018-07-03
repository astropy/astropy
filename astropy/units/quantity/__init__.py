# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines the `Quantity` object, which represents a number with some
associated units. `Quantity` objects support operations like ordinary numbers,
but will deal with unit conversions internally.
"""

from .quantity import _UNIT_NOT_INITIALISED, _unquantify_allclose_arguments, conf
from .quantity import *
