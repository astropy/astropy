# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for celestial coordinates
of astronomical objects. It also contains a framework for conversions
between coordinate systems.
"""

from .core import *
from .errors import *
from .utils import *

# isort: split
from . import core as _core
from . import errors as _errors
from . import utils as _utils

__all__ = []
__all__ += _core.__all__
__all__ += _errors.__all__
__all__ += _utils.__all__
