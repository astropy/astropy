# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains helper functions and classes for handling metadata."""

from . import core as _core
from . import exceptions as _exceptions
from . import merge as _merge
from . import utils as _utils
from .core import *
from .exceptions import *
from .merge import *
from .utils import *

__all__ = []
__all__ += _core.__all__
__all__ += _merge.__all__
__all__ += _exceptions.__all__
__all__ += _utils.__all__
