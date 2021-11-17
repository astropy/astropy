# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Unified I/O Registry.
"""

from . import base, compat, core, entry_points, interface
from .base import *
from .compat import *
from .compat import _identifiers, _readers, _writers  # for backwards compat
from .core import *
from .entry_points import *
from .interface import *

__all__ = core.__all__ + interface.__all__ + compat.__all__ + base.__all__ + entry_points.__all__


