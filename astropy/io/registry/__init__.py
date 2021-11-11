# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Unified I/O Registry.
"""

from . import base, compat, core, interface
from .base import *
from .compat import *  # for backwards compat
from .compat import _identifiers, _readers, _writers
from .core import *
from .interface import *

__all__ = core.__all__ + interface.__all__ + compat.__all__ + base.__all__
