# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Unified I/O Registry.
"""

from . import compat, core
from .compat import *
from .core import *

__all__ = core.__all__ + compat.__all__
