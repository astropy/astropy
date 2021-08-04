# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Unified I/O Registry.
"""

from . import core, default
from .core import *
from .default import *
from .default import _identifiers, _readers, _writers

__all__ = core.__all__ + default.__all__
