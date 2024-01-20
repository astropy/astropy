"""Coordinate transformations.

This module contains a general framework for defining graphs of transformations
between coordinates, suitable for either spatial coordinates or more generalized
coordinate systems.

The fundamental idea is that each class is a node in the transformation graph,
and transitions from one node to another are defined as functions (or methods)
wrapped in transformation objects.

This module also includes more specific transformation classes for
celestial/spatial coordinate frames, generally focused around matrix-style
transformations that are typically how the algorithms are defined.
"""

from . import affine, base, composite, function, graph
from .affine import *
from .base import *
from .composite import *
from .function import *
from .graph import *

__all__: list[str] = []
__all__ += graph.__all__
__all__ += base.__all__
__all__ += composite.__all__
__all__ += affine.__all__
__all__ += function.__all__
