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

from .affine import (
    AffineTransform,
    BaseAffineTransform,
    DynamicMatrixTransform,
    StaticMatrixTransform,
)
from .base import CoordinateTransform
from .composite import CompositeTransform
from .function import FunctionTransform, FunctionTransformWithFiniteDifference
from .graph import TransformGraph

__all__ = [
    "TransformGraph",
    # transformations
    "CoordinateTransform",
    "FunctionTransform",
    "BaseAffineTransform",
    "AffineTransform",
    "StaticMatrixTransform",
    "DynamicMatrixTransform",
    "FunctionTransformWithFiniteDifference",
    "CompositeTransform",
]
