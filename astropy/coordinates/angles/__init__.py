# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for celestial coordinates
of astronomical objects. It also contains a framework for conversions
between coordinate systems.
"""

from .core import Angle, Latitude, Longitude
from .errors import (
    BoundsError,
    IllegalHourError,
    IllegalHourWarning,
    IllegalMinuteError,
    IllegalMinuteWarning,
    IllegalSecondError,
    IllegalSecondWarning,
    RangeError,
)
from .utils import (
    angular_separation,
    golden_spiral_grid,
    offset_by,
    position_angle,
    uniform_spherical_random_surface,
    uniform_spherical_random_volume,
)

__all__ = [
    "Angle",
    "Latitude",
    "Longitude",
    # utilities
    "angular_separation",
    "position_angle",
    "offset_by",
    "golden_spiral_grid",
    "uniform_spherical_random_surface",
    "uniform_spherical_random_volume",
    # errors
    "RangeError",
    "BoundsError",
    "IllegalHourError",
    "IllegalMinuteError",
    "IllegalSecondError",
    "IllegalHourWarning",
    "IllegalMinuteWarning",
    "IllegalSecondWarning",
]
