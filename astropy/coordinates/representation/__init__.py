"""
In this module, we define the coordinate representation classes, which are
used to represent low-level cartesian, spherical, cylindrical, and other
coordinates.
"""

from .base import BaseDifferential, BaseRepresentation, BaseRepresentationOrDifferential
from .cartesian import CartesianDifferential, CartesianRepresentation
from .cylindrical import CylindricalDifferential, CylindricalRepresentation
from .spherical import (
    BaseSphericalCosLatDifferential,
    BaseSphericalDifferential,
    PhysicsSphericalDifferential,
    PhysicsSphericalRepresentation,
    RadialDifferential,
    RadialRepresentation,
    SphericalCosLatDifferential,
    SphericalDifferential,
    SphericalRepresentation,
    UnitSphericalCosLatDifferential,
    UnitSphericalDifferential,
    UnitSphericalRepresentation,
)

# The following imports are included for backwards compatibility.
# isort: split
from .base import DIFFERENTIAL_CLASSES  # noqa: F401
from .base import DUPLICATE_REPRESENTATIONS  # noqa: F401
from .base import REPRESENTATION_CLASSES  # noqa: F401
from .base import _array2string  # noqa: F401
from .base import _invalidate_reprdiff_cls_hash  # noqa: F401
from .base import get_reprdiff_cls_hash  # noqa: F401

__all__ = [
    "BaseRepresentationOrDifferential",
    "BaseRepresentation",
    "CartesianRepresentation",
    "SphericalRepresentation",
    "UnitSphericalRepresentation",
    "RadialRepresentation",
    "PhysicsSphericalRepresentation",
    "CylindricalRepresentation",
    "BaseDifferential",
    "CartesianDifferential",
    "BaseSphericalDifferential",
    "BaseSphericalCosLatDifferential",
    "SphericalDifferential",
    "SphericalCosLatDifferential",
    "UnitSphericalDifferential",
    "UnitSphericalCosLatDifferential",
    "RadialDifferential",
    "CylindricalDifferential",
    "PhysicsSphericalDifferential",
]
