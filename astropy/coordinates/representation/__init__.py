"""
In this module, we define the coordinate representation classes, which are
used to represent low-level cartesian, spherical, cylindrical, and other
coordinates.
"""
from .base import BaseRepresentationOrDifferential, BaseRepresentation, BaseDifferential
from .cartesian import CartesianRepresentation, CartesianDifferential
from .cylindrical import CylindricalRepresentation, CylindricalDifferential
from .geodetic import (
    BaseGeodeticRepresentation,
    BaseBodycentricRepresentation,
    WGS84GeodeticRepresentation,
    WGS72GeodeticRepresentation,
    GRS80GeodeticRepresentation,
)
from .spherical import (
    SphericalRepresentation,
    UnitSphericalRepresentation,
    RadialRepresentation,
    PhysicsSphericalRepresentation,
    BaseSphericalDifferential,
    BaseSphericalCosLatDifferential,
    SphericalDifferential,
    SphericalCosLatDifferential,
    UnitSphericalDifferential,
    UnitSphericalCosLatDifferential,
    RadialDifferential,
    PhysicsSphericalDifferential,
)

# The following imports are included for backwards compatibility.
# isort: split
from .base import DIFFERENTIAL_CLASSES
from .base import DUPLICATE_REPRESENTATIONS
from .base import REPRESENTATION_CLASSES
from .base import get_reprdiff_cls_hash


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
    "BaseGeodeticRepresentation",
    "BaseBodycentricRepresentation",
    "WGS84GeodeticRepresentation",
    "WGS72GeodeticRepresentation",
    "GRS80GeodeticRepresentation",
]
