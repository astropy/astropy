"""
In this module, we define the coordinate representation classes, which are
used to represent low-level cartesian, spherical, cylindrical, and other
coordinates.
"""

from .base import BaseDifferential, BaseRepresentation, BaseRepresentationOrDifferential
from .cartesian import CartesianDifferential, CartesianRepresentation
from .cylindrical import CylindricalDifferential, CylindricalRepresentation
from .geodetic import (
    BaseBodycentricRepresentation,
    BaseGeodeticRepresentation,
    GRS80GeodeticRepresentation,
    WGS72GeodeticRepresentation,
    WGS84GeodeticRepresentation,
)
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
from .base import (
    DIFFERENTIAL_CLASSES,
    DUPLICATE_REPRESENTATIONS,
    REPRESENTATION_CLASSES,
    get_reprdiff_cls_hash,
)

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
