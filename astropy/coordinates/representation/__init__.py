"""
In this module, we define the coordinate representation classes, which are
used to represent low-level cartesian, spherical, cylindrical, and other
coordinates.
"""

from .base import (
    BaseDifferential,
    BasePhysicalDifferential,
    BaseRepresentation,
    BaseRepresentationOrDifferential,
)
from .cartesian import CartesianDifferential, CartesianRepresentation
from .cylindrical import (
    CylindricalDifferential,
    CylindricalPhysicalDifferential,
    CylindricalRepresentation,
)
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
    PhysicsSphericalPhysicalDifferential,
    PhysicsSphericalRepresentation,
    RadialDifferential,
    RadialRepresentation,
    SphericalCosLatDifferential,
    SphericalDifferential,
    SphericalPhysicalDifferential,
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
    "BaseBodycentricRepresentation",
    "BaseDifferential",
    "BaseGeodeticRepresentation",
    "BasePhysicalDifferential",
    "BaseRepresentation",
    "BaseRepresentationOrDifferential",
    "BaseSphericalCosLatDifferential",
    "BaseSphericalDifferential",
    "CartesianDifferential",
    "CartesianRepresentation",
    "CylindricalDifferential",
    "CylindricalPhysicalDifferential",
    "CylindricalRepresentation",
    "GRS80GeodeticRepresentation",
    "PhysicsSphericalDifferential",
    "PhysicsSphericalPhysicalDifferential",
    "PhysicsSphericalRepresentation",
    "RadialDifferential",
    "RadialRepresentation",
    "SphericalCosLatDifferential",
    "SphericalDifferential",
    "SphericalPhysicalDifferential",
    "SphericalRepresentation",
    "UnitSphericalCosLatDifferential",
    "UnitSphericalDifferential",
    "UnitSphericalRepresentation",
    "WGS72GeodeticRepresentation",
    "WGS84GeodeticRepresentation",
]
