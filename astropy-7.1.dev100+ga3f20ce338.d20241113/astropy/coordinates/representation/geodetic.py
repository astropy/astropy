# Licensed under a 3-clause BSD style license - see LICENSE.rst

import erfa
import numpy as np

from astropy import units as u
from astropy.coordinates.angles import Latitude, Longitude
from astropy.utils.decorators import format_doc

from .base import BaseRepresentation
from .cartesian import CartesianRepresentation

ELLIPSOIDS = {}
"""Available ellipsoids (defined in erfam.h, with numbers exposed in erfa)."""
# Note: they get filled by the creation of the geodetic classes.


geodetic_base_doc = """{__doc__}

    Parameters
    ----------
    lon, lat : angle-like
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle` and either
        `~astropy.coordinates.Longitude` not `~astropy.coordinates.Latitude`,
        depending on the parameter.

    height : `~astropy.units.Quantity` ['length']
        The height to the point(s).

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
"""


@format_doc(geodetic_base_doc)
class BaseGeodeticRepresentation(BaseRepresentation):
    """
    Base class for geodetic representations.

    Subclasses need to set attributes ``_equatorial_radius`` and ``_flattening``
    to quantities holding correct values (with units of length and dimensionless,
    respectively), or alternatively an ``_ellipsoid`` attribute to the relevant ERFA
    index (as passed in to `erfa.eform`). The geodetic latitude is defined by the
    angle between the vertical to the surface at a specific point of the spheroid and
    its projection onto the equatorial plane.
    """

    attr_classes = {"lon": Longitude, "lat": Latitude, "height": u.Quantity}

    def __init_subclass__(cls, **kwargs):
        if "_ellipsoid" in cls.__dict__:
            equatorial_radius, flattening = erfa.eform(getattr(erfa, cls._ellipsoid))
            cls._equatorial_radius = equatorial_radius * u.m
            cls._flattening = flattening * u.dimensionless_unscaled
            ELLIPSOIDS[cls._ellipsoid] = cls
        elif (
            "_equatorial_radius" not in cls.__dict__
            or "_flattening" not in cls.__dict__
        ):
            raise AttributeError(
                f"{cls.__name__} requires '_ellipsoid' or '_equatorial_radius' and '_flattening'."
            )
        super().__init_subclass__(**kwargs)

    def __init__(self, lon, lat=None, height=None, copy=True):
        if height is None and not isinstance(lon, self.__class__):
            height = 0 << u.m

        super().__init__(lon, lat, height, copy=copy)
        if not self.height.unit.is_equivalent(u.m):
            raise u.UnitTypeError(
                f"{self.__class__.__name__} requires height with units of length."
            )

    def to_cartesian(self):
        """
        Converts geodetic coordinates to 3D rectangular (geocentric)
        cartesian coordinates.
        """
        xyz = erfa.gd2gce(
            self._equatorial_radius,
            self._flattening,
            self.lon,
            self.lat,
            self.height,
        )
        return CartesianRepresentation(xyz, xyz_axis=-1, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates (assumed geocentric) to
        geodetic coordinates.
        """
        # Compute geodetic/planetodetic angles
        lon, lat, height = erfa.gc2gde(
            cls._equatorial_radius, cls._flattening, cart.get_xyz(xyz_axis=-1)
        )
        return cls(lon, lat, height, copy=False)


@format_doc(geodetic_base_doc)
class BaseBodycentricRepresentation(BaseRepresentation):
    """Representation of points in bodycentric 3D coordinates.

    Subclasses need to set attributes ``_equatorial_radius`` and ``_flattening``
    to quantities holding correct values (with units of length and dimensionless,
    respectively). the bodycentric latitude and longitude are spherical latitude
    and longitude relative to the barycenter of the body.
    """

    attr_classes = {"lon": Longitude, "lat": Latitude, "height": u.Quantity}

    def __init_subclass__(cls, **kwargs):
        if (
            "_equatorial_radius" not in cls.__dict__
            or "_flattening" not in cls.__dict__
        ):
            raise AttributeError(
                f"{cls.__name__} requires '_equatorial_radius' and '_flattening'."
            )
        super().__init_subclass__(**kwargs)

    def __init__(self, lon, lat=None, height=None, copy=True):
        if height is None and not isinstance(lon, self.__class__):
            height = 0 << u.m

        super().__init__(lon, lat, height, copy=copy)
        if not self.height.unit.is_equivalent(u.m):
            raise u.UnitTypeError(
                f"{self.__class__.__name__} requires height with units of length."
            )

    def to_cartesian(self):
        """
        Converts bodycentric coordinates to 3D rectangular (geocentric)
        cartesian coordinates.
        """
        coslat = np.cos(self.lat)
        sinlat = np.sin(self.lat)
        coslon = np.cos(self.lon)
        sinlon = np.sin(self.lon)
        r = (
            self._equatorial_radius * np.hypot(coslat, (1 - self._flattening) * sinlat)
            + self.height
        )
        x = r * coslon * coslat
        y = r * sinlon * coslat
        z = r * sinlat
        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates (assumed geocentric) to
        bodycentric coordinates.
        """
        # Compute bodycentric latitude
        p = np.hypot(cart.x, cart.y)
        d = np.hypot(p, cart.z)
        lat = np.arctan2(cart.z, p)
        p_spheroid = cls._equatorial_radius * np.cos(lat)
        z_spheroid = cls._equatorial_radius * (1 - cls._flattening) * np.sin(lat)
        r_spheroid = np.hypot(p_spheroid, z_spheroid)
        height = d - r_spheroid
        lon = np.arctan2(cart.y, cart.x)
        return cls(lon, lat, height, copy=False)


@format_doc(geodetic_base_doc)
class WGS84GeodeticRepresentation(BaseGeodeticRepresentation):
    """Representation of points in WGS84 3D geodetic coordinates."""

    _ellipsoid = "WGS84"


@format_doc(geodetic_base_doc)
class WGS72GeodeticRepresentation(BaseGeodeticRepresentation):
    """Representation of points in WGS72 3D geodetic coordinates."""

    _ellipsoid = "WGS72"


@format_doc(geodetic_base_doc)
class GRS80GeodeticRepresentation(BaseGeodeticRepresentation):
    """Representation of points in GRS80 3D geodetic coordinates."""

    _ellipsoid = "GRS80"
