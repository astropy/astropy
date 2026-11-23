# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Cylindrical representations and differentials."""

import operator
from typing import ClassVar

import numpy as np

import astropy.units as u
from astropy.coordinates.angles import Angle
from astropy.utils import classproperty

from .base import BaseDifferential, BasePhysicalDifferential, BaseRepresentation
from .cartesian import CartesianRepresentation
from .spherical import PhysicsSphericalRepresentation, _spherical_op_funcs


class CylindricalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cylindrical coordinates.

    Parameters
    ----------
    rho : `~astropy.units.Quantity`
        The distance from the z axis to the point(s).

    phi : `~astropy.units.Quantity` or str
        The azimuth of the point(s), in angular units, which will be wrapped
        to an angle between 0 and 360 degrees. This can also be instances of
        `~astropy.coordinates.Angle`,

    z : `~astropy.units.Quantity`
        The z coordinate(s) of the point(s)

    differentials : dict, `~astropy.coordinates.CylindricalDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `~astropy.coordinates.CylindricalDifferential` instance, or a dictionary of differential
        instances with keys set to a string representation of the SI unit with
        which the differential (derivative) is taken. For example, for a
        velocity differential on a positional representation, the key would be
        ``'s'`` for seconds, indicating that the derivative is a time
        derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    attr_classes = {"rho": u.Quantity, "phi": Angle, "z": u.Quantity}

    def __init__(self, rho, phi=None, z=None, differentials=None, copy=True):
        super().__init__(rho, phi, z, copy=copy, differentials=differentials)

        if not self._rho.unit.is_equivalent(self._z.unit):
            raise u.UnitsError("rho and z should have matching physical types")

    @property
    def rho(self):
        """
        The distance of the point(s) from the z-axis.
        """
        return self._rho

    @property
    def phi(self):
        """
        The azimuth of the point(s).
        """
        return self._phi

    @property
    def z(self):
        """
        The height of the point(s).
        """
        return self._z

    @classproperty
    def _compatible_differentials(cls) -> list[type[BaseDifferential]]:
        return [*super()._compatible_differentials, CylindricalPhysicalDifferential]

    def unit_vectors(self):
        sinphi, cosphi = np.sin(self.phi), np.cos(self.phi)
        l = np.broadcast_to(1.0, self.shape)
        return {
            "rho": CartesianRepresentation(cosphi, sinphi, 0, copy=None),
            "phi": CartesianRepresentation(-sinphi, cosphi, 0, copy=None),
            "z": CartesianRepresentation(0, 0, l, unit=u.one, copy=None),
        }

    def scale_factors(self):
        rho = self.rho / u.radian
        l = np.broadcast_to(1.0 * u.one, self.shape, subok=True)
        return {"rho": l, "phi": rho, "z": l}

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to cylindrical polar
        coordinates.
        """
        rho = np.hypot(cart.x, cart.y)
        phi = np.arctan2(cart.y, cart.x)
        z = cart.z

        return cls(rho=rho, phi=phi, z=z, copy=False)

    def to_cartesian(self):
        """
        Converts cylindrical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        z = self.z

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    def _scale_operation(self, op, *args):
        if any(
            differential.base_representation is not self.__class__
            for differential in self.differentials.values()
        ):
            return super()._scale_operation(op, *args)

        phi_op, _, rho_op = _spherical_op_funcs(op, *args)
        z_op = lambda x: op(x, *args)

        # For a negative scale, the point is reflected to (phi + 180°, -z), which flips
        # the rho and phi unit vectors. For a typical differential, d_rho then scales by
        # |scale|, d_phi is unchanged, and d_z scales by the scale. For a physical
        # differential, d_phi instead scales by |scale| (because it is actually v_phi =
        # d_phi/d_t * rho). See SphericalRepresentation._scale_operation for details.
        diff_ops = (rho_op, operator.pos, z_op)
        phys_ops = (rho_op, rho_op, z_op)

        result = self.__class__(
            rho_op(self.rho), phi_op(self.phi), z_op(self.z), copy=None
        )
        for key, differential in self.differentials.items():
            new_comps = (
                op(getattr(differential, comp))
                for op, comp in zip(
                    phys_ops
                    if isinstance(differential, BasePhysicalDifferential)
                    else diff_ops,
                    differential.components,
                )
            )
            result.differentials[key] = differential.__class__(*new_comps, copy=False)
        return result

    def represent_as(self, other_class, differential_class=None):
        if isinstance(other_class, type):
            if issubclass(other_class, PhysicsSphericalRepresentation):
                diffs = self._re_represent_differentials(
                    other_class, differential_class
                )
                r = np.hypot(self.rho, self.z)
                return other_class(
                    r=r,
                    theta=np.arctan2(self.rho, self.z),
                    phi=self.phi,
                    differentials=diffs,
                )

        return super().represent_as(other_class, differential_class)


class CylindricalDifferential(BaseDifferential):
    """Differential(s) of points in cylindrical coordinates.

    Parameters
    ----------
    d_rho : `~astropy.units.Quantity` ['speed']
        The differential cylindrical radius.
    d_phi : `~astropy.units.Quantity` ['angular speed']
        The differential azimuth.
    d_z : `~astropy.units.Quantity` ['speed']
        The differential height.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = CylindricalRepresentation

    def __init__(self, d_rho, d_phi=None, d_z=None, copy=True):
        super().__init__(d_rho, d_phi, d_z, copy=copy)
        if not self._d_rho.unit.is_equivalent(self._d_z.unit):
            raise u.UnitsError("d_rho and d_z should have equivalent units.")

    def represent_as(self, other_class, base=None):
        from .spherical import (
            SphericalCosLatDifferential,
            SphericalDifferential,
            SphericalPhysicalDifferential,
            SphericalRepresentation,
        )

        if issubclass(other_class, CylindricalPhysicalDifferential):
            self._check_base(base)
            scale = base.rho / u.radian
            return other_class(self.d_rho, self.d_phi * scale, self.d_z)
        elif issubclass(other_class, SphericalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_r = self.d_rho * coslat + self.d_z * sinlat
            v_lat = -self.d_rho * sinlat + self.d_z * coslat
            scale = base_sph.distance / u.radian
            return other_class(self.d_phi, v_lat / scale, v_r)
        elif issubclass(other_class, SphericalCosLatDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_r = self.d_rho * coslat + self.d_z * sinlat
            v_lat = -self.d_rho * sinlat + self.d_z * coslat
            scale = base_sph.distance / u.radian
            d_lon_coslat = self.d_phi * coslat
            return other_class(d_lon_coslat, v_lat / scale, v_r)
        elif issubclass(other_class, SphericalPhysicalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_r = self.d_rho * coslat + self.d_z * sinlat
            v_lat = -self.d_rho * sinlat + self.d_z * coslat
            scale = base_sph.distance / u.radian
            d_lon_phys = (self.d_phi * coslat) * scale
            return other_class(d_lon_phys, v_lat, v_r)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        from .spherical import (
            SphericalCosLatDifferential,
            SphericalDifferential,
            SphericalPhysicalDifferential,
            SphericalRepresentation,
        )

        if isinstance(representation, CylindricalPhysicalDifferential):
            cls._check_base(base)
            scale = base.rho / u.radian
            return cls(representation.d_rho, representation.d_phi_phys / scale, representation.d_z)
        elif isinstance(representation, SphericalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            scale = base_sph.distance / u.radian
            v_lat = scale * representation.d_lat
            v_r = representation.d_distance
            d_rho = v_r * coslat - v_lat * sinlat
            d_z = v_r * sinlat + v_lat * coslat
            return cls(d_rho, representation.d_lon, d_z)
        elif isinstance(representation, SphericalCosLatDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            scale = base_sph.distance / u.radian
            v_lat = scale * representation.d_lat
            v_r = representation.d_distance
            d_rho = v_r * coslat - v_lat * sinlat
            d_z = v_r * sinlat + v_lat * coslat
            d_phi = representation.d_lon_coslat / coslat
            return cls(d_rho, d_phi, d_z)
        elif isinstance(representation, SphericalPhysicalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_lat = representation.d_lat_phys
            v_r = representation.d_distance
            d_rho = v_r * coslat - v_lat * sinlat
            d_z = v_r * sinlat + v_lat * coslat
            scale = base_sph.distance / u.radian
            d_phi = (representation.d_lon_phys / scale) / coslat
            return cls(d_rho, d_phi, d_z)

        return super().from_representation(representation, base)


class CylindricalPhysicalDifferential(BasePhysicalDifferential):
    """Differential(s) of points in cylindrical coordinates, in physical units.

    Components are along the local unit vectors, so that, e.g., for a
    velocity, ``d_phi_phys = rho * dphi/dt``, in the same units as ``d_rho`` and
    ``d_z`` (e.g., km/s).

    Parameters
    ----------
    d_rho, d_phi_phys, d_z : `~astropy.units.Quantity`
        The differentials along the cylindrical radius, azimuth and height
        unit vectors, all with equivalent units.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation: ClassVar = CylindricalRepresentation
    attr_classes: ClassVar = {
        "d_rho": u.Quantity,
        "d_phi_phys": u.Quantity,
        "d_z": u.Quantity,
    }

    def __init__(
        self,
        d_rho: u.Quantity,
        d_phi_phys: u.Quantity | None = None,
        d_z: u.Quantity | None = None,
        copy: bool = True,
    ) -> None:
        super().__init__(d_rho, d_phi_phys, d_z, copy=copy)

    def represent_as(self, other_class, base=None):
        from .spherical import (
            SphericalCosLatDifferential,
            SphericalDifferential,
            SphericalPhysicalDifferential,
            SphericalRepresentation,
        )

        if issubclass(other_class, CylindricalDifferential):
            self._check_base(base)
            scale = base.rho / u.radian
            return other_class(self.d_rho, self.d_phi_phys / scale, self.d_z)
        elif issubclass(other_class, SphericalPhysicalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_r = self.d_rho * coslat + self.d_z * sinlat
            v_lat = -self.d_rho * sinlat + self.d_z * coslat
            return other_class(self.d_phi_phys, v_lat, v_r)
        elif issubclass(other_class, SphericalCosLatDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_r = self.d_rho * coslat + self.d_z * sinlat
            v_lat = -self.d_rho * sinlat + self.d_z * coslat
            scale = base_sph.distance / u.radian
            d_lon_coslat = self.d_phi_phys / scale
            return other_class(d_lon_coslat, v_lat / scale, v_r)
        elif issubclass(other_class, SphericalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_r = self.d_rho * coslat + self.d_z * sinlat
            v_lat = -self.d_rho * sinlat + self.d_z * coslat
            scale = base_sph.distance / u.radian
            d_lon = (self.d_phi_phys / scale) / coslat
            return other_class(d_lon, v_lat / scale, v_r)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        from .spherical import (
            SphericalCosLatDifferential,
            SphericalDifferential,
            SphericalPhysicalDifferential,
            SphericalRepresentation,
        )

        if isinstance(representation, CylindricalDifferential):
            cls._check_base(base)
            scale = base.rho / u.radian
            return cls(representation.d_rho, representation.d_phi * scale, representation.d_z)
        elif isinstance(representation, SphericalPhysicalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            v_lat = representation.d_lat_phys
            v_r = representation.d_distance
            d_rho = v_r * coslat - v_lat * sinlat
            d_z = v_r * sinlat + v_lat * coslat
            return cls(d_rho, representation.d_lon_phys, d_z)
        elif isinstance(representation, SphericalCosLatDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            scale = base_sph.distance / u.radian
            v_lat = scale * representation.d_lat
            v_r = representation.d_distance
            d_rho = v_r * coslat - v_lat * sinlat
            d_z = v_r * sinlat + v_lat * coslat
            d_phi_phys = representation.d_lon_coslat * scale
            return cls(d_rho, d_phi_phys, d_z)
        elif isinstance(representation, SphericalDifferential):
            base_sph = base.represent_as(SphericalRepresentation)
            coslat = np.cos(base_sph.lat)
            sinlat = np.sin(base_sph.lat)
            scale = base_sph.distance / u.radian
            v_lat = scale * representation.d_lat
            v_r = representation.d_distance
            d_rho = v_r * coslat - v_lat * sinlat
            d_z = v_r * sinlat + v_lat * coslat
            d_phi_phys = (representation.d_lon * coslat) * scale
            return cls(d_rho, d_phi_phys, d_z)

        return super().from_representation(representation, base)
