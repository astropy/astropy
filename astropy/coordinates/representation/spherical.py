# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Spherical representations and differentials."""
import operator

import numpy as np
from erfa import ufunc as erfa_ufunc

import astropy.units as u
from astropy.coordinates.angles import Angle, Latitude, Longitude
from astropy.coordinates.distances import Distance
from astropy.coordinates.matrix_utilities import is_O3
from astropy.utils import classproperty

from .base import BaseDifferential, BaseRepresentation
from .cartesian import CartesianRepresentation


class UnitSphericalRepresentation(BaseRepresentation):
    """
    Representation of points on a unit sphere.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity` ['angle'] or str
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    differentials : dict, `~astropy.coordinates.BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `~astropy.coordinates.BaseDifferential`
        instance (see `._compatible_differentials` for valid types), or a
        dictionary of of differential instances with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    attr_classes = {"lon": Longitude, "lat": Latitude}

    @classproperty
    def _dimensional_representation(cls):
        return SphericalRepresentation

    def __init__(self, lon, lat=None, differentials=None, copy=True):
        super().__init__(lon, lat, differentials=differentials, copy=copy)

    @classproperty
    def _compatible_differentials(cls):
        return [
            UnitSphericalDifferential,
            UnitSphericalCosLatDifferential,
            SphericalDifferential,
            SphericalCosLatDifferential,
            RadialDifferential,
        ]

    # Could let the metaclass define these automatically, but good to have
    # a bit clearer docstrings.
    @property
    def lon(self):
        """
        The longitude of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The latitude of the point(s).
        """
        return self._lat

    def unit_vectors(self):
        sinlon, coslon = np.sin(self.lon), np.cos(self.lon)
        sinlat, coslat = np.sin(self.lat), np.cos(self.lat)
        return {
            "lon": CartesianRepresentation(-sinlon, coslon, 0.0, copy=False),
            "lat": CartesianRepresentation(
                -sinlat * coslon, -sinlat * sinlon, coslat, copy=False
            ),
        }

    def scale_factors(self, omit_coslat=False):
        sf_lat = np.broadcast_to(1.0 / u.radian, self.shape, subok=True)
        sf_lon = sf_lat if omit_coslat else np.cos(self.lat) / u.radian
        return {"lon": sf_lon, "lat": sf_lat}

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        # erfa s2c: Convert [unit]spherical coordinates to Cartesian.
        p = erfa_ufunc.s2c(self.lon, self.lat)
        return CartesianRepresentation(p, xyz_axis=-1, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """
        p = cart.get_xyz(xyz_axis=-1)
        # erfa c2s: P-vector to [unit]spherical coordinates.
        return cls(*erfa_ufunc.c2s(p), copy=False)

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation
        # TODO! for differential_class. This cannot (currently) be implemented
        # like in the other Representations since `_re_represent_differentials`
        # keeps differentials' unit keys, but this can result in a mismatch
        # between the UnitSpherical expected key (e.g. "s") and that expected
        # in the other class (here "s / m"). For more info, see PR #11467
        if isinstance(other_class, type) and not differential_class:
            if issubclass(other_class, PhysicsSphericalRepresentation):
                return other_class(
                    phi=self.lon, theta=90 * u.deg - self.lat, r=1.0, copy=False
                )
            elif issubclass(other_class, SphericalRepresentation):
                return other_class(lon=self.lon, lat=self.lat, distance=1.0, copy=False)

        return super().represent_as(other_class, differential_class)

    def transform(self, matrix):
        r"""Transform the unit-spherical coordinates using a 3x3 matrix.

        This returns a new representation and does not modify the original one.
        Any differentials attached to this representation will also be
        transformed.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 matrix, such as a rotation matrix (or a stack of matrices).

        Returns
        -------
        `~astropy.coordinates.UnitSphericalRepresentation` or `~astropy.coordinates.SphericalRepresentation`
            If ``matrix`` is O(3) -- :math:`M \dot M^T = I` -- like a rotation,
            then the result is a `~astropy.coordinates.UnitSphericalRepresentation`.
            All other matrices will change the distance, so the dimensional
            representation is used instead.

        """
        # the transformation matrix does not need to be a rotation matrix,
        # so the unit-distance is not guaranteed. For speed, we check if the
        # matrix is in O(3) and preserves lengths.
        if np.all(is_O3(matrix)):  # remain in unit-rep
            xyz = erfa_ufunc.s2c(self.lon, self.lat)
            p = erfa_ufunc.rxp(matrix, xyz)
            lon, lat = erfa_ufunc.c2s(p)
            rep = self.__class__(lon=lon, lat=lat)
            # handle differentials
            new_diffs = {
                k: d.transform(matrix, self, rep) for k, d in self.differentials.items()
            }
            rep = rep.with_differentials(new_diffs)

        else:  # switch to dimensional representation
            rep = self._dimensional_representation(
                lon=self.lon, lat=self.lat, distance=1, differentials=self.differentials
            ).transform(matrix)

        return rep

    def _scale_operation(self, op, *args):
        return self._dimensional_representation(
            lon=self.lon, lat=self.lat, distance=1.0, differentials=self.differentials
        )._scale_operation(op, *args)

    def __neg__(self):
        if any(
            differential.base_representation is not self.__class__
            for differential in self.differentials.values()
        ):
            return super().__neg__()

        result = self.__class__(self.lon + 180.0 * u.deg, -self.lat, copy=False)
        for key, differential in self.differentials.items():
            new_comps = (
                op(getattr(differential, comp))
                for op, comp in zip(
                    (operator.pos, operator.neg), differential.components
                )
            )
            result.differentials[key] = differential.__class__(*new_comps, copy=False)
        return result

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units, which is
        always unity for vectors on the unit sphere.

        Returns
        -------
        norm : `~astropy.units.Quantity` ['dimensionless']
            Dimensionless ones, with the same shape as the representation.
        """
        return u.Quantity(np.ones(self.shape), u.dimensionless_unscaled, copy=False)

    def _combine_operation(self, op, other, reverse=False):
        self._raise_if_has_differentials(op.__name__)

        result = self.to_cartesian()._combine_operation(op, other, reverse)
        if result is NotImplemented:
            return NotImplemented
        else:
            return self._dimensional_representation.from_cartesian(result)

    def mean(self, *args, **kwargs):
        """Vector mean.

        The representation is converted to cartesian, the means of the x, y,
        and z components are calculated, and the result is converted to a
        `~astropy.coordinates.SphericalRepresentation`.

        Refer to `~numpy.mean` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials("mean")
        return self._dimensional_representation.from_cartesian(
            self.to_cartesian().mean(*args, **kwargs)
        )

    def sum(self, *args, **kwargs):
        """Vector sum.

        The representation is converted to cartesian, the sums of the x, y,
        and z components are calculated, and the result is converted to a
        `~astropy.coordinates.SphericalRepresentation`.

        Refer to `~numpy.sum` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials("sum")
        return self._dimensional_representation.from_cartesian(
            self.to_cartesian().sum(*args, **kwargs)
        )

    def cross(self, other):
        """Cross product of two representations.

        The calculation is done by converting both ``self`` and ``other``
        to `~astropy.coordinates.CartesianRepresentation`, and converting the
        result back to `~astropy.coordinates.SphericalRepresentation`.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            The representation to take the cross product with.

        Returns
        -------
        cross_product : `~astropy.coordinates.SphericalRepresentation`
            With vectors perpendicular to both ``self`` and ``other``.
        """
        self._raise_if_has_differentials("cross")
        return self._dimensional_representation.from_cartesian(
            self.to_cartesian().cross(other)
        )


class RadialRepresentation(BaseRepresentation):
    """
    Representation of the distance of points from the origin.

    Note that this is mostly intended as an internal helper representation.
    It can do little else but being used as a scale in multiplication.

    Parameters
    ----------
    distance : `~astropy.units.Quantity` ['length']
        The distance of the point(s) from the origin.

    differentials : dict, `~astropy.coordinates.BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `~astropy.coordinates.BaseDifferential`
        instance (see `._compatible_differentials` for valid types), or a
        dictionary of of differential instances with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    attr_classes = {"distance": u.Quantity}

    def __init__(self, distance, differentials=None, copy=True):
        super().__init__(distance, differentials=differentials, copy=copy)

    @property
    def distance(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    def unit_vectors(self):
        """Cartesian unit vectors are undefined for radial representation."""
        raise NotImplementedError(
            f"Cartesian unit vectors are undefined for {self.__class__} instances"
        )

    def scale_factors(self):
        l = np.broadcast_to(1.0 * u.one, self.shape, subok=True)
        return {"distance": l}

    def to_cartesian(self):
        """Cannot convert radial representation to cartesian."""
        raise NotImplementedError(
            f"cannot convert {self.__class__} instance to cartesian."
        )

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to radial coordinate.
        """
        return cls(distance=cart.norm(), copy=False)

    def __mul__(self, other):
        if isinstance(other, BaseRepresentation):
            return self.distance * other
        else:
            return super().__mul__(other)

    def norm(self):
        """Vector norm.

        Just the distance itself.

        Returns
        -------
        norm : `~astropy.units.Quantity` ['dimensionless']
            Dimensionless ones, with the same shape as the representation.
        """
        return self.distance

    def _combine_operation(self, op, other, reverse=False):
        return NotImplemented

    def transform(self, matrix):
        """Radial representations cannot be transformed by a Cartesian matrix.

        Parameters
        ----------
        matrix : array-like
            The transformation matrix in a Cartesian basis.
            Must be a multiplication: a diagonal matrix with identical elements.
            Must have shape (..., 3, 3), where the last 2 indices are for the
            matrix on each other axis. Make sure that the matrix shape is
            compatible with the shape of this representation.

        Raises
        ------
        ValueError
            If the matrix is not a multiplication.

        """
        scl = matrix[..., 0, 0]
        # check that the matrix is a scaled identity matrix on the last 2 axes.
        if np.any(matrix != scl[..., np.newaxis, np.newaxis] * np.identity(3)):
            raise ValueError(
                "Radial representations can only be "
                "transformed by a scaled identity matrix"
            )

        return self * scl


def _spherical_op_funcs(op, *args):
    """For given operator, return functions that adjust lon, lat, distance."""
    if op is operator.neg:
        return lambda x: x + 180 * u.deg, operator.neg, operator.pos

    try:
        scale_sign = np.sign(args[0])
    except Exception:
        # This should always work, even if perhaps we get a negative distance.
        return operator.pos, operator.pos, lambda x: op(x, *args)

    scale = abs(args[0])
    return (
        lambda x: x + 180 * u.deg * np.signbit(scale_sign),
        lambda x: x * scale_sign,
        lambda x: op(x, scale),
    )


class SphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity` ['angle']
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    distance : `~astropy.units.Quantity` ['length']
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, otherwise
        it is passed to the :class:`~astropy.units.Quantity` class.

    differentials : dict, `~astropy.coordinates.BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `~astropy.coordinates.BaseDifferential`
        instance (see `._compatible_differentials` for valid types), or a
        dictionary of of differential instances with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    attr_classes = {"lon": Longitude, "lat": Latitude, "distance": u.Quantity}
    _unit_representation = UnitSphericalRepresentation

    def __init__(self, lon, lat=None, distance=None, differentials=None, copy=True):
        super().__init__(lon, lat, distance, copy=copy, differentials=differentials)
        if (
            not isinstance(self._distance, Distance)
            and self._distance.unit.physical_type == "length"
        ):
            try:
                self._distance = Distance(self._distance, copy=False)
            except ValueError as e:
                if e.args[0].startswith("distance must be >= 0"):
                    raise ValueError(
                        "Distance must be >= 0. To allow negative distance values, you"
                        " must explicitly pass in a `Distance` object with the "
                        "argument 'allow_negative=True'."
                    ) from e
                else:
                    raise

    @classproperty
    def _compatible_differentials(cls):
        return [
            UnitSphericalDifferential,
            UnitSphericalCosLatDifferential,
            SphericalDifferential,
            SphericalCosLatDifferential,
            RadialDifferential,
        ]

    @property
    def lon(self):
        """
        The longitude of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The latitude of the point(s).
        """
        return self._lat

    @property
    def distance(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    def unit_vectors(self):
        sinlon, coslon = np.sin(self.lon), np.cos(self.lon)
        sinlat, coslat = np.sin(self.lat), np.cos(self.lat)
        return {
            "lon": CartesianRepresentation(-sinlon, coslon, 0.0, copy=False),
            "lat": CartesianRepresentation(
                -sinlat * coslon, -sinlat * sinlon, coslat, copy=False
            ),
            "distance": CartesianRepresentation(
                coslat * coslon, coslat * sinlon, sinlat, copy=False
            ),
        }

    def scale_factors(self, omit_coslat=False):
        sf_lat = self.distance / u.radian
        sf_lon = sf_lat if omit_coslat else sf_lat * np.cos(self.lat)
        sf_distance = np.broadcast_to(1.0 * u.one, self.shape, subok=True)
        return {"lon": sf_lon, "lat": sf_lat, "distance": sf_distance}

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation

        if isinstance(other_class, type):
            if issubclass(other_class, PhysicsSphericalRepresentation):
                diffs = self._re_represent_differentials(
                    other_class, differential_class
                )
                return other_class(
                    phi=self.lon,
                    theta=90 * u.deg - self.lat,
                    r=self.distance,
                    differentials=diffs,
                    copy=False,
                )

            elif issubclass(other_class, UnitSphericalRepresentation):
                diffs = self._re_represent_differentials(
                    other_class, differential_class
                )
                return other_class(
                    lon=self.lon, lat=self.lat, differentials=diffs, copy=False
                )

        return super().represent_as(other_class, differential_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        # We need to convert Distance to Quantity to allow negative values.
        if isinstance(self.distance, Distance):
            d = self.distance.view(u.Quantity)
        else:
            d = self.distance

        # erfa s2p: Convert spherical polar coordinates to p-vector.
        p = erfa_ufunc.s2p(self.lon, self.lat, d)

        return CartesianRepresentation(p, xyz_axis=-1, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """
        p = cart.get_xyz(xyz_axis=-1)
        # erfa p2s: P-vector to spherical polar coordinates.
        return cls(*erfa_ufunc.p2s(p), copy=False)

    def transform(self, matrix):
        """Transform the spherical coordinates using a 3x3 matrix.

        This returns a new representation and does not modify the original one.
        Any differentials attached to this representation will also be
        transformed.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 matrix, such as a rotation matrix (or a stack of matrices).

        """
        xyz = erfa_ufunc.s2c(self.lon, self.lat)
        p = erfa_ufunc.rxp(matrix, xyz)
        lon, lat, ur = erfa_ufunc.p2s(p)
        rep = self.__class__(lon=lon, lat=lat, distance=self.distance * ur)

        # handle differentials
        new_diffs = {
            k: d.transform(matrix, self, rep) for k, d in self.differentials.items()
        }
        return rep.with_differentials(new_diffs)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.  For
        spherical coordinates, this is just the absolute value of the distance.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return np.abs(self.distance)

    def _scale_operation(self, op, *args):
        # TODO: expand special-casing to UnitSpherical and RadialDifferential.
        if any(
            differential.base_representation is not self.__class__
            for differential in self.differentials.values()
        ):
            return super()._scale_operation(op, *args)

        lon_op, lat_op, distance_op = _spherical_op_funcs(op, *args)

        result = self.__class__(
            lon_op(self.lon), lat_op(self.lat), distance_op(self.distance), copy=False
        )
        for key, differential in self.differentials.items():
            new_comps = (
                op(getattr(differential, comp))
                for op, comp in zip(
                    (operator.pos, lat_op, distance_op), differential.components
                )
            )
            result.differentials[key] = differential.__class__(*new_comps, copy=False)
        return result


class PhysicsSphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates (using the physics
    convention of using ``phi`` and ``theta`` for azimuth and inclination
    from the pole).

    Parameters
    ----------
    phi, theta : `~astropy.units.Quantity` or str
        The azimuth and inclination of the point(s), in angular units. The
        inclination should be between 0 and 180 degrees, and the azimuth will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`.  If ``copy`` is False, `phi`
        will be changed inplace if it is not between 0 and 360 degrees.

    r : `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, otherwise
        it is passed to the :class:`~astropy.units.Quantity` class.

    differentials : dict, `~astropy.coordinates.PhysicsSphericalDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `~astropy.coordinates.PhysicsSphericalDifferential` instance, or a dictionary of of
        differential instances with keys set to a string representation of the
        SI unit with which the differential (derivative) is taken. For example,
        for a velocity differential on a positional representation, the key
        would be ``'s'`` for seconds, indicating that the derivative is a time
        derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    attr_classes = {"phi": Angle, "theta": Angle, "r": u.Quantity}

    def __init__(self, phi, theta=None, r=None, differentials=None, copy=True):
        super().__init__(phi, theta, r, copy=copy, differentials=differentials)

        # Wrap/validate phi/theta
        # Note that _phi already holds our own copy if copy=True.
        self._phi.wrap_at(360 * u.deg, inplace=True)

        if np.any(self._theta < 0.0 * u.deg) or np.any(self._theta > 180.0 * u.deg):
            raise ValueError(
                "Inclination angle(s) must be within 0 deg <= angle <= 180 deg, "
                f"got {theta.to(u.degree)}"
            )

        if self._r.unit.physical_type == "length":
            self._r = self._r.view(Distance)

    @property
    def phi(self):
        """
        The azimuth of the point(s).
        """
        return self._phi

    @property
    def theta(self):
        """
        The elevation of the point(s).
        """
        return self._theta

    @property
    def r(self):
        """
        The distance from the origin to the point(s).
        """
        return self._r

    def unit_vectors(self):
        sinphi, cosphi = np.sin(self.phi), np.cos(self.phi)
        sintheta, costheta = np.sin(self.theta), np.cos(self.theta)
        return {
            "phi": CartesianRepresentation(-sinphi, cosphi, 0.0, copy=False),
            "theta": CartesianRepresentation(
                costheta * cosphi, costheta * sinphi, -sintheta, copy=False
            ),
            "r": CartesianRepresentation(
                sintheta * cosphi, sintheta * sinphi, costheta, copy=False
            ),
        }

    def scale_factors(self):
        r = self.r / u.radian
        sintheta = np.sin(self.theta)
        l = np.broadcast_to(1.0 * u.one, self.shape, subok=True)
        return {"phi": r * sintheta, "theta": r, "r": l}

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation

        if isinstance(other_class, type):
            if issubclass(other_class, SphericalRepresentation):
                diffs = self._re_represent_differentials(
                    other_class, differential_class
                )
                return other_class(
                    lon=self.phi,
                    lat=90 * u.deg - self.theta,
                    distance=self.r,
                    differentials=diffs,
                    copy=False,
                )
            elif issubclass(other_class, UnitSphericalRepresentation):
                diffs = self._re_represent_differentials(
                    other_class, differential_class
                )
                return other_class(
                    lon=self.phi,
                    lat=90 * u.deg - self.theta,
                    differentials=diffs,
                    copy=False,
                )

        return super().represent_as(other_class, differential_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        # We need to convert Distance to Quantity to allow negative values.
        if isinstance(self.r, Distance):
            d = self.r.view(u.Quantity)
        else:
            d = self.r

        x = d * np.sin(self.theta) * np.cos(self.phi)
        y = d * np.sin(self.theta) * np.sin(self.phi)
        z = d * np.cos(self.theta)

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """
        s = np.hypot(cart.x, cart.y)
        r = np.hypot(s, cart.z)

        phi = np.arctan2(cart.y, cart.x)
        theta = np.arctan2(s, cart.z)

        return cls(phi=phi, theta=theta, r=r, copy=False)

    def transform(self, matrix):
        """Transform the spherical coordinates using a 3x3 matrix.

        This returns a new representation and does not modify the original one.
        Any differentials attached to this representation will also be
        transformed.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 matrix, such as a rotation matrix (or a stack of matrices).

        """
        # apply transformation in unit-spherical coordinates
        xyz = erfa_ufunc.s2c(self.phi, 90 * u.deg - self.theta)
        p = erfa_ufunc.rxp(matrix, xyz)
        lon, lat, ur = erfa_ufunc.p2s(p)  # `ur` is transformed unit-`r`
        # create transformed physics-spherical representation,
        # reapplying the distance scaling
        rep = self.__class__(phi=lon, theta=90 * u.deg - lat, r=self.r * ur)

        new_diffs = {
            k: d.transform(matrix, self, rep) for k, d in self.differentials.items()
        }
        return rep.with_differentials(new_diffs)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.  For
        spherical coordinates, this is just the absolute value of the radius.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return np.abs(self.r)

    def _scale_operation(self, op, *args):
        if any(
            differential.base_representation is not self.__class__
            for differential in self.differentials.values()
        ):
            return super()._scale_operation(op, *args)

        phi_op, adjust_theta_sign, r_op = _spherical_op_funcs(op, *args)
        # Also run phi_op on theta to ensure theta remains between 0 and 180:
        # any time the scale is negative, we do -theta + 180 degrees.
        result = self.__class__(
            phi_op(self.phi),
            phi_op(adjust_theta_sign(self.theta)),
            r_op(self.r),
            copy=False,
        )
        for key, differential in self.differentials.items():
            new_comps = (
                op(getattr(differential, comp))
                for op, comp in zip(
                    (operator.pos, adjust_theta_sign, r_op), differential.components
                )
            )
            result.differentials[key] = differential.__class__(*new_comps, copy=False)
        return result


class BaseSphericalDifferential(BaseDifferential):
    def _d_lon_coslat(self, base):
        """Convert longitude differential d_lon to d_lon_coslat.

        Parameters
        ----------
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        self._check_base(base)
        return self.d_lon * np.cos(base.lat)

    @classmethod
    def _get_d_lon(cls, d_lon_coslat, base):
        """Convert longitude differential d_lon_coslat to d_lon.

        Parameters
        ----------
        d_lon_coslat : `~astropy.units.Quantity`
            Longitude differential that includes ``cos(lat)``.
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        cls._check_base(base)
        return d_lon_coslat / np.cos(base.lat)

    def _combine_operation(self, op, other, reverse=False):
        """Combine two differentials, or a differential with a representation.

        If ``other`` is of the same differential type as ``self``, the
        components will simply be combined.  If both are different parts of
        a `~astropy.coordinates.SphericalDifferential` (e.g., a
        `~astropy.coordinates.UnitSphericalDifferential` and a
        `~astropy.coordinates.RadialDifferential`), they will combined
        appropriately.

        If ``other`` is a representation, it will be used as a base for which
        to evaluate the differential, and the result is a new representation.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            The other differential or representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        if (
            isinstance(other, BaseSphericalDifferential)
            and not isinstance(self, type(other))
            or isinstance(other, RadialDifferential)
        ):
            all_components = set(self.components) | set(other.components)
            first, second = (self, other) if not reverse else (other, self)
            result_args = {
                c: op(getattr(first, c, 0.0), getattr(second, c, 0.0))
                for c in all_components
            }
            return SphericalDifferential(**result_args)

        return super()._combine_operation(op, other, reverse)


class UnitSphericalDifferential(BaseSphericalDifferential):
    """Differential(s) of points on a unit sphere.

    Parameters
    ----------
    d_lon, d_lat : `~astropy.units.Quantity`
        The longitude and latitude of the differentials.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = UnitSphericalRepresentation

    @classproperty
    def _dimensional_differential(cls):
        return SphericalDifferential

    def __init__(self, d_lon, d_lat=None, copy=True):
        super().__init__(d_lon, d_lat, copy=copy)
        if not self._d_lon.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError("d_lon and d_lat should have equivalent units.")

    @classmethod
    def from_cartesian(cls, other, base):
        # Go via the dimensional equivalent, so that the longitude and latitude
        # differentials correctly take into account the norm of the base.
        dimensional = cls._dimensional_differential.from_cartesian(other, base)
        return dimensional.represent_as(cls)

    def to_cartesian(self, base):
        if isinstance(base, SphericalRepresentation):
            scale = base.distance
        elif isinstance(base, PhysicsSphericalRepresentation):
            scale = base.r
        else:
            return super().to_cartesian(base)

        base = base.represent_as(UnitSphericalRepresentation)
        return scale * super().to_cartesian(base)

    def represent_as(self, other_class, base=None):
        # Only have enough information to represent other unit-spherical.
        if issubclass(other_class, UnitSphericalCosLatDifferential):
            return other_class(self._d_lon_coslat(base), self.d_lat)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude.
        if isinstance(representation, SphericalDifferential):
            return cls(representation.d_lon, representation.d_lat)
        elif isinstance(
            representation,
            (SphericalCosLatDifferential, UnitSphericalCosLatDifferential),
        ):
            d_lon = cls._get_d_lon(representation.d_lon_coslat, base)
            return cls(d_lon, representation.d_lat)
        elif isinstance(representation, PhysicsSphericalDifferential):
            return cls(representation.d_phi, -representation.d_theta)

        return super().from_representation(representation, base)

    def transform(self, matrix, base, transformed_base):
        """Transform differential using a 3x3 matrix in a Cartesian basis.

        This returns a new differential and does not modify the original one.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 (or stack thereof) matrix, such as a rotation matrix.
        base : instance of ``cls.base_representation``
            Base relative to which the differentials are defined.  If the other
            class is a differential representation, the base will be converted
            to its ``base_representation``.
        transformed_base : instance of ``cls.base_representation``
            Base relative to which the transformed differentials are defined.
            If the other class is a differential representation, the base will
            be converted to its ``base_representation``.
        """
        # the transformation matrix does not need to be a rotation matrix,
        # so the unit-distance is not guaranteed. For speed, we check if the
        # matrix is in O(3) and preserves lengths.
        if np.all(is_O3(matrix)):  # remain in unit-rep
            # TODO! implement without Cartesian intermediate step.
            # some of this can be moved to the parent class.
            diff = super().transform(matrix, base, transformed_base)

        else:  # switch to dimensional representation
            du = self.d_lon.unit / base.lon.unit  # derivative unit
            diff = self._dimensional_differential(
                d_lon=self.d_lon, d_lat=self.d_lat, d_distance=0 * du
            ).transform(matrix, base, transformed_base)

        return diff

    def _scale_operation(self, op, *args, scaled_base=False):
        if scaled_base:
            return self.copy()
        else:
            return super()._scale_operation(op, *args)


class SphericalDifferential(BaseSphericalDifferential):
    """Differential(s) of points in 3D spherical coordinates.

    Parameters
    ----------
    d_lon, d_lat : `~astropy.units.Quantity`
        The differential longitude and latitude.
    d_distance : `~astropy.units.Quantity`
        The differential distance.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = SphericalRepresentation
    _unit_differential = UnitSphericalDifferential

    def __init__(self, d_lon, d_lat=None, d_distance=None, copy=True):
        super().__init__(d_lon, d_lat, d_distance, copy=copy)
        if not self._d_lon.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError("d_lon and d_lat should have equivalent units.")

    def represent_as(self, other_class, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude.
        if issubclass(other_class, UnitSphericalDifferential):
            return other_class(self.d_lon, self.d_lat)
        elif issubclass(other_class, RadialDifferential):
            return other_class(self.d_distance)
        elif issubclass(other_class, SphericalCosLatDifferential):
            return other_class(self._d_lon_coslat(base), self.d_lat, self.d_distance)
        elif issubclass(other_class, UnitSphericalCosLatDifferential):
            return other_class(self._d_lon_coslat(base), self.d_lat)
        elif issubclass(other_class, PhysicsSphericalDifferential):
            return other_class(self.d_lon, -self.d_lat, self.d_distance)
        else:
            return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # Other spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude.
        if isinstance(representation, SphericalCosLatDifferential):
            d_lon = cls._get_d_lon(representation.d_lon_coslat, base)
            return cls(d_lon, representation.d_lat, representation.d_distance)
        elif isinstance(representation, PhysicsSphericalDifferential):
            return cls(
                representation.d_phi, -representation.d_theta, representation.d_r
            )

        return super().from_representation(representation, base)

    def _scale_operation(self, op, *args, scaled_base=False):
        if scaled_base:
            return self.__class__(self.d_lon, self.d_lat, op(self.d_distance, *args))
        else:
            return super()._scale_operation(op, *args)


class BaseSphericalCosLatDifferential(BaseDifferential):
    """Differentials from points on a spherical base representation.

    With cos(lat) assumed to be included in the longitude differential.
    """

    @classmethod
    def _get_base_vectors(cls, base):
        """Get unit vectors and scale factors from (unit)spherical base.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            The points for which the unit vectors and scale factors should be
            retrieved.

        Returns
        -------
        unit_vectors : dict of `~astropy.coordinates.CartesianRepresentation`
            In the directions of the coordinates of base.
        scale_factors : dict of `~astropy.units.Quantity`
            Scale factors for each of the coordinates.  The scale factor for
            longitude does not include the cos(lat) factor.

        Raises
        ------
        TypeError : if the base is not of the correct type
        """
        cls._check_base(base)
        return base.unit_vectors(), base.scale_factors(omit_coslat=True)

    def _d_lon(self, base):
        """Convert longitude differential with cos(lat) to one without.

        Parameters
        ----------
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        self._check_base(base)
        return self.d_lon_coslat / np.cos(base.lat)

    @classmethod
    def _get_d_lon_coslat(cls, d_lon, base):
        """Convert longitude differential d_lon to d_lon_coslat.

        Parameters
        ----------
        d_lon : `~astropy.units.Quantity`
            Value of the longitude differential without ``cos(lat)``.
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        cls._check_base(base)
        return d_lon * np.cos(base.lat)

    def _combine_operation(self, op, other, reverse=False):
        """Combine two differentials, or a differential with a representation.

        If ``other`` is of the same differential type as ``self``, the
        components will simply be combined.  If both are different parts of
        a `~astropy.coordinates.SphericalDifferential` (e.g., a
        `~astropy.coordinates.UnitSphericalDifferential` and a
        `~astropy.coordinates.RadialDifferential`), they will combined
        appropriately.

        If ``other`` is a representation, it will be used as a base for which
        to evaluate the differential, and the result is a new representation.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            The other differential or representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        if (
            isinstance(other, BaseSphericalCosLatDifferential)
            and not isinstance(self, type(other))
            or isinstance(other, RadialDifferential)
        ):
            all_components = set(self.components) | set(other.components)
            first, second = (self, other) if not reverse else (other, self)
            result_args = {
                c: op(getattr(first, c, 0.0), getattr(second, c, 0.0))
                for c in all_components
            }
            return SphericalCosLatDifferential(**result_args)

        return super()._combine_operation(op, other, reverse)


class UnitSphericalCosLatDifferential(BaseSphericalCosLatDifferential):
    """Differential(s) of points on a unit sphere.

    Parameters
    ----------
    d_lon_coslat, d_lat : `~astropy.units.Quantity`
        The longitude and latitude of the differentials.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = UnitSphericalRepresentation
    attr_classes = {"d_lon_coslat": u.Quantity, "d_lat": u.Quantity}

    @classproperty
    def _dimensional_differential(cls):
        return SphericalCosLatDifferential

    def __init__(self, d_lon_coslat, d_lat=None, copy=True):
        super().__init__(d_lon_coslat, d_lat, copy=copy)
        if not self._d_lon_coslat.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError("d_lon_coslat and d_lat should have equivalent units.")

    @classmethod
    def from_cartesian(cls, other, base):
        # Go via the dimensional equivalent, so that the longitude and latitude
        # differentials correctly take into account the norm of the base.
        dimensional = cls._dimensional_differential.from_cartesian(other, base)
        return dimensional.represent_as(cls)

    def to_cartesian(self, base):
        if isinstance(base, SphericalRepresentation):
            scale = base.distance
        elif isinstance(base, PhysicsSphericalRepresentation):
            scale = base.r
        else:
            return super().to_cartesian(base)

        base = base.represent_as(UnitSphericalRepresentation)
        return scale * super().to_cartesian(base)

    def represent_as(self, other_class, base=None):
        # Only have enough information to represent other unit-spherical.
        if issubclass(other_class, UnitSphericalDifferential):
            return other_class(self._d_lon(base), self.d_lat)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though w/o CosLat needs base for the latitude.
        if isinstance(representation, SphericalCosLatDifferential):
            return cls(representation.d_lon_coslat, representation.d_lat)
        elif isinstance(
            representation, (SphericalDifferential, UnitSphericalDifferential)
        ):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_lon, base)
            return cls(d_lon_coslat, representation.d_lat)
        elif isinstance(representation, PhysicsSphericalDifferential):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_phi, base)
            return cls(d_lon_coslat, -representation.d_theta)

        return super().from_representation(representation, base)

    def transform(self, matrix, base, transformed_base):
        """Transform differential using a 3x3 matrix in a Cartesian basis.

        This returns a new differential and does not modify the original one.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 (or stack thereof) matrix, such as a rotation matrix.
        base : instance of ``cls.base_representation``
            Base relative to which the differentials are defined.  If the other
            class is a differential representation, the base will be converted
            to its ``base_representation``.
        transformed_base : instance of ``cls.base_representation``
            Base relative to which the transformed differentials are defined.
            If the other class is a differential representation, the base will
            be converted to its ``base_representation``.
        """
        # the transformation matrix does not need to be a rotation matrix,
        # so the unit-distance is not guaranteed. For speed, we check if the
        # matrix is in O(3) and preserves lengths.
        if np.all(is_O3(matrix)):  # remain in unit-rep
            # TODO! implement without Cartesian intermediate step.
            diff = super().transform(matrix, base, transformed_base)

        else:  # switch to dimensional representation
            du = self.d_lat.unit / base.lat.unit  # derivative unit
            diff = self._dimensional_differential(
                d_lon_coslat=self.d_lon_coslat, d_lat=self.d_lat, d_distance=0 * du
            ).transform(matrix, base, transformed_base)

        return diff

    def _scale_operation(self, op, *args, scaled_base=False):
        if scaled_base:
            return self.copy()
        else:
            return super()._scale_operation(op, *args)


class SphericalCosLatDifferential(BaseSphericalCosLatDifferential):
    """Differential(s) of points in 3D spherical coordinates.

    Parameters
    ----------
    d_lon_coslat, d_lat : `~astropy.units.Quantity`
        The differential longitude (with cos(lat) included) and latitude.
    d_distance : `~astropy.units.Quantity`
        The differential distance.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = SphericalRepresentation
    _unit_differential = UnitSphericalCosLatDifferential
    attr_classes = {
        "d_lon_coslat": u.Quantity,
        "d_lat": u.Quantity,
        "d_distance": u.Quantity,
    }

    def __init__(self, d_lon_coslat, d_lat=None, d_distance=None, copy=True):
        super().__init__(d_lon_coslat, d_lat, d_distance, copy=copy)
        if not self._d_lon_coslat.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError("d_lon_coslat and d_lat should have equivalent units.")

    def represent_as(self, other_class, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though some need base for the latitude to remove cos(lat).
        if issubclass(other_class, UnitSphericalCosLatDifferential):
            return other_class(self.d_lon_coslat, self.d_lat)
        elif issubclass(other_class, RadialDifferential):
            return other_class(self.d_distance)
        elif issubclass(other_class, SphericalDifferential):
            return other_class(self._d_lon(base), self.d_lat, self.d_distance)
        elif issubclass(other_class, UnitSphericalDifferential):
            return other_class(self._d_lon(base), self.d_lat)
        elif issubclass(other_class, PhysicsSphericalDifferential):
            return other_class(self._d_lon(base), -self.d_lat, self.d_distance)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # Other spherical differentials can be done without going to Cartesian,
        # though we need base for the latitude to remove coslat.
        if isinstance(representation, SphericalDifferential):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_lon, base)
            return cls(d_lon_coslat, representation.d_lat, representation.d_distance)
        elif isinstance(representation, PhysicsSphericalDifferential):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_phi, base)
            return cls(d_lon_coslat, -representation.d_theta, representation.d_r)

        return super().from_representation(representation, base)

    def _scale_operation(self, op, *args, scaled_base=False):
        if scaled_base:
            return self.__class__(
                self.d_lon_coslat, self.d_lat, op(self.d_distance, *args)
            )
        else:
            return super()._scale_operation(op, *args)


class RadialDifferential(BaseDifferential):
    """Differential(s) of radial distances.

    Parameters
    ----------
    d_distance : `~astropy.units.Quantity`
        The differential distance.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = RadialRepresentation

    def to_cartesian(self, base):
        unit_vec = base.represent_as(UnitSphericalRepresentation).to_cartesian()
        return self.d_distance * unit_vec

    def norm(self, base=None):
        return self.d_distance

    @classmethod
    def from_cartesian(cls, other, base):
        return cls(
            other.dot(base.represent_as(UnitSphericalRepresentation)), copy=False
        )

    @classmethod
    def from_representation(cls, representation, base=None):
        if isinstance(
            representation, (SphericalDifferential, SphericalCosLatDifferential)
        ):
            return cls(representation.d_distance)
        elif isinstance(representation, PhysicsSphericalDifferential):
            return cls(representation.d_r)
        else:
            return super().from_representation(representation, base)

    def _combine_operation(self, op, other, reverse=False):
        if isinstance(other, self.base_representation):
            if reverse:
                first, second = other.distance, self.d_distance
            else:
                first, second = self.d_distance, other.distance
            return other.__class__(op(first, second), copy=False)
        elif isinstance(
            other, (BaseSphericalDifferential, BaseSphericalCosLatDifferential)
        ):
            all_components = set(self.components) | set(other.components)
            first, second = (self, other) if not reverse else (other, self)
            result_args = {
                c: op(getattr(first, c, 0.0), getattr(second, c, 0.0))
                for c in all_components
            }
            return SphericalDifferential(**result_args)

        else:
            return super()._combine_operation(op, other, reverse)


class PhysicsSphericalDifferential(BaseDifferential):
    """Differential(s) of 3D spherical coordinates using physics convention.

    Parameters
    ----------
    d_phi, d_theta : `~astropy.units.Quantity`
        The differential azimuth and inclination.
    d_r : `~astropy.units.Quantity`
        The differential radial distance.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = PhysicsSphericalRepresentation

    def __init__(self, d_phi, d_theta=None, d_r=None, copy=True):
        super().__init__(d_phi, d_theta, d_r, copy=copy)
        if not self._d_phi.unit.is_equivalent(self._d_theta.unit):
            raise u.UnitsError("d_phi and d_theta should have equivalent units.")

    def represent_as(self, other_class, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude. For those, explicitly
        # do the equivalent of self._d_lon_coslat in SphericalDifferential.
        if issubclass(other_class, SphericalDifferential):
            return other_class(self.d_phi, -self.d_theta, self.d_r)
        elif issubclass(other_class, UnitSphericalDifferential):
            return other_class(self.d_phi, -self.d_theta)
        elif issubclass(other_class, SphericalCosLatDifferential):
            self._check_base(base)
            d_lon_coslat = self.d_phi * np.sin(base.theta)
            return other_class(d_lon_coslat, -self.d_theta, self.d_r)
        elif issubclass(other_class, UnitSphericalCosLatDifferential):
            self._check_base(base)
            d_lon_coslat = self.d_phi * np.sin(base.theta)
            return other_class(d_lon_coslat, -self.d_theta)
        elif issubclass(other_class, RadialDifferential):
            return other_class(self.d_r)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # Other spherical differentials can be done without going to Cartesian,
        # though we need base for the latitude to remove coslat. For that case,
        # do the equivalent of cls._d_lon in SphericalDifferential.
        if isinstance(representation, SphericalDifferential):
            return cls(
                representation.d_lon, -representation.d_lat, representation.d_distance
            )
        elif isinstance(representation, SphericalCosLatDifferential):
            cls._check_base(base)
            d_phi = representation.d_lon_coslat / np.sin(base.theta)
            return cls(d_phi, -representation.d_lat, representation.d_distance)

        return super().from_representation(representation, base)

    def _scale_operation(self, op, *args, scaled_base=False):
        if scaled_base:
            return self.__class__(self.d_phi, self.d_theta, op(self.d_r, *args))
        else:
            return super()._scale_operation(op, *args)
