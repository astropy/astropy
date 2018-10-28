# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Implements rotations, including spherical rotations as defined in WCS Paper II
[1]_

`RotateNative2Celestial` and `RotateCelestial2Native` follow the convention in
WCS Paper II to rotate to/from a native sphere and the celestial sphere.

The implementation uses `EulerAngleRotation`. The model parameters are
three angles: the longitude (``lon``) and latitude (``lat``) of the fiducial point
in the celestial system (``CRVAL`` keywords in FITS), and the longitude of the celestial
pole in the native system (``lon_pole``). The Euler angles are ``lon+90``, ``90-lat``
and ``-(lon_pole-90)``.


References
----------
.. [1] Calabretta, M.R., Greisen, E.W., 2002, A&A, 395, 1077 (Paper II)
"""

import math

import numpy as np

from .core import Model
from .parameters import Parameter
from ..coordinates.matrix_utilities import rotation_matrix, matrix_product
from .. import units as u
from ..utils.decorators import deprecated
from .utils import _to_radian, _to_orig_unit

__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial', 'Rotation2D',
           'EulerAngleRotation']


class _EulerRotation:
    """
    Base class which does the actual computation.
    """

    _separable = False

    def _create_matrix(self, phi, theta, psi, axes_order):
        matrices = []
        for angle, axis in zip([phi, theta, psi], axes_order):
            if isinstance(angle, u.Quantity):
                angle = angle.value
            angle = angle.item()
            matrices.append(rotation_matrix(angle, axis, unit=u.rad))
        result = matrix_product(*matrices[::-1])
        return result

    @staticmethod
    def spherical2cartesian(alpha, delta):
        alpha = np.deg2rad(alpha)
        delta = np.deg2rad(delta)
        x = np.cos(alpha) * np.cos(delta)
        y = np.cos(delta) * np.sin(alpha)
        z = np.sin(delta)
        return np.array([x, y, z])

    @staticmethod
    def cartesian2spherical(x, y, z):
        h = np.hypot(x, y)
        alpha = np.rad2deg(np.arctan2(y, x))
        delta = np.rad2deg(np.arctan2(z, h))
        return alpha, delta

    @deprecated(2.0)
    @staticmethod
    def rotation_matrix_from_angle(angle):
        """
        Clockwise rotation matrix.

        Parameters
        ----------
        angle : float
            Rotation angle in radians.
        """
        return np.array([[math.cos(angle), math.sin(angle)],
                         [-math.sin(angle), math.cos(angle)]])

    def evaluate(self, alpha, delta, phi, theta, psi, axes_order):
        shape = None
        if isinstance(alpha, np.ndarray) and alpha.ndim == 2:
            alpha = alpha.flatten()
            delta = delta.flatten()
            shape = alpha.shape
        inp = self.spherical2cartesian(alpha, delta)
        matrix = self._create_matrix(phi, theta, psi, axes_order)
        result = np.dot(matrix, inp)
        a, b = self.cartesian2spherical(*result)
        if shape is not None:
            a.shape = shape
            b.shape = shape
        return a, b

    _input_units_strict = True

    _input_units_allow_dimensionless = True

    @property
    def input_units(self):
        """ Input units. """
        return {'alpha': u.deg, 'delta': u.deg}

    @property
    def return_units(self):
        """ Output units. """
        return {'alpha': u.deg, 'delta': u.deg}


class EulerAngleRotation(_EulerRotation, Model):
    """
    Implements Euler angle intrinsic rotations.

    Rotates one coordinate system into another (fixed) coordinate system.
    All coordinate systems are right-handed. The sign of the angles is
    determined by the right-hand rule..

    Parameters
    ----------
    phi, theta, psi : float or `~astropy.units.Quantity`
        "proper" Euler angles in deg.
        If floats, they should be in deg.
    axes_order : str
        A 3 character string, a combination of 'x', 'y' and 'z',
        where each character denotes an axis in 3D space.
    """

    inputs = ('alpha', 'delta')
    outputs = ('alpha', 'delta')

    phi = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian)
    theta = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian)
    psi = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian)

    def __init__(self, phi, theta, psi, axes_order, **kwargs):
        self.axes = ['x', 'y', 'z']
        if len(axes_order) != 3:
            raise TypeError(
                "Expected axes_order to be a character sequence of length 3,"
                "got {0}".format(axes_order))
        unrecognized = set(axes_order).difference(self.axes)
        if unrecognized:
            raise ValueError("Unrecognized axis label {0}; "
                             "should be one of {1} ".format(unrecognized, self.axes))
        self.axes_order = axes_order
        qs = [isinstance(par, u.Quantity) for par in [phi, theta, psi]]
        if any(qs) and not all(qs):
            raise TypeError("All parameters should be of the same type - float or Quantity.")

        super().__init__(phi=phi, theta=theta, psi=psi, **kwargs)

    def inverse(self):
        return self.__class__(phi=-self.psi,
                              theta=-self.theta,
                              psi=-self.phi,
                              axes_order=self.axes_order[::-1])

    def evaluate(self, alpha, delta, phi, theta, psi):
        a, b = super().evaluate(alpha, delta, phi, theta, psi, self.axes_order)
        return a, b


class _SkyRotation(_EulerRotation, Model):
    """
    Base class for RotateNative2Celestial and RotateCelestial2Native.
    """

    lon = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian)
    lat = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian)
    lon_pole = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian)

    def __init__(self, lon, lat, lon_pole, **kwargs):
        qs = [isinstance(par, u.Quantity) for par in [lon, lat, lon_pole]]
        if any(qs) and not all(qs):
            raise TypeError("All parameters should be of the same type - float or Quantity.")
        super().__init__(lon, lat, lon_pole, **kwargs)
        self.axes_order = 'zxz'

    def _evaluate(self, phi, theta, lon, lat, lon_pole):
        alpha, delta = super().evaluate(phi, theta, lon, lat, lon_pole,
                                        self.axes_order)
        mask = alpha < 0
        if isinstance(mask, np.ndarray):
            alpha[mask] += 360
        else:
            alpha += 360
        return alpha, delta


class RotateNative2Celestial(_SkyRotation):
    """
    Transform from Native to Celestial Spherical Coordinates.

    Parameters
    ----------
    lon : float or or `~astropy.units.Quantity`
        Celestial longitude of the fiducial point.
    lat : float or or `~astropy.units.Quantity`
        Celestial latitude of the fiducial point.
    lon_pole : float or or `~astropy.units.Quantity`
        Longitude of the celestial pole in the native system.

    Notes
    -----
    If ``lon``, ``lat`` and ``lon_pole`` are numerical values they should be in units of deg.
    """

    #: Inputs are angles on the native sphere
    inputs = ('phi_N', 'theta_N')

    #: Outputs are angles on the celestial sphere
    outputs = ('alpha_C', 'delta_C')

    @property
    def input_units(self):
        """ Input units. """
        return {'phi_N': u.deg, 'theta_N': u.deg}

    @property
    def return_units(self):
        """ Output units. """
        return {'alpha_C': u.deg, 'delta_C': u.deg}

    def __init__(self, lon, lat, lon_pole, **kwargs):
        super().__init__(lon, lat, lon_pole, **kwargs)

    def evaluate(self, phi_N, theta_N, lon, lat, lon_pole):
        """
        Parameters
        ----------
        phi_N, theta_N : float (deg) or `~astropy.units.Quantity`
            Angles in the Native coordinate system.
        lon, lat, lon_pole : float (in deg) or `~astropy.units.Quantity`
            Parameter values when the model was initialized.

        Returns
        -------
        alpha_C, delta_C : float (deg) or `~astropy.units.Quantity`
            Angles on the Celestial sphere.
        """
        # The values are in radians since they have already been through the setter.
        if isinstance(lon, u.Quantity):
            lon = lon.value
            lat = lat.value
            lon_pole = lon_pole.value
        # Convert to Euler angles
        phi = lon_pole - np.pi / 2
        theta = - (np.pi / 2 - lat)
        psi = -(np.pi / 2 + lon)
        alpha_C, delta_C = super()._evaluate(phi_N, theta_N, phi, theta, psi)
        return alpha_C, delta_C

    @property
    def inverse(self):
        # convert to angles on the celestial sphere
        return RotateCelestial2Native(self.lon, self.lat, self.lon_pole)


class RotateCelestial2Native(_SkyRotation):
    """
    Transform from Celestial to Native Spherical Coordinates.

    Parameters
    ----------
    lon : float or or `~astropy.units.Quantity`
        Celestial longitude of the fiducial point.
    lat : float or or `~astropy.units.Quantity`
        Celestial latitude of the fiducial point.
    lon_pole : float or or `~astropy.units.Quantity`
        Longitude of the celestial pole in the native system.

    Notes
    -----
    If ``lon``, ``lat`` and ``lon_pole`` are numerical values they should be in units of deg.
    """

    #: Inputs are angles on the celestial sphere
    inputs = ('alpha_C', 'delta_C')

    #: Outputs are angles on the native sphere
    outputs = ('phi_N', 'theta_N')

    @property
    def input_units(self):
        """ Input units. """
        return {'alpha_C': u.deg, 'delta_C': u.deg}

    @property
    def return_units(self):
        """ Output units. """
        return {'phi_N': u.deg, 'theta_N': u.deg}

    def __init__(self, lon, lat, lon_pole, **kwargs):
        super().__init__(lon, lat, lon_pole, **kwargs)

    def evaluate(self, alpha_C, delta_C, lon, lat, lon_pole):
        """
        Parameters
        ----------
        alpha_C, delta_C : float (deg) or `~astropy.units.Quantity`
            Angles in the Celestial coordinate frame.
        lon, lat, lon_pole : float (deg) or `~astropy.units.Quantity`
            Parameter values when the model was initialized.

        Returns
        -------
        phi_N, theta_N : float (deg) or `~astropy.units.Quantity`
            Angles on the Native sphere.

        """
        if isinstance(lon, u.Quantity):
            lon = lon.value
            lat = lat.value
            lon_pole = lon_pole.value
        # Convert to Euler angles
        phi = (np.pi / 2 + lon)
        theta = (np.pi / 2 - lat)
        psi = -(lon_pole - np.pi / 2)
        phi_N, theta_N = super()._evaluate(alpha_C, delta_C, phi, theta, psi)

        return phi_N, theta_N

    @property
    def inverse(self):
        return RotateNative2Celestial(self.lon, self.lat, self.lon_pole)


class Rotation2D(Model):
    """
    Perform a 2D rotation given an angle.

    Positive angles represent a counter-clockwise rotation and vice-versa.

    Parameters
    ----------
    angle : float or `~astropy.units.Quantity`
        Angle of rotation (if float it should be in deg).
    """

    inputs = ('x', 'y')
    outputs = ('x', 'y')
    _separable = False

    angle = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)

    @property
    def inverse(self):
        """Inverse rotation."""

        return self.__class__(angle=-self.angle)

    @classmethod
    def evaluate(cls, x, y, angle):
        """
        Rotate (x, y) about ``angle``.

        Parameters
        ----------
        x, y : ndarray-like
            Input quantities
        angle : float (deg) or `~astropy.units.Quantity`
            Angle of rotations.

        """

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        # If one argument has units, enforce they both have units and they are compatible.
        x_unit = getattr(x, 'unit', None)
        y_unit = getattr(y, 'unit', None)
        has_units = x_unit is not None and y_unit is not None
        if x_unit != y_unit:
            if has_units and y_unit.is_equivalent(x_unit):
                y = y.to(x_unit)
                y_unit = x_unit
            else:
                raise u.UnitsError("x and y must have compatible units")

        # Note: If the original shape was () (an array scalar) convert to a
        # 1-element 1-D array on output for consistency with most other models
        orig_shape = x.shape or (1,)
        inarr = np.array([x.flatten(), y.flatten()])
        if isinstance(angle, u.Quantity):
            angle = angle.to_value(u.rad)
        result = np.dot(cls._compute_matrix(angle), inarr)
        x, y = result[0], result[1]
        x.shape = y.shape = orig_shape
        if has_units:
            return u.Quantity(x, unit=x_unit), u.Quantity(y, unit=y_unit)
        else:
            return x, y

    @staticmethod
    def _compute_matrix(angle):
        return np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]],
                        dtype=np.float64)
