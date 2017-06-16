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

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import math

import numpy as np

from .core import Model
from .parameters import Parameter
from ..extern.six.moves import zip
from ..coordinates.matrix_utilities import rotation_matrix, matrix_product
from .. import units as u
from ..utils.decorators import deprecated


__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial', 'Rotation2D',
           'EulerAngleRotation']


class _EulerRotation(object):
    """
    Base class which does the actual computation.
    """

    # We allow values without units to be passed when evaluating the model, and
    # in this case the input (alpha, delta) values are assumed to be deg.
    input_units_allow_dimensionless = True

    def _create_matrix(self, phi, theta, psi, axes_order):
        matrices = []
        for angle, axis in zip([phi, theta, psi], axes_order):
            angle = np.asscalar(angle)
            matrices.append(rotation_matrix(angle, axis, unit=u.rad))
        return matrix_product(*matrices[::-1])

    @staticmethod
    def spherical2cartesian(alpha, delta):
        if isinstance(alpha, u.Quantity):
            alpha = alpha.to(u.rad)
        else:
            alpha = np.deg2rad(alpha)
        if isinstance(delta, u.Quantity):
            delta = delta.to(u.rad)
        else:
            delta = np.deg2rad(delta)
        x = np.cos(alpha) * np.cos(delta)
        y = np.cos(delta) * np.sin(alpha)
        z = np.sin(delta)
        return np.array([x, y, z])

    @staticmethod
    def cartesian2spherical(x, y, z):
        h = np.hypot(x, y)
        alpha  = np.rad2deg(np.arctan2(y, x))
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


class EulerAngleRotation(_EulerRotation, Model):
    """
    Implements Euler angle intrinsic rotations.

    Notes
    -----
    Rotates one coordinate system into another (fixed) coordinate system.
    All coordinate systems are right-handed. The sign of the angles is
    determined by the right-hand rule.
    The result has the same units as the initializing Euler angles.
    If the initializing Euler angles are numbers they should be in units of deg,
    in which case the output units are degrees.

    Parameters
    ----------
    phi, theta, psi : float or `~astropy.units.Quantity`
        "proper" Euler angles
    axes_order : str
        A 3 character string, a combination of 'x', 'y' and 'z',
        where each character denotes an axis in 3D space.

    Returns
    -------
    alpha, delta : float or `~astropy.units.Quantity`
        Rotated angles in the units of the Euler angles.
        Default units are degrees.

    """

    inputs = ('x', 'y')
    outputs = ('alpha', 'delta')

    phi = Parameter(default=0)
    theta = Parameter(default=0)
    psi = Parameter(default=0)

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

        if not isinstance(phi, u.Quantity):
            phi = np.deg2rad(phi)
        else:
            phi = phi.to(u.rad)
        if not isinstance(psi, u.Quantity):
            psi = np.deg2rad(psi)
        else:
            psi = psi.to(u.rad)
        if not isinstance(theta, u.Quantity):
            theta = np.deg2rad(theta)
        else:
            theta = theta.to(u.rad)
        super(EulerAngleRotation, self).__init__(phi=phi, theta=theta, psi=psi, **kwargs)

    def inverse(self):
        return self.__class__(phi=-self.psi,
                              theta=-self.theta,
                              psi=-self.phi,
                              axes_order=self.axes_order[::-1])

    def evaluate(self, alpha, delta, phi, theta, psi):
        a, b = super(EulerAngleRotation, self).evaluate(alpha, delta, phi, theta, psi, self.axes_order)
        if isinstance(alpha, u.Quantity):
            if isinstance(a, u.Quantity):
                a = a.to(alpha.unit)
                b = b.to(alpha.unit)
            else:
                # it's in deg
                a = u.Quantity(a, unit=u.deg).to(alpha.unit)
                b = u.Quantity(b, unit=u.deg).to(alpha.unit)
        return a, b


class _SkyRotation(_EulerRotation, Model):
    """
    Base class for RotateNative2Celestial and RotateCelestial2Native.
    """

    lon = Parameter(default=0)
    lat = Parameter(default=0)
    lon_pole = Parameter(default=0)

    def __init__(self, lon, lat, lon_pole, **kwargs):

        self.axes_order = 'zxz'

        qs = [isinstance(par, u.Quantity) for par in [lon, lat, lon_pole]]
        if any(qs) and not all(qs):
            raise TypeError("All parameters should be of the same type - float or Quantity.")

        if not isinstance(lon, u.Quantity):
            lon = np.deg2rad(lon)
        else:
            lon = lon.to(u.rad)
        if not isinstance(lon_pole, u.Quantity):
            lon_pole = np.deg2rad(lon_pole)
        else:
            lon_pole = lon_pole.to(u.rad)
        if not isinstance(lat, u.Quantity):
            lat = np.deg2rad(lat)
        else:
            lat = lat.to(u.rad)

        super(_SkyRotation, self).__init__(lon, lat, lon_pole, **kwargs)

    def _evaluate(self, phi, theta, lon, lat, lon_pole):
        alpha, delta = super(_SkyRotation, self).evaluate(phi, theta, lon, lat,
                                                          lon_pole, self.axes_order)
        if isinstance(alpha, u.Quantity):
            avalue = alpha.value
        else:
            avalue = alpha
        mask = avalue < 0
        if isinstance(mask, np.ndarray):
            avalue[mask] += 360
        else:
            avalue += 360
        if isinstance(alpha, u.Quantity):
            alpha = u.Quantity(avalue, alpha.unit)
        return alpha, delta

    def _convert_parameters(self):
        if self.lon.unit is None:
            lon = np.rad2deg(self.lon)
        else:
            lon = self.lon
        if self.lat.unit is None:
            lat = np.rad2deg(self.lat)
        else:
            lat = self.lat
        if self.lon_pole.unit is None:
            lon_pole = np.rad2deg(self.lon_pole)
        else:
            lon_pole = self.lon_pole
        return lon, lat, lon_pole

class RotateNative2Celestial(_SkyRotation):
    """
    Transform from Native to Celestial Spherical Coordinates.

    Parameters
    ----------
    lon : float
        Celestial longitude of the fiducial point.
    lat : float
        Celestial latitude of the fiducial point.
    lon_pole : float
        Longitude of the celestial pole in the native system.
    """

    # angles in degrees on the native sphere
    inputs = ('phi_N', 'theta_N')
    outputs = ('alpha_C', 'delta_C')

    def __init__(self, lon, lat, lon_pole, **kwargs):
        super(RotateNative2Celestial, self).__init__(lon, lat, lon_pole, **kwargs)

    def evaluate(self, phi_N, theta_N, lon, lat, lon_pole):
        # Convert to Euler angles
        if isinstance(lon, u.Quantity):
            phi = lon_pole.value - np.pi / 2
            theta = - (np.pi / 2 - lat.value)
            psi = -(np.pi / 2 + lon.value)
        else:
            phi = lon_pole - np.pi / 2
            theta = - (np.pi / 2 - lat)
            psi = -(np.pi / 2 + lon)
        alpha_C, delta_C = super(RotateNative2Celestial, self)._evaluate(phi_N, theta_N,
                                                                         phi, theta, psi)
        return alpha_C, delta_C

    @property
    def inverse(self):
        # convert to deg
        lon, lat, lon_pole  = self._convert_parameters()
        return RotateCelestial2Native(lon, lat, lon_pole)


class RotateCelestial2Native(_SkyRotation):
    """
    Transform from Celestial to Native Spherical Coordinates.

    Parameters
    ----------
    lon : float
        Celestial longitude of the fiducial point.
    lat : float
        Celestial latitude of the fiducial point.
    lon_pole : float
        Longitude of the celestial pole in the native system.
    """

    # angles in degrees on the celestial sphere
    inputs = ('alpha_C', 'delta_C')
    # angles in degrees on the native sphere
    outputs = ('phi_N', 'theta_N')


    def __init__(self, lon, lat, lon_pole, **kwargs):
        super(RotateCelestial2Native, self).__init__(lon, lat, lon_pole, **kwargs)

    def evaluate(self, alpha_C, delta_C, lon, lat, lon_pole):
        # Convert to Euler angles
        if isinstance(lon, u.Quantity):
            phi = (np.pi / 2 + lon.value)
            theta =  (np.pi / 2 - lat.value)
            psi = -(lon_pole.value - np.pi / 2)
        else:
            phi = (np.pi / 2 + lon)
            theta =  (np.pi / 2 - lat)
            psi = -(lon_pole - np.pi / 2)
        phi_N, theta_N = super(RotateCelestial2Native, self)._evaluate(alpha_C, delta_C,
                                                                       phi, theta, psi)
        return phi_N, theta_N

    @property
    def inverse(self):
        # convert to deg
        lon, lat, lon_pole  = self._convert_parameters()
        return RotateNative2Celestial(lon, lat, lon_pole)


class Rotation2D(Model):
    """
    Perform a 2D rotation given an angle in degrees.

    Positive angles represent a counter-clockwise rotation and vice-versa.

    Parameters
    ----------
    angle : float
        angle of rotation in deg
    """

    inputs = ('x', 'y')
    outputs = ('x', 'y')

    angle = Parameter(default=0.0)

    # We allow values without units to be passed when evaluating the model, and
    # in this case the input (x, y) values are assumed to be deg.
    input_units_allow_dimensionless = True

    def __init__(self, angle, **kwargs):
        if isinstance(angle, u.Quantity):
            angle = angle.to(u.rad)
        else:
            angle = np.deg2rad(angle)

        super(Rotation2D, self).__init__(angle=angle, **kwargs)

    @property
    def inverse(self):
        """Inverse rotation."""
        if self.angle.unit is not None:
            return self.__class__(angle=-self.angle)
        else:
            return self.__class__(angle=-np.rad2deg(self.angle))

    @classmethod
    def evaluate(cls, x, y, angle):
        """
        Apply the rotation to a set of 2D Cartesian coordinates given as two
        lists--one for the x coordinates and one for a y coordinates--or a
        single coordinate pair.
        """
        unit = None
        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        # Note: If the original shape was () (an array scalar) convert to a
        # 1-element 1-D array on output for consistency with most other models
        orig_shape = x.shape or (1,)
        if isinstance(x, u.Quantity):
            unit = x.unit
        inarr = np.array([x.flatten(), y.flatten()])
        if isinstance(angle, u.Quantity):
            angle = angle.value

        result = np.dot(cls._compute_matrix(angle), inarr)
        x, y = result[0], result[1]
        x.shape = y.shape = orig_shape
        if unit is not None:
            return u.Quantity(x, unit=unit), u.Quantity(y, unit=unit)
        else:
            return x, y

    @staticmethod
    def _compute_matrix(angle):
        return np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]],
                        dtype=np.float64)
