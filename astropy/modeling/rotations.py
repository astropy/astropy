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


__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial', 'Rotation2D',
           'EulerAngleRotation']


class _EulerRotation(object):
    """
    Base class which does the actual computation.
    """
    def _create_matrix(self, phi, theta, psi, axes_order):
        matrices = []
        for angle, axis in zip([phi, theta, psi], axes_order):
            matrix = np.zeros((3, 3), dtype=np.float)
            if axis == 'x':
                mat = self.rotation_matrix_from_angle(angle)
                matrix[0, 0] = 1
                matrix[1:, 1:] = mat
            elif axis == 'y':
                mat = self.rotation_matrix_from_angle(-angle)
                matrix[1, 1] = 1
                matrix[0, 0] = mat[0, 0]
                matrix[0, 2] = mat[0, 1]
                matrix[2, 0] = mat[1, 0]
                matrix[2, 2] = mat[1, 1]
            elif axis == 'z':
                mat = self.rotation_matrix_from_angle(angle)
                matrix[2, 2] = 1
                matrix[:2, :2] = mat
            else:
                raise ValueError("Expected axes_order to be a combination of characters"
                                 "'x', 'y' and 'z', got {0}".format(
                                     set(axes_order).difference(self.axes)))
            matrices.append(matrix)
        return np.dot(matrices[2], np.dot(matrices[1], matrices[0]))

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
        alpha  = np.rad2deg(np.arctan2(y, x))
        delta = np.rad2deg(np.arctan2(z, h))
        return alpha, delta

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

    Rotates one coordinate system into another (fixed) coordinate system.
    All coordinate systems are right-handed. The sign of the angles is
    determined by the right-hand rule..

    Parameters
    ----------
    phi, theta, psi : float
        "proper" Euler angles in deg
    axes_order : str
        A 3 character string, a combination of 'x', 'y' and 'z',
        where each character denotes an axis in 3D space.
    """

    inputs = ('alpha', 'delta')
    outputs = ('alpha', 'delta')

    phi = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)
    theta = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)
    psi = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)

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
        super(EulerAngleRotation, self).__init__(phi=phi, theta=theta, psi=psi, **kwargs)

    def inverse(self):
        return self.__class__(phi=-self.psi,
                              theta=-self.theta,
                              psi=-self.phi,
                              axes_order=self.axes_order[::-1])

    def evaluate(self, alpha, delta, phi, theta, psi):
        shape = None
        if isinstance(alpha, np.ndarray) and alpha.ndim == 2:
            alpha = alpha.flatten()
            delta = delta.flatten()
            shape = alpha.shape
        inp = self.spherical2cartesian(alpha, delta)
        matrix = self._create_matrix(phi, theta, psi, self.axes_order)
        result = np.dot(matrix, inp)
        a, b = self.cartesian2spherical(*result)
        if shape is not None:
            a.shape = shape
            b.shape = shape
        return a, b


class _SkyRotation(_EulerRotation, Model):
    """
    Base class for RotateNative2Celestial and RotateCelestial2Native.
    """

    lon = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)
    lat = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)
    lon_pole = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, lon, lat, lon_pole, **kwargs):
        super(_SkyRotation, self).__init__(lon, lat, lon_pole, **kwargs)
        self.axes_order = 'zxz'

    def _evaluate(self, phi, theta, lon, lat, lon_pole):
        alpha, delta = super(_SkyRotation, self).evaluate(phi, theta, lon, lat,
                                                          lon_pole, self.axes_order)
        mask = alpha < 0
        if isinstance(mask, np.ndarray):
            alpha[mask] +=360
        else:
            alpha +=360
        return alpha, delta


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
        phi = lon_pole - np.pi / 2
        theta = - (np.pi / 2 - lat)
        psi = -(np.pi / 2 + lon)
        alpha_C, delta_C = super(RotateNative2Celestial, self)._evaluate(phi_N, theta_N,
                                                                         phi, theta, psi)
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
        phi = (np.pi / 2 + lon)
        theta =  (np.pi / 2 - lat)
        psi = -(lon_pole - np.pi / 2)
        phi_N, theta_N = super(RotateCelestial2Native, self)._evaluate(alpha_C, delta_C,
                                                                       phi, theta, psi)
        return phi_N, theta_N

    @property
    def inverse(self):
        return RotateNative2Celestial(self.lon, self.lat, self.lon_pole)


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

    angle = Parameter(default=0.0, getter=np.rad2deg, setter=np.deg2rad)

    @property
    def inverse(self):
        """Inverse rotation."""

        return self.__class__(angle=-self.angle)

    @classmethod
    def evaluate(cls, x, y, angle):
        """
        Apply the rotation to a set of 2D Cartesian coordinates given as two
        lists--one for the x coordinates and one for a y coordinates--or a
        single coordinate pair.
        """

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        # Note: If the original shape was () (an array scalar) convert to a
        # 1-element 1-D array on output for consistency with most other models
        orig_shape = x.shape or (1,)

        inarr = np.array([x.flatten(), y.flatten()])
        result = np.dot(cls._compute_matrix(angle), inarr)

        x, y = result[0], result[1]
        x.shape = y.shape = orig_shape

        return x, y

    @staticmethod
    def _compute_matrix(angle):
        return np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]],
                        dtype=np.float64)
