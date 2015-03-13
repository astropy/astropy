# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Implements rotations, including spherical rotations as defined in WCS Paper II
[1]_

`RotateNative2Celestial` and `RotateCelestial2Native` follow the convention in
WCS Paper II to rotate to/from a native sphere and the celestial sphere.

The user interface sets and displays angles in degrees but the values are
stored internally in radians.  This is managed through the parameter
setters/getters.

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


__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial', 'Rotation2D',
           'EulerAngleRotation']


class EulerAngleRotation(Model):
    """
    Implements Euler angle intrinsic rotations.

    Rotates one coordinate system into another (fixed) coordinate system.
    All coordinate systems are right-handed. All intrinsic rotations are
    performed clockwise.

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

    def __init__(self, phi, theta, psi, axes_order):
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
        super(EulerAngleRotation, self).__init__(phi=phi, theta=theta, psi=psi)

    def _create_matrix(self, phi, theta, psi, axes_order):
        matrices = []
        for angle, axis in zip([phi, theta, psi], axes_order):
            matrix = np.zeros((3, 3), dtype=np.float)
            if axis == 'x':
                mat = self._rotation_matrix_from_angle(angle)
                matrix[0, 0] = 1
                matrix[1:, 1:] = mat
            elif axis == 'y':
                mat = self._rotation_matrix_from_angle(-angle)
                matrix[1, 1] = 1
                matrix[0, 0] = mat[0, 0]
                matrix[0, 2] = mat[0, 1]
                matrix[2, 0] = mat[1, 0]
                matrix[2, 2] = mat[1, 1]
            elif axis == 'z':
                mat = self._rotation_matrix_from_angle(angle)
                matrix[2, 2] = 1
                matrix[:2, :2] = mat
            else:
                raise ValueError("Expected axes_order to be a combination of characters"
                                 "'x', 'y' and 'z', got {0}".format(
                                     set(axes_order).difference(self.axes)))
            matrices.append(matrix)
        return np.dot(matrices[2], np.dot(matrices[1], matrices[0]))


    def _rotation_matrix_from_angle(self, angle):
        """
        Clockwise rotation matrix.
        """
        return np.array([[math.cos(angle), math.sin(angle)],
                         [-math.sin(angle), math.cos(angle)]])

    @staticmethod
    def spherical2cartesian(alpha, delta):
        alpha = np.deg2rad(alpha)
        delta = np.deg2rad(delta)

        x = np.cos(alpha) * np.cos(delta)
        y = np.cos(delta) * np.sin(alpha)
        z = np.sin(delta)
        return np.array([x, y, z])

    def inverse(self):
        return self.__class__(phi=-self.psi,
                              theta=-self.theta,
                              psi=-self.phi,
                              axes_order=self.axes_order[::-1])

    def evaluate(self, alpha, delta, phi, theta, psi):
        inp = self.spherical2cartesian(alpha, delta)
        matrix = self._create_matrix(phi, theta, psi, self.axes_order)
        result = np.dot(matrix, inp)
        return (np.rad2deg(np.arctan2(result[1], result[0])),
                np.rad2deg(np.arcsin(result[2])))


class _SkyRotation(Model):
    """
    Base class for FITS WCS sky rotations.

    Parameters
    ----------
    lon : float
        Celestial longitude of the fiducial point.
    lat : float
        Celestial latitude of the fiducial point.
    lon_pole : float
        Longitude of the celestial pole in the native system.
    """

    lon = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)
    lat = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)
    lon_pole = Parameter(default=0, getter=np.rad2deg, setter=np.deg2rad)

    @staticmethod
    def _rotate_zxz(phi_i, theta_i, lon, lat, lon_pole):
        """
        Defines a ZXZ rotation from initial coordinates phi_i, theta_i.

        All inputs and outputs are in radians.
        """

        cos_theta_i = np.cos(theta_i)
        sin_theta_i = np.sin(theta_i)
        cos_lat = np.cos(lat)
        sin_lat = np.sin(lat)
        delta = phi_i - lon_pole
        cos_delta = np.cos(delta)

        phi_f = lon + np.arctan2(-cos_theta_i * np.sin(delta),
                                 sin_theta_i * cos_lat -
                                 cos_theta_i * sin_lat * cos_delta)

        theta_f = np.arcsin(sin_theta_i * sin_lat +
                            cos_theta_i * cos_lat * cos_delta)

        return phi_f, theta_f


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

    inputs = ('phi_N', 'theta_N')
    outputs = ('alpha_C', 'delta_C')

    @property
    def inverse(self):
        return RotateCelestial2Native(self.lon, self.lat, self.lon_pole)

    @classmethod
    def evaluate(cls, phi_N, theta_N, lon, lat, lon_pole):
        """
        Rotate native spherical coordinates into celestial coordinates.
        """

        phi_N = np.deg2rad(phi_N)
        theta_N = np.deg2rad(theta_N)

        alpha_C, delta_C = cls._rotate_zxz(phi_N, theta_N, lon, lat, lon_pole)

        alpha_C = np.rad2deg(alpha_C)
        delta_C = np.rad2deg(delta_C)

        mask = alpha_C < 0
        if isinstance(mask, np.ndarray):
            alpha_C[mask] += 360
        elif mask:
            alpha_C += 360

        return alpha_C, delta_C


class RotateCelestial2Native(_SkyRotation):
    """
    Transform from Celestial to Native to Spherical Coordinates.

    Parameters
    ----------
    lon : float
        Celestial longitude of the fiducial point.
    lat : float
        Celestial latitude of the fiducial point.
    lon_pole : float
        Longitude of the celestial pole in the native system.
    """

    inputs = ('alpha_C', 'delta_C')
    outputs = ('phi_N', 'theta_N')

    @property
    def inverse(self):
        return RotateNative2Celestial(self.lon, self.lat, self.lon_pole)

    @classmethod
    def evaluate(cls, alpha_C, delta_C, lon, lat, lon_pole):
        """
        Rotate celestial coordinates into native spherical coordinates.

        This is the inverse transformation of RotateNative2Celestial.
        """

        alpha_C = np.deg2rad(alpha_C)
        delta_C = np.deg2rad(delta_C)

        phi_N, theta_N = cls._rotate_zxz(alpha_C, delta_C, lon_pole, lat, lon)

        phi_N = np.rad2deg(phi_N)
        theta_N = np.rad2deg(theta_N)

        mask = phi_N > 180
        if isinstance(mask, np.ndarray):
            phi_N[mask] -= 360
        elif mask:
            phi_N -= 360

        return phi_N, theta_N


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
