# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Implements roations, including spherical rotations as defined in WCS Paper II
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
from .parameters import Parameter, InputParameterError


__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial', 'Rotation2D']


class EulerAngleRotation(Model):
    """
    Base class for Euler angle rotations.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in deg
    """

    phi = Parameter(getter=np.rad2deg, setter=np.deg2rad)
    theta = Parameter(getter=np.rad2deg, setter=np.deg2rad)
    psi = Parameter(getter=np.rad2deg, setter=np.deg2rad)

    @staticmethod
    def _rotate_zxz(phi_i, theta_i, phi, theta, psi):
        """
        Defines a ZXZ rotation from initial coordinates phi_i, theta_i.

        All inputs and outputs are in radians.
        """

        cos_theta_i = np.cos(theta_i)
        sin_theta_i = np.sin(theta_i)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        delta = phi_i - psi
        cos_delta = np.cos(delta)

        phi_f = phi + np.arctan2(-cos_theta_i * np.sin(delta),
                                 sin_theta_i * cos_theta -
                                 cos_theta_i * sin_theta * cos_delta)

        theta_f = np.arcsin(sin_theta_i * sin_theta +
                            cos_theta_i * cos_theta * cos_delta)

        return phi_f, theta_f


class RotateNative2Celestial(EulerAngleRotation):
    """
    Transformation from Native to Celestial Spherical Coordinates.

    Defines a ZXZ rotation.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in deg
    """

    inputs = ('phi_N', 'theta_N')
    outputs = ('alpha_C', 'delta_C')

    @property
    def inverse(self):
        return RotateCelestial2Native(self.phi, self.theta, self.psi)

    @classmethod
    def evaluate(cls, phi_N, theta_N, phi, theta, psi):
        """
        Evaluate ZXZ rotation into celestial coordinates.

        Note currently the parameter values are passed in as degrees so we need
        to convert all inputs to radians and all outputs back to degrees.
        """

        # TODO: Unfortunately right now this superfluously converts from
        # radians to degrees back to radians again--this will be addressed in a
        # future change
        phi_N = np.deg2rad(phi_N)
        theta_N = np.deg2rad(theta_N)
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        psi = np.deg2rad(psi)

        alpha_C, delta_C = cls._rotate_zxz(phi_N, theta_N, phi, theta, psi)

        alpha_C = np.rad2deg(alpha_C)
        delta_C = np.rad2deg(delta_C)

        mask = alpha_C < 0
        if isinstance(mask, np.ndarray):
            alpha_C[mask] += 360
        elif mask:
            alpha_C += 360

        return alpha_C, delta_C


class RotateCelestial2Native(EulerAngleRotation):
    """
    Transformation from Celestial to Native to Spherical Coordinates.

    Defines a ZXZ rotation.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in deg
    """

    inputs = ('alpha_C', 'delta_C')
    outputs = ('phi_N', 'theta_N')

    @property
    def inverse(self):
        return RotateNative2Celestial(self.phi, self.theta, self.psi)

    @classmethod
    def evaluate(cls, alpha_C, delta_C, phi, theta, psi):
        """
        Evaluate ZXZ rotation into native coordinates.

        This is like RotateNative2Celestial.evaluate except phi and psi are
        swapped in ZXZ rotation.

        Note currently the parameter values are passed in as degrees so we need
        to convert all inputs to radians and all outputs back to degrees.
        """

        alpha_C = np.deg2rad(alpha_C)
        delta_C = np.deg2rad(delta_C)
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        psi = np.deg2rad(psi)

        phi_N, theta_N = cls._rotate_zxz(alpha_C, delta_C, psi, theta, phi)

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

    def __init__(self, angle=angle.default):
        super(Rotation2D, self).__init__(angle)
        self._matrix = self._compute_matrix(np.deg2rad(angle))

    @property
    def inverse(self):
        """Inverse rotation."""

        return self.__class__(angle=-self.angle)

    # TODO: This needs to be refactored so that it can work as a
    # static/classmethod
    def evaluate(self, x, y, angle):
        """
        Apply the rotation to a set of 2D Cartesian coordinates given as two
        lists--one for the x coordinates and one for a y coordinates--or a
        single coordinate pair.

        Parameters
        ----------
        x, y : array, float
            x and y coordinates
        """

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        shape = x.shape
        inarr = np.array([x.flatten(), y.flatten()], dtype=np.float64)
        if inarr.shape[0] != 2 or inarr.ndim != 2:
            raise ValueError("Incompatible input shapes")
        result = np.dot(self._matrix, inarr)
        x, y = result[0], result[1]
        if x.shape != shape:
            x.shape = shape
            y.shape = shape

        return x, y

    @staticmethod
    def _compute_matrix(angle):
        return np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]],
                        dtype=np.float64)
