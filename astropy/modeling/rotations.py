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

    n_inputs = 2
    n_outputs = 2
    phi = Parameter(getter=np.rad2deg, setter=np.deg2rad)
    theta = Parameter(getter=np.rad2deg, setter=np.deg2rad)
    psi = Parameter(getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, phi, theta, psi):
        super(EulerAngleRotation, self).__init__(phi, theta, psi)


class RotateNative2Celestial(EulerAngleRotation):
    """
    Transformation from Native to Celestial Spherical Coordinates.

    Defines a ZXZ rotation.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in deg
    """

    def inverse(self):
        return RotateCelestial2Native(self.phi, self.theta, self.psi)

    def __call__(self, nphi, ntheta):
        nphi = np.deg2rad(nphi)
        ntheta = np.deg2rad(ntheta)
        # TODO: Unfortunately right now this superfluously converts from
        # radians to degrees back to radians again--this will be addressed in a
        # future change
        phi = np.deg2rad(self.phi)
        psi = np.deg2rad(self.psi)
        theta = np.deg2rad(self.theta)

        calpha = np.rad2deg(
            phi +
            np.arctan2(-np.cos(ntheta) * np.sin(nphi - psi),
                       np.sin(ntheta) * np.cos(theta) -
                       np.cos(ntheta) * np.sin(theta) * np.cos(nphi - psi)))

        cdelta = np.rad2deg(
            np.arcsin(np.sin(ntheta) * np.sin(theta) +
                      np.cos(ntheta) * np.cos(theta) * np.cos(nphi - psi)))

        ind = calpha < 0
        if isinstance(ind, np.ndarray):
            calpha[ind] += 360
        elif ind:
            calpha += 360

        return calpha, cdelta


class RotateCelestial2Native(EulerAngleRotation):
    """
    Transformation from Celestial to Native to Spherical Coordinates.

    Defines a ZXZ rotation.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in deg
    """

    def inverse(self):
        return RotateNative2Celestial(self.phi, self.theta, self.psi)

    def __call__(self, calpha, cdelta):
        calpha = np.deg2rad(calpha)
        cdelta = np.deg2rad(cdelta)

        # TODO: Unfortunately right now this superfluously converts from
        # radians to degrees back to radians again--this will be addressed in a
        # future change
        phi = np.deg2rad(self.phi)
        psi = np.deg2rad(self.psi)
        theta = np.deg2rad(self.theta)

        nphi = np.rad2deg(
            psi +
            np.arctan2(-np.cos(cdelta) * np.sin(calpha - phi),
                       np.sin(cdelta) * np.cos(theta) -
                       np.cos(cdelta) * np.sin(theta) * np.cos(calpha - phi)))

        ntheta = np.rad2deg(
            np.arcsin(np.sin(cdelta) * np.sin(theta) +
                      np.cos(cdelta) * np.cos(theta) * np.cos(calpha - phi)))

        ind = nphi > 180
        if isinstance(ind, np.ndarray):
            nphi[ind] -= 360
        elif ind:
            nphi -= 360

        return nphi, ntheta


class Rotation2D(Model):
    """
    Perform a 2D rotation given an angle in degrees.

    Positive angles represent a counter-clockwise rotation and vice-versa.

    Parameters
    ----------
    angle : float
        angle of rotation in deg
    """

    n_inputs = 2
    n_outputs = 2

    angle = Parameter(default=0.0, getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, angle=angle.default):
        super(Rotation2D, self).__init__(angle)
        self._matrix = self._compute_matrix(np.deg2rad(angle))

    def inverse(self):
        """Inverse rotation."""

        return self.__class__(angle=-self.angle)

    def __call__(self, x, y):
        """
        Apply the rotation to a set of 2D Cartesian coordinates given as two
        lists--one for the x coordinates and one for a y coordinates--or a
        single coordinate pair.

        Parameters
        ----------
        x, y : array, float
            x and y coordinates
        """

        x = np.asarray(x)
        y = np.asarray(y)
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
