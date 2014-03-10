# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Implements spherical rotations, defined in WCS Paper II [1]_

RotateNative2Celestial and RotateCelestial2Native follow the convention
in WCS paper II to rotate to/from a native sphere and the celestial sphere.

The user interface uses angles in degrees but the values are stored internally
in radians.  This is managed through the parameter setters/getters.

References
----------
.. [1] Calabretta, M.R., Greisen, E.W., 2002, A&A, 395, 1077 (Paper II)
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import math
import numbers

import numpy as np

from .core import Model
from .parameters import Parameter, InputParameterError


__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial',
           'RotateByAngle2D']


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
    phi = Parameter('phi', getter=np.rad2deg, setter=np.deg2rad)
    theta = Parameter('theta', getter=np.rad2deg, setter=np.deg2rad)
    psi = Parameter('psi', getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, phi, theta, psi):
        super(EulerAngleRotation, self).__init__(phi=phi, theta=theta, psi=psi)
        self.phi = phi
        self.theta = theta
        self.psi = psi


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
        phi = self._phi
        psi = self._psi
        theta = self._theta
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
        psi = self._psi
        phi = self._phi
        theta = self._theta

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


class RotateByAngle2D(Model):
    """
    Perform clockwise rotation through an angle.

    Parameters
    ----------
    angle : float
        angle of rotation in deg
    """

    def _validate_angle(angle):
        """Validates that an input angle is a number and converts it from
        degrees to radians.
        """

        if not isinstance(angle, numbers.Number):
            raise TypeError("Expected angle to be a number")

        return np.deg2rad(angle)

    n_inputs = 2
    n_outputs = 2

    angle = Parameter('angle', getter=np.rad2deg, setter=_validate_angle)

    def __init__(self, angle, param_dim=1):
        super(RotateByAngle2D, self).__init__(angle=angle, param_dim=param_dim)
        self.angle = angle
        self._matrix = self._compute_matrix(self._angle)

    def _compute_matrix(self, angle):
        return np.array([[math.cos(angle), math.sin(angle)],
                         [-math.sin(angle), math.cos(angle)]],
                        dtype=np.float64)

    def inverse(self):
        return RotateByAngle2D(angle=-np.rad2deg(self._angle))

    def __call__(self, x, y):
        """
        Parameters
        ----------
        x, y : 1D array or list
              x and y coordinates
        """

        x = np.asarray(x)
        y = np.asarray(y)
        assert x.shape == y.shape
        shape = x.shape
        inarr = np.array([x.flatten(), y.flatten()], dtype=np.float64)
        assert inarr.shape[0] == 2 and inarr.ndim == 2, \
            "Incompatible shape of input in RotateByAngle2D, expected (2, n)"
        result = np.dot(self._matrix, inarr)
        x, y = result[0], result[1]
        if x.shape != shape:
            x.shape = shape
            y.shape = shape
        return x, y
