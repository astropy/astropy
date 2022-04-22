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
# pylint: disable=invalid-name, too-many-arguments, no-member

import math

import numpy as np

from astropy import units as u
from astropy.coordinates.matrix_utilities import matrix_product, rotation_matrix

from .core import Model
from .parameters import Parameter
from .utils import _to_orig_unit, _to_radian

__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial', 'Rotation2D',
           'EulerAngleRotation', 'RotationSequence3D', 'SphericalRotationSequence']


def _create_matrix(angles, axes_order):
    matrices = []
    for angle, axis in zip(angles, axes_order):
        if isinstance(angle, u.Quantity):
            angle = angle.value
        angle = angle.item()
        matrices.append(rotation_matrix(angle, axis, unit=u.rad))
    result = matrix_product(*matrices[::-1])
    return result


def spherical2cartesian(alpha, delta):
    alpha = np.deg2rad(alpha)
    delta = np.deg2rad(delta)
    x = np.cos(alpha) * np.cos(delta)
    y = np.cos(delta) * np.sin(alpha)
    z = np.sin(delta)
    return np.array([x, y, z])


def cartesian2spherical(x, y, z):
    h = np.hypot(x, y)
    alpha = np.rad2deg(np.arctan2(y, x))
    delta = np.rad2deg(np.arctan2(z, h))
    return alpha, delta


class RotationSequence3D(Model):
    """
    Perform a series of rotations about different axis in 3D space.

    Positive angles represent a counter-clockwise rotation.

    Parameters
    ----------
    angles : array-like
        Angles of rotation in deg in the order of axes_order.
    axes_order : str
        A sequence of 'x', 'y', 'z' corresponding to axis of rotation.

    Examples
    --------
    >>> model = RotationSequence3D([1.1, 2.1, 3.1, 4.1], axes_order='xyzx')

    """
    standard_broadcasting = False
    _separable = False
    n_inputs = 3
    n_outputs = 3

    angles = Parameter(default=[], getter=_to_orig_unit, setter=_to_radian,
                       description="Angles of rotation in deg in the order of axes_order")

    def __init__(self, angles, axes_order, name=None):
        self.axes = ['x', 'y', 'z']
        unrecognized = set(axes_order).difference(self.axes)
        if unrecognized:
            raise ValueError(f"Unrecognized axis label {unrecognized}; "
                             f"should be one of {self.axes} ")
        self.axes_order = axes_order
        if len(angles) != len(axes_order):
            raise ValueError(f"The number of angles {len(angles)} should match "
                             f"the number of axes {len(axes_order)}.")
        super().__init__(angles, name=name)
        self._inputs = ('x', 'y', 'z')
        self._outputs = ('x', 'y', 'z')

    @property
    def inverse(self):
        """Inverse rotation."""
        angles = self.angles.value[::-1] * -1
        return self.__class__(angles, axes_order=self.axes_order[::-1])

    def evaluate(self, x, y, z, angles):
        """
        Apply the rotation to a set of 3D Cartesian coordinates.
        """
        if x.shape != y.shape or x.shape != z.shape:
            raise ValueError("Expected input arrays to have the same shape")
        # Note: If the original shape was () (an array scalar) convert to a
        # 1-element 1-D array on output for consistency with most other models
        orig_shape = x.shape or (1,)
        inarr = np.array([x.flatten(), y.flatten(), z.flatten()])
        result = np.dot(_create_matrix(angles[0], self.axes_order), inarr)
        x, y, z = result[0], result[1], result[2]
        x.shape = y.shape = z.shape = orig_shape
        return x, y, z


class SphericalRotationSequence(RotationSequence3D):
    """
    Perform a sequence of rotations about arbitrary number of axes
    in spherical coordinates.

    Parameters
    ----------
    angles : list
        A sequence of angles (in deg).
    axes_order : str
        A sequence of characters ('x', 'y', or 'z') corresponding to the
        axis of rotation and matching the order in ``angles``.

    """
    def __init__(self, angles, axes_order, name=None, **kwargs):
        self._n_inputs = 2
        self._n_outputs = 2
        super().__init__(angles, axes_order=axes_order, name=name, **kwargs)
        self._inputs = ("lon", "lat")
        self._outputs = ("lon", "lat")

    @property
    def n_inputs(self):
        return self._n_inputs

    @property
    def n_outputs(self):
        return self._n_outputs

    def evaluate(self, lon, lat, angles):
        x, y, z = spherical2cartesian(lon, lat)
        x1, y1, z1 = super().evaluate(x, y, z, angles)
        lon, lat = cartesian2spherical(x1, y1, z1)
        return lon, lat


class _EulerRotation:
    """
    Base class which does the actual computation.
    """

    _separable = False

    def evaluate(self, alpha, delta, phi, theta, psi, axes_order):
        shape = None
        if isinstance(alpha, np.ndarray):
            alpha = alpha.flatten()
            delta = delta.flatten()
            shape = alpha.shape
        inp = spherical2cartesian(alpha, delta)
        matrix = _create_matrix([phi, theta, psi], axes_order)
        result = np.dot(matrix, inp)
        a, b = cartesian2spherical(*result)
        if shape is not None:
            a.shape = shape
            b.shape = shape
        return a, b

    _input_units_strict = True

    _input_units_allow_dimensionless = True

    @property
    def input_units(self):
        """ Input units. """
        return {self.inputs[0]: u.deg,
                self.inputs[1]: u.deg}

    @property
    def return_units(self):
        """ Output units. """
        return {self.outputs[0]: u.deg,
                self.outputs[1]: u.deg}


class EulerAngleRotation(_EulerRotation, Model):
    """
    Implements Euler angle intrinsic rotations.

    Rotates one coordinate system into another (fixed) coordinate system.
    All coordinate systems are right-handed. The sign of the angles is
    determined by the right-hand rule..

    Parameters
    ----------
    phi, theta, psi : float or `~astropy.units.Quantity` ['angle']
        "proper" Euler angles in deg.
        If floats, they should be in deg.
    axes_order : str
        A 3 character string, a combination of 'x', 'y' and 'z',
        where each character denotes an axis in 3D space.
    """

    n_inputs = 2
    n_outputs = 2

    phi = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian,
                    description="1st Euler angle (Quantity or value in deg)")
    theta = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian,
                      description="2nd Euler angle (Quantity or value in deg)")
    psi = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian,
                    description="3rd Euler angle (Quantity or value in deg)")

    def __init__(self, phi, theta, psi, axes_order, **kwargs):
        self.axes = ['x', 'y', 'z']
        if len(axes_order) != 3:
            raise TypeError(
                "Expected axes_order to be a character sequence of length 3, "
                f"got {axes_order}")
        unrecognized = set(axes_order).difference(self.axes)
        if unrecognized:
            raise ValueError(f"Unrecognized axis label {unrecognized}; "
                             f"should be one of {self.axes}")
        self.axes_order = axes_order
        qs = [isinstance(par, u.Quantity) for par in [phi, theta, psi]]
        if any(qs) and not all(qs):
            raise TypeError("All parameters should be of the same type - float or Quantity.")

        super().__init__(phi=phi, theta=theta, psi=psi, **kwargs)
        self._inputs = ('alpha', 'delta')
        self._outputs = ('alpha', 'delta')

    @property
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

    lon = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian,
                    description="Latitude")
    lat = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian,
                    description="Longtitude")
    lon_pole = Parameter(default=0, getter=_to_orig_unit, setter=_to_radian,
                         description="Longitude of a pole")

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
    lon : float or `~astropy.units.Quantity` ['angle']
        Celestial longitude of the fiducial point.
    lat : float or `~astropy.units.Quantity` ['angle']
        Celestial latitude of the fiducial point.
    lon_pole : float or `~astropy.units.Quantity` ['angle']
        Longitude of the celestial pole in the native system.

    Notes
    -----
    If ``lon``, ``lat`` and ``lon_pole`` are numerical values they
    should be in units of deg. Inputs are angles on the native sphere.
    Outputs are angles on the celestial sphere.
    """

    n_inputs = 2
    n_outputs = 2

    @property
    def input_units(self):
        """ Input units. """
        return {self.inputs[0]: u.deg,
                self.inputs[1]: u.deg}

    @property
    def return_units(self):
        """ Output units. """
        return {self.outputs[0]: u.deg, self.outputs[1]: u.deg}

    def __init__(self, lon, lat, lon_pole, **kwargs):
        super().__init__(lon, lat, lon_pole, **kwargs)
        self.inputs = ('phi_N', 'theta_N')
        self.outputs = ('alpha_C', 'delta_C')

    def evaluate(self, phi_N, theta_N, lon, lat, lon_pole):
        """
        Parameters
        ----------
        phi_N, theta_N : float or `~astropy.units.Quantity` ['angle']
            Angles in the Native coordinate system.
            it is assumed that numerical only inputs are in degrees.
            If float, assumed in degrees.
        lon, lat, lon_pole : float or `~astropy.units.Quantity` ['angle']
            Parameter values when the model was initialized.
            If float, assumed in degrees.

        Returns
        -------
        alpha_C, delta_C : float or `~astropy.units.Quantity` ['angle']
            Angles on the Celestial sphere.
            If float, in degrees.
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
    lon : float or `~astropy.units.Quantity` ['angle']
        Celestial longitude of the fiducial point.
    lat : float or `~astropy.units.Quantity` ['angle']
        Celestial latitude of the fiducial point.
    lon_pole : float or `~astropy.units.Quantity` ['angle']
        Longitude of the celestial pole in the native system.

    Notes
    -----
    If ``lon``, ``lat`` and ``lon_pole`` are numerical values they should be
    in units of deg. Inputs are angles on the celestial sphere.
    Outputs are angles on the native sphere.
    """
    n_inputs = 2
    n_outputs = 2

    @property
    def input_units(self):
        """ Input units. """
        return {self.inputs[0]: u.deg,
                self.inputs[1]: u.deg}

    @property
    def return_units(self):
        """ Output units. """
        return {self.outputs[0]: u.deg,
                self.outputs[1]: u.deg}

    def __init__(self, lon, lat, lon_pole, **kwargs):
        super().__init__(lon, lat, lon_pole, **kwargs)

        # Inputs are angles on the celestial sphere
        self.inputs = ('alpha_C', 'delta_C')
        # Outputs are angles on the native sphere
        self.outputs = ('phi_N', 'theta_N')

    def evaluate(self, alpha_C, delta_C, lon, lat, lon_pole):
        """
        Parameters
        ----------
        alpha_C, delta_C : float or `~astropy.units.Quantity` ['angle']
            Angles in the Celestial coordinate frame.
            If float, assumed in degrees.
        lon, lat, lon_pole : float or `~astropy.units.Quantity` ['angle']
            Parameter values when the model was initialized.
            If float, assumed in degrees.

        Returns
        -------
        phi_N, theta_N : float or `~astropy.units.Quantity` ['angle']
            Angles on the Native sphere.
            If float, in degrees.

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
    angle : float or `~astropy.units.Quantity` ['angle']
        Angle of rotation (if float it should be in deg).
    """
    n_inputs = 2
    n_outputs = 2

    _separable = False

    angle = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian,
                      description="Angle of rotation (Quantity or value in deg)")

    def __init__(self, angle=angle, **kwargs):
        super().__init__(angle=angle, **kwargs)
        self._inputs = ("x", "y")
        self._outputs = ("x", "y")

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
        x, y : array-like
            Input quantities
        angle : float or `~astropy.units.Quantity` ['angle']
            Angle of rotations.
            If float, assumed in degrees.

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
        return x, y

    @staticmethod
    def _compute_matrix(angle):
        return np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]],
                        dtype=np.float64)
