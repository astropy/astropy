# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

"""
Implements projections--particularly sky projections defined in WCS Paper II
[1]_

All angles are set and and displayed in degrees but internally computations are
performed in radians.

References
----------
.. [1] Calabretta, M.R., Greisen, E.W., 2002, A&A, 395, 1077 (Paper II)
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import abc

import numpy as np

from .core import Model
from .parameters import Parameter, InputParameterError

from ..utils.compat import ignored


projcodes = ['TAN', 'AZP', 'SZP', 'STG', 'SIN', 'ARC', 'ZPN', 'ZEA', 'AIR',
             'CYP', 'CEA', 'MER']


__all__ = ['Projection', 'Pix2SkyProjection', 'Sky2PixProjection',
           'Zenithal', 'Cylindrical',
           'Pix2Sky_AZP', 'Sky2Pix_AZP', 'Pix2Sky_CAR', 'Sky2Pix_CAR',
           'Pix2Sky_CEA', 'Sky2Pix_CEA', 'Pix2Sky_CYP', 'Sky2Pix_CYP',
           'Pix2Sky_MER', 'Sky2Pix_MER',
           'Pix2Sky_SIN', 'Sky2Pix_SIN', 'Pix2Sky_STG', 'Sky2Pix_STG',
           'Pix2Sky_TAN', 'Sky2Pix_TAN',
           'AffineTransformation2D']


class Projection(Model):
    """Base class for all sky projections."""

    # the radius of the projection sphere, by which x,y are scaled
    r0 = 180 / np.pi

    @abc.abstractproperty
    def inverse(self):
        """
        Inverse projection--all projection models must provide an inverse.
        """


class Pix2SkyProjection(Projection):
    """Base class for all Pix2Sky projections."""

    inputs = ('x', 'y')
    outputs = ('phi', 'theta')


class Sky2PixProjection(Projection):
    """Base class for all Sky2Pix projections."""

    inputs = ('phi', 'theta')
    outputs = ('x', 'y')


class Zenithal(Projection):
    r"""Base class for all Zenithal projections.

    Zenithal (or azimuthal) projections map the sphere directly onto a
    plane.  All zenithal projections are specified by defining the
    radius as a function of native latitude, :math:`R_\theta`.

    The pixel-to-sky transformation is defined as:

    .. math::
        \phi &= \arg(-y, x) \\
        R_\theta &= \sqrt{x^2 + y^2}

    and the inverse (sky-to-pixel) is defined as:

    .. math::
        x &= R_\theta \sin \phi \\
        y &= R_\theta \cos \phi
    """


class Pix2Sky_AZP(Pix2SkyProjection, Zenithal):
    r"""
    AZP : Zenithal perspective projection - pixel to sky.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = \arg(\rho, 1) - \sin^{-1}\left(\frac{\rho \mu}{\sqrt{\rho^2 + 1}}\right)

    where:

    .. math::
        \rho = \frac{\pi}{180^{\circ}}\frac{R_\theta}{\mu + 1}

    Parameters
    --------------
    mu : float
        distance from point of projection to center of sphere
        in spherical radii, default is 0.
    gamma : float
        look angle in deg, default is 0.
    """


    def _validate_mu(mu):
        if np.asarray(mu == -1).any():
            raise ValueError("AZP projection is not defined for mu=-1")
        return mu

    mu = Parameter(default=0.0, setter=_validate_mu)
    gamma = Parameter(default=0.0, getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, mu=mu.default, gamma=gamma.default, **kwargs):
        self.check_mu(mu)
        # units : mu - in spherical radii, gamma - in deg
        # TODO: Support quantity objects here and in similar contexts
        super(Pix2Sky_AZP, self).__init__(mu, gamma, **kwargs)

    def check_mu(self, val):
        if np.asarray(val == -1).any():
            raise ValueError("AZP projection is not defined for mu=-1")

    @property
    def inverse(self):
        return Sky2Pix_AZP(self.mu.value, self.gamma.value)

    @classmethod
    def evaluate(cls, x, y, mu, gamma):
        phi = np.arctan2(x / np.cos(gamma), -y)
        r = cls._compute_r_theta(x, y, gamma)
        pho = r / (cls.r0 * (mu + 1) + y * np.sin(gamma))
        psi = np.arctan2(1, pho)
        omega = np.arcsin((pho * mu) / np.sqrt(pho ** 2 + 1))

        theta1 = np.rad2deg(psi - omega)
        theta2 = np.rad2deg(psi + omega) + 180

        if np.abs(mu) < 1:
            if theta1 < 90 and theta1 > -90:
                theta = theta1
            else:
                theta = theta2
        else:
            # theta1dif = 90 - theta1
            # theta2dif = 90 - theta2
            if theta1 < theta2:
                theta = theta1
            else:
                theta = theta2

        phi = np.rad2deg(phi)
        return phi, theta

    @staticmethod
    def _compute_r_theta(x, y, gamma):
        return np.sqrt(x ** 2 + y ** 2 * (np.cos(gamma)) ** 2)


class Sky2Pix_AZP(Sky2PixProjection, Zenithal):
    r"""
    AZP : Zenital perspective projection - sky to pixel.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\frac{(\mu + 1) \cos \theta}{\mu + \sin \theta}

    Parameters
    ----------
    mu : float
        distance from point of projection to center of sphere
        in spherical radii, default is 0.

    gamma : float
        look angle in deg, default is 0.
    """

    def _validate_mu(mu):
        if np.asarray(mu == -1).any():
            raise ValueError("AZP projection is not defined for mu=-1")
        return mu

    mu = Parameter(default=0.0, setter=_validate_mu)
    gamma = Parameter(default=0.0, getter=np.rad2deg, setter=np.deg2rad)

    def check_mu(self, val):
        if np.asarray(val == -1).any():
            raise ValueError("AZP projection is not defined for mu=-1")

    @property
    def inverse(self):
        return Pix2Sky_AZP(self.mu.value, self.gamma.value)

    @classmethod
    def evaluate(cls, phi, theta, mu, gamma):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)

        r = cls._compute_r_theta(phi, theta, mu, gamma)
        x = r * np.sin(phi)
        y = (-r * np.cos(phi)) / np.cos(gamma)

        return x, y

    @classmethod
    def _compute_r_theta(cls, phi, theta, mu, gamma):
        return ((cls.r0 * (mu + 1) * np.cos(theta)) /
                (mu + np.sin(theta) +
                 np.cos(theta) * np.cos(phi) * np.tan(gamma)))


class Pix2Sky_TAN(Pix2SkyProjection, Zenithal):
    r"""
    TAN : Gnomonic projection - pixel to sky.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = \tan^{-1}\left(\frac{180^{\circ}}{\pi R_\theta}\right)
    """

    @property
    def inverse(self):
        return Sky2Pix_TAN()

    @classmethod
    def evaluate(cls, x, y):
        phi = np.rad2deg(np.arctan2(x, -y))
        r_theta = cls._compute_r_theta(x, y)
        theta = np.rad2deg(np.arctan2(cls.r0, r_theta))

        return phi, theta

    @staticmethod
    def _compute_r_theta(x, y):
        return np.sqrt(x ** 2 + y ** 2)


class Sky2Pix_TAN(Sky2PixProjection, Zenithal):
    r"""
    TAN : Gnomonic Projection - sky to pixel.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\cot \theta
    """

    @property
    def inverse(self):
        return Pix2Sky_TAN()

    @classmethod
    def evaluate(cls, phi, theta):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)

        r_theta = cls._compute_r_theta(theta)
        x = np.rad2deg(r_theta * np.sin(phi))
        y = -np.rad2deg(r_theta * np.cos(phi))

        return x, y

    @staticmethod
    def _compute_r_theta(theta):
        return 1 / np.tan(theta)


class Pix2Sky_STG(Pix2SkyProjection, Zenithal):
    r"""
    STG : Stereographic Projection - pixel to sky.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = 90^{\circ} - 2 \tan^{-1}\left(\frac{\pi R_\theta}{360^{\circ}}\right)
    """

    @property
    def inverse(self):
        return Sky2Pix_STG()

    @classmethod
    def evaluate(cls, x, y):
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = cls._compute_r_theta(x, y)
        theta = 90 - np.rad2deg(2 * np.arctan(rtheta / (2 * cls.r0)))

        return phi, theta

    @staticmethod
    def _compute_r_theta(x, y):
        return np.sqrt(x ** 2 + y ** 2)


class Sky2Pix_STG(Sky2PixProjection, Zenithal):
    r"""
    STG : Stereographic Projection - sky to pixel.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\frac{2 \cos \theta}{1 + \sin \theta}
    """

    @property
    def inverse(self):
        return Pix2Sky_STG()

    @classmethod
    def evaluate(cls, phi, theta):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)

        r_theta = cls._compute_r_theta(theta)
        x = r_theta * np.sin(phi)
        y = -r_theta * np.cos(phi)

        return x, y

    @classmethod
    def _compute_r_theta(cls, theta):
        return (cls.r0 * 2 * np.cos(theta)) / (1 + np.sin(theta))


class Pix2Sky_SIN(Pix2SkyProjection, Zenithal):
    r"""
    SIN : Slant orthographic projection - pixel to sky.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = \cos^{-1}\left(\frac{\pi}{180^{\circ}}R_\theta\right)
    """

    @property
    def inverse(self):
        return Sky2Pix_SIN()

    @classmethod
    def evaluate(cls, x, y):
        r_theta = cls._compute_r_theta(x, y)
        phi = np.rad2deg(np.arctan2(x, -y))
        theta = np.rad2deg(np.arccos(r_theta / cls.r0))

        return phi, theta

    @staticmethod
    def _compute_r_theta(x, y):
        return np.sqrt(x ** 2 + y ** 2)


class Sky2Pix_SIN(Sky2PixProjection, Zenithal):
    r"""
    SIN : Slant orthographic projection - sky to pixel.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\cos \theta
    """

    @property
    def inverse(self):
        return Pix2Sky_SIN()

    @classmethod
    def evaluate(cls, phi, theta):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        r_theta = cls._compute_r_theta(theta)
        x = r_theta * np.sin(phi)
        y = -r_theta * np.cos(phi)

        return x, y

    @classmethod
    def _compute_r_theta(cls, theta):
        return cls.r0 * np.cos(theta)


class Cylindrical(Projection):
    r"""Base class for Cylindrical projections.

    Cylindrical projections are so-named because the surface of
    projection is a cylinder.
    """


class Pix2Sky_CYP(Pix2SkyProjection, Cylindrical):
    r"""
    CYP : Cylindrical perspective - pixel to sky.

    .. math::
        \phi &= \frac{x}{\lambda} \\
        \theta &= \arg(1, \eta) + \sin{-1}\left(\frac{\eta \mu}{\sqrt{\eta^2 + 1}}\right)

    where:

    .. math::
        \eta = \frac{\pi}{180^{\circ}}\frac{y}{\mu + \lambda}

    Parameters
    ----------
    mu : float
        distance from center of sphere in the direction opposite the
        projected surface, in spherical radii, default is 0.

    lam : float
        radius of the cylinder in spherical radii, default is 0.
    """

    def _validate_mu(mu, model):
        with ignored(AttributeError):
            # An attribute error can occur if model.lam has not been set yet
            if np.asarray(mu == -model.lam).any():
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return mu

    def _validate_lam(lam, model):
        with ignored(AttributeError):
            # An attribute error can occur if model.lam has not been set yet
            if np.asarray(lam == -model.mu).any():
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return lam

    mu = Parameter(default=0, setter=_validate_mu)
    lam = Parameter(default=0, setter=_validate_lam)

    @property
    def inverse(self):
        return Sky2Pix_CYP(self.mu.value, self.lam.value)

    @classmethod
    def evaluate(cls, x, y, mu, lam):
        phi = x / lam
        eta = y / (cls.r0 * (mu + lam))
        theta = (np.arctan2(eta, 1) +
                 np.arcsin(eta * mu / np.sqrt(eta ** 2 + 1)))

        return phi, np.rad2deg(theta)


class Sky2Pix_CYP(Sky2PixProjection, Cylindrical):
    r"""
    CYP : Cylindrical Perspective - sky to pixel.

    .. math::
        x &= \lambda \phi
        y &= \frac{180^{\circ}}{\pi}\left(\frac{\mu + \lambda}{\mu + \cos \theta}\right)\sin \theta

    Parameters
    ----------
    mu : float
        distance from center of sphere in the direction opposite the
        projected surface, in spherical radii, default is 0.

    lam : float
        radius of the cylinder in spherical radii, default is 0.
    """

    # TODO: Eliminate duplication on these
    def _validate_mu(mu, model):
        with ignored(AttributeError):
            if np.asarray(mu == -model.lam).any():
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return mu

    def _validate_lam(lam, model):
        with ignored(AttributeError):
            if np.asarray(lam == -model.mu).any():
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return lam

    mu = Parameter(default=0, setter=_validate_mu)
    lam = Parameter(default=0, setter=_validate_lam)

    @property
    def inverse(self):
        return Pix2Sky_CYP(self.mu, self.lam)

    @classmethod
    def evaluate(cls, phi, theta, mu, lam):
        theta = np.deg2rad(theta)
        x = lam * phi
        y = (cls.r0 * (mu + lam) / (mu + np.cos(theta))) * np.sin(theta)

        return x, y


class Pix2Sky_CEA(Pix2SkyProjection, Cylindrical):
    r"""
    CEA : Cylindrical equal area projection - pixel to sky.

    .. math::
        \phi &= x \\
        \theta &= \sin^{-1}\left(\frac{\pi}{180^{\circ}}\lambda y\right)

    Parameters
    ----------
    lam : float
        radius of the cylinder in spherical radii, default is 0.
    """

    lam = Parameter(default=1)

    @property
    def inverse(self):
        return Sky2Pix_CEA(self.lam)

    @classmethod
    def evaluate(cls, x, y, lam):
        phi = x.copy()
        theta = np.rad2deg(np.arcsin(1 / cls.r0 * lam * y))

        return phi, theta


class Sky2Pix_CEA(Sky2PixProjection, Cylindrical):
    r"""
    CEA: Cylindrical equal area projection - sky to pixel.

    .. math::
        x &= \phi \\
        y &= \frac{180^{\circ}}{\pi}\frac{\sin \theta}{\lambda}

    Parameters
    ----------
    lam : float
        radius of the cylinder in spherical radii, default is 0.
    """

    lam = Parameter(default=1)

    @property
    def inverse(self):
        return Pix2Sky_CEA(self.lam)

    @classmethod
    def evaluate(cls, phi, theta, lam):
        x = phi.copy()
        theta = np.deg2rad(theta)
        y = cls.r0 * np.sin(theta) / lam

        return x, y


class Pix2Sky_CAR(Pix2SkyProjection, Cylindrical):
    r"""
    CAR: Plate carrée projection - pixel to sky.

    .. math::
        \phi &= x \\
        \theta &= y
    """

    @property
    def inverse(self):
        return Sky2Pix_CAR()

    @staticmethod
    def evaluate(x, y):
        # The intermediate variables are only used here for clarity
        phi = x.copy()
        theta = y.copy()

        return phi, theta


class Sky2Pix_CAR(Sky2PixProjection, Cylindrical):
    r"""
    CAR: Plate carrée projection - sky to pixel.

    .. math::
        x &= \phi \\
        y &= \theta
    """

    @property
    def inverse(self):
        return Pix2Sky_CAR()

    @staticmethod
    def evaluate(phi, theta):
        # The intermediate variables are only used here for clarity
        x = phi.copy()
        y = theta.copy()

        return x, y


class Pix2Sky_MER(Pix2SkyProjection, Cylindrical):
    r"""
    MER: Mercator - pixel to sky.

    .. math::
        \phi &= x \\
        \theta &= 2 \tan^{-1}\left(e^{y \pi / 180^{\circ}}\right)-90^{\circ}
    """

    @property
    def inverse(self):
        return Sky2Pix_MER()

    @classmethod
    def evaluate(cls, x, y):
        phi = x.copy()
        theta = np.rad2deg(2 * np.arctan(np.exp(y / cls.r0))) - 90

        return phi, theta


class Sky2Pix_MER(Sky2PixProjection, Cylindrical):
    r"""
    MER: Mercator - sky to pixel.

    .. math::
        x &= \phi \\
        y &= \frac{180^{\circ}}{\pi}\ln \tan \left(\frac{90^{\circ} + \theta}{2}\right)
    """

    @property
    def inverse(self):
        return Pix2Sky_MER()

    @classmethod
    def evaluate(cls, phi, theta):
        x = phi.copy()
        theta = np.deg2rad(theta)
        y = cls.r0 * np.log(np.tan((np.pi / 2 + theta) / 2))

        return x, y


class AffineTransformation2D(Model):
    """
    Perform an affine transformation in 2 dimensions.

    Parameters
    ----------
    matrix : array
        A 2x2 matrix specifying the linear transformation to apply to the
        inputs

    translation : array
        A 2D vector (given as either a 2x1 or 1x2 array) specifying a
        translation to apply to the inputs
    """

    inputs = ('x', 'y')
    outputs = ('x', 'y')

    standard_broadcasting = False

    matrix = Parameter(
        setter=lambda m: AffineTransformation2D._validate_matrix(m),
        default=[[1.0, 0.0], [0.0, 1.0]])
    translation = Parameter(
        setter=lambda t: AffineTransformation2D._validate_vector(t),
        default=[0.0, 0.0])

    @property
    def inverse(self):
        """
        Inverse transformation.

        Raises `~astropy.modeling.InputParameterError` if the transformation cannot be inverted.
        """

        det = np.linalg.det(self.matrix.value)

        if det == 0:
            raise InputParameterError(
                "Transformation matrix is singular; {0} model does not "
                "have an inverse".format(self.__class__.__name__))

        matrix = np.linalg.inv(self.matrix.value)
        translation = -np.dot(matrix, self.translation.value)

        return self.__class__(matrix=matrix, translation=translation)

    @classmethod
    def evaluate(cls, x, y, matrix, translation):
        """
        Apply the transformation to a set of 2D Cartesian coordinates given as
        two lists--one for the x coordinates and one for a y coordinates--or a
        single coordinate pair.

        Parameters
        ----------
        x, y : array, float
              x and y coordinates
        """

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        shape = x.shape or (1,)
        inarr = np.vstack([x.flatten(), y.flatten(), np.ones(x.size)])

        if inarr.shape[0] != 3 or inarr.ndim != 2:
            raise ValueError("Incompatible input shapes")

        augmented_matrix = cls._create_augmented_matrix(matrix, translation)
        result = np.dot(augmented_matrix, inarr)

        x, y = result[0], result[1]
        x.shape = y.shape = shape

        return x, y

    @staticmethod
    def _validate_matrix(matrix):
        """Validates that the input matrix is a 2x2 2D array."""

        matrix = np.array(matrix)
        if matrix.shape != (2, 2):
            raise ValueError(
                "Expected transformation matrix to be a 2x2 array")
        return matrix

    @staticmethod
    def _validate_vector(vector):
        """
        Validates that the translation vector is a 2D vector.  This allows
        either a "row" vector or a "column" vector where in the latter case the
        resultant Numpy array has ``ndim=2`` but the shape is ``(1, 2)``.
        """

        vector = np.array(vector)

        if vector.ndim == 1:
            vector = vector[:, np.newaxis]

        if vector.shape != (2, 1):
            raise ValueError(
                "Expected translation vector to be a 2 element row or column "
                "vector array")

        return vector

    @staticmethod
    def _create_augmented_matrix(matrix, translation):
        augmented_matrix = np.empty((3, 3), dtype=np.float)
        augmented_matrix[0:2, 0:2] = matrix
        augmented_matrix[0:2, 2:].flat = translation
        augmented_matrix[2] = [0, 0, 1]
        return augmented_matrix
