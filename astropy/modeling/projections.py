# Licensed under a 3-clause BSD style license - see LICENSE.rst

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

import numpy as np

from .core import (Model, format_input)
from .parameters import Parameter, InputParameterError

from ..utils.compat import ignored


projcodes = ['TAN', 'AZP', 'SZP', 'STG', 'SIN', 'ARC', 'ZPN', 'ZEA', 'AIR',
             'CYP', 'CEA', 'MER']


__all__ = ['Pix2Sky_AZP', 'Sky2Pix_AZP', 'Pix2Sky_CAR', 'Sky2Pix_CAR',
           'Pix2Sky_CEA', 'Sky2Pix_CEA', 'Pix2Sky_CYP', 'Sky2Pix_CYP',
           'Pix2Sky_MER', 'Sky2Pix_MER',
           'Pix2Sky_SIN', 'Sky2Pix_SIN', 'Pix2Sky_STG', 'Sky2Pix_STG',
           'Pix2Sky_TAN', 'Sky2Pix_TAN',
           'AffineTransformation2D']


class Projection(Model):
    """
    Base class for all sky projections.

    """

    n_inputs = 2
    n_outputs = 2

    # the radius of the projection sphere, by which x,y are scaled
    r0 = 180 / np.pi


class Zenithal(Projection):
    """
    Base class for all Zenithal projections.

    """

    def _compute_rtheta(self):
        # Subclasses must implement this method
        raise NotImplementedError("Subclasses should implement this")

    def __call__(self, x, y):
        raise NotImplementedError("Subclasses should implement this")


class Pix2Sky_AZP(Zenithal):
    """
    AZP : Zenital perspective projection - pixel to sky.

    Parameters
    --------------
    mu : float
        distance from point of projection to center of sphere
        in spherical radii, default is 0.
    gamma : float
        look angle in deg, default is 0.
    """


    def _validate_mu(mu):
        if mu == -1:
            raise ValueError("AZP projection is not defined for mu=-1")
        return mu

    mu = Parameter(default=0.0, setter=_validate_mu)
    gamma = Parameter(default=0.0, getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, mu=mu.default, gamma=gamma.default):
        self.check_mu(mu)
        # units : mu - in spherical radii, gamma - in deg
        # TODO: Support quantity objects here and in similar contexts
        super(Pix2Sky_AZP, self).__init__(mu, gamma)

    def check_mu(self, val):
        if val == -1:
            raise ValueError("AZP projection is not defined for mu=-1")

    def inverse(self):
        return Sky2Pix_AZP(self.mu.value, self.gamma.value)

    # TODO: This contains many superfluous conversions of gamma from radians to
    # degrees back to radians again--this is going to be reworked in a future
    # changeset
    def _compute_rtheta(self, x, y):
        gamma = np.deg2rad(self.gamma)
        return np.sqrt(x ** 2 + y ** 2 * (np.cos(gamma)) ** 2)

    @format_input
    def __call__(self, x, y, model_set_axis=None):
        gamma = np.deg2rad(self.gamma)
        phi = np.rad2deg(np.arctan2(x / np.cos(gamma), -y))
        r = self._compute_rtheta(x, y)
        pho = r / (self.r0 * (self.mu + 1) +
                   y * np.sin(gamma))
        psi = np.arctan2(1, pho)
        omega = np.arcsin((pho * self.mu) / (np.sqrt(pho ** 2 + 1)))
        theta1 = np.rad2deg(psi - omega)
        theta2 = np.rad2deg(psi + omega) + 180
        if np.abs(self.mu) < 1:
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
        return phi, theta


class Sky2Pix_AZP(Zenithal):
    """
    AZP : Zenital perspective projection - sky to pixel.

    Parameters
    --------------
    mu : float
        distance from point of projection to center of sphere
        in spherical radii, default is 0.
    gamma : float
        look angle in deg, default is 0.
    """

    def _validate_mu(mu):
        if mu == -1:
            raise ValueError("AZP projection is not defined for mu=-1")
        return mu

    mu = Parameter(default=0.0, setter=_validate_mu)
    gamma = Parameter(default=0.0, getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, mu=mu.default, gamma=gamma.default):
        super(Sky2Pix_AZP, self).__init__(mu, gamma)

    def check_mu(self, val):
        if val == -1:
            raise ValueError("AZP projection is not defined for mu=-1")

    def inverse(self):
        return Pix2Sky_AZP(self.mu.value, self.gamma.value)

    def _compute_rtheta(self, phi, theta):
        gamma = np.deg2rad(self.gamma)
        rtheta = (self.r0 * (self.mu + 1) *
                  np.cos(theta)) / ((self.mu + np.sin(theta)) +
                                    np.cos(theta) * np.cos(phi) *
                                    np.tan(gamma))
        return rtheta

    @format_input
    def __call__(self, phi, theta, model_set_axis=None):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        gamma = np.deg2rad(self.gamma)
        r = self._compute_rtheta(phi, theta)
        x = r * np.sin(phi)
        y = (-r * np.cos(phi)) / np.cos(gamma)
        return x, y


class Pix2Sky_TAN(Zenithal):
    """
    TAN : Gnomonic projection - pixel to sky.
    """

    def _compute_rtheta(self, x, y):
        return np.sqrt(x ** 2 + y ** 2)

    def inverse(self):
        return Sky2Pix_TAN()

    @format_input
    def __call__(self, x, y, model_set_axis=None):
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = self._compute_rtheta(x, y)
        theta = np.rad2deg(np.arctan2(self.r0, rtheta))
        return phi, theta


class Sky2Pix_TAN(Zenithal):
    """
    TAN : Gnomonic Projection - sky to pixel.
    """

    def _compute_rtheta(self, theta):
        return 1 / np.tan(theta)

    def inverse(self):
        return Pix2Sky_TAN()

    @format_input
    def __call__(self, phi, theta, model_set_axis=None):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        rtheta = self._compute_rtheta(theta)
        x = np.rad2deg(rtheta * np.sin(phi))
        y = - np.rad2deg(rtheta * np.cos(phi))
        return x, y


class Pix2Sky_STG(Zenithal):
    """
    STG : Stereographic Projection - pixel to sky.
    """

    def _compute_rtheta(self, x, y):
        return np.sqrt(x ** 2 + y ** 2)

    def inverse(self):
        return Sky2Pix_STG()

    def __call__(self, x, y):
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = self._compute_rtheta(x, y)
        theta = 90 - np.rad2deg(2 * np.arctan(rtheta / (2 * self.r0)))
        return phi, theta


class Sky2Pix_STG(Zenithal):
    """
    STG : Stereographic Projection - sky to pixel.
    """

    def inverse(self):
        return Pix2Sky_STG()

    def _compute_rtheta(self, theta):
        return (self.r0 * 2 * np.cos(theta)) / (1 + np.sin(theta))

    @format_input
    def __call__(self, phi, theta, model_set_axis=None):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        rtheta = self._compute_rtheta(theta)
        x = rtheta * np.sin(phi)
        y = - rtheta * np.cos(phi)
        return x, y


class Pix2Sky_SIN(Zenithal):
    """
    SIN : Slant orthographic projection - pixel to sky.
    """

    def inverse(self):
        return Sky2Pix_SIN()

    def _compute_rtheta(self, x, y):
        return np.sqrt(x ** 2 + y ** 2)

    @format_input
    def __call__(self, x, y, model_set_axis=None):
        rtheta = self._compute_rtheta(x, y)
        theta = np.rad2deg(np.arccos(rtheta / self.r0))
        phi = np.rad2deg(np.arctan2(x, -y))
        return phi, theta


class Sky2Pix_SIN(Zenithal):
    """
    SIN : Slant othographic projection - sky to pixel.
    """

    def inverse(self):
        return Pix2Sky_SIN()

    def _compute_rtheta(self, theta):
        return self.r0 * np.cos(theta)

    def __call__(self, phi, theta):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        rtheta = self._compute_rtheta(theta)
        x = rtheta * np.sin(phi)
        y = - rtheta * np.cos(phi)
        return x, y


class Cylindrical(Projection):
    """
    Base class for Cylindrical projections.
    """

    def inverse(self):
        raise NotImplementedError()

    def __call__(self, x, y):
        raise NotImplementedError()


class Pix2Sky_CYP(Cylindrical):
    """
    CYP : Cylindrical perspective - pixel to sky.
    """

    def _validate_mu(mu, model):
        with ignored(AttributeError):
            # An attribute error can occur if model.lam has not been set yet
            if mu == -model.lam:
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return mu

    def _validate_lam(lam, model):
        with ignored(AttributeError):
            # An attribute error can occur if model.lam has not been set yet
            if lam == -model.mu:
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return lam

    mu = Parameter(setter=_validate_mu)
    lam = Parameter(setter=_validate_lam)

    def __init__(self, mu, lam):
        super(Pix2Sky_CYP, self).__init__(mu, lam)

    def inverse(self):
        return Sky2Pix_CYP(self.mu.value, self.lam.value)

    @format_input
    def __call__(self, x, y, model_set_axis=None):
        phi = x / self.lam
        eta = y / (self.r0 * (self.mu + self.lam))
        theta = np.arctan2(eta, 1) + np.arcsin(eta * self.mu /
                                               (np.sqrt(eta ** 2 + 1)))
        return phi, np.rad2deg(theta)


class Sky2Pix_CYP(Cylindrical):
    """
    CYP : Cylindrical Perspective - sky to pixel.
    """

    # TODO: Eliminate duplication on these
    def _validate_mu(mu, model):
        with ignored(AttributeError):
            if mu == -model.lam:
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return mu

    def _validate_lam(lam, model):
        with ignored(AttributeError):
            if lam == -model.mu:
                raise ValueError(
                    "CYP projection is not defined for mu=-lambda")
        return lam

    mu = Parameter(setter=_validate_mu)
    lam = Parameter(setter=_validate_lam)

    def __init__(self, mu, lam):
        super(Sky2Pix_CYP, self).__init__(mu, lam)

    def inverse(self):
        return Pix2Sky_CYP(self.mu, self.lam)

    @format_input
    def __call__(self, phi, theta, model_set_axis=None):
        theta = np.deg2rad(theta)
        x = self.lam * phi
        y = (self.r0 * ((self.mu + self.lam) /
                        (self.mu + np.cos(theta))) * np.sin(theta))
        return x, y


class Pix2Sky_CEA(Cylindrical):
    """
    CEA : Cylindrical equal area projection - pixel to sky.
    """

    lam = Parameter(default=1)

    def __init__(self, lam=lam.default):
        super(Pix2Sky_CEA, self).__init__(lam)

    def inverse(self):
        return Sky2Pix_CEA(self.lam)

    def __call__(self, x, y):
        phi = np.asarray(x)
        theta = np.rad2deg(np.arcsin(1 / self.r0 * self.lam * y))
        return phi, theta


class Sky2Pix_CEA(Cylindrical):
    """
    CEA: Cylindrical equal area projection - sky to pixel.
    """

    lam = Parameter(default=1)

    def __init__(self, lam=lam.default):
        super(Sky2Pix_CEA, self).__init__(lam)

    def inverse(self):
        return Pix2Sky_CEA(self.lam)

    @format_input
    def __call__(self, phi, theta, model_set_axis=None):
        x = phi.copy()
        theta = np.deg2rad(theta)
        y = self.r0 * np.sin(theta) / self.lam
        return x, y


class Pix2Sky_CAR(Cylindrical):
    """
    CAR: Plate carree projection - pixel to sky.
    """

    def inverse(self):
        return Sky2Pix_CAR()

    @format_input
    def __call__(self, x, y, model_set_axis=None):
        return x.copy(), y.copy()

class Sky2Pix_CAR(Cylindrical):
    """
    CAR: Plate carree projection - sky to pixel.
    """

    def inverse(self):
        return Pix2Sky_CAR()

    @format_input
    def __call__(self, phi, theta, model_set_axis=None):
        return phi.copy(), theta.copy()

class Pix2Sky_MER(Cylindrical):
    """
    MER: Mercator - pixel to sky.
    """

    def inverse(self):
        return Sky2Pix_MER()

    @format_input
    def __call__(self, x, y, model_set_axis=None):
        phi = x.copy()
        theta = np.rad2deg(2 * np.arctan(np.e ** (y * np.pi / 180.))) - 90.
        return phi, theta


class Sky2Pix_MER(Cylindrical):
    """
    MER: Mercator - sky to pixel.
    """

    def inverse(self):
        return Pix2Sky_MER()

    @format_input
    def __call__(self, phi, theta, model_set_axis=None):
        x = phi.copy()
        theta = np.deg2rad(theta)
        y = self.r0 * np.log(np.tan((np.pi / 2 + theta) / 2))
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

    n_inputs = 2
    n_outputs = 2
    standard_broadcasting = False

    matrix = Parameter(
        setter=lambda m: AffineTransformation2D._validate_matrix(m),
        default=[[1.0, 0.0], [0.0, 1.0]])
    translation = Parameter(
        setter=lambda t: AffineTransformation2D._validate_vector(t),
        default=[0.0, 0.0])

    def __init__(self, matrix=matrix.default,
                 translation=translation.default):
        super(AffineTransformation2D, self).__init__(matrix, translation)
        self._augmented_matrix = self._create_augmented_matrix(
            self.matrix.value, self.translation.value)

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

    @format_input
    def __call__(self, x, y, model_set_axis=None):
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

        shape = x.shape
        inarr = np.vstack([x.flatten(), y.flatten(), np.ones(x.size)])

        if inarr.shape[0] != 3 or inarr.ndim != 2:
            raise ValueError("Incompatible input shapes")

        result = np.dot(self._augmented_matrix, inarr)

        x, y = result[0], result[1]

        if x.shape != shape:
            x.shape = shape
            y.shape = shape

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
