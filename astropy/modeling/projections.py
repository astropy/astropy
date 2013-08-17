# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Implements sky projections defined in WCS Paper II [1]_

All angles are set and reported in deg but internally the code works and keeps
all angles in radians. For this to work, the mechanism of setting Model's
properties is bypassed by passing an empty parlist to
`~astropy.modeling.core.Model.__init__`.  This also has the effect of not
creating Projection.parameters.  Projection.param_names is created within the
Projection class.  For an example see the AZP classes.

References
----------
.. [1] Calabretta, M.R., Greisen, E.W., 2002, A&A, 395, 1077 (Paper II)
"""

from __future__ import division

import numpy as np

from .core import Model
from . import parameters


projcodes = ['TAN', 'AZP', 'SZP', 'STG', 'SIN', 'ARC', 'ZPN', 'ZEA', 'AIR',
                    'CYP', 'CEA', 'MER']


__all__ = ['Pix2Sky_AZP', 'Sky2Pix_AZP', 'Pix2Sky_CAR', 'Sky2Pix_CAR',
           'Pix2Sky_CEA', 'Sky2Pix_CEA', 'Pix2Sky_CYP', 'Sky2Pix_CYP',
           'Pix2Sky_MER', 'Sky2Pix_MER',
           'Pix2Sky_SIN', 'Sky2Pix_SIN', 'Pix2Sky_STG', 'Sky2Pix_STG',
           'Pix2Sky_TAN', 'Sky2Pix_TAN']


class Projection(Model):
    """
    Base class for all sky projections.

    Parameters
    ----------
    param_names : list of strings
        parameter names
    """

    def __init__(self, param_names):
        super(Projection, self).__init__(param_names, n_inputs=2, n_outputs=2)
        self._pdim = 1
        # the radius of the projection sphere, by which x,y are scaled
        # not sure if it's necessary to make this a user parameter
        self.r0 = 180 / np.pi

    @property
    def n_inputs(self):
        return self._n_inputs

    @property
    def n_outputs(self):
        return self._n_outputs

    @property
    def pdim(self):
        return self._pdim


class Zenithal(Projection):
    """
    Base class for all Zenithal projections.

    Parameters
    ----------
    param_names : list of strings
        parameter names
    """

    def __init__(self, param_names):
        self.phi0 = 0.
        self.theta0 = 90.
        super(Zenithal, self).__init__(param_names)

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

    param_names = ['mu', 'gamma']

    def __init__(self, mu=0., gamma=0.):
        if mu == -1:
            raise ValueError("AZP projection is not defined for mu = -1")
        # units : mu - in spherical radii, gamma - in deg
        self._mu = parameters.Parameter('mu', mu, self, 1)
        self._gamma = parameters.Parameter('gamma', np.deg2rad(gamma), self, 1)
        # bypass initialization of self.parameters
        self._parcheck = {'mu': self.check_mu}
        super(Pix2Sky_AZP, self).__init__(param_names=[])

    @property
    def mu(self):
        return self._mu.value

    @mu.setter
    def mu(self, val):
        self.check_mu(val)
        self._mu.value = val

    @property
    def gamma(self):
        return np.rad2deg(self._gamma.value)

    @gamma.setter
    def gamma(self, val):
        self._gamma.value = np.deg2rad(val)

    def check_mu(self, val):
        if val == -1:
            raise ValueError("AZP projection is not defined for mu=-1")

    def inverse(self):
        return Sky2Pix_AZP(self.mu, self.gamma)

    def _compute_rtheta(self, x, y):
        return np.sqrt(x ** 2 + y ** 2 * (np.cos(self._gamma.value)) ** 2)

    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(np.arctan2(x / np.cos(self._gamma.value), -y))
        r = self._compute_rtheta(x, y)
        pho = r / (self.r0 * (self.mu + 1) + y * np.sin(self._gamma.value))
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

    param_names = ['mu', 'gamma']

    def __init__(self, mu=0., gamma=0.):
        if mu == -1:
            raise ValueError("AZP projections is not defined for mu = -1")
        self._mu = parameters.Parameter('mu', mu, self, 1)
        self._gamma = parameters.Parameter('gamma', np.deg2rad(gamma), self, 1)
        self._parcheck = {'mu': self.check_mu}
        super(Sky2Pix_AZP, self).__init__(param_names=[])

    @property
    def mu(self):
        return self._mu.value

    @mu.setter
    def mu(self, val):
        self.check_mu(val)
        self._mu.value = val

    @property
    def gamma(self):
        return np.rad2deg(self._gamma.value)

    @gamma.setter
    def gamma(self, val):
        self._gamma.value = np.deg2rad(val)

    def check_mu(self, val):
        if val == -1:
            raise ValueError("AZP projection is not defined for mu=-1")

    def inverse(self):
        return Pix2Sky_AZP(self.mu, self.gamma)

    def _compute_rtheta(self, phi, theta):
        rtheta = (self.r0 * (self.mu + 1) *
                  np.cos(theta)) / ((self.mu + np.sin(theta)) +
                                    np.cos(theta) * np.cos(phi) *
                                    np.tan(self._gamma.value))
        return rtheta

    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        r = self._compute_rtheta(phi, theta)
        x = r * np.sin(phi)
        y = (-r * np.cos(phi)) / np.cos(self._gamma.value)
        return x, y


class Pix2Sky_TAN(Zenithal):
    """
    TAN : Gnomonic projection - pixel to sky.
    """

    def __init__(self):
        super(Pix2Sky_TAN, self).__init__(param_names=[])

    def _compute_rtheta(self, x, y):
        return np.sqrt(x ** 2 + y ** 2)

    def inverse(self):
        return Sky2Pix_TAN()

    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = self._compute_rtheta(x, y)
        theta = np.rad2deg(np.arctan2(self.r0, rtheta))
        return phi, theta


class Sky2Pix_TAN(Zenithal):
    """
    TAN : Gnomonic Projection - sky to pixel.
    """

    def __init__(self):
        super(Sky2Pix_TAN, self).__init__(param_names=[])

    def _compute_rtheta(self, theta):
        return 1 / np.tan(theta)

    def inverse(self):
        return Pix2Sky_TAN()

    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        rtheta = self._compute_rtheta(theta)
        x = np.rad2deg(rtheta * np.sin(phi))
        y = - np.rad2deg(rtheta * np.cos(phi))
        return x, y


class Pix2Sky_STG(Zenithal):
    """
    STG : Stereographic Projection - pixel to sky.
    """

    def __init__(self):
        super(Pix2Sky_STG, self).__init__(param_names=[])

    def _compute_rtheta(self, x, y):
        return np.sqrt(x ** 2 + y ** 2)

    def inverse(self):
        return Sky2Pix_STG()

    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = self._compute_rtheta(x, y)
        theta = 90 - np.rad2deg(2 * np.arctan(rtheta / (2 * self.r0)))
        return phi, theta


class Sky2Pix_STG(Zenithal):
    """
    STG : Stereographic Projection - sky to pixel.
    """

    def __init__(self):
        super(Sky2Pix_STG, self).__init__(param_names=[])

    def inverse(self):
        return Pix2Sky_STG()

    def _compute_rtheta(self, theta):
        return (self.r0 * 2 * np.cos(theta)) / (1 + np.sin(theta))

    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        rtheta = self._compute_rtheta(theta)
        x = rtheta * np.sin(phi)
        y = - rtheta * np.cos(phi)
        return x, y


class Pix2Sky_SIN(Zenithal):
    """
    SIN : Slant orthographic projection - pixel to sky.
    """

    def __init__(self):
        super(Pix2Sky_SIN, self).__init__([])

    def inverse(self):
        return Sky2Pix_SIN()

    def _compute_rtheta(self, x, y):
        return np.sqrt(x ** 2 + y ** 2)

    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        rtheta = self._compute_rtheta(x, y)
        theta = np.rad2deg(np.arccos(rtheta / self.r0))
        phi = np.rad2deg(np.arctan2(x, -y))
        return phi, theta


class Sky2Pix_SIN(Zenithal):
    """
    SIN : Slant othographic projection - sky to pixel.
    """

    def __init__(self):
        super(Sky2Pix_SIN, self).__init__([])

    def inverse(self):
        return Pix2Sky_SIN()

    def _compute_rtheta(self, theta):
        return self.r0 * np.cos(theta)

    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        rtheta = self._compute_rtheta(theta)
        x = rtheta * np.sin(phi)
        y = - rtheta * np.cos(phi)
        return x, y


class Cylindrical(Projection):
    """
    Base class for Cylindrical projections.
    """

    # TODO: define param_names
    def __init__(self, param_names):
        self.phi0 = 0.0
        self.theta0 = 0.0
        super(Cylindrical, self).__init__(param_names)

    def inverse(self):
        raise NotImplementedError()

    def __call__(self, x, y):
        raise NotImplementedError()


class Pix2Sky_CYP(Cylindrical):
    """
    CYP : Cylindrical perspective - pixel to sky.
    """

    param_names = ['mu', 'lam']

    def __init__(self, mu, lam):
        if mu == -lam:
            raise ValueError("CYP projection is not defined for mu=-lambda")
        self._mu = parameters.Parameter('mu', mu, self, 1)
        self._lam = parameters.Parameter('lam', lam, self, 1)
        self._parcheck = {'mu': self.check_mu, 'lam': self.check_lam}
        super(Pix2Sky_CYP, self).__init__(param_names=self.param_names)

    def check_mu(self, val):
        if val == -self.lam:
            raise ValueError("CYP projection is not defined for mu = -lambda")

    def check_lam(self, val):
        if val == -self.mu:
            raise ValueError("CYP projection is not defined for mu = -lambda")

    @property
    def mu(self):
        return self._mu.value

    @mu.setter
    def mu(self, val):
        self.check_mu(val)
        self._mu.value = val

    @property
    def lam(self):
        return self._lam.value

    @lam.setter
    def lam(self, val):
        self.check_lam(val)
        self._lam._value = val

    def inverse(self):
        return Sky2Pix_CYP(self.mu.value, self.lam.value)

    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = x / self.lam
        eta = y / (self.r0 * (self.mu + self.lam))
        theta = np.arctan2(eta, 1) + np.arcsin(eta * self.mu /
                                               (np.sqrt(eta ** 2 + 1)))
        return phi, np.rad2deg(theta)


class Sky2Pix_CYP(Cylindrical):
    """
    CYP : Cylindrical Perspective - sky to pixel.
    """

    param_names = ['mu', 'lam']

    def __init__(self, mu, lam):
        if mu == -lam:
            ValueError("CYP projection is not defined for mu = -lambda")
        self._mu = parameters.Parameter('mu', mu, self, 1)
        self._lam = parameters.Parameter('lam', lam, self, 1)
        self._parcheck = {'mu': self.check_mu, 'lam': self.check_lam}
        super(Sky2Pix_CYP, self).__init__(param_names=self.param_names)

    def check_mu(self, val):
        if val == -self.lam:
            raise ValueError("CYP projection is not defined for mu = -lambda")

    def check_lam(self, val):
        if val == -self.mu:
            raise ValueError("CYP projection is not defined for mu = -lambda")

    @property
    def mu(self):
        return self._mu.value

    @mu.setter
    def mu(self, val):
        self.check_mu(val)
        self._mu.value = val

    @property
    def lam(self):
        return self._lam.value

    @lam.setter
    def lam(self, val):
        self.check_lam(val)
        self._lam._value = val

    def inverse(self):
        return Pix2Sky_CYP(self.mu, self.lam)

    def __call__(self, phi, theta):
        theta = np.asarray(np.deg2rad(theta))
        x = self._lam * phi
        y = self.r0 * ((self.mu + self.lam) / (self.mu + np.cos(theta))) * \
            np.sin(theta)
        return x, y


class Pix2Sky_CEA(Cylindrical):
    """
    CEA : Cylindrical equal area projection - pixel to sky.
    """

    param_names = ['lam']

    def __init__(self, lam=1):
        self._lam = parameters.Parameter('lam', lam, self, 1)
        super(Pix2Sky_CEA, self).__init__(param_names=self.param_names)

    @property
    def lam(self):
        return self._lam.value

    @lam.setter
    def lam(self, val):
        self._lam.value = val

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

    param_names = ['lam']

    def __init__(self, lam=1):
        self._lam = parameters.Parameter('lam', lam, self, 1)
        super(Sky2Pix_CEA, self).__init__(param_names=self.param_names)

    def inverse(self):
        return Pix2Sky_CEA(self.lam)

    def __call__(self, phi, theta):
        x = np.asarray(phi, dtype=np.float)
        theta = np.asarray(np.deg2rad(theta), dtype=np.float)
        y = self.r0 * np.sin(theta) / self._lam
        return x, y


class Pix2Sky_CAR(Cylindrical):
    """
    CAR: Plate carree projection - pixel to sky.
    """

    def __init__(self):
        super(Pix2Sky_CAR, self).__init__([])

    def inverse(self):
        return Sky2Pix_CAR()

    def __call__(self, x, y):
        phi = np.asarray(x, dtype=np.float)
        theta = np.asarray(y, dtype=np.float)
        return phi, theta


class Sky2Pix_CAR(Cylindrical):
    """
    CAR: Plate carree projection - sky to pixel.
    """

    def __init__(self):
        super(Sky2Pix_CAR, self).__init__([])

    def inverse(self):
        return Pix2Sky_CAR()

    def __call__(self, phi, theta):
        x = np.asarray(phi, dtype=np.float)
        y = np.asarray(theta, dtype=np.float)
        return x, y


class Pix2Sky_MER(Cylindrical):
    """
    MER: Mercator - pixel to sky.
    """

    def __init__(self):
        super(Pix2Sky_MER, self).__init__([])

    def inverse(self):
        return Sky2Pix_MER()

    def __call__(self, x, y):
        phi = np.asarray(x, dtype=np.float)
        theta = np.asarray(np.rad2deg(2 * np.arctan(np.e ** (y * np.pi / 180.))) - 90.)
        return phi, theta


class Sky2Pix_MER(Cylindrical):
    """
    MER: Mercator - sky to pixel.
    """

    def __init__(self):
        super(Sky2Pix_MER, self).__init__([])

    def inverse(self):
        return Pix2Sky_MER()

    def __call__(self, phi, theta):
        x = np.asarray(phi, dtype=np.float)
        theta = np.asarray(np.deg2rad(theta))
        y = self.r0 * np.log(np.tan((np.pi / 2 + theta) / 2))
        return x, y
