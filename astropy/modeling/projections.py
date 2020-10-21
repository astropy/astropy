# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
Implements projections--particularly sky projections defined in WCS Paper II
[1]_.

All angles are set and and displayed in degrees but internally computations are
performed in radians. All functions expect inputs and outputs degrees.

References
----------
.. [1] Calabretta, M.R., Greisen, E.W., 2002, A&A, 395, 1077 (Paper II)
"""


import abc

import numpy as np

from astropy import units as u

from .core import Model
from .parameters import Parameter, InputParameterError
from . import _projections
from .utils import _to_radian, _to_orig_unit


projcodes = [
    'AZP', 'SZP', 'TAN', 'STG', 'SIN', 'ARC', 'ZEA', 'AIR', 'CYP',
    'CEA', 'CAR', 'MER', 'SFL', 'PAR', 'MOL', 'AIT', 'COP', 'COE',
    'COD', 'COO', 'BON', 'PCO', 'TSC', 'CSC', 'QSC', 'HPX', 'XPH'
]


__all__ = ['Projection', 'Pix2SkyProjection', 'Sky2PixProjection',
           'Zenithal', 'Cylindrical', 'PseudoCylindrical', 'Conic',
           'PseudoConic', 'QuadCube', 'HEALPix',
           'AffineTransformation2D',
           'projcodes',

           'Pix2Sky_ZenithalPerspective', 'Sky2Pix_ZenithalPerspective',
           'Pix2Sky_SlantZenithalPerspective', 'Sky2Pix_SlantZenithalPerspective',
           'Pix2Sky_Gnomonic', 'Sky2Pix_Gnomonic',
           'Pix2Sky_Stereographic', 'Sky2Pix_Stereographic',
           'Pix2Sky_SlantOrthographic', 'Sky2Pix_SlantOrthographic',
           'Pix2Sky_ZenithalEquidistant', 'Sky2Pix_ZenithalEquidistant',
           'Pix2Sky_ZenithalEqualArea', 'Sky2Pix_ZenithalEqualArea',
           'Pix2Sky_Airy', 'Sky2Pix_Airy',
           'Pix2Sky_CylindricalPerspective', 'Sky2Pix_CylindricalPerspective',
           'Pix2Sky_CylindricalEqualArea', 'Sky2Pix_CylindricalEqualArea',
           'Pix2Sky_PlateCarree', 'Sky2Pix_PlateCarree',
           'Pix2Sky_Mercator', 'Sky2Pix_Mercator',
           'Pix2Sky_SansonFlamsteed', 'Sky2Pix_SansonFlamsteed',
           'Pix2Sky_Parabolic', 'Sky2Pix_Parabolic',
           'Pix2Sky_Molleweide', 'Sky2Pix_Molleweide',
           'Pix2Sky_HammerAitoff', 'Sky2Pix_HammerAitoff',
           'Pix2Sky_ConicPerspective', 'Sky2Pix_ConicPerspective',
           'Pix2Sky_ConicEqualArea', 'Sky2Pix_ConicEqualArea',
           'Pix2Sky_ConicEquidistant', 'Sky2Pix_ConicEquidistant',
           'Pix2Sky_ConicOrthomorphic', 'Sky2Pix_ConicOrthomorphic',
           'Pix2Sky_BonneEqualArea', 'Sky2Pix_BonneEqualArea',
           'Pix2Sky_Polyconic', 'Sky2Pix_Polyconic',
           'Pix2Sky_TangentialSphericalCube', 'Sky2Pix_TangentialSphericalCube',
           'Pix2Sky_COBEQuadSphericalCube', 'Sky2Pix_COBEQuadSphericalCube',
           'Pix2Sky_QuadSphericalCube', 'Sky2Pix_QuadSphericalCube',
           'Pix2Sky_HEALPix', 'Sky2Pix_HEALPix',
           'Pix2Sky_HEALPixPolar', 'Sky2Pix_HEALPixPolar',

           # The following are short FITS WCS aliases
           'Pix2Sky_AZP', 'Sky2Pix_AZP',
           'Pix2Sky_SZP', 'Sky2Pix_SZP',
           'Pix2Sky_TAN', 'Sky2Pix_TAN',
           'Pix2Sky_STG', 'Sky2Pix_STG',
           'Pix2Sky_SIN', 'Sky2Pix_SIN',
           'Pix2Sky_ARC', 'Sky2Pix_ARC',
           'Pix2Sky_ZEA', 'Sky2Pix_ZEA',
           'Pix2Sky_AIR', 'Sky2Pix_AIR',
           'Pix2Sky_CYP', 'Sky2Pix_CYP',
           'Pix2Sky_CEA', 'Sky2Pix_CEA',
           'Pix2Sky_CAR', 'Sky2Pix_CAR',
           'Pix2Sky_MER', 'Sky2Pix_MER',
           'Pix2Sky_SFL', 'Sky2Pix_SFL',
           'Pix2Sky_PAR', 'Sky2Pix_PAR',
           'Pix2Sky_MOL', 'Sky2Pix_MOL',
           'Pix2Sky_AIT', 'Sky2Pix_AIT',
           'Pix2Sky_COP', 'Sky2Pix_COP',
           'Pix2Sky_COE', 'Sky2Pix_COE',
           'Pix2Sky_COD', 'Sky2Pix_COD',
           'Pix2Sky_COO', 'Sky2Pix_COO',
           'Pix2Sky_BON', 'Sky2Pix_BON',
           'Pix2Sky_PCO', 'Sky2Pix_PCO',
           'Pix2Sky_TSC', 'Sky2Pix_TSC',
           'Pix2Sky_CSC', 'Sky2Pix_CSC',
           'Pix2Sky_QSC', 'Sky2Pix_QSC',
           'Pix2Sky_HPX', 'Sky2Pix_HPX',
           'Pix2Sky_XPH', 'Sky2Pix_XPH'
           ]


class Projection(Model):
    """Base class for all sky projections."""

    # Radius of the generating sphere.
    # This sets the circumference to 360 deg so that arc length is measured in deg.
    r0 = 180 * u.deg / np.pi

    _separable = False

    @property
    @abc.abstractmethod
    def inverse(self):
        """
        Inverse projection--all projection models must provide an inverse.
        """


class Pix2SkyProjection(Projection):
    """Base class for all Pix2Sky projections."""

    n_inputs = 2
    n_outputs = 2

    _input_units_strict = True
    _input_units_allow_dimensionless = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.inputs = ('x', 'y')
        self.outputs = ('phi', 'theta')

    @property
    def input_units(self):
        return {self.inputs[0]: u.deg,
                self.inputs[1]: u.deg}

    @property
    def return_units(self):
        return {self.outputs[0]: u.deg,
                self.outputs[1]: u.deg}


class Sky2PixProjection(Projection):
    """Base class for all Sky2Pix projections."""

    n_inputs = 2
    n_outputs = 2

    _input_units_strict = True
    _input_units_allow_dimensionless = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.inputs = ('phi', 'theta')
        self.outputs = ('x', 'y')

    @property
    def input_units(self):
        return {self.inputs[0]: u.deg,
                self.inputs[1]: u.deg}

    @property
    def return_units(self):
        return {self.outputs[0]: u.deg,
                self.outputs[1]: u.deg}


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

    _separable = False


class Pix2Sky_ZenithalPerspective(Pix2SkyProjection, Zenithal):
    r"""
    Zenithal perspective projection - pixel to sky.

    Corresponds to the ``AZP`` projection in FITS WCS.

    .. math::
        \phi &= \arg(-y \cos \gamma, x) \\
        \theta &= \left\{\genfrac{}{}{0pt}{}{\psi - \omega}{\psi + \omega + 180^{\circ}}\right.

    where:

    .. math::
        \psi &= \arg(\rho, 1) \\
        \omega &= \sin^{-1}\left(\frac{\rho \mu}{\sqrt{\rho^2 + 1}}\right) \\
        \rho &= \frac{R}{\frac{180^{\circ}}{\pi}(\mu + 1) + y \sin \gamma} \\
        R &= \sqrt{x^2 + y^2 \cos^2 \gamma}

    Parameters
    --------------
    mu : float
        Distance from point of projection to center of sphere
        in spherical radii, μ.  Default is 0.

    gamma : float
        Look angle γ in degrees.  Default is 0°.
    """

    mu = Parameter(default=0.0)
    gamma = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)

    def __init__(self, mu=mu.default, gamma=gamma.default, **kwargs):
        # units : mu - in spherical radii, gamma - in deg
        super().__init__(mu, gamma, **kwargs)

    @mu.validator
    def mu(self, value):
        if np.any(value == -1):
            raise InputParameterError(
                "Zenithal perspective projection is not defined for mu = -1")

    @property
    def inverse(self):
        return Sky2Pix_ZenithalPerspective(self.mu.value, self.gamma.value)

    @classmethod
    def evaluate(cls, x, y, mu, gamma):
        return _projections.azpx2s(x, y, mu, _to_orig_unit(gamma))


Pix2Sky_AZP = Pix2Sky_ZenithalPerspective


class Sky2Pix_ZenithalPerspective(Sky2PixProjection, Zenithal):
    r"""
    Zenithal perspective projection - sky to pixel.

    Corresponds to the ``AZP`` projection in FITS WCS.

    .. math::
        x &= R \sin \phi \\
        y &= -R \sec \gamma \cos \theta

    where:

    .. math::
        R = \frac{180^{\circ}}{\pi} \frac{(\mu + 1) \cos \theta}{(\mu + \sin \theta) + \cos \theta \cos \phi \tan \gamma}

    Parameters
    ----------
    mu : float
        Distance from point of projection to center of sphere
        in spherical radii, μ. Default is 0.

    gamma : float
        Look angle γ in degrees. Default is 0°.
    """

    mu = Parameter(default=0.0)
    gamma = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)

    @mu.validator
    def mu(self, value):
        if np.any(value == -1):
            raise InputParameterError(
                "Zenithal perspective projection is not defined for mu = -1")

    @property
    def inverse(self):
        return Pix2Sky_AZP(self.mu.value, self.gamma.value)

    @classmethod
    def evaluate(cls, phi, theta, mu, gamma):
        return _projections.azps2x(
            phi, theta, mu, _to_orig_unit(gamma))


Sky2Pix_AZP = Sky2Pix_ZenithalPerspective


class Pix2Sky_SlantZenithalPerspective(Pix2SkyProjection, Zenithal):
    r"""
    Slant zenithal perspective projection - pixel to sky.

    Corresponds to the ``SZP`` projection in FITS WCS.

    Parameters
    --------------
    mu : float
        Distance from point of projection to center of sphere
        in spherical radii, μ.  Default is 0.

    phi0 : float
        The longitude φ₀ of the reference point, in degrees.  Default
        is 0°.

    theta0 : float
        The latitude θ₀ of the reference point, in degrees.  Default
        is 90°.
    """

    def _validate_mu(mu):
        if np.asarray(mu == -1).any():
            raise ValueError(
                "Zenithal perspective projection is not defined for mu=-1")
        return mu

    mu = Parameter(default=0.0, setter=_validate_mu)
    phi0 = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)
    theta0 = Parameter(default=90.0, getter=_to_orig_unit, setter=_to_radian)

    @property
    def inverse(self):
        return Sky2Pix_SlantZenithalPerspective(
            self.mu.value, self.phi0.value, self.theta0.value)

    @classmethod
    def evaluate(cls, x, y, mu, phi0, theta0):
        return _projections.szpx2s(
            x, y, mu, _to_orig_unit(phi0), _to_orig_unit(theta0))


Pix2Sky_SZP = Pix2Sky_SlantZenithalPerspective


class Sky2Pix_SlantZenithalPerspective(Sky2PixProjection, Zenithal):
    r"""
    Zenithal perspective projection - sky to pixel.

    Corresponds to the ``SZP`` projection in FITS WCS.

    Parameters
    ----------
    mu : float
        distance from point of projection to center of sphere
        in spherical radii, μ.  Default is 0.

    phi0 : float
        The longitude φ₀ of the reference point, in degrees.  Default
        is 0°.

    theta0 : float
        The latitude θ₀ of the reference point, in degrees.  Default
        is 90°.
    """

    def _validate_mu(mu):
        if np.asarray(mu == -1).any():
            raise ValueError("Zenithal perspective projection is not defined for mu=-1")
        return mu

    mu = Parameter(default=0.0, setter=_validate_mu)
    phi0 = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)
    theta0 = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)

    @property
    def inverse(self):
        return Pix2Sky_SlantZenithalPerspective(
            self.mu.value, self.phi0.value, self.theta0.value)

    @classmethod
    def evaluate(cls, phi, theta, mu, phi0, theta0):
        return _projections.szps2x(
            phi, theta, mu, _to_orig_unit(phi0), _to_orig_unit(theta0))


Sky2Pix_SZP = Sky2Pix_SlantZenithalPerspective


class Pix2Sky_Gnomonic(Pix2SkyProjection, Zenithal):
    r"""
    Gnomonic projection - pixel to sky.

    Corresponds to the ``TAN`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = \tan^{-1}\left(\frac{180^{\circ}}{\pi R_\theta}\right)
    """

    @property
    def inverse(self):
        return Sky2Pix_Gnomonic()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.tanx2s(x, y)


Pix2Sky_TAN = Pix2Sky_Gnomonic


class Sky2Pix_Gnomonic(Sky2PixProjection, Zenithal):
    r"""
    Gnomonic Projection - sky to pixel.

    Corresponds to the ``TAN`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\cot \theta
    """

    @property
    def inverse(self):
        return Pix2Sky_Gnomonic()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.tans2x(phi, theta)


Sky2Pix_TAN = Sky2Pix_Gnomonic


class Pix2Sky_Stereographic(Pix2SkyProjection, Zenithal):
    r"""
    Stereographic Projection - pixel to sky.

    Corresponds to the ``STG`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = 90^{\circ} - 2 \tan^{-1}\left(\frac{\pi R_\theta}{360^{\circ}}\right)
    """

    @property
    def inverse(self):
        return Sky2Pix_Stereographic()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.stgx2s(x, y)


Pix2Sky_STG = Pix2Sky_Stereographic


class Sky2Pix_Stereographic(Sky2PixProjection, Zenithal):
    r"""
    Stereographic Projection - sky to pixel.

    Corresponds to the ``STG`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\frac{2 \cos \theta}{1 + \sin \theta}
    """

    @property
    def inverse(self):
        return Pix2Sky_Stereographic()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.stgs2x(phi, theta)


Sky2Pix_STG = Sky2Pix_Stereographic


class Pix2Sky_SlantOrthographic(Pix2SkyProjection, Zenithal):
    r"""
    Slant orthographic projection - pixel to sky.

    Corresponds to the ``SIN`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    The following transformation applies when :math:`\xi` and
    :math:`\eta` are both zero.

    .. math::
        \theta = \cos^{-1}\left(\frac{\pi}{180^{\circ}}R_\theta\right)

    The parameters :math:`\xi` and :math:`\eta` are defined from the
    reference point :math:`(\phi_c, \theta_c)` as:

    .. math::
        \xi &= \cot \theta_c \sin \phi_c \\
        \eta &= - \cot \theta_c \cos \phi_c

    Parameters
    ----------
    xi : float
        Obliqueness parameter, ξ.  Default is 0.0.

    eta : float
        Obliqueness parameter, η.  Default is 0.0.
    """

    xi = Parameter(default=0.0)
    eta = Parameter(default=0.0)

    @property
    def inverse(self):
        return Sky2Pix_SlantOrthographic(self.xi.value, self.eta.value)

    @classmethod
    def evaluate(cls, x, y, xi, eta):
        return _projections.sinx2s(x, y, xi, eta)


Pix2Sky_SIN = Pix2Sky_SlantOrthographic


class Sky2Pix_SlantOrthographic(Sky2PixProjection, Zenithal):
    r"""
    Slant orthographic projection - sky to pixel.

    Corresponds to the ``SIN`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    The following transformation applies when :math:`\xi` and
    :math:`\eta` are both zero.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\cos \theta

    But more specifically are:

    .. math::
        x &= \frac{180^\circ}{\pi}[\cos \theta \sin \phi + \xi(1 - \sin \theta)] \\
        y &= \frac{180^\circ}{\pi}[\cos \theta \cos \phi + \eta(1 - \sin \theta)]
    """

    xi = Parameter(default=0.0)
    eta = Parameter(default=0.0)

    @property
    def inverse(self):
        return Pix2Sky_SlantOrthographic(self.xi.value, self.eta.value)

    @classmethod
    def evaluate(cls, phi, theta, xi, eta):
        return _projections.sins2x(phi, theta, xi, eta)


Sky2Pix_SIN = Sky2Pix_SlantOrthographic


class Pix2Sky_ZenithalEquidistant(Pix2SkyProjection, Zenithal):
    r"""
    Zenithal equidistant projection - pixel to sky.

    Corresponds to the ``ARC`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = 90^\circ - R_\theta
    """
    @property
    def inverse(self):
        return Sky2Pix_ZenithalEquidistant()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.arcx2s(x, y)


Pix2Sky_ARC = Pix2Sky_ZenithalEquidistant


class Sky2Pix_ZenithalEquidistant(Sky2PixProjection, Zenithal):
    r"""
    Zenithal equidistant projection - sky to pixel.

    Corresponds to the ``ARC`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = 90^\circ - \theta
    """
    @property
    def inverse(self):
        return Pix2Sky_ZenithalEquidistant()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.arcs2x(phi, theta)


Sky2Pix_ARC = Sky2Pix_ZenithalEquidistant


class Pix2Sky_ZenithalEqualArea(Pix2SkyProjection, Zenithal):
    r"""
    Zenithal equidistant projection - pixel to sky.

    Corresponds to the ``ZEA`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = 90^\circ - 2 \sin^{-1} \left(\frac{\pi R_\theta}{360^\circ}\right)
    """
    @property
    def inverse(self):
        return Sky2Pix_ZenithalEqualArea()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.zeax2s(x, y)


Pix2Sky_ZEA = Pix2Sky_ZenithalEqualArea


class Sky2Pix_ZenithalEqualArea(Sky2PixProjection, Zenithal):
    r"""
    Zenithal equidistant projection - sky to pixel.

    Corresponds to the ``ZEA`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta &= \frac{180^\circ}{\pi} \sqrt{2(1 - \sin\theta)} \\
                 &= \frac{360^\circ}{\pi} \sin\left(\frac{90^\circ - \theta}{2}\right)
    """
    @property
    def inverse(self):
        return Pix2Sky_ZenithalEqualArea()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.zeas2x(phi, theta)


Sky2Pix_ZEA = Sky2Pix_ZenithalEqualArea


class Pix2Sky_Airy(Pix2SkyProjection, Zenithal):
    r"""
    Airy projection - pixel to sky.

    Corresponds to the ``AIR`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    Parameters
    ----------
    theta_b : float
        The latitude :math:`\theta_b` at which to minimize the error,
        in degrees.  Default is 90°.
    """
    theta_b = Parameter(default=90.0)

    @property
    def inverse(self):
        return Sky2Pix_Airy(self.theta_b.value)

    @classmethod
    def evaluate(cls, x, y, theta_b):
        return _projections.airx2s(x, y, theta_b)


Pix2Sky_AIR = Pix2Sky_Airy


class Sky2Pix_Airy(Sky2PixProjection, Zenithal):
    r"""
    Airy - sky to pixel.

    Corresponds to the ``AIR`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = -2 \frac{180^\circ}{\pi}\left(\frac{\ln(\cos \xi)}{\tan \xi} + \frac{\ln(\cos \xi_b)}{\tan^2 \xi_b} \tan \xi \right)

    where:

    .. math::
        \xi &= \frac{90^\circ - \theta}{2} \\
        \xi_b &= \frac{90^\circ - \theta_b}{2}

    Parameters
    ----------
    theta_b : float
        The latitude :math:`\theta_b` at which to minimize the error,
        in degrees.  Default is 90°.
    """
    theta_b = Parameter(default=90.0)

    @property
    def inverse(self):
        return Pix2Sky_Airy(self.theta_b.value)

    @classmethod
    def evaluate(cls, phi, theta, theta_b):
        return _projections.airs2x(phi, theta, theta_b)


Sky2Pix_AIR = Sky2Pix_Airy


class Cylindrical(Projection):
    r"""Base class for Cylindrical projections.

    Cylindrical projections are so-named because the surface of
    projection is a cylinder.
    """
    _separable = True


class Pix2Sky_CylindricalPerspective(Pix2SkyProjection, Cylindrical):
    r"""
    Cylindrical perspective - pixel to sky.

    Corresponds to the ``CYP`` projection in FITS WCS.

    .. math::
        \phi &= \frac{x}{\lambda} \\
        \theta &= \arg(1, \eta) + \sin{-1}\left(\frac{\eta \mu}{\sqrt{\eta^2 + 1}}\right)

    where:

    .. math::
        \eta = \frac{\pi}{180^{\circ}}\frac{y}{\mu + \lambda}

    Parameters
    ----------
    mu : float
        Distance from center of sphere in the direction opposite the
        projected surface, in spherical radii, μ. Default is 1.

    lam : float
        Radius of the cylinder in spherical radii, λ. Default is 1.
    """

    mu = Parameter(default=1.0)
    lam = Parameter(default=1.0)

    @mu.validator
    def mu(self, value):
        if np.any(value == -self.lam):
            raise InputParameterError(
                "CYP projection is not defined for mu = -lambda")

    @lam.validator
    def lam(self, value):
        if np.any(value == -self.mu):
            raise InputParameterError(
                "CYP projection is not defined for lambda = -mu")

    @property
    def inverse(self):
        return Sky2Pix_CylindricalPerspective(self.mu.value, self.lam.value)

    @classmethod
    def evaluate(cls, x, y, mu, lam):
        return _projections.cypx2s(x, y, mu, lam)


Pix2Sky_CYP = Pix2Sky_CylindricalPerspective


class Sky2Pix_CylindricalPerspective(Sky2PixProjection, Cylindrical):
    r"""
    Cylindrical Perspective - sky to pixel.

    Corresponds to the ``CYP`` projection in FITS WCS.

    .. math::
        x &= \lambda \phi \\
        y &= \frac{180^{\circ}}{\pi}\left(\frac{\mu + \lambda}{\mu + \cos \theta}\right)\sin \theta

    Parameters
    ----------
    mu : float
        Distance from center of sphere in the direction opposite the
        projected surface, in spherical radii, μ.  Default is 0.

    lam : float
        Radius of the cylinder in spherical radii, λ.  Default is 0.
    """

    mu = Parameter(default=1.0)
    lam = Parameter(default=1.0)

    @mu.validator
    def mu(self, value):
        if np.any(value == -self.lam):
            raise InputParameterError(
                "CYP projection is not defined for mu = -lambda")

    @lam.validator
    def lam(self, value):
        if np.any(value == -self.mu):
            raise InputParameterError(
                "CYP projection is not defined for lambda = -mu")

    @property
    def inverse(self):
        return Pix2Sky_CylindricalPerspective(self.mu, self.lam)

    @classmethod
    def evaluate(cls, phi, theta, mu, lam):
        return _projections.cyps2x(phi, theta, mu, lam)


Sky2Pix_CYP = Sky2Pix_CylindricalPerspective


class Pix2Sky_CylindricalEqualArea(Pix2SkyProjection, Cylindrical):
    r"""
    Cylindrical equal area projection - pixel to sky.

    Corresponds to the ``CEA`` projection in FITS WCS.

    .. math::
        \phi &= x \\
        \theta &= \sin^{-1}\left(\frac{\pi}{180^{\circ}}\lambda y\right)

    Parameters
    ----------
    lam : float
        Radius of the cylinder in spherical radii, λ.  Default is 1.
    """

    lam = Parameter(default=1)

    @property
    def inverse(self):
        return Sky2Pix_CylindricalEqualArea(self.lam)

    @classmethod
    def evaluate(cls, x, y, lam):
        return _projections.ceax2s(x, y, lam)


Pix2Sky_CEA = Pix2Sky_CylindricalEqualArea


class Sky2Pix_CylindricalEqualArea(Sky2PixProjection, Cylindrical):
    r"""
    Cylindrical equal area projection - sky to pixel.

    Corresponds to the ``CEA`` projection in FITS WCS.

    .. math::
        x &= \phi \\
        y &= \frac{180^{\circ}}{\pi}\frac{\sin \theta}{\lambda}

    Parameters
    ----------
    lam : float
        Radius of the cylinder in spherical radii, λ.  Default is 0.
    """

    lam = Parameter(default=1)

    @property
    def inverse(self):
        return Pix2Sky_CylindricalEqualArea(self.lam)

    @classmethod
    def evaluate(cls, phi, theta, lam):
        return _projections.ceas2x(phi, theta, lam)


Sky2Pix_CEA = Sky2Pix_CylindricalEqualArea


class Pix2Sky_PlateCarree(Pix2SkyProjection, Cylindrical):
    r"""
    Plate carrée projection - pixel to sky.

    Corresponds to the ``CAR`` projection in FITS WCS.

    .. math::
        \phi &= x \\
        \theta &= y
    """

    @property
    def inverse(self):
        return Sky2Pix_PlateCarree()

    @staticmethod
    def evaluate(x, y):
        # The intermediate variables are only used here for clarity
        phi = np.array(x, copy=True)
        theta = np.array(y, copy=True)

        return phi, theta


Pix2Sky_CAR = Pix2Sky_PlateCarree


class Sky2Pix_PlateCarree(Sky2PixProjection, Cylindrical):
    r"""
    Plate carrée projection - sky to pixel.

    Corresponds to the ``CAR`` projection in FITS WCS.

    .. math::
        x &= \phi \\
        y &= \theta
    """

    @property
    def inverse(self):
        return Pix2Sky_PlateCarree()

    @staticmethod
    def evaluate(phi, theta):
        # The intermediate variables are only used here for clarity
        x = np.array(phi, copy=True)
        y = np.array(theta, copy=True)

        return x, y


Sky2Pix_CAR = Sky2Pix_PlateCarree


class Pix2Sky_Mercator(Pix2SkyProjection, Cylindrical):
    r"""
    Mercator - pixel to sky.

    Corresponds to the ``MER`` projection in FITS WCS.

    .. math::
        \phi &= x \\
        \theta &= 2 \tan^{-1}\left(e^{y \pi / 180^{\circ}}\right)-90^{\circ}
    """

    @property
    def inverse(self):
        return Sky2Pix_Mercator()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.merx2s(x, y)


Pix2Sky_MER = Pix2Sky_Mercator


class Sky2Pix_Mercator(Sky2PixProjection, Cylindrical):
    r"""
    Mercator - sky to pixel.

    Corresponds to the ``MER`` projection in FITS WCS.

    .. math::
        x &= \phi \\
        y &= \frac{180^{\circ}}{\pi}\ln \tan \left(\frac{90^{\circ} + \theta}{2}\right)
    """

    @property
    def inverse(self):
        return Pix2Sky_Mercator()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.mers2x(phi, theta)


Sky2Pix_MER = Sky2Pix_Mercator


class PseudoCylindrical(Projection):
    r"""Base class for pseudocylindrical projections.

    Pseudocylindrical projections are like cylindrical projections
    except the parallels of latitude are projected at diminishing
    lengths toward the polar regions in order to reduce lateral
    distortion there.  Consequently, the meridians are curved.
    """

    _separable = True


class Pix2Sky_SansonFlamsteed(Pix2SkyProjection, PseudoCylindrical):
    r"""
    Sanson-Flamsteed projection - pixel to sky.

    Corresponds to the ``SFL`` projection in FITS WCS.

    .. math::
        \phi &= \frac{x}{\cos y} \\
        \theta &= y
    """

    @property
    def inverse(self):
        return Sky2Pix_SansonFlamsteed()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.sflx2s(x, y)


Pix2Sky_SFL = Pix2Sky_SansonFlamsteed


class Sky2Pix_SansonFlamsteed(Sky2PixProjection, PseudoCylindrical):
    r"""
    Sanson-Flamsteed projection - sky to pixel.

    Corresponds to the ``SFL`` projection in FITS WCS.

    .. math::
        x &= \phi \cos \theta \\
        y &= \theta
    """

    @property
    def inverse(self):
        return Pix2Sky_SansonFlamsteed()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.sfls2x(phi, theta)


Sky2Pix_SFL = Sky2Pix_SansonFlamsteed


class Pix2Sky_Parabolic(Pix2SkyProjection, PseudoCylindrical):
    r"""
    Parabolic projection - pixel to sky.

    Corresponds to the ``PAR`` projection in FITS WCS.

    .. math::
        \phi &= \frac{180^\circ}{\pi} \frac{x}{1 - 4(y / 180^\circ)^2} \\
        \theta &= 3 \sin^{-1}\left(\frac{y}{180^\circ}\right)
    """

    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_Parabolic()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.parx2s(x, y)


Pix2Sky_PAR = Pix2Sky_Parabolic


class Sky2Pix_Parabolic(Sky2PixProjection, PseudoCylindrical):
    r"""
    Parabolic projection - sky to pixel.

    Corresponds to the ``PAR`` projection in FITS WCS.

    .. math::
        x &= \phi \left(2\cos\frac{2\theta}{3} - 1\right) \\
        y &= 180^\circ \sin \frac{\theta}{3}
    """

    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_Parabolic()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.pars2x(phi, theta)


Sky2Pix_PAR = Sky2Pix_Parabolic


class Pix2Sky_Molleweide(Pix2SkyProjection, PseudoCylindrical):
    r"""
    Molleweide's projection - pixel to sky.

    Corresponds to the ``MOL`` projection in FITS WCS.

    .. math::
        \phi &= \frac{\pi x}{2 \sqrt{2 - \left(\frac{\pi}{180^\circ}y\right)^2}} \\
        \theta &= \sin^{-1}\left(\frac{1}{90^\circ}\sin^{-1}\left(\frac{\pi}{180^\circ}\frac{y}{\sqrt{2}}\right) + \frac{y}{180^\circ}\sqrt{2 - \left(\frac{\pi}{180^\circ}y\right)^2}\right)
    """

    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_Molleweide()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.molx2s(x, y)


Pix2Sky_MOL = Pix2Sky_Molleweide


class Sky2Pix_Molleweide(Sky2PixProjection, PseudoCylindrical):
    r"""
    Molleweide's projection - sky to pixel.

    Corresponds to the ``MOL`` projection in FITS WCS.

    .. math::
        x &= \frac{2 \sqrt{2}}{\pi} \phi \cos \gamma \\
        y &= \sqrt{2} \frac{180^\circ}{\pi} \sin \gamma

    where :math:`\gamma` is defined as the solution of the
    transcendental equation:

    .. math::

        \sin \theta = \frac{\gamma}{90^\circ} + \frac{\sin 2 \gamma}{\pi}
    """

    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_Molleweide()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.mols2x(phi, theta)


Sky2Pix_MOL = Sky2Pix_Molleweide


class Pix2Sky_HammerAitoff(Pix2SkyProjection, PseudoCylindrical):
    r"""
    Hammer-Aitoff projection - pixel to sky.

    Corresponds to the ``AIT`` projection in FITS WCS.

    .. math::
        \phi &= 2 \arg \left(2Z^2 - 1, \frac{\pi}{180^\circ} \frac{Z}{2}x\right) \\
        \theta &= \sin^{-1}\left(\frac{\pi}{180^\circ}yZ\right)
    """

    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_HammerAitoff()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.aitx2s(x, y)


Pix2Sky_AIT = Pix2Sky_HammerAitoff


class Sky2Pix_HammerAitoff(Sky2PixProjection, PseudoCylindrical):
    r"""
    Hammer-Aitoff projection - sky to pixel.

    Corresponds to the ``AIT`` projection in FITS WCS.

    .. math::
        x &= 2 \gamma \cos \theta \sin \frac{\phi}{2} \\
        y &= \gamma \sin \theta

    where:

    .. math::
        \gamma = \frac{180^\circ}{\pi} \sqrt{\frac{2}{1 + \cos \theta \cos(\phi / 2)}}
    """

    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_HammerAitoff()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.aits2x(phi, theta)


Sky2Pix_AIT = Sky2Pix_HammerAitoff


class Conic(Projection):
    r"""Base class for conic projections.

    In conic projections, the sphere is thought to be projected onto
    the surface of a cone which is then opened out.

    In a general sense, the pixel-to-sky transformation is defined as:

    .. math::

        \phi &= \arg\left(\frac{Y_0 - y}{R_\theta}, \frac{x}{R_\theta}\right) / C \\
        R_\theta &= \mathrm{sign} \theta_a \sqrt{x^2 + (Y_0 - y)^2}

    and the inverse (sky-to-pixel) is defined as:

    .. math::
        x &= R_\theta \sin (C \phi) \\
        y &= R_\theta \cos (C \phi) + Y_0

    where :math:`C` is the "constant of the cone":

    .. math::
        C = \frac{180^\circ \cos \theta}{\pi R_\theta}
    """
    sigma = Parameter(default=90.0, getter=_to_orig_unit, setter=_to_radian)
    delta = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)

    _separable = False


class Pix2Sky_ConicPerspective(Pix2SkyProjection, Conic):
    r"""
    Colles' conic perspective projection - pixel to sky.

    Corresponds to the ``COP`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::
        C &= \sin \theta_a \\
        R_\theta &= \frac{180^\circ}{\pi} \cos \eta [ \cot \theta_a - \tan(\theta - \theta_a)] \\
        Y_0 &= \frac{180^\circ}{\pi} \cos \eta \cot \theta_a

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Sky2Pix_ConicPerspective(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, x, y, sigma, delta):
        return _projections.copx2s(x, y, _to_orig_unit(sigma), _to_orig_unit(delta))


Pix2Sky_COP = Pix2Sky_ConicPerspective


class Sky2Pix_ConicPerspective(Sky2PixProjection, Conic):
    r"""
    Colles' conic perspective projection - sky to pixel.

    Corresponds to the ``COP`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::
        C &= \sin \theta_a \\
        R_\theta &= \frac{180^\circ}{\pi} \cos \eta [ \cot \theta_a - \tan(\theta - \theta_a)] \\
        Y_0 &= \frac{180^\circ}{\pi} \cos \eta \cot \theta_a

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Pix2Sky_ConicPerspective(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, phi, theta, sigma, delta):
        return _projections.cops2x(phi, theta,
                                   _to_orig_unit(sigma), _to_orig_unit(delta))


Sky2Pix_COP = Sky2Pix_ConicPerspective


class Pix2Sky_ConicEqualArea(Pix2SkyProjection, Conic):
    r"""
    Alber's conic equal area projection - pixel to sky.

    Corresponds to the ``COE`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::
        C &= \gamma / 2 \\
        R_\theta &= \frac{180^\circ}{\pi} \frac{2}{\gamma} \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin \theta} \\
        Y_0 &= \frac{180^\circ}{\pi} \frac{2}{\gamma} \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin((\theta_1 + \theta_2)/2)}

    where:

    .. math::
        \gamma = \sin \theta_1 + \sin \theta_2

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Sky2Pix_ConicEqualArea(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, x, y, sigma, delta):
        return _projections.coex2s(x, y, _to_orig_unit(sigma), _to_orig_unit(delta))


Pix2Sky_COE = Pix2Sky_ConicEqualArea


class Sky2Pix_ConicEqualArea(Sky2PixProjection, Conic):
    r"""
    Alber's conic equal area projection - sky to pixel.

    Corresponds to the ``COE`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::
        C &= \gamma / 2 \\
        R_\theta &= \frac{180^\circ}{\pi} \frac{2}{\gamma} \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin \theta} \\
        Y_0 &= \frac{180^\circ}{\pi} \frac{2}{\gamma} \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin((\theta_1 + \theta_2)/2)}

    where:

    .. math::
        \gamma = \sin \theta_1 + \sin \theta_2

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Pix2Sky_ConicEqualArea(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, phi, theta, sigma, delta):
        return _projections.coes2x(phi, theta,
                                   _to_orig_unit(sigma), _to_orig_unit(delta))


Sky2Pix_COE = Sky2Pix_ConicEqualArea


class Pix2Sky_ConicEquidistant(Pix2SkyProjection, Conic):
    r"""
    Conic equidistant projection - pixel to sky.

    Corresponds to the ``COD`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::

        C &= \frac{180^\circ}{\pi} \frac{\sin\theta_a\sin\eta}{\eta} \\
        R_\theta &= \theta_a - \theta + \eta\cot\eta\cot\theta_a \\
        Y_0 = \eta\cot\eta\cot\theta_a

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Sky2Pix_ConicEquidistant(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, x, y, sigma, delta):
        return _projections.codx2s(x, y, _to_orig_unit(sigma), _to_orig_unit(delta))


Pix2Sky_COD = Pix2Sky_ConicEquidistant


class Sky2Pix_ConicEquidistant(Sky2PixProjection, Conic):
    r"""
    Conic equidistant projection - sky to pixel.

    Corresponds to the ``COD`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::

        C &= \frac{180^\circ}{\pi} \frac{\sin\theta_a\sin\eta}{\eta} \\
        R_\theta &= \theta_a - \theta + \eta\cot\eta\cot\theta_a \\
        Y_0 = \eta\cot\eta\cot\theta_a

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Pix2Sky_ConicEquidistant(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, phi, theta, sigma, delta):
        return _projections.cods2x(phi, theta,
                                   _to_orig_unit(sigma), _to_orig_unit(delta))


Sky2Pix_COD = Sky2Pix_ConicEquidistant


class Pix2Sky_ConicOrthomorphic(Pix2SkyProjection, Conic):
    r"""
    Conic orthomorphic projection - pixel to sky.

    Corresponds to the ``COO`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::

        C &= \frac{\ln \left( \frac{\cos\theta_2}{\cos\theta_1} \right)}
                  {\ln \left[ \frac{\tan\left(\frac{90^\circ-\theta_2}{2}\right)}
                                   {\tan\left(\frac{90^\circ-\theta_1}{2}\right)} \right] } \\
        R_\theta &= \psi \left[ \tan \left( \frac{90^\circ - \theta}{2} \right) \right]^C \\
        Y_0 &= \psi \left[ \tan \left( \frac{90^\circ - \theta_a}{2} \right) \right]^C

    where:

    .. math::

        \psi = \frac{180^\circ}{\pi} \frac{\cos \theta}
               {C\left[\tan\left(\frac{90^\circ-\theta}{2}\right)\right]^C}

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Sky2Pix_ConicOrthomorphic(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, x, y, sigma, delta):
        return _projections.coox2s(x, y, _to_orig_unit(sigma), _to_orig_unit(delta))


Pix2Sky_COO = Pix2Sky_ConicOrthomorphic


class Sky2Pix_ConicOrthomorphic(Sky2PixProjection, Conic):
    r"""
    Conic orthomorphic projection - sky to pixel.

    Corresponds to the ``COO`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::

        C &= \frac{\ln \left( \frac{\cos\theta_2}{\cos\theta_1} \right)}
                  {\ln \left[ \frac{\tan\left(\frac{90^\circ-\theta_2}{2}\right)}
                                   {\tan\left(\frac{90^\circ-\theta_1}{2}\right)} \right] } \\
        R_\theta &= \psi \left[ \tan \left( \frac{90^\circ - \theta}{2} \right) \right]^C \\
        Y_0 &= \psi \left[ \tan \left( \frac{90^\circ - \theta_a}{2} \right) \right]^C

    where:

    .. math::

        \psi = \frac{180^\circ}{\pi} \frac{\cos \theta}
               {C\left[\tan\left(\frac{90^\circ-\theta}{2}\right)\right]^C}

    Parameters
    ----------
    sigma : float
        :math:`(\theta_1 + \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 90.

    delta : float
        :math:`(\theta_1 - \theta_2) / 2`, where :math:`\theta_1` and
        :math:`\theta_2` are the latitudes of the standard parallels,
        in degrees.  Default is 0.
    """
    @property
    def inverse(self):
        return Pix2Sky_ConicOrthomorphic(self.sigma.value, self.delta.value)

    @classmethod
    def evaluate(cls, phi, theta, sigma, delta):
        return _projections.coos2x(phi, theta,
                                   _to_orig_unit(sigma), _to_orig_unit(delta))


Sky2Pix_COO = Sky2Pix_ConicOrthomorphic


class PseudoConic(Projection):
    r"""Base class for pseudoconic projections.

    Pseudoconics are a subclass of conics with concentric parallels.
    """


class Pix2Sky_BonneEqualArea(Pix2SkyProjection, PseudoConic):
    r"""
    Bonne's equal area pseudoconic projection - pixel to sky.

    Corresponds to the ``BON`` projection in FITS WCS.

    .. math::

        \phi &= \frac{\pi}{180^\circ} A_\phi R_\theta / \cos \theta \\
        \theta &= Y_0 - R_\theta

    where:

    .. math::

        R_\theta &= \mathrm{sign} \theta_1 \sqrt{x^2 + (Y_0 - y)^2} \\
        A_\phi &= \arg\left(\frac{Y_0 - y}{R_\theta}, \frac{x}{R_\theta}\right)

    Parameters
    ----------
    theta1 : float
        Bonne conformal latitude, in degrees.
    """
    theta1 = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)
    _separable = True

    @property
    def inverse(self):
        return Sky2Pix_BonneEqualArea(self.theta1.value)

    @classmethod
    def evaluate(cls, x, y, theta1):
        return _projections.bonx2s(x, y, _to_orig_unit(theta1))


Pix2Sky_BON = Pix2Sky_BonneEqualArea


class Sky2Pix_BonneEqualArea(Sky2PixProjection, PseudoConic):
    r"""
    Bonne's equal area pseudoconic projection - sky to pixel.

    Corresponds to the ``BON`` projection in FITS WCS.

    .. math::
        x &= R_\theta \sin A_\phi \\
        y &= -R_\theta \cos A_\phi + Y_0

    where:

    .. math::
        A_\phi &= \frac{180^\circ}{\pi R_\theta} \phi \cos \theta \\
        R_\theta &= Y_0 - \theta \\
        Y_0 &= \frac{180^\circ}{\pi} \cot \theta_1 + \theta_1

    Parameters
    ----------
    theta1 : float
        Bonne conformal latitude, in degrees.
    """
    theta1 = Parameter(default=0.0, getter=_to_orig_unit, setter=_to_radian)
    _separable = True

    @property
    def inverse(self):
        return Pix2Sky_BonneEqualArea(self.theta1.value)

    @classmethod
    def evaluate(cls, phi, theta, theta1):
        return _projections.bons2x(phi, theta,
                                   _to_orig_unit(theta1))


Sky2Pix_BON = Sky2Pix_BonneEqualArea


class Pix2Sky_Polyconic(Pix2SkyProjection, PseudoConic):
    r"""
    Polyconic projection - pixel to sky.

    Corresponds to the ``PCO`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_Polyconic()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.pcox2s(x, y)


Pix2Sky_PCO = Pix2Sky_Polyconic


class Sky2Pix_Polyconic(Sky2PixProjection, PseudoConic):
    r"""
    Polyconic projection - sky to pixel.

    Corresponds to the ``PCO`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_Polyconic()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.pcos2x(phi, theta)


Sky2Pix_PCO = Sky2Pix_Polyconic


class QuadCube(Projection):
    r"""Base class for quad cube projections.

    Quadrilateralized spherical cube (quad-cube) projections belong to
    the class of polyhedral projections in which the sphere is
    projected onto the surface of an enclosing polyhedron.

    The six faces of the quad-cube projections are numbered and laid
    out as::

              0
        4 3 2 1 4 3 2
              5

    """


class Pix2Sky_TangentialSphericalCube(Pix2SkyProjection, QuadCube):
    r"""
    Tangential spherical cube projection - pixel to sky.

    Corresponds to the ``TSC`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_TangentialSphericalCube()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.tscx2s(x, y)


Pix2Sky_TSC = Pix2Sky_TangentialSphericalCube


class Sky2Pix_TangentialSphericalCube(Sky2PixProjection, QuadCube):
    r"""
    Tangential spherical cube projection - sky to pixel.

    Corresponds to the ``PCO`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_TangentialSphericalCube()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.tscs2x(phi, theta)


Sky2Pix_TSC = Sky2Pix_TangentialSphericalCube


class Pix2Sky_COBEQuadSphericalCube(Pix2SkyProjection, QuadCube):
    r"""
    COBE quadrilateralized spherical cube projection - pixel to sky.

    Corresponds to the ``CSC`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_COBEQuadSphericalCube()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.cscx2s(x, y)


Pix2Sky_CSC = Pix2Sky_COBEQuadSphericalCube


class Sky2Pix_COBEQuadSphericalCube(Sky2PixProjection, QuadCube):
    r"""
    COBE quadrilateralized spherical cube projection - sky to pixel.

    Corresponds to the ``CSC`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_COBEQuadSphericalCube()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.cscs2x(phi, theta)


Sky2Pix_CSC = Sky2Pix_COBEQuadSphericalCube


class Pix2Sky_QuadSphericalCube(Pix2SkyProjection, QuadCube):
    r"""
    Quadrilateralized spherical cube projection - pixel to sky.

    Corresponds to the ``QSC`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_QuadSphericalCube()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.qscx2s(x, y)


Pix2Sky_QSC = Pix2Sky_QuadSphericalCube


class Sky2Pix_QuadSphericalCube(Sky2PixProjection, QuadCube):
    r"""
    Quadrilateralized spherical cube projection - sky to pixel.

    Corresponds to the ``QSC`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_QuadSphericalCube()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.qscs2x(phi, theta)


Sky2Pix_QSC = Sky2Pix_QuadSphericalCube


class HEALPix(Projection):
    r"""Base class for HEALPix projections.
    """


class Pix2Sky_HEALPix(Pix2SkyProjection, HEALPix):
    r"""
    HEALPix - pixel to sky.

    Corresponds to the ``HPX`` projection in FITS WCS.

    Parameters
    ----------
    H : float
        The number of facets in longitude direction.

    X : float
        The number of facets in latitude direction.
    """
    _separable = True

    H = Parameter(default=4.0)
    X = Parameter(default=3.0)

    @property
    def inverse(self):
        return Sky2Pix_HEALPix(self.H.value, self.X.value)

    @classmethod
    def evaluate(cls, x, y, H, X):
        return _projections.hpxx2s(x, y, H, X)


Pix2Sky_HPX = Pix2Sky_HEALPix


class Sky2Pix_HEALPix(Sky2PixProjection, HEALPix):
    r"""
    HEALPix projection - sky to pixel.

    Corresponds to the ``HPX`` projection in FITS WCS.

    Parameters
    ----------
    H : float
        The number of facets in longitude direction.

    X : float
        The number of facets in latitude direction.
    """
    _separable = True

    H = Parameter(default=4.0)
    X = Parameter(default=3.0)

    @property
    def inverse(self):
        return Pix2Sky_HEALPix(self.H.value, self.X.value)

    @classmethod
    def evaluate(cls, phi, theta, H, X):
        return _projections.hpxs2x(phi, theta, H, X)


Sky2Pix_HPX = Sky2Pix_HEALPix


class Pix2Sky_HEALPixPolar(Pix2SkyProjection, HEALPix):
    r"""
    HEALPix polar, aka "butterfly" projection - pixel to sky.

    Corresponds to the ``XPH`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Sky2Pix_HEALPix()

    @classmethod
    def evaluate(cls, x, y):
        return _projections.xphx2s(x, y)


Pix2Sky_XPH = Pix2Sky_HEALPixPolar


class Sky2Pix_HEALPixPolar(Sky2PixProjection, HEALPix):
    r"""
    HEALPix polar, aka "butterfly" projection - pixel to sky.

    Corresponds to the ``XPH`` projection in FITS WCS.
    """
    _separable = False

    @property
    def inverse(self):
        return Pix2Sky_HEALPix()

    @classmethod
    def evaluate(cls, phi, theta):
        return _projections.hpxs2x(phi, theta)


Sky2Pix_XPH = Sky2Pix_HEALPixPolar


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

    _separable = False

    matrix = Parameter(default=[[1.0, 0.0], [0.0, 1.0]])
    translation = Parameter(default=[0.0, 0.0])

    @matrix.validator
    def matrix(self, value):
        """Validates that the input matrix is a 2x2 2D array."""

        if np.shape(value) != (2, 2):
            raise InputParameterError(
                "Expected transformation matrix to be a 2x2 array")

    @translation.validator
    def translation(self, value):
        """
        Validates that the translation vector is a 2D vector.  This allows
        either a "row" vector or a "column" vector where in the latter case the
        resultant Numpy array has ``ndim=2`` but the shape is ``(1, 2)``.
        """

        if not ((np.ndim(value) == 1 and np.shape(value) == (2,)) or
                (np.ndim(value) == 2 and np.shape(value) == (1, 2))):
            raise InputParameterError(
                "Expected translation vector to be a 2 element row or column "
                "vector array")

    def __init__(self, matrix=matrix, translation=translation, **kwargs):
        super().__init__(matrix=matrix, translation=translation, **kwargs)
        self.inputs = ("x", "y")
        self.outputs = ("x", "y")

    @property
    def inverse(self):
        """
        Inverse transformation.

        Raises `~astropy.modeling.InputParameterError` if the transformation cannot be inverted.
        """

        det = np.linalg.det(self.matrix.value)

        if det == 0:
            raise InputParameterError(
                "Transformation matrix is singular; {} model does not "
                "have an inverse".format(self.__class__.__name__))

        matrix = np.linalg.inv(self.matrix.value)
        if self.matrix.unit is not None:
            matrix = matrix * self.matrix.unit
        # If matrix has unit then translation has unit, so no need to assign it.
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
        # Use asarray to ensure loose the units.
        inarr = np.vstack([np.asarray(x).ravel(),
                           np.asarray(y).ravel(),
                           np.ones(x.size, x.dtype)])

        if inarr.shape[0] != 3 or inarr.ndim != 2:
            raise ValueError("Incompatible input shapes")

        augmented_matrix = cls._create_augmented_matrix(matrix, translation)
        result = np.dot(augmented_matrix, inarr)
        x, y = result[0], result[1]
        x.shape = y.shape = shape

        return x, y

    @staticmethod
    def _create_augmented_matrix(matrix, translation):
        unit = None
        if any([hasattr(translation, 'unit'), hasattr(matrix, 'unit')]):
            if not all([hasattr(translation, 'unit'), hasattr(matrix, 'unit')]):
                raise ValueError("To use AffineTransformation with quantities, "
                                 "both matrix and unit need to be quantities.")
            unit = translation.unit
            # matrix should have the same units as translation
            if not (matrix.unit / translation.unit) == u.dimensionless_unscaled:
                raise ValueError("matrix and translation must have the same units.")

        augmented_matrix = np.empty((3, 3), dtype=float)
        augmented_matrix[0:2, 0:2] = matrix
        augmented_matrix[0:2, 2:].flat = translation
        augmented_matrix[2] = [0, 0, 1]
        if unit is not None:
            return augmented_matrix * unit
        return augmented_matrix

    @property
    def input_units(self):
        if self.translation.unit is None and self.matrix.unit is None:
            return None
        elif self.translation.unit is not None:
            return dict(zip(self.inputs, [self.translation.unit] * 2))
        else:
            return dict(zip(self.inputs, [self.matrix.unit] * 2))
