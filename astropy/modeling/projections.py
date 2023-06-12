# Licensed under a 3-clause BSD style license - see LICENSE.rst
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
from itertools import chain, product

import numpy as np

from astropy import units as u
from astropy import wcs

from .core import Model
from .parameters import InputParameterError, Parameter
from .utils import _to_orig_unit, _to_radian

# List of tuples of the form
# (long class name without suffix, short WCSLIB projection code):
_PROJ_NAME_CODE = [
    ("ZenithalPerspective", "AZP"),
    ("SlantZenithalPerspective", "SZP"),
    ("Gnomonic", "TAN"),
    ("Stereographic", "STG"),
    ("SlantOrthographic", "SIN"),
    ("ZenithalEquidistant", "ARC"),
    ("ZenithalEqualArea", "ZEA"),
    ("Airy", "AIR"),
    ("CylindricalPerspective", "CYP"),
    ("CylindricalEqualArea", "CEA"),
    ("PlateCarree", "CAR"),
    ("Mercator", "MER"),
    ("SansonFlamsteed", "SFL"),
    ("Parabolic", "PAR"),
    ("Molleweide", "MOL"),
    ("HammerAitoff", "AIT"),
    ("ConicPerspective", "COP"),
    ("ConicEqualArea", "COE"),
    ("ConicEquidistant", "COD"),
    ("ConicOrthomorphic", "COO"),
    ("BonneEqualArea", "BON"),
    ("Polyconic", "PCO"),
    ("TangentialSphericalCube", "TSC"),
    ("COBEQuadSphericalCube", "CSC"),
    ("QuadSphericalCube", "QSC"),
    ("HEALPix", "HPX"),
    ("HEALPixPolar", "XPH"),
]

_NOT_SUPPORTED_PROJ_CODES = ["ZPN"]

_PROJ_NAME_CODE_MAP = dict(_PROJ_NAME_CODE)

projcodes = [code for _, code in _PROJ_NAME_CODE]


__all__ = [
    "Projection",
    "Pix2SkyProjection",
    "Sky2PixProjection",
    "Zenithal",
    "Cylindrical",
    "PseudoCylindrical",
    "Conic",
    "PseudoConic",
    "QuadCube",
    "HEALPix",
    "AffineTransformation2D",
    "projcodes",
] + list(map("_".join, product(["Pix2Sky", "Sky2Pix"], chain(*_PROJ_NAME_CODE))))


class _ParameterDS(Parameter):
    """
    Same as `Parameter` but can indicate its modified status via the ``dirty``
    property. This flag also gets set automatically when a parameter is
    modified.

    This ability to track parameter's modified status is needed for automatic
    update of WCSLIB's prjprm structure (which may be a more-time intensive
    operation) *only as required*.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dirty = True

    def validate(self, value):
        super().validate(value)
        self.dirty = True


class Projection(Model):
    """Base class for all sky projections."""

    # Radius of the generating sphere.
    # This sets the circumference to 360 deg so that arc length is measured in deg.
    r0 = 180 * u.deg / np.pi

    _separable = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._prj = wcs.Prjprm()

    @property
    @abc.abstractmethod
    def inverse(self):
        """
        Inverse projection--all projection models must provide an inverse.
        """

    @property
    def prjprm(self):
        """WCSLIB ``prjprm`` structure."""
        self._update_prj()
        return self._prj

    def _update_prj(self):
        """
        A default updater for projection's pv.

        .. warning::
            This method assumes that PV0 is never modified. If a projection
            that uses PV0 is ever implemented in this module, that projection
            class should override this method.

        .. warning::
            This method assumes that the order in which PVi values (i>0)
            are to be assigned is identical to the order of model parameters
            in ``param_names``. That is, pv[1] = model.parameters[0], ...

        """
        if not self.param_names:
            return

        pv = []
        dirty = False

        for p in self.param_names:
            param = getattr(self, p)
            pv.append(float(param.value))
            dirty |= param.dirty
            param.dirty = False

        if dirty:
            self._prj.pv = None, *pv
            self._prj.set()

    def __getstate__(self):
        return {
            "p": self.parameters,
            "fixed": self.fixed,
            "tied": self.tied,
            "bounds": self.bounds,
        }

    def __setstate__(self, state):
        params = state.pop("p")
        return self.__init__(*params, **state)


class Pix2SkyProjection(Projection):
    """Base class for all Pix2Sky projections."""

    n_inputs = 2
    n_outputs = 2

    _input_units_strict = True
    _input_units_allow_dimensionless = True

    def __new__(cls, *args, **kwargs):
        long_name = cls.name.split("_")[1]
        cls.prj_code = _PROJ_NAME_CODE_MAP[long_name]
        return super().__new__(cls)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._prj.code = self.prj_code
        self._update_prj()
        if not self.param_names:
            # force initial call to Prjprm.set() for projections
            # with no parameters:
            self._prj.set()

        self.inputs = ("x", "y")
        self.outputs = ("phi", "theta")

    @property
    def input_units(self):
        return {self.inputs[0]: u.deg, self.inputs[1]: u.deg}

    @property
    def return_units(self):
        return {self.outputs[0]: u.deg, self.outputs[1]: u.deg}

    def evaluate(self, x, y, *args, **kwargs):
        self._update_prj()
        return self._prj.prjx2s(x, y)

    @property
    def inverse(self):
        pv = [getattr(self, param).value for param in self.param_names]
        return self._inv_cls(*pv)


class Sky2PixProjection(Projection):
    """Base class for all Sky2Pix projections."""

    n_inputs = 2
    n_outputs = 2

    _input_units_strict = True
    _input_units_allow_dimensionless = True

    def __new__(cls, *args, **kwargs):
        long_name = cls.name.split("_")[1]
        cls.prj_code = _PROJ_NAME_CODE_MAP[long_name]
        return super().__new__(cls)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._prj.code = self.prj_code
        self._update_prj()
        if not self.param_names:
            # force initial call to Prjprm.set() for projections
            # without parameters:
            self._prj.set()

        self.inputs = ("phi", "theta")
        self.outputs = ("x", "y")

    @property
    def input_units(self):
        return {self.inputs[0]: u.deg, self.inputs[1]: u.deg}

    @property
    def return_units(self):
        return {self.outputs[0]: u.deg, self.outputs[1]: u.deg}

    def evaluate(self, phi, theta, *args, **kwargs):
        self._update_prj()
        return self._prj.prjs2x(phi, theta)

    @property
    def inverse(self):
        pv = [getattr(self, param).value for param in self.param_names]
        return self._inv_cls(*pv)


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
    ----------
    mu : float
        Distance from point of projection to center of sphere
        in spherical radii, μ.  Default is 0.

    gamma : float
        Look angle γ in degrees.  Default is 0°.

    """

    mu = _ParameterDS(
        default=0.0, description="Distance from point of projection to center of sphere"
    )
    gamma = _ParameterDS(
        default=0.0,
        getter=_to_orig_unit,
        setter=_to_radian,
        description="Look angle γ in degrees (Default = 0°)",
    )

    def _mu_validator(self, value):
        if np.any(np.equal(value, -1.0)):
            raise InputParameterError(
                "Zenithal perspective projection is not defined for mu = -1"
            )

    mu._validator = _mu_validator


class Sky2Pix_ZenithalPerspective(Sky2PixProjection, Zenithal):
    r"""
    Zenithal perspective projection - sky to pixel.

    Corresponds to the ``AZP`` projection in FITS WCS.

    .. math::
        x &= R \sin \phi \\
        y &= -R \sec \gamma \cos \theta

    where:

    .. math::
        R = \frac{180^{\circ}}{\pi} \frac{(\mu + 1) \cos \theta}
            {(\mu + \sin \theta) + \cos \theta \cos \phi \tan \gamma}

    Parameters
    ----------
    mu : float
        Distance from point of projection to center of sphere
        in spherical radii, μ. Default is 0.

    gamma : float
        Look angle γ in degrees. Default is 0°.

    """

    mu = _ParameterDS(
        default=0.0, description="Distance from point of projection to center of sphere"
    )
    gamma = _ParameterDS(
        default=0.0,
        getter=_to_orig_unit,
        setter=_to_radian,
        description="Look angle γ in degrees (Default=0°)",
    )

    def _mu_validator(self, value):
        if np.any(np.equal(value, -1.0)):
            raise InputParameterError(
                "Zenithal perspective projection is not defined for mu = -1"
            )

    mu._validator = _mu_validator


class Pix2Sky_SlantZenithalPerspective(Pix2SkyProjection, Zenithal):
    r"""
    Slant zenithal perspective projection - pixel to sky.

    Corresponds to the ``SZP`` projection in FITS WCS.

    Parameters
    ----------
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

    mu = _ParameterDS(
        default=0.0, description="Distance from point of projection to center of sphere"
    )
    phi0 = _ParameterDS(
        default=0.0,
        getter=_to_orig_unit,
        setter=_to_radian,
        description="The longitude φ₀ of the reference point in degrees (Default=0°)",
    )
    theta0 = _ParameterDS(
        default=90.0,
        getter=_to_orig_unit,
        setter=_to_radian,
        description="The latitude θ₀ of the reference point, in degrees (Default=0°)",
    )

    def _mu_validator(self, value):
        if np.any(np.equal(value, -1.0)):
            raise InputParameterError(
                "Zenithal perspective projection is not defined for mu = -1"
            )

    mu._validator = _mu_validator


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

    mu = _ParameterDS(
        default=0.0, description="Distance from point of projection to center of sphere"
    )
    phi0 = _ParameterDS(
        default=0.0,
        getter=_to_orig_unit,
        setter=_to_radian,
        description="The longitude φ₀ of the reference point in degrees",
    )
    theta0 = _ParameterDS(
        default=0.0,
        getter=_to_orig_unit,
        setter=_to_radian,
        description="The latitude θ₀ of the reference point, in degrees",
    )

    def _mu_validator(self, value):
        if np.any(np.equal(value, -1.0)):
            raise InputParameterError(
                "Zenithal perspective projection is not defined for mu = -1"
            )

    mu._validator = _mu_validator


class Pix2Sky_Gnomonic(Pix2SkyProjection, Zenithal):
    r"""
    Gnomonic projection - pixel to sky.

    Corresponds to the ``TAN`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = \tan^{-1}\left(\frac{180^{\circ}}{\pi R_\theta}\right)
    """


class Sky2Pix_Gnomonic(Sky2PixProjection, Zenithal):
    r"""
    Gnomonic Projection - sky to pixel.

    Corresponds to the ``TAN`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\cot \theta
    """


class Pix2Sky_Stereographic(Pix2SkyProjection, Zenithal):
    r"""
    Stereographic Projection - pixel to sky.

    Corresponds to the ``STG`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = 90^{\circ} - 2 \tan^{-1}\left(\frac{\pi R_\theta}{360^{\circ}}\right)
    """


class Sky2Pix_Stereographic(Sky2PixProjection, Zenithal):
    r"""
    Stereographic Projection - sky to pixel.

    Corresponds to the ``STG`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = \frac{180^{\circ}}{\pi}\frac{2 \cos \theta}{1 + \sin \theta}
    """


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

    xi = _ParameterDS(default=0.0, description="Obliqueness parameter")
    eta = _ParameterDS(default=0.0, description="Obliqueness parameter")


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

    xi = _ParameterDS(default=0.0)
    eta = _ParameterDS(default=0.0)


class Pix2Sky_ZenithalEquidistant(Pix2SkyProjection, Zenithal):
    r"""
    Zenithal equidistant projection - pixel to sky.

    Corresponds to the ``ARC`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = 90^\circ - R_\theta
    """


class Sky2Pix_ZenithalEquidistant(Sky2PixProjection, Zenithal):
    r"""
    Zenithal equidistant projection - sky to pixel.

    Corresponds to the ``ARC`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = 90^\circ - \theta
    """


class Pix2Sky_ZenithalEqualArea(Pix2SkyProjection, Zenithal):
    r"""
    Zenithal equidistant projection - pixel to sky.

    Corresponds to the ``ZEA`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        \theta = 90^\circ - 2 \sin^{-1} \left(\frac{\pi R_\theta}{360^\circ}\right)
    """


class Sky2Pix_ZenithalEqualArea(Sky2PixProjection, Zenithal):
    r"""
    Zenithal equidistant projection - sky to pixel.

    Corresponds to the ``ZEA`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta &= \frac{180^\circ}{\pi} \sqrt{2(1 - \sin\theta)} \\
                 &= \frac{360^\circ}{\pi} \sin\left(\frac{90^\circ - \theta}{2}\right)
    """


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

    theta_b = _ParameterDS(default=90.0)


class Sky2Pix_Airy(Sky2PixProjection, Zenithal):
    r"""
    Airy - sky to pixel.

    Corresponds to the ``AIR`` projection in FITS WCS.

    See `Zenithal` for a definition of the full transformation.

    .. math::
        R_\theta = -2 \frac{180^\circ}{\pi}\left(\frac{\ln(\cos \xi)}{\tan \xi} +
            \frac{\ln(\cos \xi_b)}{\tan^2 \xi_b} \tan \xi \right)

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

    theta_b = _ParameterDS(
        default=90.0,
        description="The latitude at which to minimize the error,in degrees",
    )


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

    mu = _ParameterDS(default=1.0)
    lam = _ParameterDS(default=1.0)

    def _mu_validator(self, value):
        if np.any(value == -self.lam):
            raise InputParameterError("CYP projection is not defined for mu = -lambda")

    mu._validator = _mu_validator

    def _lam_validator(self, value):
        if np.any(value == -self.mu):
            raise InputParameterError("CYP projection is not defined for lambda = -mu")

    lam._validator = _lam_validator


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

    mu = _ParameterDS(
        default=1.0, description="Distance from center of sphere in spherical radii"
    )
    lam = _ParameterDS(
        default=1.0, description="Radius of the cylinder in spherical radii"
    )

    def _mu_validator(self, value):
        if np.any(value == -self.lam):
            raise InputParameterError("CYP projection is not defined for mu = -lambda")

    mu._validator = _mu_validator

    def _lam_validator(self, value):
        if np.any(value == -self.mu):
            raise InputParameterError("CYP projection is not defined for lambda = -mu")

    lam._validator = _lam_validator


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

    lam = _ParameterDS(default=1)


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

    lam = _ParameterDS(default=1)


class Pix2Sky_PlateCarree(Pix2SkyProjection, Cylindrical):
    r"""
    Plate carrée projection - pixel to sky.

    Corresponds to the ``CAR`` projection in FITS WCS.

    .. math::
        \phi &= x \\
        \theta &= y
    """

    @staticmethod
    def evaluate(x, y):
        # The intermediate variables are only used here for clarity
        phi = np.array(x)
        theta = np.array(y)
        return phi, theta


class Sky2Pix_PlateCarree(Sky2PixProjection, Cylindrical):
    r"""
    Plate carrée projection - sky to pixel.

    Corresponds to the ``CAR`` projection in FITS WCS.

    .. math::
        x &= \phi \\
        y &= \theta
    """

    @staticmethod
    def evaluate(phi, theta):
        # The intermediate variables are only used here for clarity
        x = np.array(phi)
        y = np.array(theta)
        return x, y


class Pix2Sky_Mercator(Pix2SkyProjection, Cylindrical):
    r"""
    Mercator - pixel to sky.

    Corresponds to the ``MER`` projection in FITS WCS.

    .. math::
        \phi &= x \\
        \theta &= 2 \tan^{-1}\left(e^{y \pi / 180^{\circ}}\right)-90^{\circ}
    """


class Sky2Pix_Mercator(Sky2PixProjection, Cylindrical):
    r"""
    Mercator - sky to pixel.

    Corresponds to the ``MER`` projection in FITS WCS.

    .. math::
        x &= \phi \\
        y &= \frac{180^{\circ}}{\pi}\ln \tan \left(\frac{90^{\circ} + \theta}{2}\right)
    """


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


class Sky2Pix_SansonFlamsteed(Sky2PixProjection, PseudoCylindrical):
    r"""
    Sanson-Flamsteed projection - sky to pixel.

    Corresponds to the ``SFL`` projection in FITS WCS.

    .. math::
        x &= \phi \cos \theta \\
        y &= \theta
    """


class Pix2Sky_Parabolic(Pix2SkyProjection, PseudoCylindrical):
    r"""
    Parabolic projection - pixel to sky.

    Corresponds to the ``PAR`` projection in FITS WCS.

    .. math::
        \phi &= \frac{180^\circ}{\pi} \frac{x}{1 - 4(y / 180^\circ)^2} \\
        \theta &= 3 \sin^{-1}\left(\frac{y}{180^\circ}\right)
    """


class Sky2Pix_Parabolic(Sky2PixProjection, PseudoCylindrical):
    r"""
    Parabolic projection - sky to pixel.

    Corresponds to the ``PAR`` projection in FITS WCS.

    .. math::
        x &= \phi \left(2\cos\frac{2\theta}{3} - 1\right) \\
        y &= 180^\circ \sin \frac{\theta}{3}
    """


class Pix2Sky_Molleweide(Pix2SkyProjection, PseudoCylindrical):
    r"""
    Molleweide's projection - pixel to sky.

    Corresponds to the ``MOL`` projection in FITS WCS.

    .. math::
        \phi &= \frac{\pi x}{2 \sqrt{2 - \left(\frac{\pi}{180^\circ}y\right)^2}} \\
        \theta &= \sin^{-1}\left(
                \frac{1}{90^\circ}\sin^{-1}\left(\frac{\pi}{180^\circ}\frac{y}{\sqrt{2}}\right)
                + \frac{y}{180^\circ}\sqrt{2 - \left(\frac{\pi}{180^\circ}y\right)^2}
            \right)
    """


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


class Pix2Sky_HammerAitoff(Pix2SkyProjection, PseudoCylindrical):
    r"""
    Hammer-Aitoff projection - pixel to sky.

    Corresponds to the ``AIT`` projection in FITS WCS.

    .. math::
        \phi &= 2 \arg \left(2Z^2 - 1, \frac{\pi}{180^\circ} \frac{Z}{2}x\right) \\
        \theta &= \sin^{-1}\left(\frac{\pi}{180^\circ}yZ\right)
    """


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

    sigma = _ParameterDS(default=90.0, getter=_to_orig_unit, setter=_to_radian)
    delta = _ParameterDS(default=0.0, getter=_to_orig_unit, setter=_to_radian)


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


class Pix2Sky_ConicEqualArea(Pix2SkyProjection, Conic):
    r"""
    Alber's conic equal area projection - pixel to sky.

    Corresponds to the ``COE`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::
        C &= \gamma / 2 \\
        R_\theta &= \frac{180^\circ}{\pi} \frac{2}{\gamma}
            \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin \theta} \\
        Y_0 &= \frac{180^\circ}{\pi} \frac{2}{\gamma}
            \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin((\theta_1 + \theta_2)/2)}

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


class Sky2Pix_ConicEqualArea(Sky2PixProjection, Conic):
    r"""
    Alber's conic equal area projection - sky to pixel.

    Corresponds to the ``COE`` projection in FITS WCS.

    See `Conic` for a description of the entire equation.

    The projection formulae are:

    .. math::
        C &= \gamma / 2 \\
        R_\theta &= \frac{180^\circ}{\pi} \frac{2}{\gamma}
            \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin \theta} \\
        Y_0 &= \frac{180^\circ}{\pi} \frac{2}{\gamma}
            \sqrt{1 + \sin \theta_1 \sin \theta_2 - \gamma \sin((\theta_1 + \theta_2)/2)}

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

    _separable = True

    theta1 = _ParameterDS(default=0.0, getter=_to_orig_unit, setter=_to_radian)


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

    _separable = True

    theta1 = _ParameterDS(
        default=0.0,
        getter=_to_orig_unit,
        setter=_to_radian,
        description="Bonne conformal latitude, in degrees",
    )


class Pix2Sky_Polyconic(Pix2SkyProjection, PseudoConic):
    r"""
    Polyconic projection - pixel to sky.

    Corresponds to the ``PCO`` projection in FITS WCS.
    """


class Sky2Pix_Polyconic(Sky2PixProjection, PseudoConic):
    r"""
    Polyconic projection - sky to pixel.

    Corresponds to the ``PCO`` projection in FITS WCS.
    """


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


class Sky2Pix_TangentialSphericalCube(Sky2PixProjection, QuadCube):
    r"""
    Tangential spherical cube projection - sky to pixel.

    Corresponds to the ``TSC`` projection in FITS WCS.
    """


class Pix2Sky_COBEQuadSphericalCube(Pix2SkyProjection, QuadCube):
    r"""
    COBE quadrilateralized spherical cube projection - pixel to sky.

    Corresponds to the ``CSC`` projection in FITS WCS.
    """


class Sky2Pix_COBEQuadSphericalCube(Sky2PixProjection, QuadCube):
    r"""
    COBE quadrilateralized spherical cube projection - sky to pixel.

    Corresponds to the ``CSC`` projection in FITS WCS.
    """


class Pix2Sky_QuadSphericalCube(Pix2SkyProjection, QuadCube):
    r"""
    Quadrilateralized spherical cube projection - pixel to sky.

    Corresponds to the ``QSC`` projection in FITS WCS.
    """


class Sky2Pix_QuadSphericalCube(Sky2PixProjection, QuadCube):
    r"""
    Quadrilateralized spherical cube projection - sky to pixel.

    Corresponds to the ``QSC`` projection in FITS WCS.
    """


class HEALPix(Projection):
    r"""Base class for HEALPix projections."""


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

    H = _ParameterDS(
        default=4.0, description="The number of facets in longitude direction."
    )
    X = _ParameterDS(
        default=3.0, description="The number of facets in latitude direction."
    )


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

    H = _ParameterDS(
        default=4.0, description="The number of facets in longitude direction."
    )
    X = _ParameterDS(
        default=3.0, description="The number of facets in latitude direction."
    )


class Pix2Sky_HEALPixPolar(Pix2SkyProjection, HEALPix):
    r"""
    HEALPix polar, aka "butterfly" projection - pixel to sky.

    Corresponds to the ``XPH`` projection in FITS WCS.
    """


class Sky2Pix_HEALPixPolar(Sky2PixProjection, HEALPix):
    r"""
    HEALPix polar, aka "butterfly" projection - pixel to sky.

    Corresponds to the ``XPH`` projection in FITS WCS.
    """


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

    def _matrix_validator(self, value):
        """Validates that the input matrix is a 2x2 2D array."""
        if np.shape(value) != (2, 2):
            raise InputParameterError(
                "Expected transformation matrix to be a 2x2 array"
            )

    matrix._validator = _matrix_validator

    def _translation_validator(self, value):
        """
        Validates that the translation vector is a 2D vector.  This allows
        either a "row" vector or a "column" vector where in the latter case the
        resultant Numpy array has ``ndim=2`` but the shape is ``(1, 2)``.
        """
        if not (
            (np.ndim(value) == 1 and np.shape(value) == (2,))
            or (np.ndim(value) == 2 and np.shape(value) == (1, 2))
        ):
            raise InputParameterError(
                "Expected translation vector to be a 2 element row or column "
                "vector array"
            )

    translation._validator = _translation_validator

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
                f"Transformation matrix is singular; {self.__class__.__name__} model"
                " does not have an inverse"
            )

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
        inarr = np.vstack(
            [np.asarray(x).ravel(), np.asarray(y).ravel(), np.ones(x.size, x.dtype)]
        )

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
        if any([hasattr(translation, "unit"), hasattr(matrix, "unit")]):
            if not all([hasattr(translation, "unit"), hasattr(matrix, "unit")]):
                raise ValueError(
                    "To use AffineTransformation with quantities, "
                    "both matrix and unit need to be quantities."
                )
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
        translation_unit = self.translation.input_unit
        matrix_unit = self.matrix.input_unit

        if translation_unit is None and matrix_unit is None:
            return None
        elif translation_unit is not None:
            return dict(zip(self.inputs, [translation_unit] * 2))
        else:
            return dict(zip(self.inputs, [matrix_unit] * 2))


for long_name, short_name in _PROJ_NAME_CODE:
    # define short-name projection equivalent classes:
    globals()["Pix2Sky_" + short_name] = globals()["Pix2Sky_" + long_name]
    globals()["Sky2Pix_" + short_name] = globals()["Sky2Pix_" + long_name]
    # set inverse classes:
    globals()["Pix2Sky_" + long_name]._inv_cls = globals()["Sky2Pix_" + long_name]
    globals()["Sky2Pix_" + long_name]._inv_cls = globals()["Pix2Sky_" + long_name]
