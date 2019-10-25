# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Models that have physical origins.
"""

import warnings

import numpy as np

from .core import Fittable1DModel, Model
from .parameters import Parameter, InputParameterError

from astropy import constants as const
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning


__all__ = ["BlackBody", "Drude1D", "WavelengthFromGratingEquation",
           "AnglesFromGratingEquation3D"]


class BlackBody(Fittable1DModel):
    """
    Blackbody model using the Planck function.

    Parameters
    ----------
    temperature : :class:`~astropy.units.Quantity`
        Blackbody temperature.

    scale : float or :class:`~astropy.units.Quantity`
        Scale factor

    Notes
    -----

    Model formula:

        .. math:: B_{\\nu}(T) = A \\frac{2 h \\nu^{3} / c^{2}}{exp(h \\nu / k T) - 1}

    Examples
    --------
    >>> from astropy.modeling import models
    >>> from astropy import units as u
    >>> bb = models.BlackBody(temperature=5000*u.K)
    >>> bb(6000 * u.AA)  # doctest: +FLOAT_CMP
    <Quantity 1.53254685e-05 erg / (cm2 Hz s sr)>

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import BlackBody
        from astropy import units as u
        from astropy.visualization import quantity_support

        bb = BlackBody(temperature=5778*u.K)
        wav = np.arange(1000, 110000) * u.AA
        flux = bb(wav)

        with quantity_support():
            plt.figure()
            plt.semilogx(wav, flux)
            plt.axvline(bb.nu_max.to(u.AA, equivalencies=u.spectral()).value, ls='--')
            plt.show()
    """

    # We parametrize this model with a temperature and a scale.
    temperature = Parameter(default=5000.0, min=0, unit=u.K)
    scale = Parameter(default=1.0, min=0)

    # We allow values without units to be passed when evaluating the model, and
    # in this case the input x values are assumed to be frequencies in Hz.
    _input_units_allow_dimensionless = True

    # We enable the spectral equivalency by default for the spectral axis
    input_units_equivalencies = {"x": u.spectral()}

    def evaluate(self, x, temperature, scale):
        """Evaluate the model.

        Parameters
        ----------
        x : float, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Frequency at which to compute the blackbody. If no units are given,
            this defaults to Hz.

        temperature : float, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Temperature of the blackbody. If no units are given, this defaults
            to Kelvin.

        scale : float, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Desired scale for the blackbody.

        Returns
        -------
        y : number or ndarray
            Blackbody spectrum. The units are determined from the units of
            ``scale``.

        .. note::

            Use `numpy.errstate` to suppress Numpy warnings, if desired.

        .. warning::

            Output values might contain ``nan`` and ``inf``.

        Raises
        ------
        ValueError
            Invalid temperature.

        ZeroDivisionError
            Wavelength is zero (when converting to frequency).
        """
        if not isinstance(temperature, u.Quantity):
            in_temp = u.Quantity(temperature, u.K)
        else:
            in_temp = temperature

        # Convert to units for calculations, also force double precision
        with u.add_enabled_equivalencies(u.spectral() + u.temperature()):
            freq = u.Quantity(x, u.Hz, dtype=np.float64)
            temp = u.Quantity(in_temp, u.K)

        # check the units of scale and setup the output units
        bb_unit = u.erg / (u.cm ** 2 * u.s * u.Hz * u.sr)  # default unit
        # use the scale that was used at initialization for determining the units to return
        # to support returning the right units when fitting where units are stripped
        if hasattr(self.scale, "unit") and self.scale.unit is not None:
            # check that the units on scale are covertable to surface brightness units
            if not self.scale.unit.is_equivalent(bb_unit, u.spectral_density(x)):
                raise ValueError(
                    f"scale units not surface brightness: {self.scale.unit}"
                )
            # use the scale passed to get the value for scaling
            if hasattr(scale, "unit"):
                mult_scale = scale.value
            else:
                mult_scale = scale
            bb_unit = self.scale.unit
        else:
            mult_scale = scale

        # Check if input values are physically possible
        if np.any(temp < 0):
            raise ValueError(f"Temperature should be positive: {temp}")
        if not np.all(np.isfinite(freq)) or np.any(freq <= 0):
            warnings.warn(
                "Input contains invalid wavelength/frequency value(s)",
                AstropyUserWarning,
            )

        log_boltz = const.h * freq / (const.k_B * temp)
        boltzm1 = np.expm1(log_boltz)

        # Calculate blackbody flux
        bb_nu = 2.0 * const.h * freq ** 3 / (const.c ** 2 * boltzm1) / u.sr

        y = mult_scale * bb_nu.to(bb_unit, u.spectral_density(freq))

        # If the temperature parameter has no unit, we should return a unitless
        # value. This occurs for instance during fitting, since we drop the
        # units temporarily.
        if hasattr(temperature, "unit"):
            return y
        else:
            return y.value

    @property
    def input_units(self):
        # The input units are those of the 'x' value, which should always be
        # Hz. Because we do this, and because input_units_allow_dimensionless
        # is set to True, dimensionless values are assumed to be in Hz.
        return {"x": u.Hz}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {"temperature": u.K}

    @property
    def bolometric_flux(self):
        """Bolometric flux."""
        # bolometric flux in the native units of the planck function
        native_bolflux = (
            self.scale.value * const.sigma_sb * self.temperature ** 4 / np.pi
        )
        # return in more "astro" units
        return native_bolflux.to(u.erg / (u.cm ** 2 * u.s))

    @property
    def lambda_max(self):
        """Peak wavelength when the curve is expressed as power density."""
        return const.b_wien / self.temperature

    @property
    def nu_max(self):
        """Peak frequency when the curve is expressed as power density."""
        return 2.8214391 * const.k_B * self.temperature / const.h


class Drude1D(Fittable1DModel):
    """
    Drude model based one the behavior of electons in materials (esp. metals).

    Parameters
    ----------
    amplitude : float
        Peak value
    x_0 : float
        Position of the peak
    fwhm : float
        Full width at half maximum

    Model formula:

        .. math:: f(x) = A \\frac{(fwhm/x_0)^2}{((x/x_0 - x_0/x)^2 + (fwhm/x_0)^2}

    Examples
    --------

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Drude1D

        fig, ax = plt.subplots()

        # generate the curves and plot them
        x = np.arange(7.5 , 12.5 , 0.1)

        dmodel = Drude1D(amplitude=1.0, fwhm=1.0, x_0=10.0)
        ax.plot(x, dmodel(x))

        ax.set_xlabel('x')
        ax.set_ylabel('F(x)')

        ax.legend(loc='best')
        plt.show()
    """

    amplitude = Parameter(default=1.0)
    x_0 = Parameter(default=1.0)
    fwhm = Parameter(default=1.0)

    @staticmethod
    def evaluate(x, amplitude, x_0, fwhm):
        """
        One dimensional Drude model function
        """
        return (
            amplitude
            * ((fwhm / x_0) ** 2)
            / ((x / x_0 - x_0 / x) ** 2 + (fwhm / x_0) ** 2)
        )

    @staticmethod
    def fit_deriv(x, amplitude, x_0, fwhm):
        """
        Drude1D model function derivatives.
        """
        d_amplitude = (fwhm / x_0) ** 2 / ((x / x_0 - x_0 / x) ** 2 + (fwhm / x_0) ** 2)
        d_x_0 = (
            -2
            * amplitude
            * d_amplitude
            * (
                (1 / x_0)
                + d_amplitude
                * (x_0 ** 2 / fwhm ** 2)
                * (
                    (-x / x_0 - 1 / x) * (x / x_0 - x_0 / x)
                    - (2 * fwhm ** 2 / x_0 ** 3)
                )
            )
        )
        d_fwhm = (2 * amplitude * d_amplitude / fwhm) * (1 - d_amplitude)
        return [d_amplitude, d_x_0, d_fwhm]

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {"x": self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_0": inputs_unit["x"],
            "fwhm": inputs_unit["x"],
            "amplitude": outputs_unit["y"],
        }

    @property
    def return_units(self):
        if self.amplitude.unit is None:
            return None
        else:
            return {'y': self.amplitude.unit}

    @x_0.validator
    def x_0(self, val):
        if val == 0:
            raise InputParameterError("0 is not an allowed value for x_0")

    def bounding_box(self, factor=50):
        """Tuple defining the default ``bounding_box`` limits,
        ``(x_low, x_high)``.

        Parameters
        ----------
        factor : float
            The multiple of FWHM used to define the limits.
        """
        x0 = self.x_0
        dx = factor * self.fwhm

        return (x0 - dx, x0 + dx)

class WavelengthFromGratingEquation(Model):
    r""" Solve the Grating Dispersion Law for the wavelength.

    .. Note:: This form of the equation can be used for paraxial
    (small angle approximation) as well as oblique incident angles.
    With paraxial systems the inputs are sin of the angles and it
    transforms to :math:`\sin(alpha_in) + \sin(alpha_out) / m * d` .
    With oblique angles the inputs are the direction cosines of the
    angles.

    Parameters
    ----------
    groove_density : int
        Grating ruling density in units of 1/length.
    spectral_order : int
        Spectral order.

    Examples
    --------
    >>> from astropy.modeling.models import math
    >>> model = WavelengthFromGratingEquation(groove_density=20000*1/u.m, spectral_order=-1)
    >>> alpha_in = (math.Deg2radUfunc() | math.SinUfunc())(.0001 * u.deg)
    >>> alpha_out = (math.Deg2radUfunc() | math.SinUfunc())(.0001 * u.deg)
    >>> lam = model(alpha_in, alpha_out)
    >>> print(lam)
    1.745241985530844e-06 m
    """

    _separable = False

    linear = False

    n_inputs = 2
    n_outputs = 1

    groove_density = Parameter(default=1)
    """ Grating ruling density in units of 1/m."""
    spectral_order = Parameter(default=1)
    """ Spectral order."""

    def __init__(self, groove_density, spectral_order, **kwargs):
        super().__init__(groove_density=groove_density,
                         spectral_order=spectral_order, **kwargs)
        self.inputs = ("alpha_in", "alpha_out")
        """ Sine function of the angles or the direction cosines."""
        self.outputs = ("wavelength",)
        """ Wavelength."""

    def evaluate(self, alpha_in, alpha_out, groove_density, spectral_order):
        return (alpha_in + alpha_out) / (groove_density * spectral_order)

    @property
    def return_units(self):
        if self.groove_density.unit is None:
            return None
        else:
            return {'wavelength': u.Unit(1 / self.groove_density.unit)}


class AnglesFromGratingEquation3D(Model):
    """
    Solve the 3D Grating Dispersion Law for the refracted angle.

    Parameters
    ----------
    groove_density : int
        Grating ruling density in units of 1/m.
    order : int
        Spectral order.


    Examples
    --------
    >>> from astropy.modeling.models import math
    >>> model = AnglesFromGratingEquation3D(groove_density=20000*1/u.m, spectral_order=-1)
    >>> alpha_in = (math.Deg2radUfunc() | math.SinUfunc())(.0001 * u.deg)
    >>> beta_in = (math.Deg2radUfunc() | math.SinUfunc())(.0001 * u.deg)
    >>> lam = 2e-6 * u.m
    >>> alpha_out, beta_out = model(lam, alpha_in, beta_in)
    >>> print(alpha_out, beta_out)
    (<Quantity -0.03999825>, <Quantity -1.74532925e-06>)
    """

    _separable = False

    linear = False

    n_inputs = 3
    n_outputs = 2

    groove_density = Parameter(default=1)
    """ Grating ruling density in units 1/ length."""

    spectral_order = Parameter(default=1)
    """ Spectral order."""

    def __init__(self, groove_density, spectral_order, **kwargs):
        super().__init__(groove_density=groove_density,
                         spectral_order=spectral_order, **kwargs)
        self.inputs = ("wavelength", "alpha_in", "beta_in")
        """ Wavelength and 2 angle coordinates going into the grating."""

        self.outputs = ("alpha_out", "beta_out")
        """ Two angles coming out of the grating. """

    def evaluate(self, wavelength, alpha_in, beta_in,
                 groove_density, spectral_order):
        if alpha_in.shape != beta_in.shape != gamma_in.shape:
            raise ValueError("Expected input arrays to have the same shape.")
        alpha_out = groove_density * spectral_order * wavelength + alpha_in
        beta_out = - beta_in
        return alpha_out, beta_out

    @property
    def return_units(self):
        if self.groove_density.unit is None:
            return None
        else:
            return {'alpha_out': u.Unit(1),
                    'beta_out': u.Unit(1)}
