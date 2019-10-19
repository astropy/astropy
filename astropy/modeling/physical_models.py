# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Models that have physical origins.
"""

import warnings

import numpy as np

from .core import Fittable1DModel
from .parameters import Parameter
from astropy import constants as const
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["BlackBody"]


class BlackBody(Fittable1DModel):
    """
    One dimensional blackbody model.

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
        native_bolflux = self.scale.value * const.sigma_sb * self.temperature ** 4 / np.pi
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
