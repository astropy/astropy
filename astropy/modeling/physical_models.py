# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Models that have physical origins.
"""

import warnings
from collections import OrderedDict

import numpy as np

from .core import Fittable1DModel
from .parameters import Parameter
from astropy import constants as const
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["BlackBody"]

# Units
FNU = u.erg / (u.cm ** 2 * u.s * u.Hz)
FLAM = u.erg / (u.cm ** 2 * u.s * u.AA)
SNU = u.erg / (u.cm ** 2 * u.s * u.Hz * u.sr)
SLAM = u.erg / (u.cm ** 2 * u.s * u.AA * u.sr)

# Some platform implementations of expm1() are buggy and Numpy uses
# them anyways--the bug is that on certain large inputs it returns
# NaN instead of INF like it should (it should only return NaN on a
# NaN input
# See https://github.com/astropy/astropy/issues/4171
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    _has_buggy_expm1 = np.isnan(np.expm1(1000)) or np.isnan(np.expm1(1e10))


class BlackBody(Fittable1DModel):
    """
    One dimensional blackbody model.

    Parameters
    ----------
    temperature : :class:`~astropy.units.Quantity`
        Blackbody temperature.
    scale : :class:`~astropy.units.Quantity`
        Scale factor

    Notes
    -----

    Model formula:

        .. math:: B_{\\nu}(T) = A \\frac{2 h \\nu^{3} / c^{2}}{exp(h \\nu / k T) - 1}

    Examples
    --------
    >>> from astropy.modeling import models
    >>> from astropy import units as u
    >>> bb = models.BlackBody()
    >>> bb(6000 * u.AA)  # doctest: +FLOAT_CMP
    <Quantity 1.3585381201978953e-15 erg / (cm2 Hz s)>

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import BlackBody
        from astropy.modeling.physical_models import SLAM
        from astropy import units as u
        from astropy.visualization import quantity_support

        bb = BlackBody(temperature=5778*u.K)
        wav = np.arange(1000, 110000) * u.AA
        flux = bb(wav)

        with quantity_support():
            plt.figure()
            plt.semilogx(wav, flux)
            plt.axvline(bb.lambda_max.to(u.AA).value, ls='--')
            plt.show()
    """

    # We parametrize this model with a temperature and a scale.
    temperature = Parameter(default=5000, min=0, unit=u.K)
    scale = Parameter(default=1, min=0)

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
        # Convert to units for calculations, also force double precision
        with u.add_enabled_equivalencies(u.spectral() + u.temperature()):
            freq = u.Quantity(x, u.Hz, dtype=np.float64)
            temp = u.Quantity(temperature, u.K, dtype=np.float64)

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

        if _has_buggy_expm1:
            # Replace incorrect nan results with infs--any result of 'nan' is
            # incorrect unless the input (in log_boltz) happened to be nan to begin
            # with.  (As noted in #4393 ideally this would be replaced by a version
            # of expm1 that doesn't have this bug, rather than fixing incorrect
            # results after the fact...)
            boltzm1_nans = np.isnan(boltzm1)
            if np.any(boltzm1_nans):
                if boltzm1.isscalar and not np.isnan(log_boltz):
                    boltzm1 = np.inf
                else:
                    boltzm1[np.where(~np.isnan(log_boltz) & boltzm1_nans)] = np.inf

        # Calculate blackbody flux
        bb_nu = 2.0 * const.h * freq ** 3 / (const.c ** 2 * boltzm1)
        flux = bb_nu.to(FNU, u.spectral_density(freq))

        # Add per steradian to output flux unit
        fnu = scale * flux / u.sr

        # If the bolometric_flux parameter has no unit, we should drop the /Hz
        # and return a unitless value. This occurs for instance during fitting,
        # since we drop the units temporarily.
        if hasattr(scale, "unit"):
            return fnu
        else:
            return fnu

    @property
    def input_units(self):
        # The input units are those of the 'x' value, which should always be
        # Hz. Because we do this, and because input_units_allow_dimensionless
        # is set to True, dimensionless values are assumed to be in Hz.
        return {"x": u.Hz}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict(
            [("temperature", u.K), ("bolometric_flux", outputs_unit["y"] * u.Hz)]
        )

    @property
    def lambda_max(self):
        """Peak wavelength when the curve is expressed as power density."""
        return const.b_wien / self.temperature
