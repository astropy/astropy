# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Power law model variants.
"""

# pylint: disable=invalid-name
import numpy as np

from astropy.units import Magnitude, Quantity, UnitsError, dimensionless_unscaled, mag

from .core import Fittable1DModel
from .parameters import InputParameterError, Parameter

__all__ = [
    "BrokenPowerLaw1D",
    "ExponentialCutoffPowerLaw1D",
    "LogParabola1D",
    "PowerLaw1D",
    "Schechter1D",
    "SmoothlyBrokenPowerLaw1D",
]


class PowerLaw1D(Fittable1DModel):
    """
    One dimensional power law model.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the reference point
    x_0 : float
        Reference point
    alpha : float
        Power law index

    See Also
    --------
    BrokenPowerLaw1D, ExponentialCutoffPowerLaw1D, LogParabola1D

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude`` and :math:`\\alpha` for ``alpha``):

        .. math:: f(x) = A (x / x_0) ^ {-\\alpha}

    """

    amplitude = Parameter(default=1, description="Peak value at the reference point")
    x_0 = Parameter(default=1, description="Reference point")
    alpha = Parameter(default=1, description="Power law index")

    @staticmethod
    def evaluate(x, amplitude, x_0, alpha):
        """One dimensional power law model function."""
        xx = x / x_0
        return amplitude * xx ** (-alpha)

    @staticmethod
    def fit_deriv(x, amplitude, x_0, alpha):
        """One dimensional power law derivative with respect to parameters."""
        xx = x / x_0

        d_amplitude = xx ** (-alpha)
        d_x_0 = amplitude * alpha * d_amplitude / x_0
        d_alpha = -amplitude * d_amplitude * np.log(xx)

        return [d_amplitude, d_x_0, d_alpha]

    @property
    def input_units(self):
        if self.x_0.input_unit is None:
            return None
        return {self.inputs[0]: self.x_0.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }


class BrokenPowerLaw1D(Fittable1DModel):
    """
    One dimensional power law model with a break.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the break point.
    x_break : float
        Break point.
    alpha_1 : float
        Power law index for x < x_break.
    alpha_2 : float
        Power law index for x > x_break.

    See Also
    --------
    PowerLaw1D, ExponentialCutoffPowerLaw1D, LogParabola1D

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude`` and :math:`\\alpha_1`
    for ``alpha_1`` and :math:`\\alpha_2` for ``alpha_2``):

        .. math::

            f(x) = \\left \\{
                     \\begin{array}{ll}
                       A (x / x_{break}) ^ {-\\alpha_1} & : x < x_{break} \\\\
                       A (x / x_{break}) ^ {-\\alpha_2} & :  x > x_{break} \\\\
                     \\end{array}
                   \\right.
    """

    amplitude = Parameter(default=1, description="Peak value at break point")
    x_break = Parameter(default=1, description="Break point")
    alpha_1 = Parameter(default=1, description="Power law index before break point")
    alpha_2 = Parameter(default=1, description="Power law index after break point")

    @staticmethod
    def evaluate(x, amplitude, x_break, alpha_1, alpha_2):
        """One dimensional broken power law model function."""
        alpha = np.where(x < x_break, alpha_1, alpha_2)
        xx = x / x_break
        return amplitude * xx ** (-alpha)

    @staticmethod
    def fit_deriv(x, amplitude, x_break, alpha_1, alpha_2):
        """One dimensional broken power law derivative with respect to parameters."""
        alpha = np.where(x < x_break, alpha_1, alpha_2)
        xx = x / x_break

        d_amplitude = xx ** (-alpha)
        d_x_break = amplitude * alpha * d_amplitude / x_break
        d_alpha = -amplitude * d_amplitude * np.log(xx)
        d_alpha_1 = np.where(x < x_break, d_alpha, 0)
        d_alpha_2 = np.where(x >= x_break, d_alpha, 0)

        return [d_amplitude, d_x_break, d_alpha_1, d_alpha_2]

    @property
    def input_units(self):
        if self.x_break.input_unit is None:
            return None
        return {self.inputs[0]: self.x_break.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_break": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }


class SmoothlyBrokenPowerLaw1D(Fittable1DModel):
    """One dimensional smoothly broken power law model.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the break point.
    x_break : float
        Break point.
    alpha_1 : float
        Power law index for ``x << x_break``.
    alpha_2 : float
        Power law index for ``x >> x_break``.
    delta : float
        Smoothness parameter.

    See Also
    --------
    BrokenPowerLaw1D

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude``, :math:`x_b` for
    ``x_break``, :math:`\\alpha_1` for ``alpha_1``,
    :math:`\\alpha_2` for ``alpha_2`` and :math:`\\Delta` for
    ``delta``):

        .. math::

            f(x) = A \\left( \\frac{x}{x_b} \\right) ^ {-\\alpha_1}
                   \\left\\{
                      \\frac{1}{2}
                      \\left[
                        1 + \\left( \\frac{x}{x_b}\\right)^{1 / \\Delta}
                      \\right]
                   \\right\\}^{(\\alpha_1 - \\alpha_2) \\Delta}


    The change of slope occurs between the values :math:`x_1`
    and :math:`x_2` such that:

        .. math::
            \\log_{10} \\frac{x_2}{x_b} = \\log_{10} \\frac{x_b}{x_1}
            \\sim \\Delta


    At values :math:`x \\lesssim x_1` and :math:`x \\gtrsim x_2` the
    model is approximately a simple power law with index
    :math:`\\alpha_1` and :math:`\\alpha_2` respectively.  The two
    power laws are smoothly joined at values :math:`x_1 < x < x_2`,
    hence the :math:`\\Delta` parameter sets the "smoothness" of the
    slope change.

    The ``delta`` parameter is bounded to values greater than 1e-3
    (corresponding to :math:`x_2 / x_1 \\gtrsim 1.002`) to avoid
    overflow errors.

    The ``amplitude`` parameter is bounded to positive values since
    this model is typically used to represent positive quantities.


    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.modeling import models

        x = np.logspace(0.7, 2.3, 500)
        f = models.SmoothlyBrokenPowerLaw1D(amplitude=1, x_break=20,
                                            alpha_1=-2, alpha_2=2)

        plt.figure()
        plt.title("amplitude=1, x_break=20, alpha_1=-2, alpha_2=2")

        f.delta = 0.5
        plt.loglog(x, f(x), '--', label='delta=0.5')

        f.delta = 0.3
        plt.loglog(x, f(x), '-.', label='delta=0.3')

        f.delta = 0.1
        plt.loglog(x, f(x), label='delta=0.1')

        plt.axis([x.min(), x.max(), 0.1, 1.1])
        plt.legend(loc='lower center')
        plt.grid(True)
        plt.show()

    """

    amplitude = Parameter(
        default=1, min=0, description="Peak value at break point", mag=True
    )
    x_break = Parameter(default=1, description="Break point")
    alpha_1 = Parameter(default=-2, description="Power law index before break point")
    alpha_2 = Parameter(default=2, description="Power law index after break point")
    delta = Parameter(default=1, min=1.0e-3, description="Smoothness Parameter")

    def _amplitude_validator(self, value):
        if np.any(value <= 0):
            raise InputParameterError("amplitude parameter must be > 0")

    amplitude._validator = _amplitude_validator

    def _delta_validator(self, value):
        if np.any(value < 0.001):
            raise InputParameterError("delta parameter must be >= 0.001")

    delta._validator = _delta_validator

    @staticmethod
    def evaluate(x, amplitude, x_break, alpha_1, alpha_2, delta):
        """One dimensional smoothly broken power law model function."""
        # Pre-calculate `x/x_b`
        xx = x / x_break

        # Initialize the return value
        f = np.zeros_like(xx, subok=False)

        if isinstance(amplitude, Quantity):
            return_unit = amplitude.unit
            amplitude = amplitude.value
        else:
            return_unit = None

        # The quantity `t = (x / x_b)^(1 / delta)` can become quite
        # large.  To avoid overflow errors we will start by calculating
        # its natural logarithm:
        logt = np.log(xx) / delta

        # When `t >> 1` or `t << 1` we don't actually need to compute
        # the `t` value since the main formula (see docstring) can be
        # significantly simplified by neglecting `1` or `t`
        # respectively.  In the following we will check whether `t` is
        # much greater, much smaller, or comparable to 1 by comparing
        # the `logt` value with an appropriate threshold.
        threshold = 30  # corresponding to exp(30) ~ 1e13
        i = logt > threshold
        if i.max():
            # In this case the main formula reduces to a simple power
            # law with index `alpha_2`.
            f[i] = (
                amplitude * xx[i] ** (-alpha_2) / (2.0 ** ((alpha_1 - alpha_2) * delta))
            )

        i = logt < -threshold
        if i.max():
            # In this case the main formula reduces to a simple power
            # law with index `alpha_1`.
            f[i] = (
                amplitude * xx[i] ** (-alpha_1) / (2.0 ** ((alpha_1 - alpha_2) * delta))
            )

        i = np.abs(logt) <= threshold
        if i.max():
            # In this case the `t` value is "comparable" to 1, hence we
            # we will evaluate the whole formula.
            t = np.exp(logt[i])
            r = (1.0 + t) / 2.0
            f[i] = amplitude * xx[i] ** (-alpha_1) * r ** ((alpha_1 - alpha_2) * delta)

        if return_unit:
            return Quantity(f, unit=return_unit, copy=False, subok=True)
        return f

    @staticmethod
    def fit_deriv(x, amplitude, x_break, alpha_1, alpha_2, delta):
        """One dimensional smoothly broken power law derivative with respect
        to parameters.
        """
        # Pre-calculate `x_b` and `x/x_b` and `logt` (see comments in
        # SmoothlyBrokenPowerLaw1D.evaluate)
        xx = x / x_break
        logt = np.log(xx) / delta

        # Initialize the return values
        f = np.zeros_like(xx)
        d_amplitude = np.zeros_like(xx)
        d_x_break = np.zeros_like(xx)
        d_alpha_1 = np.zeros_like(xx)
        d_alpha_2 = np.zeros_like(xx)
        d_delta = np.zeros_like(xx)

        threshold = 30  # (see comments in SmoothlyBrokenPowerLaw1D.evaluate)
        i = logt > threshold
        if i.max():
            f[i] = (
                amplitude * xx[i] ** (-alpha_2) / (2.0 ** ((alpha_1 - alpha_2) * delta))
            )

            d_amplitude[i] = f[i] / amplitude
            d_x_break[i] = f[i] * alpha_2 / x_break
            d_alpha_1[i] = f[i] * (-delta * np.log(2))
            d_alpha_2[i] = f[i] * (-np.log(xx[i]) + delta * np.log(2))
            d_delta[i] = f[i] * (-(alpha_1 - alpha_2) * np.log(2))

        i = logt < -threshold
        if i.max():
            f[i] = (
                amplitude * xx[i] ** (-alpha_1) / (2.0 ** ((alpha_1 - alpha_2) * delta))
            )

            d_amplitude[i] = f[i] / amplitude
            d_x_break[i] = f[i] * alpha_1 / x_break
            d_alpha_1[i] = f[i] * (-np.log(xx[i]) - delta * np.log(2))
            d_alpha_2[i] = f[i] * delta * np.log(2)
            d_delta[i] = f[i] * (-(alpha_1 - alpha_2) * np.log(2))

        i = np.abs(logt) <= threshold
        if i.max():
            t = np.exp(logt[i])
            r = (1.0 + t) / 2.0
            f[i] = amplitude * xx[i] ** (-alpha_1) * r ** ((alpha_1 - alpha_2) * delta)

            d_amplitude[i] = f[i] / amplitude
            d_x_break[i] = (
                f[i] * (alpha_1 - (alpha_1 - alpha_2) * t / 2.0 / r) / x_break
            )
            d_alpha_1[i] = f[i] * (-np.log(xx[i]) + delta * np.log(r))
            d_alpha_2[i] = f[i] * (-delta * np.log(r))
            d_delta[i] = (
                f[i]
                * (alpha_1 - alpha_2)
                * (np.log(r) - t / (1.0 + t) / delta * np.log(xx[i]))
            )

        return [d_amplitude, d_x_break, d_alpha_1, d_alpha_2, d_delta]

    @property
    def input_units(self):
        if self.x_break.input_unit is None:
            return None
        return {self.inputs[0]: self.x_break.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_break": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }


class ExponentialCutoffPowerLaw1D(Fittable1DModel):
    """
    One dimensional power law model with an exponential cutoff.

    Parameters
    ----------
    amplitude : float
        Model amplitude
    x_0 : float
        Reference point
    alpha : float
        Power law index
    x_cutoff : float
        Cutoff point

    See Also
    --------
    PowerLaw1D, BrokenPowerLaw1D, LogParabola1D

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude`` and :math:`\\alpha` for ``alpha``):

        .. math:: f(x) = A (x / x_0) ^ {-\\alpha} \\exp (-x / x_{cutoff})

    """

    amplitude = Parameter(default=1, description="Peak value of model")
    x_0 = Parameter(default=1, description="Reference point")
    alpha = Parameter(default=1, description="Power law index")
    x_cutoff = Parameter(default=1, description="Cutoff point")

    @staticmethod
    def evaluate(x, amplitude, x_0, alpha, x_cutoff):
        """One dimensional exponential cutoff power law model function."""
        xx = x / x_0
        return amplitude * xx ** (-alpha) * np.exp(-x / x_cutoff)

    @staticmethod
    def fit_deriv(x, amplitude, x_0, alpha, x_cutoff):
        """
        One dimensional exponential cutoff power law derivative with respect to parameters.
        """
        xx = x / x_0
        xc = x / x_cutoff

        d_amplitude = xx ** (-alpha) * np.exp(-xc)
        d_x_0 = alpha * amplitude * d_amplitude / x_0
        d_alpha = -amplitude * d_amplitude * np.log(xx)
        d_x_cutoff = amplitude * x * d_amplitude / x_cutoff**2

        return [d_amplitude, d_x_0, d_alpha, d_x_cutoff]

    @property
    def input_units(self):
        if self.x_0.input_unit is None:
            return None
        return {self.inputs[0]: self.x_0.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "x_cutoff": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }


class LogParabola1D(Fittable1DModel):
    """
    One dimensional log parabola model (sometimes called curved power law).

    Parameters
    ----------
    amplitude : float
        Model amplitude
    x_0 : float
        Reference point
    alpha : float
        Power law index
    beta : float
        Power law curvature

    See Also
    --------
    PowerLaw1D, BrokenPowerLaw1D, ExponentialCutoffPowerLaw1D

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude`` and
    :math:`\\alpha` for ``alpha`` and :math:`\\beta` for ``beta``):

        .. math:: f(x) = A \\left(
                \\frac{x}{x_{0}}\\right)^{- \\alpha - \\beta \\log{\\left (\\frac{x}{x_{0}}
            \\right )}}

    """

    amplitude = Parameter(default=1, description="Peak value of model")
    x_0 = Parameter(default=1, description="Reference point")
    alpha = Parameter(default=1, description="Power law index")
    beta = Parameter(default=0, description="Power law curvature")

    @staticmethod
    def evaluate(x, amplitude, x_0, alpha, beta):
        """One dimensional log parabola model function."""
        xx = x / x_0
        exponent = -alpha - beta * np.log(xx)
        return amplitude * xx**exponent

    @staticmethod
    def fit_deriv(x, amplitude, x_0, alpha, beta):
        """One dimensional log parabola derivative with respect to parameters."""
        xx = x / x_0
        log_xx = np.log(xx)
        exponent = -alpha - beta * log_xx

        d_amplitude = xx**exponent
        d_beta = -amplitude * d_amplitude * log_xx**2
        d_x_0 = amplitude * d_amplitude * (beta * log_xx / x_0 - exponent / x_0)
        d_alpha = -amplitude * d_amplitude * log_xx
        return [d_amplitude, d_x_0, d_alpha, d_beta]

    @property
    def input_units(self):
        if self.x_0.input_unit is None:
            return None
        return {self.inputs[0]: self.x_0.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }


class Schechter1D(Fittable1DModel):
    r"""
    Schechter luminosity function (`Schechter 1976
    <https://ui.adsabs.harvard.edu/abs/1976ApJ...203..297S/abstract>`_),
    parameterized in terms of magnitudes.

    Parameters
    ----------
    phi_star : float
        The normalization factor in units of number density.

    m_star : float
        The characteristic magnitude where the power-law form of the
        function cuts off.

    alpha : float
        The power law index, also known as the faint-end slope. Must not
        have units.

    See Also
    --------
    PowerLaw1D, ExponentialCutoffPowerLaw1D, BrokenPowerLaw1D

    Notes
    -----
    Model formula (with :math:`\phi^{*}` for ``phi_star``, :math:`M^{*}`
    for ``m_star``, and :math:`\alpha` for ``alpha``):

    .. math::

        n(M) \ dM = (0.4 \ln 10) \ \phi^{*} \
            [{10^{0.4 (M^{*} - M)}}]^{\alpha + 1} \
            \exp{[-10^{0.4 (M^{*} - M)}]} \ dM

    ``phi_star`` is the normalization factor in units of number density.
    ``m_star`` is the characteristic magnitude where the power-law form
    of the function cuts off into the exponential form. ``alpha`` is
    the power-law index, defining the faint-end slope of the luminosity
    function.

    Examples
    --------
    .. plot::
        :include-source:

        from astropy.modeling.models import Schechter1D
        import astropy.units as u
        import matplotlib.pyplot as plt
        import numpy as np

        phi_star = 4.3e-4 * (u.Mpc ** -3)
        m_star = -20.26
        alpha = -1.98
        model = Schechter1D(phi_star, m_star, alpha)
        mag = np.linspace(-25, -17)

        fig, ax = plt.subplots()
        ax.plot(mag, model(mag))
        ax.set_yscale('log')
        ax.set_xlim(-22.6, -17)
        ax.set_ylim(1.e-7, 1.e-2)
        ax.set_xlabel('$M_{UV}$')
        ax.set_ylabel(r'$\phi$ [mag$^{-1}$ Mpc$^{-3}]$')

    References
    ----------
    .. [1] Schechter 1976; ApJ 203, 297
           (https://ui.adsabs.harvard.edu/abs/1976ApJ...203..297S/abstract)

    .. [2] `Luminosity function <https://en.wikipedia.org/wiki/Luminosity_function_(astronomy)>`_
    """

    phi_star = Parameter(
        default=1.0, description="Normalization factor in units of number density"
    )
    m_star = Parameter(default=-20.0, description="Characteristic magnitude", mag=True)
    alpha = Parameter(default=-1.0, description="Faint-end slope")

    @staticmethod
    def _factor(magnitude, m_star):
        factor_exp = magnitude - m_star

        if isinstance(factor_exp, Quantity):
            if factor_exp.unit == mag:
                factor_exp = Magnitude(factor_exp.value, unit=mag)

                return factor_exp.to(dimensionless_unscaled)
            else:
                raise UnitsError(
                    "The units of magnitude and m_star must be a magnitude"
                )
        else:
            return 10 ** (-0.4 * factor_exp)

    def evaluate(self, mag, phi_star, m_star, alpha):
        """Schechter luminosity function model function."""
        factor = self._factor(mag, m_star)

        return 0.4 * np.log(10) * phi_star * factor ** (alpha + 1) * np.exp(-factor)

    def fit_deriv(self, mag, phi_star, m_star, alpha):
        """
        Schechter luminosity function derivative with respect to
        parameters.
        """
        factor = self._factor(mag, m_star)

        d_phi_star = 0.4 * np.log(10) * factor ** (alpha + 1) * np.exp(-factor)
        func = phi_star * d_phi_star
        d_m_star = (alpha + 1) * 0.4 * np.log(10) * func - (
            0.4 * np.log(10) * func * factor
        )
        d_alpha = func * np.log(factor)

        return [d_phi_star, d_m_star, d_alpha]

    @property
    def input_units(self):
        if self.m_star.input_unit is None:
            return None
        return {self.inputs[0]: self.m_star.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "m_star": inputs_unit[self.inputs[0]],
            "phi_star": outputs_unit[self.outputs[0]],
        }
