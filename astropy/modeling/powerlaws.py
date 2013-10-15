# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Power law model variants
"""
from __future__ import division
import numpy as np
from .core import Parametric1DModel

__all__ = sorted(['PowerLaw1DModel', 'BrokenPowerLaw1DModel',
                  'ExponentialCutoffPowerLaw1DModel'])


class PowerLaw1DModel(Parametric1DModel):
    """
    One dimensional power law model.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the reference point
    x_0 : float
        Reference point for model amplitude
    alpha : float
        Power law index

    See Also
    --------
    BrokenPowerLaw1DModel, ExponentialCutoffPowerLaw1DModel

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude`` and :math:`\\alpha` for ``alpha``):

        .. math:: f(x) = A (x / x_0) ^ {-\\alpha}

    """
    param_names = ['amplitude', 'x_0', 'alpha']

    def __init__(self, amplitude, x_0, alpha, **constraints):
        super(PowerLaw1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, alpha):
        """
        Model function PowerLaw1D.
        """
        xx = x / x_0
        return amplitude * xx ** (-alpha)

    def deriv(self, x, amplitude, x_0, alpha):
        """
        Model derivative PowerLaw1D.
        """
        xx = x / x_0
        d_amplitude = xx ** (-alpha)
        d_x_0 = amplitude * alpha * d_amplitude / x_0
        d_alpha = -amplitude * d_amplitude * np.log(xx)
        return [d_amplitude, d_x_0, d_alpha]


class BrokenPowerLaw1DModel(Parametric1DModel):
    """
    One dimensional power law model with a break.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the break point
    x_break : float
        Break point
    alpha_1 : float
        Power law index for x < x_break
    alpha_2 : float
        Power law index for x > x_break

    See Also
    --------
    PowerLaw1DModel, ExponentialCutoffPowerLaw1DModel

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
    param_names = ['amplitude', 'x_break', 'alpha_1', 'alpha_2']

    def __init__(self, amplitude, x_break, alpha_1, alpha_2, **constraints):
        super(BrokenPowerLaw1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_break, alpha_1, alpha_2):
        """
        Model function BrokenPowerLaw1D.
        """
        alpha = np.where(x < x_break, alpha_1, alpha_2)
        xx = x / x_break
        return amplitude * xx ** (-alpha)

    def deriv(self, x, amplitude, x_break, alpha_1, alpha_2):
        """
        Model derivative BrokenPowerLaw1D.
        """
        alpha = np.where(x < x_break, alpha_1, alpha_2)
        xx = x / x_break
        d_amplitude = xx ** (-alpha)
        d_x_break = amplitude * alpha * d_amplitude / x_break
        d_alpha = -amplitude * d_amplitude * np.log(xx)
        d_alpha_1 = np.where(x < x_break, d_alpha, 0)
        d_alpha_2 = np.where(x >= x_break, d_alpha, 0)
        return [d_amplitude, d_x_break, d_alpha_1, d_alpha_2]


class ExponentialCutoffPowerLaw1DModel(Parametric1DModel):
    """
    One dimensional power law model with an exponential cutoff.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the reference point
    x_0 : float
        Reference point for model amplitude
    alpha : float
        Power law index
    x_cutoff : float
        Cutoff point

    See Also
    --------
    PowerLaw1DModel, BrokenPowerLaw1DModel

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude`` and :math:`\\alpha` for ``alpha``):

        .. math:: f(x) = A (x / x_0) ^ {-\\alpha} \\exp (-x / x_{cutoff})

    """
    param_names = ['amplitude', 'x_0', 'alpha', 'x_cutoff']

    def __init__(self, amplitude, x_0, alpha, x_cutoff, **constraints):
        super(ExponentialCutoffPowerLaw1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, alpha, x_cutoff):
        """
        Model function ExponentialCutoffPowerLaw1D.
        """
        xx = x / x_0
        return amplitude * xx ** (-alpha) * np.exp(-x / x_cutoff)

    def deriv(self, x, amplitude, x_0, alpha, x_cutoff):
        """
        Model derivative ExponentialCutoffPowerLaw1D.
        """
        xx = x / x_0
        xc = x / x_cutoff
        d_amplitude = xx ** (-alpha) * np.exp(-xc)
        d_x_0 = alpha * amplitude * d_amplitude / x_0
        d_alpha = -amplitude * d_amplitude * np.log(xx)
        d_x_cutoff = amplitude * x * d_amplitude / x_cutoff ** 2
        return [d_amplitude, d_x_0, d_alpha, d_x_cutoff]
