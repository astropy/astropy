# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Power law model variants
"""
from __future__ import division
import numpy as np
from .core import Parametric1DModel

__all__ = sorted(['PowerLaw1DModel'])


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
    ExponentialCutoffPowerLaw1DModel

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
