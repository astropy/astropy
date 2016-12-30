# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Power law model variants
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from .core import Fittable1DModel
from .parameters import Parameter


__all__ = ['PowerLaw1D', 'BrokenPowerLaw1D', 'SmoothlyBrokenPowerLaw1D', 'ExponentialCutoffPowerLaw1D',
           'LogParabola1D']


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

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=1)
    alpha = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, alpha):
        """One dimensional power law model function"""

        xx = x / x_0
        return amplitude * xx ** (-alpha)

    @staticmethod
    def fit_deriv(x, amplitude, x_0, alpha):
        """One dimensional power law derivative with respect to parameters"""

        xx = x / x_0

        d_amplitude = xx ** (-alpha)
        d_x_0 = amplitude * alpha * d_amplitude / x_0
        d_alpha = -amplitude * d_amplitude * np.log(xx)

        return [d_amplitude, d_x_0, d_alpha]


class BrokenPowerLaw1D(Fittable1DModel):
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

    amplitude = Parameter(default=1)
    x_break = Parameter(default=1)
    alpha_1 = Parameter(default=1)
    alpha_2 = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_break, alpha_1, alpha_2):
        """One dimensional broken power law model function"""

        alpha = np.where(x < x_break, alpha_1, alpha_2)
        xx = x / x_break
        return amplitude * xx ** (-alpha)

    @staticmethod
    def fit_deriv(x, amplitude, x_break, alpha_1, alpha_2):
        """One dimensional broken power law derivative with respect to parameters"""

        alpha = np.where(x < x_break, alpha_1, alpha_2)
        xx = x / x_break

        d_amplitude = xx ** (-alpha)
        d_x_break = amplitude * alpha * d_amplitude / x_break
        d_alpha = -amplitude * d_amplitude * np.log(xx)
        d_alpha_1 = np.where(x < x_break, d_alpha, 0)
        d_alpha_2 = np.where(x >= x_break, d_alpha, 0)

        return [d_amplitude, d_x_break, d_alpha_1, d_alpha_2]


class SmoothlyBrokenPowerLaw1D(Fittable1DModel):
    """One dimensional smoothly broken power law model.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the break point
    x_break : float
        Break point
    alpha_1 : float
        Power law index for x << x_break
    alpha_2 : float
        Power law index for x >> x_break
    smooth : float
        Smoothness parameter

    See Also
    --------
    BrokenPowerLaw1D

    Notes ----- Model formula (with :math:`A` for ``amplitude`` and
    :math:`\\alpha_1` for ``alpha_1`` and :math:`\\alpha_2` for
    ``alpha_2`` and :math:`C` for ``smooth`` and :math:`S` for the
    sign of ``alpha_1-alpha_2``):

        .. math::

            f(x) = A (x / x_break) ^ {-\\alpha_1}
                   \\left[
                      \\frac{1 + (x / x_{break})^{|\\alpha_1-\\alpha_2| * C}}{2}
                   \\right]}^{S / C}

    """

    amplitude = Parameter(default=1, bounds=[1.e-10, None])
    x_break = Parameter(default=1, bounds=[1.e-10, None])
    alpha_1 = Parameter(default=1, bounds=[-6, 6])
    alpha_2 = Parameter(default=1, bounds=[-6, 6])
    smooth = Parameter(default=1, bounds=[1.e-4, 1.e4])

    @staticmethod
    def evaluate(x, amplitude, x_break, alpha_1, alpha_2, smooth):
        """One dimensional smoothly broken power law model function"""

        if (amplitude <= 0):
            raise ValueError("amplitude must be positive ("+str(amplitude)+")")
        if (x_break <= 0):
            raise ValueError("x_break value must be positive ("+str(x_break)+")")
        if (smooth <= 0):
            raise ValueError("smooth value must be positive ("+str(smooth)+")")

        xx = x / x_break
        logt =  np.abs(alpha_1 - alpha_2) * smooth * np.log(xx)
        sign = np.sign(alpha_1 - alpha_2)
        if (sign == 0):
            sign = 1

        f = xx*0.

        i = np.where(logt > 30)
        if (i[0].size > 0):
            f[i] = amplitude * xx[i]**(-alpha_2) / (2.**(sign/smooth))

        i = np.where(logt < -30)
        if (i[0].size > 0):
            f[i] = amplitude * xx[i]**(-alpha_1) / (2.**(sign/smooth))

        i = np.where(np.abs(logt) <= 30)
        if (i[0].size > 0):
            t = np.exp(logt[i])
            r = (1. + t) / 2.

            f[i] = amplitude * xx[i]**(-alpha_1) * r**(sign/smooth)

        return f

    @staticmethod
    def fit_deriv(x, amplitude, x_break, alpha_1, alpha_2, smooth):
        """One dimensional smoothly broken power law derivative with respect to parameters"""

        if (amplitude <= 0):
            raise ValueError("amplitude must be positive ("+str(amplitude)+")")
        if (x_break <= 0):
            raise ValueError("x_break value must be positive ("+str(x_break)+")")
        if (smooth <= 0):
            raise ValueError("smooth value must be positive ("+str(smooth)+")")

        xx = x / x_break
        logt = np.abs(alpha_1 - alpha_2) * smooth * np.log(xx)
        sign = np.sign(alpha_1 - alpha_2)
        if (sign == 0):
            sign = 1

        f           = xx * 0.
        d_amplitude = xx * 0.
        d_x_break   = xx * 0.
        d_alpha_1   = xx * 0.
        d_alpha_2   = xx * 0.
        d_smooth    = xx * 0.

        threshold = 30
        i = np.where(logt > threshold)
        if (i[0].size > 0):
            f[i] = amplitude * xx[i]**(-alpha_2) / (2.**(sign/smooth))

            d_amplitude[i] = f[i] / amplitude
            d_x_break[i]   = f[i] * alpha_2 / x_break
            d_alpha_1[i]   = 0.
            d_alpha_2[i]   = f[i] * (-np.log(xx[i]))
            d_smooth[i]    = f[i] * np.log(2.**(sign/smooth/smooth))

        i = np.where(logt < -threshold)
        if (i[0].size > 0):
            f[i] = amplitude * xx[i]**(-alpha_1) / (2.**(sign/smooth))

            d_amplitude[i] = f[i] / amplitude
            d_x_break[i]   = f[i] * alpha_1 / x_break
            d_alpha_1[i]   = f[i] * (-np.log(xx[i]))
            d_alpha_2[i]   = 0.
            d_smooth[i]    = f[i] * np.log(2.**(sign/smooth/smooth))

        i = np.where(np.abs(logt) <= threshold)
        if (i[0].size > 0):
            t = np.exp(logt[i])
            r = (1. + t) / 2.
            f[i] = amplitude * xx[i]**(-alpha_1) * r**(sign/smooth)

            d_amplitude[i] = f[i] / amplitude
            d_x_break[i]   = f[i] * 1./x_break * (alpha_1 - (alpha_1-alpha_2) * t/2. / r)
            d_alpha_1[i]   = f[i] * np.log(xx[i]) * (-1. + t/2. / r)
            d_alpha_2[i]   = f[i] * (-np.log(xx[i]) * t/2. / r)
            d_smooth[i]    = f[i] / smooth**2. * (-sign * np.log(r) - (alpha_1-alpha_2)*smooth/(1.+t) * np.log(xx[i]) * t)

        return [d_amplitude, d_x_break, d_alpha_1, d_alpha_2, d_smooth]


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

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=1)
    alpha = Parameter(default=1)
    x_cutoff = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, alpha, x_cutoff):
        """One dimensional exponential cutoff power law model function"""

        xx = x / x_0
        return amplitude * xx ** (-alpha) * np.exp(-x / x_cutoff)

    @staticmethod
    def fit_deriv(x, amplitude, x_0, alpha, x_cutoff):
        """One dimensional exponential cutoff power law derivative with respect to parameters"""

        xx = x / x_0
        xc = x / x_cutoff

        d_amplitude = xx ** (-alpha) * np.exp(-xc)
        d_x_0 = alpha * amplitude * d_amplitude / x_0
        d_alpha = -amplitude * d_amplitude * np.log(xx)
        d_x_cutoff = amplitude * x * d_amplitude / x_cutoff ** 2

        return [d_amplitude, d_x_0, d_alpha, d_x_cutoff]


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
    Model formula (with :math:`A` for ``amplitude`` and :math:`\\alpha` for ``alpha`` and :math:`\\beta` for ``beta``):

        .. math:: f(x) = A \\left(\\frac{x}{x_{0}}\\right)^{- \\alpha - \\beta \\log{\\left (\\frac{x}{x_{0}} \\right )}}

    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=1)
    alpha = Parameter(default=1)
    beta = Parameter(default=0)

    @staticmethod
    def evaluate(x, amplitude, x_0, alpha, beta):
        """One dimensional log parabola model function"""

        xx = x / x_0
        exponent = -alpha - beta * np.log(xx)
        return amplitude * xx ** exponent

    @staticmethod
    def fit_deriv(x, amplitude, x_0, alpha, beta):
        """One dimensional log parabola derivative with respect to parameters"""

        xx = x / x_0
        log_xx = np.log(xx)
        exponent = -alpha - beta * log_xx

        d_amplitude = xx ** exponent
        d_beta = -amplitude * d_amplitude * log_xx ** 2
        d_x_0 = amplitude * d_amplitude * (beta * log_xx / x_0 - exponent / x_0)
        d_alpha = -amplitude * d_amplitude * log_xx
        return [d_amplitude, d_x_0, d_alpha, d_beta]
