# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Mathematical models."""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import collections

from textwrap import dedent

import numpy as np

from .core import (Fittable1DModel, Fittable2DModel,
                   Model, format_input, ModelDefinitionError)
from .parameters import Parameter, InputParameterError
from ..utils import find_current_module
from ..extern import six

__all__ = sorted([
    'AiryDisk2D', 'Beta1D', 'Beta2D', 'Box1D',
    'Box2D', 'Const1D', 'Const2D', 'Disk2D',
    'Gaussian1D', 'GaussianAbsorption1D', 'Gaussian2D', 'Linear1D',
    'Lorentz1D', 'MexicanHat1D', 'MexicanHat2D', 'Scale', 'Redshift', 'Shift',
    'Sine1D', 'Trapezoid1D', 'TrapezoidDisk2D', 'Ring2D',
    'custom_model_1d'
])


class Gaussian1D(Fittable1DModel):
    """
    One dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the Gaussian.
    mean : float
        Mean of the Gaussian.
    stddev : float
        Standard deviation of the Gaussian.

    Notes
    -----

    Model formula:

        .. math:: f(x) = A e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}

    Examples
    --------
    >>> from astropy.modeling import models
    >>> def tie_center(model):
    ...         mean = 50 * model.stddev
    ...         return mean
    >>> tied_parameters = {'mean': tie_center}

    Specify that 'mean' is a tied parameter in one of two ways:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                             tied=tied_parameters)

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.mean.tied
    False
    >>> g1.mean.tied = tie_center
    >>> g1.mean.tied
    <function tie_center at 0x...>

    Fixed parameters:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                        fixed={'stddev': True})
    >>> g1.stddev.fixed
    True

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.stddev.fixed
    False
    >>> g1.stddev.fixed = True
    >>> g1.stddev.fixed
    True

    See Also
    --------
    Gaussian2D, Box1D, Beta1D, Lorentz1D
    """

    amplitude = Parameter()
    mean = Parameter()
    stddev = Parameter()

    def __init__(self, amplitude, mean, stddev, **kwargs):
        super(Gaussian1D, self).__init__(
            amplitude=amplitude, mean=mean, stddev=stddev, **kwargs)

    @staticmethod
    def eval(x, amplitude, mean, stddev):
        """
        Gaussian1D model function.
        """
        return amplitude * np.exp(- 0.5 * (x - mean) ** 2 / stddev ** 2)

    @staticmethod
    def fit_deriv(x, amplitude, mean, stddev):
        """
        Gaussian1D model function derivatives.
        """

        d_amplitude = np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
        d_mean = amplitude * d_amplitude * (x - mean) / stddev ** 2
        d_stddev = amplitude * d_amplitude * (x - mean) ** 2 / stddev ** 3
        return [d_amplitude, d_mean, d_stddev]


class GaussianAbsorption1D(Fittable1DModel):
    """
    One dimensional Gaussian absorption line model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian absorption.
    mean : float
        Mean of the gaussian.
    stddev : float
        Standard deviation of the gaussian.

    Notes
    -----

    Model formula:

        .. math:: f(x) = 1 - A e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}

    See Also
    --------
    Gaussian1D
    """

    amplitude = Parameter()
    mean = Parameter()
    stddev = Parameter()

    def __init__(self, amplitude, mean, stddev, **kwargs):
        super(GaussianAbsorption1D, self).__init__(
            amplitude=amplitude, mean=mean, stddev=stddev, **kwargs)

    @staticmethod
    def eval(x, amplitude, mean, stddev):
        """
        GaussianAbsorption1D model function.
        """
        return 1.0 - Gaussian1D.eval(x, amplitude, mean, stddev)

    @staticmethod
    def fit_deriv(x, amplitude, mean, stddev):
        """
        GaussianAbsorption1D model function derivatives.
        """
        import operator
        return list(six.moves.map(
            operator.neg, Gaussian1D.fit_deriv(x, amplitude, mean, stddev)))


class Gaussian2D(Fittable2DModel):
    """
    Two dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the Gaussian.
    x_mean : float
        Mean of the Gaussian in x.
    y_mean : float
        Mean of the Gaussian in y.
    x_stddev : float
        Standard deviation of the Gaussian in x.
        ``x_stddev`` and ``y_stddev`` must be specified unless a covariance
        matrix (``cov_matrix``) is input.
    y_stddev : float
        Standard deviation of the Gaussian in y.
        ``x_stddev`` and ``y_stddev`` must be specified unless a covariance
        matrix (``cov_matrix``) is input.
    theta : float, optional
        Rotation angle in radians. The rotation angle increases
        counterclockwise.
    cov_matrix : ndarray, optional
        A 2x2 covariance matrix. If specified, overrides the ``x_stddev``,
        ``y_stddev``, and ``theta`` specification.

    Notes
    -----
    Model formula:

        .. math::

            f(x, y) = A e^{-a\\left(x - x_{0}\\right)^{2}  -b\\left(x - x_{0}\\right)
            \\left(y - y_{0}\\right)  -c\\left(y - y_{0}\\right)^{2}}

    Using the following definitions:

        .. math::
            a = \\left(\\frac{\\cos^{2}{\\left (\\theta \\right )}}{2 \\sigma_{x}^{2}} +
            \\frac{\\sin^{2}{\\left (\\theta \\right )}}{2 \\sigma_{y}^{2}}\\right)

            b = \\left(\\frac{\\sin{\\left (2 \\theta \\right )}}{2 \\sigma_{x}^{2}} -
            \\frac{\\sin{\\left (2 \\theta \\right )}}{2 \\sigma_{y}^{2}}\\right)

            c = \\left(\\frac{\\sin^{2}{\\left (\\theta \\right )}}{2 \\sigma_{x}^{2}} +
            \\frac{\\cos^{2}{\\left (\\theta \\right )}}{2 \\sigma_{y}^{2}}\\right)

    If using a ``cov_matrix``, the model is of the form:
        .. math::
            f(x, y) = A e^{-0.5 \\left(\\vec{x} - \\vec{x}_{0}\\right)^{T} \\Sigma^{-1} \\left(\\vec{x} - \\vec{x}_{0}\\right)}

    where :math:`\\vec{x} = [x, y]`, :math:`\\vec{x}_{0} = [x_{0}, y_{0}]`,
    and :math:`\\Sigma` is the covariance matrix:

        .. math::
            \\Sigma = \\left(\\begin{array}{ccc}
            \\sigma_x^2               & \\rho \\sigma_x \\sigma_y \\\\
            \\rho \\sigma_x \\sigma_y & \\sigma_y^2
            \end{array}\\right)

    :math:`\\rho` is the correlation between ``x`` and ``y``, which should
    be between -1 and +1.  Positive correlation corresponds to a
    ``theta`` in the range 0 to 90 degrees.  Negative correlation
    corresponds to a ``theta`` in the range of 0 to -90 degrees.

    See [1]_ for more details about the 2D Gaussian function.

    See Also
    --------
    Gaussian1D, Box2D, Beta2D

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Gaussian_function
    """

    amplitude = Parameter()
    x_mean = Parameter()
    y_mean = Parameter()
    x_stddev = Parameter()
    y_stddev = Parameter()
    theta = Parameter()

    def __init__(self, amplitude, x_mean, y_mean, x_stddev=None, y_stddev=None,
                 theta=0.0, cov_matrix=None, **kwargs):
        if y_stddev is None and cov_matrix is None:
            raise InputParameterError(
                "Either x/y_stddev must be specified, or a "
                "covariance matrix.")
        elif x_stddev is None and cov_matrix is None:
            raise InputParameterError(
                "Either x/y_stddev must be specified, or a "
                "covariance matrix.")
        elif cov_matrix is not None and (x_stddev is not None or
                                         y_stddev is not None):
            raise InputParameterError(
                "Cannot specify both cov_matrix and x/y_stddev")

        # Compute principle coordinate system transformation
        elif cov_matrix is not None:
            cov_matrix = np.array(cov_matrix)
            if cov_matrix.shape != (2, 2):
                raise ValueError("Covariance matrix must be 2x2")
            eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
            x_stddev, y_stddev = np.sqrt(eig_vals)
            y_vec = eig_vecs[:, 0]
            theta = np.arctan2(y_vec[1], y_vec[0])

        super(Gaussian2D, self).__init__(
            amplitude=amplitude, x_mean=x_mean, y_mean=y_mean,
            x_stddev=x_stddev, y_stddev=y_stddev, theta=theta, **kwargs)

    @staticmethod
    def eval(x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """Two dimensional Gaussian function"""

        cost2 = np.cos(theta) ** 2
        sint2 = np.sin(theta) ** 2
        sin2t = np.sin(2. * theta)
        xstd2 = x_stddev ** 2
        ystd2 = y_stddev ** 2
        xdiff = x - x_mean
        ydiff = y - y_mean
        a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2))
        b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2))
        c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2))
        return amplitude * np.exp(-((a * xdiff ** 2) + (b * xdiff * ydiff) +
                                    (c * ydiff ** 2)))

    @staticmethod
    def fit_deriv(x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """Two dimensional Gaussian function derivative with respect to parameters"""

        cost = np.cos(theta)
        sint = np.sin(theta)
        cost2 = np.cos(theta) ** 2
        sint2 = np.sin(theta) ** 2
        cos2t = np.cos(2. * theta)
        sin2t = np.sin(2. * theta)
        xstd2 = x_stddev ** 2
        ystd2 = y_stddev ** 2
        xstd3 = x_stddev ** 3
        ystd3 = y_stddev ** 3
        xdiff = x - x_mean
        ydiff = y - y_mean
        xdiff2 = xdiff ** 2
        ydiff2 = ydiff ** 2
        a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2))
        b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2))
        c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2))
        g = amplitude * np.exp(-((a * xdiff2) + (b * xdiff * ydiff) +
                                 (c * ydiff2)))
        da_dtheta = (sint * cost * ((1. / ystd2) - (1. / xstd2)))
        da_dx_stddev = -cost2 / xstd3
        da_dy_stddev = -sint2 / ystd3
        db_dtheta = (cos2t / xstd2) - (cos2t / ystd2)
        db_dx_stddev = -sin2t / xstd3
        db_dy_stddev = sin2t / ystd3
        dc_dtheta = -da_dtheta
        dc_dx_stddev = -sint2 / xstd3
        dc_dy_stddev = -cost2 / ystd3
        dg_dA = g / amplitude
        dg_dx_mean = g * ((2. * a * xdiff) + (b * ydiff))
        dg_dy_mean = g * ((b * xdiff) + (2. * c * ydiff))
        dg_dx_stddev = g * (-(da_dx_stddev * xdiff2 +
                              db_dx_stddev * xdiff * ydiff +
                              dc_dx_stddev * ydiff2))
        dg_dy_stddev = g * (-(da_dy_stddev * xdiff2 +
                              db_dy_stddev * xdiff * ydiff +
                              dc_dy_stddev * ydiff2))
        dg_dtheta = g * (-(da_dtheta * xdiff2 +
                           db_dtheta * xdiff * ydiff +
                           dc_dtheta * ydiff2))
        return [dg_dA, dg_dx_mean, dg_dy_mean, dg_dx_stddev, dg_dy_stddev,
                dg_dtheta]


class Shift(Model):
    """
    Shift a coordinate.

    Parameters
    ----------
    offsets : float or a list of floats
        offsets to be applied to a coordinate
        if a list - each value in the list is an offset to be applied to a
        column in the input coordinate array
    """

    offsets = Parameter()

    def __init__(self, offsets, **kwargs):
        super(Shift, self).__init__(offsets, **kwargs)

    # TODO: Might need to do some work to ensure that cases like this work
    # consistently.  Should iterating over self.offsets mean iterating over its
    # parameter sets?  Maybe something like this should just work seamlessly
    def inverse(self):
        inv = self.copy()
        inv.offsets *= -1
        return inv

    @format_input
    def __call__(self, x, model_set_axis=None):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        return self.offsets + x


class Scale(Model):
    """
    Multiply a model by a factor.

    Parameters
    ----------
    factors : float or a list of floats
        scale for a coordinate
    """

    factors = Parameter()
    linear = True

    def __init__(self, factors, **kwargs):
        super(Scale, self).__init__(factors, **kwargs)

    def inverse(self):
        inv = self.copy()
        inv.factors = 1 / self.factors
        return inv

    @format_input
    def __call__(self, x, model_set_axis=None):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        return self.factors * x


class Redshift(Fittable1DModel):
    """
    One dimensional redshift model.

    Parameters
    ----------
    z : float or a list of floats
        Redshift value(s).

    Notes
    -----
    Model formula:

        .. math:: \\lambda_{obs} = (1 + z) \\lambda_{rest}

    """
    z = Parameter(description='redshift')

    def __init__(self, z, **kwargs):
        super(Redshift, self).__init__(z=z, **kwargs)

    @staticmethod
    def eval(x, z):
        """One dimensional Redshift model function"""
        return (1 + z) * x

    @staticmethod
    def fit_deriv(x, z):
        """One dimensional Redshift model derivative"""
        d_z = x
        return [d_z]

    def inverse(self):
        """Inverse Redshift model"""

        inv = self.copy()
        inv.z = 1.0 / (1.0 + self.z) - 1.0
        return inv


class Sine1D(Fittable1DModel):
    """
    One dimensional Sine model.

    Parameters
    ----------
    amplitude : float
        Oscillation amplitude
    frequency : float
        Oscillation frequency

    See Also
    --------
    Const1D, Linear1D


    Notes
    -----
    Model formula:

        .. math:: f(x) = A \\sin(2 \\pi f x)
    """

    amplitude = Parameter()
    frequency = Parameter()

    def __init__(self, amplitude, frequency, **kwargs):
        super(Sine1D, self).__init__(
            amplitude=amplitude, frequency=frequency, **kwargs)

    @staticmethod
    def eval(x, amplitude, frequency):
        """One dimensional Sine model function"""

        return amplitude * np.sin(2 * np.pi * frequency * x)

    @staticmethod
    def fit_deriv(x, amplitude, frequency):
        """One dimensional Sine model derivative"""

        d_amplitude = np.sin(2 * np.pi * frequency * x)
        d_frequency = (2 * np.pi * x * amplitude *
                       np.cos(2 * np.pi * frequency * x))
        return [d_amplitude, d_frequency]


class Linear1D(Fittable1DModel):
    """
    One dimensional Line model.

    Parameters
    ----------
    slope : float
        Slope of the straight line

    intercept : float
        Intercept of the straight line

    See Also
    --------
    Const1D

    Notes
    -----
    Model formula:

        .. math:: f(x) = a x + b
    """

    slope = Parameter()
    intercept = Parameter()
    linear = True

    def __init__(self, slope, intercept, **kwargs):
        super(Linear1D, self).__init__(
            slope=slope, intercept=intercept, **kwargs)

    @staticmethod
    def eval(x, slope, intercept):
        """One dimensional Line model function"""

        return slope * x + intercept

    @staticmethod
    def fit_deriv(x, slope, intercept):
        """One dimensional Line model derivative with respect to parameters"""

        d_slope = x
        d_intercept = np.ones_like(x)
        return [d_slope, d_intercept]


class Lorentz1D(Fittable1DModel):
    """
    One dimensional Lorentzian model.

    Parameters
    ----------
    amplitude : float
        Peak value
    x_0 : float
        Position of the peak
    fwhm : float
        Full width at half maximum

    See Also
    --------
    Gaussian1D, Box1D, MexicanHat1D

    Notes
    -----
    Model formula:

    .. math::

        f(x) = \\frac{A \\gamma^{2}}{\\gamma^{2} + \\left(x - x_{0}\\right)^{2}}
    """

    amplitude = Parameter()
    x_0 = Parameter()
    fwhm = Parameter()

    def __init__(self, amplitude, x_0, fwhm, **kwargs):
        super(Lorentz1D, self).__init__(
            amplitude=amplitude, x_0=x_0, fwhm=fwhm, **kwargs)

    @staticmethod
    def eval(x, amplitude, x_0, fwhm):
        """One dimensional Lorentzian model function"""

        return (amplitude * ((fwhm / 2.) ** 2) / ((x - x_0) ** 2 +
                (fwhm / 2.) ** 2))

    @staticmethod
    def fit_deriv(x, amplitude, x_0, fwhm):
        """One dimensional Lorentzian model derivative wiht respect to parameters"""

        d_amplitude = fwhm ** 2 / (fwhm ** 2 + (x - x_0) ** 2)
        d_x_0 = (amplitude * d_amplitude * (2 * x - 2 * x_0) /
                 (fwhm ** 2 + (x - x_0) ** 2))
        d_fwhm = 2 * amplitude * d_amplitude / fwhm * (1 - d_amplitude)
        return [d_amplitude, d_x_0, d_fwhm]


class Const1D(Fittable1DModel):
    """
    One dimensional Constant model.

    Parameters
    ----------
    amplitude : float
        Value of the constant function

    See Also
    --------
    Const2D

    Notes
    -----
    Model formula:

        .. math:: f(x) = A
    """

    amplitude = Parameter()
    linear = True

    def __init__(self, amplitude, **kwargs):
        super(Const1D, self).__init__(amplitude=amplitude, **kwargs)

    @staticmethod
    def eval(x, amplitude):
        """One dimensional Constant model function"""

        if amplitude.size == 1:
            # This is slighly faster than using ones_like and multiplying
            x = np.empty_like(x)
            x.fill(amplitude.item())
        else:
            # This case is less likely but could occur if the amplitude
            # parameter is given an array-like value
            x = amplitude * np.ones_like(x)

        return x

    @staticmethod
    def fit_deriv(x, amplitude):
        """One dimensional Constant model derivative with respect to parameters"""

        d_amplitude = np.ones_like(x)
        return [d_amplitude]


class Const2D(Fittable2DModel):
    """
    Two dimensional Constant model.

    Parameters
    ----------
    amplitude : float
        Value of the constant function

    See Also
    --------
    Const1D

    Notes
    -----
    Model formula:

        .. math:: f(x, y) = A
    """

    amplitude = Parameter()
    linear = True

    def __init__(self, amplitude, **kwargs):
        super(Const2D, self).__init__(amplitude=amplitude, **kwargs)

    @staticmethod
    def eval(x, y, amplitude):
        """Two dimensional Constant model function"""

        if amplitude.size == 1:
            # This is slighly faster than using ones_like and multiplying
            x = np.empty_like(x)
            x.fill(amplitude.item())
        else:
            # This case is less likely but could occur if the amplitude
            # parameter is given an array-like value
            x = amplitude * np.ones_like(x)

        return x


class Disk2D(Fittable2DModel):
    """
    Two dimensional radial symmetric Disk model.

    Parameters
    ----------
    amplitude : float
        Value of the disk function
    x_0 : float
        x position center of the disk
    y_0 : float
        y position center of the disk
    R_0 : float
        Radius of the disk

    See Also
    --------
    Box2D, TrapezoidDisk2D

    Notes
    -----
    Model formula:

        .. math::

            f(r) = \\left \\{
                     \\begin{array}{ll}
                       A & : r \\leq R_0 \\\\
                       0 & : r > R_0
                     \\end{array}
                   \\right.
    """

    amplitude = Parameter()
    x_0 = Parameter()
    y_0 = Parameter()
    R_0 = Parameter()

    def __init__(self, amplitude, x_0, y_0, R_0, **kwargs):
        super(Disk2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, R_0=R_0, **kwargs)

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, R_0):
        """Two dimensional Disk model function"""

        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        return np.select([rr <= R_0 ** 2], [amplitude])


class Ring2D(Fittable2DModel):

    """
    Two dimensional radial symmetric Ring model.

    Parameters
    ----------
    amplitude : float
        Value of the disk function
    x_0 : float
        x position center of the disk
    y_0 : float
        y position center of the disk
    r_in : float
        Inner radius of the ring
    width : float
        Width of the ring.
    r_out : float
        Outer Radius of the ring. Can be specified instead of width.

    See Also
    --------
    Disk2D, TrapezoidDisk2D

    Notes
    -----
    Model formula:

        .. math::

            f(r) = \\left \\{
                     \\begin{array}{ll}
                       A & : r_{in} \\leq r \\leq r_{out} \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.

    Where :math:`r_{out} = r_{in} + r_{width}`.
    """

    amplitude = Parameter()
    x_0 = Parameter()
    y_0 = Parameter()
    r_in = Parameter()
    width = Parameter()

    def __init__(self, amplitude, x_0, y_0, r_in, width=None, r_out=None,
                 **kwargs):
        if r_out is not None:
            width = r_out - r_in
        if r_out is None and width is None:
            raise ModelDefinitionError("Either specify width or r_out.")

        super(Ring2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, r_in=r_in, width=width,
            **kwargs)

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, r_in, width):
        """Two dimensional Ring model function."""

        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        r_range = np.logical_and(rr >= r_in ** 2, rr <= (r_in + width) ** 2)
        return np.select([r_range], [amplitude])


class Delta1D(Fittable1DModel):
    """One dimensional Dirac delta function."""

    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Delta2D(Fittable2DModel):
    """Two dimensional Dirac delta function."""

    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Box1D(Fittable1DModel):
    """
    One dimensional Box model.

    Parameters
    ----------
    amplitude : float
        Amplitude A
    x_0 : float
        Position of the center of the box function
    width : float
        Width of the box

    See Also
    --------
    Box2D, TrapezoidDisk2D

    Notes
    -----
    Model formula:

      .. math::

            f(x) = \\left \\{
                     \\begin{array}{ll}
                       A & : x_0 - w/2 \\geq x \\geq x_0 + w/2 \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.
    """

    amplitude = Parameter()
    x_0 = Parameter()
    width = Parameter()

    def __init__(self, amplitude, x_0, width, **kwargs):
        super(Box1D, self).__init__(
            amplitude=amplitude, x_0=x_0, width=width, **kwargs)

    @staticmethod
    def eval(x, amplitude, x_0, width):
        """One dimensional Box model function"""

        return np.select([np.logical_and(x >= x_0 - width / 2.,
                                         x <= x_0 + width / 2.)],
                         [amplitude], 0)

    @classmethod
    def fit_deriv(cls, x, amplitude, x_0, width):
        """One dimensional Box model derivative with respect to parameters"""

        d_amplitude = cls.eval(x, 1, x_0, width)
        d_x_0 = np.zeros_like(x)
        d_width = np.zeros_like(x)
        return [d_amplitude, d_x_0, d_width]


class Box2D(Fittable2DModel):
    """
    Two dimensional Box model.

    Parameters
    ----------
    amplitude : float
        Amplitude A
    x_0 : float
        x position of the center of the box function
    x_width : float
        Width in x direction of the box
    y_0 : float
        y position of the center of the box function
    y_width : float
        Width in y direction of the box

    See Also
    --------
    Box1D, Gaussian2D, Beta2D

    Notes
    -----
    Model formula:

      .. math::

            f(x, y) = \\left \\{
                     \\begin{array}{ll}
                       A & : x_0 - w_x/2 \\geq x \\geq x_0 + w_x/2 \\\\
                       A & : y_0 - w_y/2 \\geq y \\geq y_0 + w_y/2 \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.

    """

    amplitude = Parameter()
    x_0 = Parameter()
    y_0 = Parameter()
    x_width = Parameter()
    y_width = Parameter()

    def __init__(self, amplitude, x_0, y_0, x_width, y_width, **kwargs):
        super(Box2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, x_width=x_width,
            y_width=y_width, **kwargs)

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, x_width, y_width):
        """Two dimensional Box model function"""
        x_range = np.logical_and(x >= x_0 - x_width / 2.,
                                 x <= x_0 + x_width / 2.)
        y_range = np.logical_and(y >= y_0 - y_width / 2.,
                                 y <= y_0 + y_width / 2.)
        return np.select([np.logical_and(x_range, y_range)], [amplitude], 0)


class Trapezoid1D(Fittable1DModel):
    """
    One dimensional Trapezoid model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the trapezoid
    x_0 : float
        Center position of the trapezoid
    width : float
        Width of the constant part of the trapezoid.
    slope : float
        Slope of the tails of the trapezoid

    See Also
    --------
    Box1D, Gaussian1D, Beta1D
    """

    amplitude = Parameter()
    x_0 = Parameter()
    width = Parameter()
    slope = Parameter()

    def __init__(self, amplitude, x_0, width, slope, **kwargs):
        super(Trapezoid1D, self).__init__(
            amplitude=amplitude, x_0=x_0, width=width, slope=slope, **kwargs)

    @staticmethod
    def eval(x, amplitude, x_0, width, slope):
        """One dimensional Trapezoid model function"""
        # Compute the four points where the trapezoid changes slope
        # x1 <= x2 <= x3 <= x4
        x2 = x_0 - width / 2.
        x3 = x_0 + width / 2.
        x1 = x2 - amplitude / slope
        x4 = x3 + amplitude / slope

        # Compute model values in pieces between the change points
        range_a = np.logical_and(x >= x1, x < x2)
        range_b = np.logical_and(x >= x2, x < x3)
        range_c = np.logical_and(x >= x3, x < x4)
        val_a = slope * (x - x1)
        val_b = amplitude
        val_c = slope * (x4 - x)
        return np.select([range_a, range_b, range_c], [val_a, val_b, val_c])


class TrapezoidDisk2D(Fittable2DModel):
    """
    Two dimensional circular Trapezoid model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the trapezoid
    x_0 : float
        x position of the center of the trapezoid
    y_0 : float
        y position of the center of the trapezoid
    R_0 : float
        Radius of the constant part of the trapezoid.
    slope : float
        Slope of the tails of the trapezoid in x direction.

    See Also
    --------
    Disk2D, Box2D
    """

    amplitude = Parameter()
    x_0 = Parameter()
    y_0 = Parameter()
    R_0 = Parameter()
    slope = Parameter()

    def __init__(self, amplitude, x_0, y_0, R_0, slope, **kwargs):
        super(TrapezoidDisk2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, R_0=R_0, slope=slope,
            **kwargs)

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, R_0, slope):
        """Two dimensional Trapezoid Disk model function"""

        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2)
        range_1 = r <= R_0
        range_2 = np.logical_and(r > R_0,  r <= R_0 + amplitude / slope)
        val_1 = amplitude
        val_2 = amplitude + slope * (R_0 - r)
        return np.select([range_1, range_2], [val_1, val_2])


class MexicanHat1D(Fittable1DModel):
    """
    One dimensional Mexican Hat model.

    Parameters
    ----------
    amplitude : float
        Amplitude
    x_0 : float
        Position of the peak
    sigma : float
        Width of the Mexican hat

    See Also
    --------
    MexicanHat2D, Box1D, Gaussian1D, Trapezoid1D

    Notes
    -----
    Model formula:

    .. math::

        f(x) = {A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}}{\\sigma^{2}}\\right)
        e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}}

    """

    amplitude = Parameter()
    x_0 = Parameter()
    sigma = Parameter()

    def __init__(self, amplitude, x_0, sigma, **kwargs):
        super(MexicanHat1D, self).__init__(
            amplitude=amplitude, x_0=x_0, sigma=sigma, **kwargs)

    @staticmethod
    def eval(x, amplitude, x_0, sigma):
        """One dimensional Mexican Hat model function"""

        xx_ww = (x - x_0) ** 2 / (2 * sigma ** 2)
        return amplitude * (1 - 2 * xx_ww) * np.exp(-xx_ww)


class MexicanHat2D(Fittable2DModel):
    """
    Two dimensional symmetric Mexican Hat model.

    Parameters
    ----------
    amplitude : float
        Amplitude
    x_0 : float
        x position of the peak
    y_0 : float
        y position of the peak
    sigma : float
        Width of the Mexican hat

    See Also
    --------
    MexicanHat1D, Gaussian2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}
        + \\left(y - y_{0}\\right)^{2}}{\\sigma^{2}}\\right)
        e^{\\frac{- \\left(x - x_{0}\\right)^{2}
        - \\left(y - y_{0}\\right)^{2}}{2 \\sigma^{2}}}
    """

    amplitude = Parameter()
    x_0 = Parameter()
    y_0 = Parameter()
    sigma = Parameter()

    def __init__(self, amplitude, x_0, y_0, sigma, **kwargs):
        super(MexicanHat2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, sigma=sigma, **kwargs)

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, sigma):
        """Two dimensional Mexican Hat model function"""

        rr_ww = ((x - x_0) ** 2 + (y - y_0) ** 2) / (2 * sigma ** 2)
        return amplitude * (1 - rr_ww) * np.exp(- rr_ww)


class AiryDisk2D(Fittable2DModel):
    """
    Two dimensional Airy disk model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the Airy function.
    x_0 : float
        x position of the maximum of the Airy function.
    y_0 : float
        y position of the maximum of the Airy function.
    radius : float
        The radius of the Airy disk (radius of the first zero).

    See Also
    --------
    Box2D, TrapezoidDisk2D, Gaussian2D

    Notes
    -----
    Model formula:

        .. math:: f(r) = A \\left[\\frac{2 J_1(\\frac{\\pi r}{R/R_z})}{\\frac{\\pi r}{R/R_z}}\\right]^2

    Where :math:`J_1` is the first order Bessel function of the first
    kind, :math:`r` is radial distance from the maximum of the Airy
    function (:math:`r = \\sqrt{(x - x_0)^2 + (y - y_0)^2}`), :math:`R`
    is the input ``radius`` parameter, and :math:`R_z =
    1.2196698912665045`).

    For an optical system, the radius of the first zero represents the
    limiting angular resolution and is approximately 1.22 * lambda / D,
    where lambda is the wavelength of the light and D is the diameter of
    the aperture.

    See [1]_ for more details about the Airy disk.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Airy_disk
    """

    amplitude = Parameter()
    x_0 = Parameter()
    y_0 = Parameter()
    radius = Parameter()
    _j1 = None

    def __init__(self, amplitude, x_0, y_0, radius, **kwargs):
        if self._j1 is None:
            try:
                from scipy.special import j1, jn_zeros
                self.__class__._j1 = j1
                self.__class__._rz = jn_zeros(1, 1)[0] / np.pi
            # add a ValueError here for python3 + scipy < 0.12
            except ValueError:
                raise ImportError("AiryDisk2D model requires scipy > 0.11.")

        super(AiryDisk2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, radius=radius, **kwargs)

    def __deepcopy__(self, memo):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    def __copy__(self):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    @classmethod
    def eval(cls, x, y, amplitude, x_0, y_0, radius):
        """Two dimensional Airy model function"""

        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2) / (radius / cls._rz)
        # Since r can be zero, we have to take care to treat that case
        # separately so as not to raise a numpy warning
        z = np.ones(r.shape)
        rt = np.pi * r[r > 0]
        z[r > 0] = (2.0 * cls._j1(rt) / rt)**2
        z *= amplitude
        return z


class Beta1D(Fittable1DModel):
    """
    One dimensional Beta model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Beta model.
    gamma : float
        Core width of the Beta model.
    alpha : float
        Power index of the beta model.

    See Also
    --------
    Gaussian1D, Box1D

    Notes
    -----
    Model formula:

    .. math::

        f(x) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter()
    x_0 = Parameter()
    gamma = Parameter()
    alpha = Parameter()

    def __init__(self, amplitude, x_0, gamma, alpha, **kwargs):
        super(Beta1D, self).__init__(
            amplitude=amplitude, x_0=x_0, gamma=gamma, alpha=alpha, **kwargs)

    @staticmethod
    def eval(x, amplitude, x_0, gamma, alpha):
        """One dimensional Beta model function"""

        return amplitude * (1 + ((x - x_0) / gamma) ** 2) ** (-alpha)

    @staticmethod
    def fit_deriv(x, amplitude, x_0, gamma, alpha):
        """One dimensional Beta model derivative with respect to parameters"""

        d_A = (1 + (x - x_0) ** 2 / gamma ** 2) ** (-alpha)
        d_x_0 = (-amplitude * alpha * d_A * (-2 * x + 2 * x_0) /
                 (gamma ** 2 * d_A ** alpha))
        d_gamma = (2 * amplitude * alpha * d_A * (x - x_0) ** 2 /
                   (gamma ** 3 * d_A ** alpha))
        d_alpha = -amplitude * d_A * np.log(1 + (x - x_0) ** 2 / gamma ** 2)
        return [d_A, d_x_0, d_gamma, d_alpha]


class Beta2D(Fittable2DModel):
    """
    Two dimensional Beta model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Beta model.
    y_0 : float
        y position of the maximum of the Beta model.
    gamma : float
        Core width of the Beta model.
    alpha : float
        Power index of the beta model.

    See Also
    --------
    Gaussian2D, Box2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2} +
        \\left(y - y_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter()
    x_0 = Parameter()
    y_0 = Parameter()
    gamma = Parameter()
    alpha = Parameter()

    def __init__(self, amplitude, x_0, y_0, gamma, alpha, **kwargs):
        super(Beta2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, gamma=gamma, alpha=alpha,
            **kwargs)

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Beta model function"""

        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        return amplitude * (1 + rr_gg) ** (-alpha)

    @staticmethod
    def fit_deriv(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Beta model derivative with respect to parameters"""

        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        d_A = (1 + rr_gg) ** (-alpha)
        d_x_0 = (-amplitude * alpha * d_A * (-2 * x + 2 * x_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_y_0 = (-amplitude * alpha * d_A * (-2 * y + 2 * y_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_alpha = -amplitude * d_A * np.log(1 + rr_gg)
        d_gamma = 2 * amplitude * alpha * d_A * (rr_gg / (gamma * (1 + rr_gg)))
        return [d_A, d_x_0, d_y_0, d_gamma, d_alpha]


def custom_model_1d(func, func_fit_deriv=None):
    """
    Create a one dimensional model from a user defined function. The
    parameters of the model will be inferred from the arguments of
    the function.

    .. note::

        All model parameters have to be defined as keyword arguments
        with default values in the model function.

    If you want to use parameter sets in the model, the parameters should be
    treated as lists or arrays.

    Parameters
    ----------
    func : function
        Function which defines the model.  It should take one positional
        argument (the independent variable in the model), and any number of
        keyword arguments (the parameters).  It must return the value of the
        model (typically as an array, but can also be a scalar for scalar
        inputs).  This corresponds to the
        `~astropy.modeling.Fittable1DModel.eval` method.
    func_fit_deriv : function, optional
        Function which defines the Jacobian derivative of the model. I.e., the
        derivive with respect to the *parameters* of the model.  It should
        have the same argument signature as ``func``, but should return a
        sequence where each element of the sequence is the derivative
        with respect to the correseponding argument. This corresponds to the
        :meth:`~astropy.modeling.FittableModel.fit_deriv` method.


    Examples
    --------
    Define a sinusoidal model function as a custom 1D model:

        >>> from astropy.modeling.models import custom_model_1d
        >>> import numpy as np
        >>> def sine_model(x, amplitude=1., frequency=1.):
        ...     return amplitude * np.sin(2 * np.pi * frequency * x)
        >>> def sine_deriv(x, amplitude=1., frequency=1.):
        ...     return 2 * np.pi * amplitude * np.cos(2 * np.pi * frequency * x)
        >>> SineModel = custom_model_1d(sine_model, func_fit_deriv=sine_deriv)

    Create an instance of the custom model and evaluate it:

        >>> model = SineModel()
        >>> model(0.25)
        1.0

    This model instance can now be used like a usual astropy model.
    """

    if not six.callable(func):
        raise ModelDefinitionError("Not callable. Must be function")

    if func_fit_deriv is not None and not six.callable(func_fit_deriv):
        raise ModelDefinitionError(
                "func_fit_deriv not callable. Must be function")

    model_name = func.__name__
    param_values = six.get_function_defaults(func)

    # Check if all parameters are keyword arguments
    nparams = len(param_values)

    if (func_fit_deriv is not None and
            len(six.get_function_defaults(func_fit_deriv)) != nparams):
        raise ModelDefinitionError("derivative function should accept "
                                   "same number of parameters as func.")

    func_code = six.get_function_code(func)
    if func_code.co_argcount == nparams + 1:
        param_names = func_code.co_varnames[1:nparams + 1]
    else:
        raise ModelDefinitionError(
            "All parameters must be keyword arguments")

    params = dict((name, Parameter(name, default=default))
                  for name, default in zip(param_names, param_values))

    arg_signature_1 = ', '.join('{0}=None'.format(name)
                                for name in param_names)
    arg_signature_2 = ', '.join('{0}={0}'.format(name)
                                for name in param_names)

    mod = find_current_module(2)
    if mod:
        filename = mod.__file__
        modname = mod.__name__
    else:
        filename = '<string>'
        modname = '__main__'

    members = {'eval': staticmethod(func)}

    eval_globals = {}

    init_code_string = dedent("""
        def __init__(self, {0}, **kwargs):
            super(self.__class__, self).__init__({1}, **kwargs)
    """).format(arg_signature_1, arg_signature_2)

    eval(compile(init_code_string, filename, 'single'), eval_globals)

    if func_fit_deriv is not None:
        members['fit_deriv'] = staticmethod(func_fit_deriv)

    members['__init__'] = eval_globals['__init__']
    members.update(params)

    cls = type(model_name, (Fittable1DModel,), members)
    cls.__module__ = modname

    return cls
