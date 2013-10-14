# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Mathematical models."""

from __future__ import division

import collections

from textwrap import dedent

import numpy as np

from .core import (ParametricModel, Parametric1DModel, Parametric2DModel,
                   Model, format_input, ModelDefinitionError)
from .parameters import Parameter, InputParameterError
from ..utils import find_current_module


__all__ = sorted([
    'AiryDisk2DModel', 'Beta1DModel', 'Beta2DModel', 'Box1DModel',
    'Box2DModel', 'Const1DModel', 'Const2DModel', 'Disk2DModel',
    'Gaussian1DModel', 'Gaussian2DModel', 'Linear1DModel', 'Lorentz1DModel',
    'MexicanHat1DModel', 'MexicanHat2DModel', 'ScaleModel', 'ShiftModel',
    'Sine1DModel', 'Trapezoid1DModel', 'TrapezoidDisk2DModel', 'Ring2DModel',
    'custom_model_1d'
])


class Gaussian1DModel(Parametric1DModel):
    """
    One dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian
    mean : float
        Mean of the gaussian
    stddev : float
        Standard deviation of the gaussian

    Notes
    -----

    Model formula:

        .. math:: f(x) = A e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}

    See Also
    --------
    Gaussian2DModel, Box1DModel, Beta1DModel, Lorentz1DModel
    """

    amplitude = Parameter('amplitude')
    mean = Parameter('mean')
    stddev = Parameter('stddev')

    def __init__(self, amplitude, mean, stddev, **constraints):
        super(Gaussian1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude, mean, stddev):
        """
        Model function Gauss1D
        """
        return amplitude * np.exp(- 0.5 * (x - mean) ** 2 / stddev ** 2)

    @staticmethod
    def deriv(x, amplitude, mean, stddev):
        """
        Model function derivatives Gauss1D
        """

        d_amplitude = np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
        d_mean = amplitude * d_amplitude * (x - mean) / stddev ** 2
        d_stddev = amplitude * d_amplitude * (x - mean) ** 2 / stddev ** 3
        return [d_amplitude, d_mean, d_stddev]


class Gaussian2DModel(Parametric2DModel):
    """
    Two dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian
    x_mean : float
        Mean of the gaussian in x
    y_mean : float
        Mean of the gaussian in y
    x_stddev : float
        Standard deviation of the gaussian in x
        Either x_fwhm or x_stddev must be specified
    y_stddev : float
        Standard deviation of the gaussian in y
        Either y_fwhm or y_stddev must be specified
    theta : float
        Rotation angle in radians. Note: increases clockwise.
    cov_matrix : ndarray
        A 2x2 covariance matrix. If specified, overrides stddev, fwhm, and
        theta specification.

    Notes
    -----
    Model formula:

        .. math::

            f(x, y) = A e^{-a\\left(x - x_{0}\\right)^{2}  -b\\left(x - x_{0}\\right)
            \\left(y - y_{0}\\right)  -c\\left(y - y_{0}\\right)^{2}}

    Using the following definitions:

        .. math::
            a = \\left(- \\frac{\\sin^{2}{\\left (\\theta \\right )}}{2 \\sigma_{y}^{2}} -
            \\frac{\\cos^{2}{\\left (\\theta \\right )}}{2 \\sigma_{x}^{2}}\\right)

            b = \\left(\\frac{\\sin{\\left (2 \\theta \\right )}}{2 \\sigma_{y}^{2}} -
            \\frac{\\sin{\\left (2 \\theta \\right )}}{2 \\sigma_{x}^{2}}\\right)

            c = \\left(\\frac{\\cos^{2}{\\left (\\theta \\right )}}{2 \\sigma_{y}^{2}} +
            \\frac{\\sin^{2}{\\left (\\theta \\right )}}{2 \\sigma_{x}^{2}}\\right)


    See Also
    --------
    Gaussian1DModel, Box2DModel, Beta2DModel
    """

    amplitude = Parameter('amplitude')
    x_mean = Parameter('x_mean')
    y_mean = Parameter('y_mean')
    x_stddev = Parameter('x_stddev')
    y_stddev = Parameter('y_stddev')
    theta = Parameter('theta')

    def __init__(self, amplitude, x_mean, y_mean, x_stddev=None, y_stddev=None,
                 theta=0.0, cov_matrix=None, **constraints):
        if y_stddev is None and cov_matrix is None:
            raise InputParameterError(
                "Either y_stddev must be specified, or a "
                "covariance matrix.")
        elif x_stddev is None and cov_matrix is None:
            raise InputParameterError(
                "Either x_stddev must be specified, or a "
                "covariance matrix.")
        elif cov_matrix is not None and (x_stddev is not None or
                                         y_stddev is not None):
            raise InputParameterError(
                "Cannot specify both cov_matrix and x/y_stddev")

        # Compute principle coordinate system transformation
        elif cov_matrix is not None:
            cov_matrix = np.array(cov_matrix)
            assert cov_matrix.shape == (2, 2), "Covariance matrix must be 2x2"
            eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
            x_stddev, y_stddev = np.sqrt(eig_vals)
            y_vec = eig_vecs[:, 0]
            theta = np.arctan2(y_vec[1], y_vec[0])

        super(Gaussian2DModel, self).__init__(locals())

    @staticmethod
    def eval(x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """Two dimensional Gaussian function"""

        a = 0.5 * ((np.cos(theta) / x_stddev) ** 2 +
                         (np.sin(theta) / y_stddev) ** 2)
        b = 0.5 * (np.cos(theta) * np.sin(theta) *
                         (1. / x_stddev ** 2 - 1. / y_stddev ** 2))
        c = 0.5 * ((np.sin(theta) / x_stddev) ** 2 +
                         (np.cos(theta) / y_stddev) ** 2)

        return amplitude * np.exp(-(a * (x - x_mean) ** 2 +
                                    b * (x - x_mean) * (y - y_mean) +
                                    c * (y - y_mean) ** 2))

    @staticmethod
    def deriv(x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """Two dimensional Gaussian function derivative"""

        # Helper quantities
        # Derivatives are not checked yet
        a = 0.5 * ((np.cos(theta) / x_stddev) ** 2 +
                         (np.sin(theta) / y_stddev) ** 2)
        b = 0.5 * (np.cos(theta) * np.sin(theta) *
                         (1. / x_stddev ** 2 - 1. / y_stddev ** 2))
        c = 0.5 * ((np.sin(theta) / x_stddev) ** 2 +
                         (np.cos(theta) / y_stddev) ** 2)

        da_dtheta = np.sin(theta) * np.cos(theta) * (1/y_stddev**2 - 1/x_stddev**2)
        da_dx_stddev = -np.cos(theta)**2/x_stddev**3
        da_dy_stddev = -np.sin(theta)**2/y_stddev**3
        db_dtheta = 0.5*np.cos(2*theta) * (1/x_stddev**2 - 1/y_stddev**2)
        db_dx_stddev = -np.cos(theta)*np.sin(theta)/x_stddev**3
        db_dy_stddev = np.cos(theta)*np.sin(theta)/y_stddev**3
        dc_dtheta = np.cos(theta)*np.sin(theta)*(1/x_stddev**2 - 1/y_stddev**2)
        dc_dx_stddev = -np.sin(theta)**2/x_stddev**3
        dc_dy_stddev = -np.cos(theta)**2/y_stddev**3

        d_A = np.exp(- a * (x - x_mean) ** 2
                     - b * (x - x_mean) * (y - y_mean)
                     - c * (y - y_mean) ** 2)
        d_theta = amplitude * (-(x - x_mean) ** 2 * da_dtheta - (x - x_mean) *
                               (y - y_mean) * db_dtheta - (y - y_mean) ** 2 * dc_dtheta) * d_A

        d_x_stddev = amplitude * (-(x - x_mean) ** 2 *da_dx_stddev -
                                  (x - x_mean) * (y - y_mean) * db_dx_stddev
                                  - (y - y_mean) ** 2 * dc_dx_stddev) * d_A
        d_y_stddev = amplitude * (-(x - x_mean) ** 2 *da_dy_stddev -
                                  (x - x_mean) * (y - y_mean) * db_dy_stddev
                                  - (y - y_mean) ** 2 * dc_dy_stddev) * d_A
        d_y_mean = amplitude * (+(x - x_mean) * b + 2 * (y - y_mean) * c) * d_A
        d_x_mean = amplitude * ((y - y_mean) * b + 2 * (x - x_mean) * a) * d_A
        return [d_A, d_x_mean, d_y_mean, d_x_stddev, d_y_stddev, d_theta]


class ShiftModel(Model):
    """
    Shift a coordinate.

    Parameters
    ----------
    offsets : float or a list of floats
        offsets to be applied to a coordinate
        if a list - each value in the list is an offset to be applied to a
        column in the input coordinate array
    """

    offsets = Parameter('offsets')

    def __init__(self, offsets, param_dim=1):
        if not isinstance(offsets, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(offsets)

        self._offsets = offsets

        super(ShiftModel, self).__init__(n_inputs=1, n_outputs=1,
                                         param_dim=param_dim)

    def inverse(self):
        if self.param_dim == 1:
            return ShiftModel(offsets=(-1) * self._offsets)
        else:
            return ShiftModel(offsets=[off * (-1) for off in self._offsets])

    @format_input
    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        return self._offsets + x


class ScaleModel(Model):
    """
    Multiply a model by a factor.

    Parameters
    ----------
    factors : float or a list of floats
        scale for a coordinate
    """

    factors = Parameter('factors')

    def __init__(self, factors, param_dim=1):
        if not isinstance(factors, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(factors)

        self._factors = factors

        super(ScaleModel, self).__init__(n_inputs=1, n_outputs=1,
                                         param_dim=param_dim)

    def inverse(self):
        if self.param_dim == 1:
            return ScaleModel(factors=1. / self._factors)
        else:
            return ScaleModel(factors=[1 / factor for factor in self._factors])

    @format_input
    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        return self._factors * x


class Sine1DModel(Parametric1DModel):
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
    Const1DModel, Linear1DModel


    Notes
    -----
    Model formula:

        .. math:: f(x) = A \\sin(2 \\pi f x)
    """

    amplitude = Parameter('amplitude')
    frequency = Parameter('frequency')

    def __init__(self, amplitude, frequency, **constraints):
        super(Sine1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude, frequency):
        """One dimensional Sine model function"""

        return amplitude * np.sin(2 * np.pi * frequency * x)

    @staticmethod
    def deriv(x, amplitude, frequency):
        """One dimensional Sine model derivative"""

        d_amplitude = np.sin(2 * np.pi * frequency * x)
        d_frequency = (2 * np.pi * x * amplitude *
                       np.cos(2 * np.pi * frequency * x))
        return [d_amplitude, d_frequency]


class Linear1DModel(Parametric1DModel):
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
    Const1DModel

    Notes
    -----
    Model formula:

        .. math:: f(x) = a x + b
    """

    slope = Parameter('slope')
    intercept = Parameter('intercept')

    def __init__(self, slope, intercept, **constraints):
        super(Linear1DModel, self).__init__(locals())
        self.linear = True

    @staticmethod
    def eval(x, slope, intercept):
        """One dimensional Line model function"""

        return slope * x + intercept

    @staticmethod
    def deriv(x, slope, intercept):
        """One dimensional Line model derivative"""

        d_slope = x
        d_intercept = np.ones_like(x)
        return [d_slope, d_intercept]


class Lorentz1DModel(Parametric1DModel):
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
    Gaussian1DModel, Box1DModel, MexicanHat1DModel

    Notes
    -----
    Model formula:

    .. math::

        f(x) = \\frac{A \\gamma^{2}}{\\gamma^{2} + \\left(x - x_{0}\\right)^{2}}
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    fwhm = Parameter('fwhm')

    def __init__(self, amplitude, x_0, fwhm, **constraints):
        super(Lorentz1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude, x_0, fwhm):
        """One dimensional Lorentzian model function"""

        return (amplitude * ((fwhm / 2.) ** 2) / ((x - x_0) ** 2 +
                (fwhm / 2.) ** 2))

    @staticmethod
    def deriv(x, amplitude, x_0, fwhm):
        """One dimensional Lorentzian model derivative"""

        d_amplitude = fwhm ** 2 / (fwhm ** 2 + (x - x_0) ** 2)
        d_x_0 = (amplitude * d_amplitude * (2 * x - 2 * x_0) /
                 (fwhm ** 2 + (x - x_0) ** 2))
        d_fwhm = 2 * amplitude * d_amplitude / fwhm * (1 - d_amplitude)
        return [d_amplitude, d_x_0, d_fwhm]


class Const1DModel(Parametric1DModel):
    """
    One dimensional Constant model.

    Parameters
    ----------
    amplitude : float
        Value of the constant function

    See Also
    --------
    Const2DModel

    Notes
    -----
    Model formula:

        .. math:: f(x) = A
    """

    amplitude = Parameter('amplitude')

    def __init__(self, amplitude, **constraints):
        super(Const1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude):
        """One dimensional Constant model function"""

        return amplitude * np.ones_like(x)

    @staticmethod
    def deriv(x, amplitude):
        """One dimensional Constant model derivative"""

        d_amplitude = np.ones_like(x)
        return [d_amplitude]


class Const2DModel(Parametric2DModel):
    """
    Two dimensional Constant model.

    Parameters
    ----------
    amplitude : float
        Value of the constant function

    See Also
    --------
    Const1DModel

    Notes
    -----
    Model formula:

        .. math:: f(x, y) = A
    """

    amplitude = Parameter('amplitude')

    def __init__(self, amplitude, **constraints):
        super(Const2DModel, self).__init__(locals())

    @staticmethod
    def eval(x, y, amplitude):
        """Two dimensional Constant model function"""

        return amplitude * np.ones_like(x)


class Disk2DModel(Parametric2DModel):
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
    Box2DModel, TrapezoidDisk2DModel

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

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    y_0 = Parameter('y_0')
    R_0 = Parameter('R_0')

    def __init__(self, amplitude, x_0, y_0, R_0, **constraints):
        super(Disk2DModel, self).__init__(locals())

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, R_0):
        """Two dimensional Disk model function"""

        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        return np.select([rr <= R_0 ** 2], [amplitude])


class Ring2DModel(Parametric2DModel):

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
    Disk2DModel, TrapezoidDisk2DModel

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

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    y_0 = Parameter('y_0')
    r_in = Parameter('r_in')
    width = Parameter('width')

    def __init__(self, amplitude, x_0, y_0, r_in, width=None, r_out=None,
                 **constraints):
        if r_out is not None:
            width = r_out - r_in
        if r_out is None and width is None:
            raise ModelDefinitionError("Either specify width or r_out.")

        super(Ring2DModel, self).__init__(locals())

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, r_in, width):
        """
        Model function Ring2D.
        """

        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        r_range = np.logical_and(rr >= r_in ** 2, rr <= (r_in + width) ** 2)
        return np.select([r_range], [amplitude])


class Delta1DModel(Parametric1DModel):
    """One dimensional Dirac delta function."""

    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Delta2DModel(Parametric2DModel):
    """Two dimensional Dirac delta function."""

    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Box1DModel(Parametric1DModel):
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
    Box2DModel, TrapezoidDisk2DModel

    Notes
    -----
    Model formula:

      .. math::

            f(x) = \\left \\{
                     \\begin{array}{ll}
                       A & : x_0 - w/2 \\geq x \\geq x_0 + w/2 \\\\
                       A/2 & :  x = x_0 + w/2 \\\\
                       A/2 & :  x = x_0 - w/2 \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    width = Parameter('width')

    def __init__(self, amplitude, x_0, width, **constraints):
        super(Box1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude, x_0, width):
        """One dimensional Box model function"""

        return np.select([np.logical_and(x > x_0 - width / 2.,
                                         x < x_0 + width / 2.),
                          np.logical_or(x == x_0 - width / 2.,
                                        x == x_0 + width / 2.)],
                         [amplitude, amplitude / 2.])

    @classmethod
    def deriv(cls, x, amplitude, x_0, width):
        """One dimensional Box model derivative"""

        d_amplitude = cls.eval(x, 1, x_0, width)
        d_x_0 = np.zeros_like(x)
        d_width = np.zeros_like(x)
        return [d_amplitude, d_x_0, d_width]


class Box2DModel(Parametric2DModel):
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
    Box1DModel, Gaussian2DModel, Beta2DModel

    Notes
    -----
    Model formula:

      .. math::

            f(x, y) = \\left \\{
                     \\begin{array}{ll}
                       A & : x_0 - w_x/2 \\geq x \\geq x_0 + w_x/2 \\\\
                       A & : y_0 - w_y/2 \\geq y \\geq y_0 + w_y/2 \\\\
                       A/2 & :  x = x_0 + w_x/2 \\\\
                       A/2 & :  x = x_0 - w_x/2 \\\\
                       A/2 & :  y = y_0 + w_y/2 \\\\
                       A/2 & :  y = y_0 - w_y/2 \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.

    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    y_0 = Parameter('y_0')
    x_width = Parameter('x_width')
    y_width = Parameter('y_width')

    def __init__(self, amplitude, x_0, y_0, x_width, y_width, **constraints):
        super(Box2DModel, self).__init__(locals())

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, x_width, y_width):
        """Two dimensional Box model function"""
        x_range = np.logical_and(x > x_0 - x_width / 2., x < x_0 + x_width / 2.)
        y_range = np.logical_and(y > y_0 - y_width / 2., y < y_0 + y_width / 2.)
        x_boundary = np.logical_or(x == x_0 - x_width / 2., x == x_0 + x_width / 2.)
        y_boundary = np.logical_or(y == y_0 - y_width / 2., y == y_0 + y_width / 2.)
        return np.select([np.logical_and(x_range, y_range),
                          np.logical_and(x_boundary, y_range),
                          np.logical_and(y_boundary, x_range),
                          np.logical_and(y_boundary, x_boundary)],
                         [amplitude, amplitude / 2., amplitude / 2., amplitude / 4.])


class Trapezoid1DModel(Parametric1DModel):
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
    Box1DModel, Gaussian1DModel, Beta1DModel
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    width = Parameter('width')
    slope = Parameter('slope')

    def __init__(self, amplitude, x_0, width, slope, **constraints):
        super(Trapezoid1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude, x_0, width, slope):
        """One dimensional Trapezoid model function"""

        range_1 = np.logical_and(x >= x_0 - width / 2. - amplitude / slope,
                                 x < x_0 - width / 2.)
        range_2 = np.logical_and(x >= x_0 - width / 2., x < x_0 + width / 2.)
        range_3 = np.logical_and(x >= x_0 + width / 2.,
                                 x < x_0 + width / 2. + amplitude / slope)
        val_1 = amplitude + slope * (x_0 + width / 2. + x)
        val_2 = amplitude
        val_3 = amplitude + slope * (x_0 + width / 2. - x)
        return np.select([range_1, range_2, range_3], [val_1, val_2, val_3])


class TrapezoidDisk2DModel(Parametric2DModel):
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
    Disk2DModel, Box2DModel
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    y_0 = Parameter('y_0')
    R_0 = Parameter('R_0')
    slope = Parameter('slope')

    def __init__(self, amplitude, x_0, y_0, R_0, slope, **constraints):
        super(TrapezoidDisk2DModel, self).__init__(locals())


    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, R_0, slope):
        """Two dimensional Trapezoid Disk model function"""

        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2)
        range_1 = r <= R_0
        range_2 = np.logical_and(r > R_0,  r <= R_0 + amplitude / slope)
        val_1 = amplitude
        val_2 = amplitude + slope * (R_0 - r)
        return np.select([range_1, range_2], [val_1, val_2])


class MexicanHat1DModel(Parametric1DModel):
    """
    One dimensional Mexican Hat model.

    Parameters
    ----------
    amplitude : float
        Amplitude
    x_0 : float
        Position of the peak
    width : float
        Width of the mexican hat

    See Also
    --------
    MexicanHat2DModel, Box1DModel, Gaussian1DModel, Trapezoid1DModel

    Notes
    -----
    Model formula:

    .. math::

        f(x) = {A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}}{\\sigma^{2}}\\right)
        e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}}

    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    sigma = Parameter('sigma')

    def __init__(self, amplitude, x_0, sigma, **constraints):
        super(MexicanHat1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude, x_0, sigma):
        """One dimensional Mexican Hat model function"""

        xx_ww = (x - x_0) ** 2 / (2 * sigma ** 2)
        return amplitude * (1 - 2 * xx_ww) * np.exp(-xx_ww)


class MexicanHat2DModel(Parametric2DModel):
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
    width : float
        Width of the mexican hat

    See Also
    --------
    MexicanHat1DModel, Gaussian2DModel

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}
        + \\left(y - y_{0}\\right)^{2}}{\\sigma^{2}}\\right)
        e^{\\frac{- \\left(x - x_{0}\\right)^{2}
        - \\left(y - y_{0}\\right)^{2}}{2 \\sigma^{2}}}
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    y_0 = Parameter('y_0')
    sigma = Parameter('sigma')

    def __init__(self, amplitude, x_0, y_0, sigma, **constraints):
        super(MexicanHat2DModel, self).__init__(locals())

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, sigma):
        """Two dimensional Mexican Hat model function"""

        rr_ww = ((x - x_0) ** 2 + (y - y_0) ** 2) / (2 * sigma ** 2)
        return amplitude * (1 - rr_ww) * np.exp(- rr_ww)


class AiryDisk2DModel(Parametric2DModel):
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
    width : float
        Width of the Airy function.

    See Also
    --------
    Box2DModel, TrapezoidDisk2DModel, Gaussian2DModel


    Notes
    -----
    Model formula:

        .. math:: f(r) = A \\frac{J_1(2 \\pi r)}{\\pi r}

    Where J1 is the first order Bessel function of first kind.
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    y_0 = Parameter('y_0')
    width = Parameter('width')

    _j1 = None

    def __init__(self, amplitude, x_0, y_0, width, **constraints):
        if self._j1 is None:
            try:
                from scipy.special import j1
                self.__class__._j1 = j1
            except ImportError:
                raise ImportError("AiryDisk2DModel requires scipy.")
        super(AiryDisk2DModel, self).__init__(locals())

    @classmethod
    def eval(cls, x, y, amplitude, x_0, y_0, width):
        """Two dimensional Airy model function"""

        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2) / width

        # Since r can be zero, we have to take care to treat that case
        # separately so as not to raise a Numpy warning
        z = np.ones(r.shape)
        z[r > 0] = (amplitude * (cls._j1(2 * np.pi * r[r > 0]) /
                    (np.pi * r[r > 0])) ** 2)

        return z


class Beta1DModel(Parametric1DModel):
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
    Gaussian1DModel, Box1DModel

    Notes
    -----
    Model formula:

    .. math::

        f(x) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    gamma = Parameter('gamma')
    alpha = Parameter('alpha')

    def __init__(self, amplitude, x_0, gamma, alpha, **constraints):
        super(Beta1DModel, self).__init__(locals())

    @staticmethod
    def eval(x, amplitude, x_0, gamma, alpha):
        """One dimensional Beta model function"""

        return amplitude * (1 + ((x - x_0) / gamma) ** 2) ** (-alpha)

    @staticmethod
    def deriv(x, amplitude, x_0, gamma, alpha):
        """One dimensional Beta model derivative"""

        d_A = (1 + (x - x_0) ** 2 / gamma ** 2) ** (-alpha)
        d_x_0 = (-amplitude * alpha * d_A * (-2 * x + 2 * x_0) /
                 (gamma ** 2 * d_A ** alpha))
        d_gamma = (2 * amplitude * alpha * d_A * (x - x_0) ** 2 /
                   (gamma ** 3 * d_A ** alpha))
        d_alpha = -amplitude * d_A * np.log(1 + (x - x_0) ** 2 / gamma ** 2)
        return [d_A, d_x_0, d_gamma, d_alpha]


class Beta2DModel(Parametric2DModel):
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
    Gaussian2DModel, Box2DModel

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2} +
        \\left(y - y_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter('amplitude')
    x_0 = Parameter('x_0')
    y_0 = Parameter('y_0')
    gamma = Parameter('gamma')
    alpha = Parameter('alpha')

    def __init__(self, amplitude, x_0, y_0, gamma, alpha, **constraints):
        super(Beta2DModel, self).__init__(locals())

    @staticmethod
    def eval(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Beta model function"""

        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        return amplitude * (1 + rr_gg) ** (-alpha)

    @staticmethod
    def deriv(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Beta model derivative"""

        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        d_A = (1 + rr_gg) ** (-alpha)
        d_x_0 = (-amplitude * alpha * d_A * (-2 * x + 2 * x_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_y_0 = (-amplitude * alpha * d_A * (-2 * y + 2 * y_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_alpha = -amplitude * d_A * np.log(1 + rr_gg)
        d_gamma = 2 * amplitude * alpha * d_A * (rr_gg / (gamma * (1 + rr_gg)))
        return [d_A, d_x_0, d_y_0, d_gamma, d_alpha]


def custom_model_1d(func):
    """
    Create one dimensional model from a user defined function.

    IMPORTANT: All model parameters have to be defined as KEYWORD ARGUMENTS
    with default values in the model function.

    If you want to work with parameter sets, the parameters have to be defined
    as lists or arrays.

    Parameters
    ----------
    func : function
        Function which defines the model
    func_deriv : function
        Function which defines the model derivatives default = None

    Examples
    --------
    Define a sinusoidal model function as a custom 1D model:

        >>> from astropy.modeling.models import custom_model_1d
        >>> import numpy as np
        >>> @custom_model_1d
        ... def SineModel(x, amplitude=1., frequency=1.):
        ...     return amplitude * np.sin(2 * np.pi * frequency * x)

    Create an instance of the custom model and evaluate it:

        >>> model = SineModel()
        >>> model(0.25)
        1.0

    This model instance can now be used like a usual astropy model.
    """

    if not callable(func):
        raise ModelDefinitionError("Not callable. Must be function")

    model_name = func.__name__
    param_values = func.func_defaults

    # Check if all parameters are keyword arguments
    nparams = len(param_values)
    if func.func_code.co_argcount == nparams + 1:
        param_names = func.func_code.co_varnames[1:nparams + 1]
    else:
        raise ModelDefinitionError(
            "All parameters must be keyword arguments")

    params = dict((name, Parameter(name)) for name in param_names)

    arg_signature = ', '.join('{0}={1!r}'.format(name, value)
                              for name, value in zip(param_names,
                                                     param_values))

    mod = find_current_module(2)
    if mod:
        filename = mod.__file__
        modname = mod.__name__
    else:
        filename = '<string>'
        modname = '__main__'

    members = {'eval': staticmethod(func)}

    eval(compile(dedent("""
        def __init__(self, {0}):
            super(self.__class__, self).__init__(locals())
    """).format(arg_signature), filename, 'single'), members)

    members.update(params)

    cls = type(model_name, (Parametric1DModel,), members)
    cls.__module__ = modname

    return cls
