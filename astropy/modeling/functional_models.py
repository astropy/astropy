# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Mathematical models.
"""

from __future__ import division
import collections
import numpy as np
from . import parameters
from .core import (Parametric1DModel, Parametric2DModel, Model,
                   _convert_input, _convert_output)
from .utils import InputParameterError, ModelDefinitionError

__all__ = sorted(['AiryDisk2DModel', 'Beta1DModel', 'Beta2DModel',
           'Box1DModel', 'Box2DModel', 'Const1DModel', 'Const2DModel',
           'Custom1DModel', 'Disk2DModel', 'Gaussian1DModel', 'Gaussian2DModel',
           'Linear1DModel', 'Lorentz1DModel', 'MexicanHat1DModel',
           'MexicanHat2DModel', 'PowerLaw1DModel', 'ScaleModel', 'ShiftModel',
           'Sine1DModel', 'Trapezoid1DModel', 'TrapezoidDisk2DModel', 'Ring2DModel'])


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

    param_names = ['amplitude', 'mean', 'stddev']

    def __init__(self, amplitude, mean, stddev, **constraints):
        super(Gaussian1DModel, self).__init__(locals())

    def eval(self, x, amplitude, mean, stddev):
        """
        Model function Gauss1D
        """
        return amplitude * np.exp(- 0.5 * (x - mean) ** 2 / stddev ** 2)

    def deriv(self, x, amplitude, mean, stddev):
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

    param_names = ['amplitude', 'x_mean', 'y_mean',
                   'x_stddev', 'y_stddev', 'theta']

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
            raise InputParameterError("Cannot specify both cov_matrix and x/y_stddev")

        # Compute principle coordinate system transformation
        elif cov_matrix is not None:
            cov_matrix = np.array(cov_matrix)
            assert cov_matrix.shape == (2, 2), "Covariance matrix must be 2D"
            eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
            x_stddev, y_stddev = np.sqrt(eig_vals)
            y_vec = eig_vecs[:, 0]
            theta = np.arctan2(y_vec[1], y_vec[0])
        super(Gaussian2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """
        Model function Gaussian2D
        """
        a = 0.5 * ((np.cos(theta) / x_stddev) ** 2 +
                         (np.sin(theta) / y_stddev) ** 2)
        b = 0.5 * (np.cos(theta) * np.sin(theta) *
                         (1. / x_stddev ** 2 - 1. / y_stddev ** 2))
        c = 0.5 * ((np.sin(theta) / x_stddev) ** 2 +
                         (np.cos(theta) / y_stddev) ** 2)

        return amplitude * np.exp(-(a * (x - x_mean) ** 2 +
                                    b * (x - x_mean) * (y - y_mean) +
                                    c * (y - y_mean) ** 2))

    def deriv(self, x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """
        Model function derivatives Gaussian2D.
        """

        # Helper quantities
        # Derivatives are not checked yet
        a = 0.5 * ((np.cos(theta) / x_stddev) ** 2 +
                         (np.sin(theta) / y_stddev) ** 2)
        b = 0.5 * (np.cos(theta) * np.sin(theta) *
                         (1. / x_stddev ** 2 - 1. / y_stddev ** 2))
        c = 0.5 * ((np.sin(theta) / x_stddev) ** 2 +
                         (np.cos(theta) / y_stddev) ** 2)

        d_a = np.sin(theta) * np.cos(theta) * (1 / x_stddev ** 2 - 1 / y_stddev ** 2)
        d_b = np.cos(2 * theta) * (1 / y_stddev ** 2 - 1 / x_stddev ** 2)
        d_c = -d_a

        d_A = np.exp(- a * (x - x_mean) ** 2
                     - b * (x - x_mean) * (y - y_mean)
                     - c * (y - y_mean) ** 2)
        d_theta = amplitude * ((x - x_mean) ** 2 * d_a - (x - x_mean) *
                               (y - y_mean) * d_b + (y - y_mean) ** 2 * d_c) * d_A
        d_y_stddev = amplitude * ((x - x_mean) ** 2 * np.sin(theta) ** 2
                    + (x - x_mean) * (y - y_mean) * np.sin(2 * theta)
                    + (y - y_mean) ** 2 * np.cos(theta) ** 2) * d_A / y_stddev**3
        d_x_stddev = amplitude * ((x - x_mean) ** 2 * np.cos(theta) ** 2
                    - (x - x_mean) * (y - y_mean) * np.sin(2 * theta)
                    + (y - y_mean) ** 2 * np.sin(theta) ** 2) * d_A / x_stddev**3
        d_y_mean = amplitude * ((x - x_mean) * b - 2 * (y - y_mean) * c) * d_A
        d_x_mean = amplitude * ((y - y_mean) * b - 2 * (x - x_mean) * a) * d_A
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
    param_names = ['offsets']

    def __init__(self, offsets):
        if not isinstance(offsets, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(offsets)
        self._offsets = parameters.Parameter('offsets', offsets, self, param_dim)
        super(ShiftModel, self).__init__(self.param_names, n_inputs=1, n_outputs=1,
                                         param_dim=param_dim)

    def inverse(self):
        if self.param_dim == 1:
            return ShiftModel(offsets=(-1) * self.offsets[0])
        else:
            return ShiftModel(offsets=[off * (-1) for off in self._offsets])

    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """
        x, fmt = _convert_input(x, self.param_dim)
        result = x + self.offsets
        return _convert_output(result, fmt)


class ScaleModel(Model):

    """

    Multiply a model by a factor.

    Parameters
    ----------
    factors : float or a list of floats
        scale for a coordinate

    """
    param_names = ['factors']

    def __init__(self, factors):
        if not isinstance(factors, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(factors)
        self._factors = parameters.Parameter('factors', factors, self, param_dim)
        super(ScaleModel, self).__init__(self.param_names, n_inputs=1, n_outputs=1,
                                         param_dim=param_dim)

    def inverse(self):
        if self.param_dim == 1:
            return ScaleModel(factors=1. / self.factors[0])
        else:
            return ScaleModel(factors=[1 / factor for factor in self._factors])

    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """
        x, fmt = _convert_input(x, self.param_dim)
        result = x * self.factors
        return _convert_output(result, fmt)


class PowerLaw1DModel(Parametric1DModel):

    """
    One dimensional Power Law model.

    The model is of the form :math:`A x^\\alpha`, where :math:`A` is
    the `scale` parameter, and :math:`\\alpha` is `alpha`.

    Parameters
    ----------
    scale : float
        Model scale
    alpha : float
        power

    See Also
    --------
    Const1DModel, Linear1DModel


    Notes
    -----
    Model formula:

        .. math:: f(x) = A x^{-\\alpha}

    """
    param_names = ['scale', 'alpha']

    def __init__(self, scale, alpha, **constraints):
        super(PowerLaw1DModel, self).__init__(locals())

    def eval(self, x, scale, alpha):
        """
        Model function PowerLaw1D.
        """
        return scale * x ** (-alpha)

    def deriv(self, x, scale, alpha):
        """
        Model derivative PowerLaw1D.
        """
        d_scale = x ** (-alpha)
        d_alpha = scale * d_scale * np.log(x)
        return [d_scale, d_alpha]


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
    param_names = ['amplitude', 'frequency']

    def __init__(self, amplitude, frequency, **constraints):
        super(Sine1DModel, self).__init__(locals())

    def eval(self, x, amplitude, frequency):
        """
        Model function Sine1D.
        """
        return amplitude * np.sin(2 * np.pi * frequency * x)

    def deriv(self, x, amplitude, frequency):
        """
        Model function Sine1D.
        """
        d_amplitude = np.sin(2 * np.pi * frequency * x)
        d_frequency = (2 * np.pi * x * amplitude
                       * np.cos(2 * np.pi * frequency * x))
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
    param_names = ['slope', 'intercept']

    def __init__(self, slope, intercept, **constraints):
        super(Linear1DModel, self).__init__(locals())
        self.linear = True

    def eval(self, x, slope, intercept):
        """
        Model function Linear1D.
        """
        return slope * x + intercept

    def deriv(self, x, slope, intercept):
        """
        Model function derivatives Linear1D.
        """
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
    param_names = ['amplitude', 'x_0', 'fwhm']

    def __init__(self, amplitude, x_0, fwhm, **constraints):
        super(Lorentz1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, fwhm):
        """
        Model function Lorentz1D.
        """
        return amplitude * ((fwhm / 2.) ** 2) / ((x - x_0) ** 2 + (fwhm / 2.) ** 2)

    def deriv(self, x, amplitude, x_0, fwhm):
        """
        Model function derivatives Lorentz1D.
        """
        d_amplitude = fwhm ** 2 / (fwhm ** 2 + (x - x_0) ** 2)
        d_x_0 = amplitude * d_amplitude * (2 * x - 2 * x_0) / (fwhm ** 2 + (x - x_0) ** 2)
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
    param_names = ['amplitude']

    def __init__(self, amplitude, **constraints):
        super(Const1DModel, self).__init__(locals())

    def eval(self, x, amplitude):
        """
        Model function Const1D
        """
        return amplitude * np.ones_like(x)

    def deriv(self, x, amplitude):
        """
        Model function derivatives Const1D
        """
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
    param_names = ['amplitude']

    def __init__(self, amplitude, **constraints):
        super(Const2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude):
        """
        Model function Const2D
        """
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
    param_names = ['amplitude', 'x_0', 'y_0', 'R_0']

    def __init__(self, amplitude, x_0, y_0, R_0, **constraints):
        super(Disk2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_0, y_0, R_0):
        """
        Model function Disk2D.
        """
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
    R_in : float
        Inner Radius of the ring
    R_out : float
        Outer Radius of the ring
    width : float
        width of the ring. Can be specified instead of R_out.

    See Also
    --------
    Disk2DModel, TrapezoidDisk2DModel

    Notes
    -----
    Model formula:

        .. math::

            f(r) = \\left \\{
                     \\begin{array}{ll}
                       A & : R_{in} \\leq r \\leq R_{out} \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.
    """
    param_names = ['amplitude', 'x_0', 'y_0', 'R_in', 'R_out']

    def __init__(self, amplitude, x_0, y_0, R_in, R_out=None, width=None, **constraints):
        if width != None:
            R_out = R_in + width
        if R_out == None:
            raise ModelDefinitionError("Either specify R_out or width.")
        super(Ring2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_0, y_0, R_in, R_out):
        """
        Model function Ring2D.
        """
        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        r_range = np.logical_and(rr >= R_in ** 2, rr <= R_out ** 2)
        return np.select([r_range], [amplitude])


class Delta1DModel(Parametric1DModel):

    """
    One dimensional Dirac delta function.
    """
    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Delta2DModel(Parametric2DModel):

    """
    Two dimensional Dirac delta function.
    """
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
    param_names = ['amplitude', 'x_0', 'width']

    def __init__(self, amplitude, x_0, width, **constraints):
        super(Box1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, width):
        """
        Model function Box1D
        """
        return np.select([np.logical_and(x > x_0 - width / 2., x < x_0 + width / 2.), 
                          np.logical_or(x == x_0 - width / 2., x == x_0 + width / 2.)],
                         [amplitude, amplitude / 2.])

    def deriv(self, x, amplitude, x_0, width):
        """
        Model function derivatives Box1D
        """
        d_amplitude = self.eval(x, 1, x_0, width)
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

    param_names = ['amplitude', 'x_0', 'y_0', 'x_width', 'y_width']

    def __init__(self, amplitude, x_0, y_0, x_width, y_width, **constraints):
        super(Box2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_0, y_0, x_width, y_width):
        """
        Model function Box2DModel.
        """
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
    param_names = ['amplitude', 'x_0', 'width', 'slope']

    def __init__(self, amplitude, x_0, width, slope, **constraints):
        super(Trapezoid1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, width, slope):
        """
        Model function Trapezoid1D.
        """
        range_1 = np.logical_and(x >= x_0 - width / 2. - amplitude / slope, x < x_0 - width / 2.)
        range_2 = np.logical_and(x >= x_0 - width / 2., x < x_0 + width / 2.)
        range_3 = np.logical_and(x >= x_0 + width / 2., x < x_0 + width / 2. + amplitude / slope)
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

    param_names = ['amplitude', 'x_0', 'y_0', 'R_0', 'slope']

    def __init__(self, amplitude, x_0, y_0, R_0, slope, **constraints):
        super(TrapezoidDisk2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_0, y_0, R_0, slope):
        """
        Model function TrapezoidDisk2D.
        """
        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2)
        range_1 = r <= R_0
        range_2 = np.logical_and(r >= R_0,  r <= R_0 + amplitude / slope)
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
    param_names = ['amplitude', 'x_0', 'sigma']

    def __init__(self, amplitude, x_0, sigma, **constraints):
        super(MexicanHat1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, sigma):
        """
        Model function MexicanHat1DModel.
        """
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
    param_names = ['amplitude', 'x_0', 'y_0', 'sigma']

    def __init__(self, amplitude, x_0, y_0, sigma, **constraints):
        super(MexicanHat2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_0, y_0, sigma):
        """
        Model function MexicanHat2DModel.
        """
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
    param_names = ['amplitude', 'x_0', 'y_0', 'width']

    def __init__(self, amplitude, x_0, y_0, width, **constraints):
        try:
            from scipy.special import j1
            self._j1 = j1
        except ImportError:
            raise ImportError("Could not import scipy.special.")
        super(AiryDisk2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_0, y_0, width):
        """
        Model function Airy2D.
        """
        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2) / width
        return np.select([r == 0], [1], amplitude * self._j1(2 * np.pi * r) / (np.pi * r))


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
    param_names = ['amplitude', 'x_0', 'gamma', 'alpha']

    def __init__(self, amplitude, x_0, gamma, alpha, **constraints):
        super(Beta1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, gamma, alpha):
        """
        Model function Beta1D.
        """
        return amplitude * (1 + ((x - x_0) / gamma) ** 2) ** (-alpha)

    def deriv(self, x, amplitude, x_0, gamma, alpha):
        """
        Model function derivatives Beta1D.
        """
        d_A = (1 + (x - x_0) ** 2 / gamma ** 2) ** (-alpha)
        d_x_0 = - amplitude * alpha * d_A * (-2 * x + 2 * x_0) / (gamma ** 2 * d_A ** alpha)
        d_gamma = 2 * amplitude * alpha * d_A * (x - x_0) ** 2 / (gamma ** 3 * d_A ** alpha)
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
    param_names = ['amplitude', 'x_0', 'y_0', 'gamma', 'alpha']

    def __init__(self, amplitude, x_0, y_0, gamma, alpha, **constraints):
        super(Beta2DModel, self).__init__(locals())

    def eval(self, x, y, amplitude, x_0, y_0, gamma, alpha):
        """
        Model function Beta2D.
        """
        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        return amplitude * (1 + rr_gg) ** (-alpha)

    def deriv(self, x, y, amplitude, x_0, y_0, gamma, alpha):
        """
        Model function derivatives Beta2D.
        """
        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        d_A = (1 + rr_gg) ** (-alpha)
        d_x_0 = - amplitude * alpha * d_A * (-2 * x + 2 * x_0) / (gamma ** 2 * (1 + rr_gg))
        d_y_0 = - amplitude * alpha * d_A * (-2 * y + 2 * y_0) / (gamma ** 2 * (1 + rr_gg))
        d_alpha = - amplitude * d_A * np.log(1 + rr_gg)
        d_gamma = 2 * amplitude * alpha * d_A * (rr_gg / (gamma * (1 + rr_gg)))
        return [d_A, d_x_0, d_y_0, d_gamma, d_alpha]


class Custom1DModel(Parametric1DModel):

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
    Define a sinusoidal model function:

        >>> from astropy.modeling.models import Custom1DModel
        >>> import numpy as np
        >>> def f(x, amplitude=1., frequency=1.):
        ...     return amplitude * np.sin(2 * np.pi * frequency * x)

    And create a custom one dimensional model from it:

        >>> sin_model = Custom1DModel(f)
        >>> sin_model(0.25)
        1.0

    This model instance can now be used like a usual astropy model.
    """

    def __init__(self, func, func_deriv=None, **constraints):
        if callable(func):
            self._func = func
        else:
            raise ModelDefinitionError("Not callable. Must be function")

        param_values = func.func_defaults

        # Check if all parameters are keyword arguments
        if func.func_code.co_argcount == len(param_values) + 1:
            self.param_names = func.func_code.co_varnames[1:1 + len(param_values)]
        else:
            raise ModelDefinitionError("All parameters must be keyword arguments")
        param_dict = dict(zip(self.param_names, param_values))
        super(Custom1DModel, self).__init__(param_dict, **constraints)

    def eval(self, x, *params):
        """
        Model function Custom1D
        """
        return self._func(x, *params)
