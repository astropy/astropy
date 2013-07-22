# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Mathematical models.
"""

from __future__ import division
import collections
import numpy as np
from . import parameters
from .core import (ParametricModel, Parametric1DModel, Parametric2DModel,
                   Model, _convert_input, _convert_output)
from .utils import InputParameterError, ModelDefinitionError

__all__ = ['Gaussian1DModel', 'Gaussian2DModel', 'ScaleModel', 'ShiftModel',
           'Custom1DModel', 'Sine1DModel', 'Linear1DModel', 'PowerLaw1DModel',
           'Const1DModel', 'Const2DModel', 'Lorentz1DModel', 'Box1DModel',
           'Box2DModel', 'Disk2DModel', 'Trapezoid1DModel', 'TrapezoidDisk2DModel',
           'MexicanHat1DModel', 'MexicanHat2DModel', 'Airy2DModel']


class Gaussian1DModel(Parametric1DModel):

    """
    Implements 1D Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian
    mean : float
        Mean of the gaussian
    stddev : float
        Standard deviation of the gaussian
    """

    param_names = ['amplitude', 'mean', 'stddev']

    def __init__(self, amplitude, mean, stddev, **constraints):
        super(Gaussian1DModel, self).__init__(locals(), **constraints)

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


class Gaussian2DModel(ParametricModel):

    """

    2D Gaussian.

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
    x_fwhm : float
        Full width at half maximum in x
        Either x_fwhm or x_stddev must be specified
    y_fwhm : float
        Full width at half maximum in y
        Either y_fwhm or y_stddev must be specified
    jacobian_func : callable or None
        if callable - a function to compute the Jacobian of
        func with derivatives across the rows.
        if None - the Jacobian will be estimated
    cov_matrix : ndarray
        A 2x2 covariance matrix. If specified, overrides stddev, fwhm, and
        theta specification.
    """

    param_names = ['amplitude', 'x_mean', 'y_mean',
                   'x_stddev', 'y_stddev', 'theta']

    def __init__(self, amplitude, x_mean, y_mean, x_stddev=None, y_stddev=None,
                 theta=0.0, x_fwhm=None, y_fwhm=None, cov_matrix=None,
                 jacobian_func=None, **constraints):

        try:
            param_dim = len(amplitude)
        except TypeError:
            param_dim = 1

        if y_stddev is None and y_fwhm is None and cov_matrix is None:
            raise InputParameterError(
                "Either y_fwhm or y_stddev must be specified, or a "
                "covariance matrix.")
        elif x_stddev is None and x_fwhm is None and cov_matrix is None:
            raise InputParameterError(
                "Either x_fwhm or x_stddev must be specified, or a "
                "covariance matrix.")
        elif x_stddev is not None and x_fwhm is not None:
            raise InputParameterError("Cannot specify both x_fwhm and x_stddev")
        elif y_stddev is not None and y_fwhm is not None:
            raise InputParameterError("Cannot specify both y_fwhm and y_stddev")
        elif cov_matrix is not None and (x_stddev is not None or
                                         y_stddev is not None or
                                         x_fwhm is not None or
                                         y_fwhm is not None):
            raise InputParameterError("Cannot specify both cov_matrix and x/y_stddev or x/y_fwhm")

        self._amplitude = parameters.Parameter('amplitude', amplitude,
                                               self, param_dim)
        if cov_matrix is None:
            if x_stddev is None:
                x_stddev = 0.42466 * x_fwhm
            if y_stddev is None:
                y_stddev = 0.42466 * y_fwhm
        else:
            cov_matrix = np.array(cov_matrix)
            assert cov_matrix.shape == (2, 2), "Covariance matrix must be 2D"
            eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
            x_stddev, y_stddev = np.sqrt(eig_vals)
            y_vec = eig_vecs[:, 0]
            theta = np.arctan2(y_vec[1], y_vec[0])

        self._x_stddev = parameters.Parameter('x_stddev', x_stddev,
                                              self, param_dim)
        self._y_stddev = parameters.Parameter('y_stddev', y_stddev,
                                              self, param_dim)
        self._x_mean = parameters.Parameter('x_mean', x_mean, self, param_dim)
        self._y_mean = parameters.Parameter('y_mean', y_mean, self, param_dim)
        self._theta = parameters.Parameter('theta', theta, self, param_dim)

        super(Gaussian2DModel, self).__init__(self.param_names, n_inputs=2, n_outputs=1,
                                              param_dim=param_dim, **constraints)
        self.linear = False
        if jacobian_func:
            self.deriv = jacobian_func
        else:
            self.deriv = None

    def eval(self, x, y, params):
        """
        Evaluate the model.

        Parameters
        ----------
        x : array like or a number
            input
        y : dummy variable - array like or a number
            input
        params : array
            parameter sets
        """
        a = 0.5 * ((np.cos(params[5]) / params[3]) ** 2 +
                   (np.sin(params[5]) / params[4]) ** 2)
        b = 0.5 * (np.cos(params[5]) * np.sin(params[5]) *
                   (1. / params[3] ** 2 - 1. / params[4] ** 2))
        c = 0.5 * ((np.sin(params[5]) / params[3]) ** 2 +
                   (np.cos(params[5]) / params[4]) ** 2)

        return params[0] * np.exp(-(a * (x - params[1]) ** 2 +
                                    b * (x - params[1]) * (y - params[2]) +
                                    c * (y - params[2]) ** 2))

    def __call__(self, x, y):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        y : array like or a number
            input

        """
        x, _ = _convert_input(x, self.param_dim)
        y, fmt = _convert_input(y, self.param_dim)
        result = self.eval(x, y, self.param_sets)
        return _convert_output(result, fmt)


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

    def __init__(self, offsets, **constraints):
        if not isinstance(offsets, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(offsets)
        self._offsets = parameters.Parameter('offsets', offsets, self, param_dim)
        super(ShiftModel, self).__init__(self.param_names, n_inputs=1, n_outputs=1,
                                         param_dim=param_dim, **constraints)

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

    def __init__(self, factors, **constraints):
        if not isinstance(factors, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(factors)
        self._factors = parameters.Parameter('factors', factors, self, param_dim)
        super(ScaleModel, self).__init__(self.param_names, n_inputs=1, n_outputs=1,
                                         param_dim=param_dim, **constraints)

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
    A power law model.

    The model is of the form :math:`A x^\\alpha`, where :math:`A` is
    the `scale` parameter, and :math:`\\alpha` is `alpha`.

    Parameters
    ----------
    scale : float
        Model scale
    alpha : float
        power

    Notes
    -----
    Model formula Python:
        f(x) = scale * x ** (-alpha)

    Model formula:

    .. math:: f(x) = A x^{-\\alpha}

    """
    param_names = ['scale', 'alpha']

    def __init__(self, scale, alpha, **constraints):
        super(PowerLaw1DModel, self).__init__(locals(), **constraints)

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
    One dimensional sine model.

    Parameters
    ----------
    amplitude : float
        Oscillation amplitude
    frequency : float
        Oscillation frequency

    Notes
    -----
    Model formula:
        f(x) = amplitude * np.sin(2 * np.pi * frequency * x)
    """
    param_names = ['amplitude', 'frequency']

    def __init__(self, amplitude, frequency, **constraints):
        super(Sine1DModel, self).__init__(locals(), **constraints)

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
    Simple one dimensional straight line model.

    Parameters
    ----------
    slope : float
        Slope of the straight line

    intercept : float
        Intercept of the straight line

    Notes
    -----
    Model formula:
        f(x) = slope * x + intercept
    """
    param_names = ['slope', 'intercept']

    def __init__(self, slope, intercept, **constraints):
        super(Linear1DModel, self).__init__(locals(), **constraints)
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
    One dimensional Lorentzian function.

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
    Model formula Python:
        f(x) = amplitude * ((fwhm / 2.) ** 2) / ((x - x_0) ** 2 + (fwhm / 2.) ** 2)

    Model formula:

    .. math::

        f(x) = \\frac{A \\gamma^{2}}{\\gamma^{2} + \\left(x - x_{0}\\right)^{2}}
    """
    param_names = ['amplitude', 'x_0', 'fwhm']

    def __init__(self, amplitude, x_0, fwhm, **constraints):
        super(Lorentz1DModel, self).__init__(locals(), **constraints)

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
    One dimensional constant function.

    Parameters
    ----------
    amplitude : float
        Value of the constant function

    Notes
    -----
    Model formula:
        f(x) = amplitude
    """
    param_names = ['amplitude']

    def __init__(self, amplitude, **constraints):
        super(Const1DModel, self).__init__(locals(), **constraints)

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
    Two dimensional constant function.
    """
    param_names = ['amplitude']

    def __init__(self, amplitude, **constraints):
        super(Const2DModel, self).__init__(locals(), **constraints)

    def eval(self, x, y, amplitude):
        """
        Model function Const2D
        """
        return amplitude * np.ones_like(x)


class Disk2DModel(Parametric2DModel):

    """
    Two dimensional radial symmetric box function.

    Parameters
    ----------
    amplitude : float
        Value of the disk function
    x_0 : float
        x position center of the disk
    y_0 : float
        y position center of the disk
    radius : float
        Radius of the disk

    Notes
    -----
    Model formula:
        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        f(x) = np.select([rr >= radius], [amplitude])
    """
    param_names = ['amplitude', 'x_0', 'y_0', 'radius']

    def __init__(self, amplitude, x_0, y_0, radius, **constraints):
        super(Disk2DModel, self).__init__(locals(), **constraints)

    def eval(self, x, y, amplitude, x_0, y_0, radius):
        """
        Model function Disk2D.
        """
        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        return np.select([rr <= radius ** 2], [amplitude])


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
    One dimensional box function.

    Parameters
    ----------
    amplitude : float
        Amplitude A
    x_0 : float
        Position of the center of the box function
    width : float
        Width of the box

    Notes
    -----
    Model function:
        f(x) = np.select([x >= x_0 - width / 2., x <= x_0 + width / 2.],
                         [amplitude, amplitude])

    Note that f(x_0 - width / 2.) = f(x_0 + width / 2.) = amplitude.
    """
    param_names = ['amplitude', 'x_0', 'width']

    def __init__(self, amplitude, x_0, width, **constraints):
        super(Box1DModel, self).__init__(locals(), **constraints)

    def eval(self, x, amplitude, x_0, width):
        """
        Model function Box1D
        """
        return np.select([np.logical_and(x >= x_0 - width / 2., x <= x_0 + width / 2.)],
                         [amplitude])

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
    Two dimensional box function.

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
    Box1DModel

    """

    param_names = ['amplitude', 'x_0', 'y_0', 'x_width', 'y_width']

    def __init__(self, amplitude, x_0, y_0, x_width, y_width, **constraints):
        super(Box2DModel, self).__init__(locals(), **constraints)

    def eval(self, x, y, amplitude, x_0, y_0, x_width, y_width):
        """
        Model function Box2DModel.
        """
        x_range = np.logical_and(x >= x_0 - x_width / 2., x <= x_0 + x_width / 2.)
        y_range = np.logical_and(y >= y_0 - y_width / 2., y <= y_0 + y_width / 2.)
        return np.select([np.logical_and(x_range, y_range)], [amplitude])


class Trapezoid1DModel(Parametric1DModel):

    """
    One dimensional trapezoid model.

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
    Box1DModel, Gaussian1DModel
    """
    param_names = ['amplitude', 'x_0', 'width', 'slope']

    def __init__(self, amplitude, x_0, width, slope, **constraints):
        super(Trapezoid1DModel, self).__init__(locals(), **constraints)

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
    Two dimensional circular trapezoid model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the trapezoid
    x_0 : float
        x position of the center of the trapezoid
    y_0 : float
        y position of the center of the trapezoid
    width : float
        Width in x direction of the constant part of the trapezoid.
    slope : float
        Slope of the tails of the trapezoid in x direction.

    See Also
    --------
    Disk2DModel, Box2DModel
    """

    param_names = ['amplitude', 'x_0', 'y_0', 'radius', 'slope']

    def __init__(self, amplitude, x_0, y_0, radius, slope, **constraints):
        super(TrapezoidDisk2DModel, self).__init__(locals(), **constraints)

    def eval(self, x, y, amplitude, x_0, y_0, radius, slope):
        """
        Model function Trapezoid2D.
        """
        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2)
        range_1 = r <= radius
        range_2 = np.logical_and(r >= radius,  r <= radius + amplitude / slope)
        val_1 = amplitude
        val_2 = amplitude + slope * (radius - r)
        return np.select([range_1, range_2], [val_1, val_2])


class MexicanHat1DModel(Parametric1DModel):

    """
    One dimensional mexican hat function.

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
    Box1DModel, Gaussian1DModel, Trapezoid1DModel

    Notes
    -----
    Model formula:

    .. math::

        f(x) = {A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}\\right)
        e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}}

    """
    param_names = ['amplitude', 'x_0', 'width']

    def __init__(self, amplitude, x_0, width, **constraints):
        super(MexicanHat1DModel, self).__init__(locals(), **constraints)

    def eval(self, x, amplitude, x_0, width):
        """
        Model function MexicanHat1DModel.
        """
        xx_ww = (x - x_0) ** 2 / (2 * width ** 2)
        return amplitude * (1 - xx_ww) * np.exp(-xx_ww)


class MexicanHat2DModel(Parametric2DModel):

    """
    Two dimensional symmetrical mexican hat function.

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

    Notes
    -----
    Model formula:
        f(x, y) = amplitude * (1 - ((x - x_0) ** 2 + (y - y_0) ** 2) / (2 * width ** 2))
                * np.exp(- 0.5 * ((x - x_0) ** 2 + (y - y_0) ** 2) / width ** 2))

    """
    param_names = ['amplitude', 'x_0', 'y_0', 'width']

    def __init__(self, amplitude, x_0, y_0, width, **constraints):
        super(MexicanHat2DModel, self).__init__(locals(), **constraints)

    def eval(self, x, y, amplitude, x_0, y_0, width):
        """
        Model function MexicanHat2DModel.
        """
        rr_ww = ((x - x_0) ** 2 + (y - y_0) ** 2) / (2 * width ** 2)
        return amplitude * (1 - rr_ww) * np.exp(- rr_ww)


class Airy2DModel(Parametric2DModel):
    """
    Two dimensional symmetric Airy model.

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

    Notes
    -----
    Model formula:
        f(r) = amplitude * j1(2 * pi * r) / (pi * r)

        Where j1 is the first order Bessel function of first kind.
    """
    param_names = ['amplitude', 'x_0', 'y_0', 'width']

    def __init__(self, amplitude, x_0, y_0, width, **constraints):
        try:
            from scipy.special import j1
            self._j1 = j1
        except ImportError:
            raise ImportError("Could not import scipy.special.")
        super(Airy2DModel, self).__init__(locals(), **constraints)

    def eval(self, x, y, amplitude, x_0, y_0, width):
        """
        Model function Airy2D.
        """
        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2)
        return np.select([r == 0], [1], amplitude * self._j1(2 * np.pi * r) / (np.pi * r))


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
