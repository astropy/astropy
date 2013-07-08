# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Mathematical models
"""
from __future__ import division
import collections
import numpy as np
from . import parameters
from .core import ParametricModel, Parametric1DModel, Model, _convert_input, _convert_output
from .utils import InputParameterError, ModelDefinitionError

__all__ = ['Gaussian1DModel', 'Gaussian2DModel', 'ScaleModel', 'ShiftModel',
           'Custom1DModel', 'Sine1DModel', 'Linear1DModel', 'PowerLaw1DModel',
           'Const1DModel', 'Lorentz1DModel', 'Box1DModel']


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

    def __init__(self, amplitude, mean, stddev, **cons):
        super(Gaussian1DModel, self).__init__(locals(), **cons)

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
        d_mean = (amplitude * np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
                            * (x - mean) / stddev ** 2)
        d_stddev = (amplitude * np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
                              * (x - mean) ** 2 / stddev ** 3)
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
                 jacobian_func=None, **cons):

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
                                              param_dim=param_dim, **cons)
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
    Model formula:
        f(x) = scale * x ** (-alpha)
    """
    param_names = ['scale', 'alpha']

    def __init__(self, scale, alpha):
        super(PowerLaw1DModel, self).__init__(locals())

    def eval(self, x, scale, alpha):
        """
        Model function PowerLaw1D
        """
        return scale * x ** (-alpha)

    def deriv(self, x, scale, alpha):
        """
        Model derivative PowerLaw1D
        """
        d_scale = x ** (-alpha)
        d_alpha = scale * ((x) ** (-alpha)) * np.log(x)
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

    def __init__(self, amplitude, frequency):
        super(Sine1DModel, self).__init__(locals())

    def eval(self, x, amplitude, frequency):
        """
        Model function Sine1D
        """
        return amplitude * np.sin(2 * np.pi * frequency * x)

    def deriv(self, x, amplitude, frequency):
        """
        Model function Sine1D
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

    def __init__(self, slope, intercept):
        super(Linear1DModel, self).__init__(locals())
        self.linear = True

    def eval(self, x, slope, intercept):
        """
        Model function Linear1D
        """
        return slope * x + intercept

    def deriv(self, x, slope, intercept):
        """
        Model function derivatives Linear1D
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

    Notes
    -----
    Model formula:
        f(x) = amplitude * ((fwhm / 2.) ** 2) / ((x - x_0) ** 2 + (fwhm / 2.) ** 2)
    """
    param_names = ['amplitude', 'x_0', 'fwhm']

    def __init__(self, amplitude, x_0, fwhm):
        super(Lorentz1DModel, self).__init__(locals())

    def eval(self, x, amplitude, x_0, fwhm):
        """
        Model function Lorentz1D
        """
        return amplitude * ((fwhm / 2.) ** 2) / ((x - x_0) ** 2 + (fwhm / 2.) ** 2)


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

    def __init__(self, amplitude):
        super(Const1DModel, self).__init__(locals())

    def eval(self, x, amplitude):
        """
        Model function Const1D
        """
        return amplitude

    def deriv(self, x, amplitude):
        """
        Model function derivatives Const1D
        """
        d_amplitude = np.ones_like(x)
        return [d_amplitude]


class Const2DModel(ParametricModel):

    """
    Two dimensional constant function.
    """
    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Disk2DModel(ParametricModel):

    """
    Two dimensional radial symmetric box function.
    """
    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Delta1DModel(Parametric1DModel):

    """
    One dimensional Dirac delta function
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

    Note that at f(x_0 - width / 2.) = f(x_0 + width / 2.) = amplitude.
    """
    param_names = ['amplitude', 'x_0', 'width']

    def __init__(self, amplitude, x_0, width):
        super(Box1DModel, self).__init__(locals())

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


class Box2DModel(ParametricModel):

    """
    Two dimensional box function.
    """
    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class MexicanHat1DModel(ParametricModel):

    """
    One dimensional mexican hat function.
    """
    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class MexicanHat2DModel(ParametricModel):

    """
    Two dimensional mexican hat function.
    """
    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Custom1DModel(Parametric1DModel):

    """
    Create one dimensional model from a user defined function.

    IMPORTANT: All model parameters have to be defined as KEYWORD ARGUMENTS
    with default values in the model function.

    If you want to work with parameter sets, the parameters have to be defined
    as lists.

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

    def __init__(self, func, func_deriv=None):
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
        super(Custom1DModel, self).__init__(param_dict)

    def eval(self, x, *params):
        """
        Model function Custom1D
        """
        return self._func(x, *params)
