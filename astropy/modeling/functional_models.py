# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Mathematical models
"""
from __future__ import division, print_function
import collections
import numpy as np
from . import parameters
from .core import ParametricModel, Model, _convert_input, _convert_output
from .utils import InputParameterError

__all__ = ['Gaussian1DModel', 'Gaussian2DModel', 'ScaleModel', 'ShiftModel']


class Gaussian1DModel(ParametricModel):

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
        Either fwhm or stddev must be specified
    fwhm : float
        Full width at half maximum
        Either fwhm or stddev must be specified
    jacobian_func : callable, 'estimated' or None
        if callable - a function to compute the Jacobian of
        func with derivatives across the rows.
        if 'estimated' - the Jacobian will be estimated
        if None - if the model has a deriv method, it will be used,
        if not the Jacobian will be estimated.

    """

    param_names = ['amplitude', 'mean', 'stddev']

    def __init__(self, amplitude, mean, stddev=None, fwhm=None,
                 jacobian_func=None, **cons):

        try:
            param_dim = len(amplitude)
        except TypeError:
            param_dim = 1

        if stddev is None and fwhm is None:
            raise InputParameterError(
                "Either fwhm or stddev must be specified")
        elif stddev is not None and fwhm is not None:
            raise InputParameterError("Cannot specify both fwhm and stddev")
        elif stddev is not None:
            stddev_val = stddev
        else:
            try:
                stddev_val = 0.42466 * fwhm
            except TypeError:
                stddev_val = [0.42466 * n for n in fwhm]

        self._stddev = parameters.Parameter('stddev', stddev_val, self, param_dim)
        self._mean = parameters.Parameter('mean', mean, self, param_dim)
        self._amplitude = parameters.Parameter('amplitude', amplitude, self, param_dim)

        super(Gaussian1DModel, self).__init__(self.param_names, n_inputs=1, n_outputs=1,
                                              param_dim=param_dim, **cons)

        self.linear = False
        if jacobian_func is 'estimated':
            self.deriv = None
        elif callable(jacobian_func):
            self.deriv = jacobian_func
        else:
            self.deriv = self.gderiv

    def eval(self, x, params):
        """
        Evaluate the model.

        Parameters
        ----------
        x : array like or a number
            input
        params : list
            a list of float parameters returned by the optimization algorithm
        """
        return params[0] * np.exp(- 0.5 * (x - params[1]) ** 2 / params[2] ** 2)

    def gderiv(self, params, x, y):
        """
        Gaussian 1D derivative.

        Parameters
        ----------
        params : list
            a list of float parameters returned by the optimization algorithm
        x : array like or a number
            input
        y : dummy variable - array like or a number
            input
        """
        amplitude, mean, stddev = params
        deriv_dict = {}
        deriv_dict['amplitude'] = np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
        deriv_dict['mean'] = (amplitude
                              * np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
                              * (x - mean) / stddev ** 2)
        deriv_dict['stddev'] = (amplitude
                                * np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
                                * (x - mean) ** 2 / stddev ** 3)
        derivval = [deriv_dict[par] for par in self.param_names]
        return np.array(derivval).T

    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        --------------
        x : array, of minimum dimensions 1

        Notes
        -----
        See the module docstring for rules for model evaluation.
        """
        x, fmt = _convert_input(x, self.param_dim)
        result = self.eval(x, self.param_sets)
        return _convert_output(result, fmt)


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
        --------------
        x, y : arrays, of min dimensions 2

        Note: See the module docstring for rules for model evaluation.
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

    def __call__(self, x):
        """
        Transforms data using this model.
        """
        x, fmt = _convert_input(x, self.param_dim)
        result = x + self.offsets
        return _convert_output(result, fmt)


class ScaleModel(Model):

    """

    Multiply a model by a factor.

    Parameters
    ---------------
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

    def __call__(self, x):
        """
        Transforms data using this model.
        """
        x, fmt = _convert_input(x, self.param_dim)
        result = x * self.factors
        return _convert_output(result, fmt)
