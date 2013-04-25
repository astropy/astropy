# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Mathematical models
"""
from __future__ import division, print_function
import collections
import numpy as np
from . import parameters
from .core import *
from .utils import InputParameterError
 
__all__ = ['Gaussian1DModel', 'Gaussian2DModel',  'ScaleModel', 'ShiftModel']
        
class Gaussian1DModel(ParametricModel):
    """
    
    Implements 1D Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian
    mu : float
        Mean of the gaussian
    fwhm : float
        FWHM
    sigma : float
        Standard deviation of the gaussian
        Either fwhm or sigma must be specified
    jacobian_func : callable, 'estimated' or None
        if callable - a function to compute the Jacobian of 
        func with derivatives across the rows.
        if 'estimated' - the Jacobian will be estimated
        if None - if the model has a deriv method, it will be used,
        if not the Jacobian will be estimated.
        
    """
    param_names = ['amplitude', 'mu', 'sigma']
    def __init__(self, amplitude, mu, fwhm=None, sigma=None,
                            jacobian_func=None, **cons):
        self._amplitude = parameters.Parameter('amplitude', amplitude, self, 1)
        if sigma is None and fwhm is None:
            raise InputParameterError(
                "Either fwhm or sigma must be specified")
        if sigma is not None:
            sigmaval = sigma
        else:
            try:
                sigmaval  = 0.42466 * fwhm
            except TypeError:
                sigmaval = [0.42466 * n for n in fwhm]
        self._sigma = parameters.Parameter('sigma', sigmaval, self, 1)
        self._mu = parameters.Parameter('mu', mu, self, 1)
        try:
            param_dim = len(self._amplitude)
            assert (len(amplitude) == len(sigmaval) == len(mu) ), \
             "Input parameters do not have the same dimension"
        except TypeError:
            param_dim = 1
        super(Gaussian1DModel, self).__init__(self.param_names, ndim=1, outdim=1,
                                                                    param_dim=param_dim, **cons)
        self.linear = False
        if jacobian_func is 'estimated':
            self.deriv = None
        elif callable(jacobian_func):
            self.deriv = jacobian_func
        else:
            self.deriv = self.gderiv
            
    def eval(self, x, params):
        return params[0] * np.exp((-(1/(params[2]**2)) * 
                                                (x-params[1])**2))
 
    def gderiv(self, p, x, y):
        amplitude, mu, sigma = p
        deriv_dict = {}
        deriv_dict['amplitude'] = np.exp((-(1/(sigma**2)) * (x-mu)**2))
        deriv_dict['mu'] = 2 * amplitude * np.exp((-(1/(sigma**2)) *
                                (x-mu)**2)) * (x-mu)/(sigma**2)
        deriv_dict['sigma'] = 2 * amplitude * np.exp((-(1/(sigma**2)) *
                                (x-mu)**2)) * ((x-mu)**2)/(sigma**3)
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
    x_mu : float
        Mean of the gaussian in x
    y_mu : float
        Mean of the gaussian in y
    fwhm : float
        Full width at half maximum
    x_sigma : float
        Standard deviation of the gaussian in x
        Either fwhm or x_sigma must be specified
    y_sigma : float
        Standard deviation of the gaussian in y
        Either y_sigma or ratio should be given
    jacobian_func : callable or None
        if callable - a function to compute the Jacobian of 
        func with derivatives across the rows.
        if None - the Jacobian will be estimated
    theta : float 
        rotation angle in radians
    """
    param_names = ['amplitude', 'x_mu', 'y_mu', 'x_sigma', 'y_sigma', 'theta']
    
    def __init__(self, amplitude, x_mu, y_mu, x_fwhm=None, y_fwhm=None, 
                 x_sigma=None, y_sigma=None, theta=0.0, cov_matrix=None,
                 jacobian_func=None, **cons):
        if y_sigma is None and y_fwhm is None and cov_matrix is None:
            raise InputParameterError(
                "Either y_fwhm or y_sigma must be specified, or a "
                "covariance matrix.")
        elif x_sigma is None and x_fwhm is None and cov_matrix is None:
            raise InputParameterError(
                "Either x_fwhm or x_sigma must be specified, or a "
                "covariance matrix.")
                
        self._amplitude = parameters.Parameter('amplitude', amplitude, self, 1)
        
        if cov_matrix is None:
            if x_sigma is None:
                x_sigma = 0.42466 * x_fwhm
            
            if y_sigma is None:
                y_sigma = 0.42466 * y_fwhm
        
        else:
            cov_matrix = np.array(cov_matrix)
            assert cov_matrix.shape == (2,2), "Covariance matrix must be 2D"
            
            eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
            x_sigma, y_sigma = np.sqrt(eig_vals)
            y_vec = eig_vecs[:,0]
            theta = np.arctan2(y_vec[1],y_vec[0])
            
        self._x_sigma = parameters.Parameter('x_sigma', x_sigma, self, 1)
        self._y_sigma = parameters.Parameter('y_sigma', y_sigma, self, 1)
        
        self._x_mu = parameters.Parameter('x_mu', x_mu, self, 1)
        self._y_mu = parameters.Parameter('y_mu', y_mu, self, 1)
        self._theta = parameters.Parameter('theta', theta, self, 1)
        
        try:
            param_dim = len(self._amplitude)
            assert (len(self._amplitude) == len(self._x_sigma) == \
                            len(self._x_mu) == len(self._y_mu) == \
                            len(self._theta) ), \
                            "Input parameters do not have the same dimension"
        except TypeError:
            param_dim = 1
        super(Gaussian2DModel, self).__init__(self.param_names, ndim=2, outdim=1,
                                              param_dim=param_dim)
        self.linear = False
        if jacobian_func:
            self.deriv = jacobian_func
        else:
            self.deriv = None
            
    def eval(self, x, y, p):
        return p[0] * np.exp(-(
            ((np.cos(p[5])/p[1])**2 + (np.sin(p[5])/p[2])**2) * ((x-p[3])**2)
            + 2*(np.cos(p[5])*np.sin(p[5])*(1/(p[1]**2)-1/(p[2]**2))) *
            (x-p[3])*(y-p[4]) + ((np.sin(p[5])/p[1])**2+(np.cos(p[5])/p[2])**2)*
            ((y-p[4])**2)))
        
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
        super(ShiftModel, self).__init__(self.param_names, ndim=1, outdim=1,
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
        super(ScaleModel, self).__init__(self.param_names, ndim=1, outdim=1,
                                                            param_dim=param_dim)
    
    def __call__(self, x):
        """
        Transforms data using this model.
        """
        x, fmt = _convert_input(x, self.param_dim)
        result = x * self.factors
        return _convert_output(result, fmt)
