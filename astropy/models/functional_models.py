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
 
__all__ = ['Gauss1DModel', 'Gauss2DModel',  'ScaleModel', 'ShiftModel']
        
class Gauss1DModel(ParametricModel):
    """
    
    Implements 1D Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian
    xcen : float
        Center of the gaussian
    fwhm : float
        FWHM
    xsigma : float
        igma of the gaussian
        Either fwhm or xsigma must be specified
    fjac : callable, 'estimated' or None
        if callable - a function to compute the Jacobian of 
        func with derivatives across the rows.
        if 'estimated' - the Jacobian will be estimated
        if None - if the model has a deriv method, it will be used,
        if not the Jacobian will be estimated.
        
    """
    parnames = ['amplitude', 'xcen', 'xsigma']
    def __init__(self, amplitude, xcen, fwhm=None, xsigma=None,
                            fjac=None, **cons):
        self._amplitude = parameters.Parameter('amplitude', amplitude, self, 1)
        if xsigma is None and fwhm is None:
            raise InputParameterError(
                "Either fwhm or xsigma must be specified")
        if xsigma is not None:
            xsigmaval = xsigma
        else:
            try:
                xsigmaval  = 0.42466 * fwhm
            except TypeError:
                xsigmaval = [0.42466 * n for n in fwhm]
        self._xsigma = parameters.Parameter('xsigma', xsigmaval, self, 1)
        self._xcen = parameters.Parameter('xcen', xcen, self, 1)
        try:
            paramdim = len(self._amplitude)
            assert (len(amplitude) == len(xsigmaval) == len(xcen) ), \
             "Input parameters do not have the same dimension"
        except TypeError:
            paramdim = 1
        super(Gauss1DModel, self).__init__(self.parnames, ndim=1, outdim=1,
                                                                    paramdim=paramdim, **cons)
        self.linear = False
        if fjac is 'estimated':
            self.deriv = None
        elif callable(fjac):
            self.deriv = fjac
        else:
            self.deriv = self.gderiv
            
    def eval(self, x, params):
        return params[0] * np.exp((-(1/(params[2]**2)) * 
                                                (x-params[1])**2))
 
    def gderiv(self, p, x, y):
        amplitude, xcen, xsigma = p
        deriv_dict = {}
        deriv_dict['amplitude'] = np.exp((-(1/(xsigma**2)) * (x-xcen)**2))
        deriv_dict['xcen'] = 2 * amplitude * np.exp((-(1/(xsigma**2)) *
                                (x-xcen)**2)) * (x-xcen)/(xsigma**2)
        deriv_dict['xsigma'] = 2 * amplitude * np.exp((-(1/(xsigma**2)) *
                                (x-xcen)**2)) * ((x-xcen)**2)/(xsigma**3)
        derivval = [deriv_dict[par] for par in self.parnames]
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
        x, fmt = _convert_input(x, self.paramdim)
        result = self.eval(x, self.psets)
        return _convert_output(result, fmt)
    
class Gauss2DModel(ParametricModel):
    """
    
    2D Gaussian.
    
    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian
    xcen : float
        Center of the gaussian in x
    ycen : float
        Center of the gaussian in y
    fwhm : float
        FWHM
    xsigma : float
        sigma of the gaussian in x
        Either fwhm or xsigma must be specified
    ysigma : float
        sigma of the gaussian in y
        Either ysigma or ratio should be given
    ratio : float
        ysigma/xsigma 
    fjac : callable or None
        if callable - a function to compute the Jacobian of 
        func with derivatives across the rows.
        if None - the Jacobian will be estimated
    theta : float 
        rotation angle in radians
    """
    parnames = ['amplitude', 'xcen', 'ycen', 'xsigma', 'ysigma', 'theta']
    
    def __init__(self, amplitude, xcen, ycen, fwhm=None, xsigma=None,
                 ysigma=None, ratio=None, theta=0.0, fjac=None, **cons):
        if ysigma is None and ratio is None:
            raise InputParameterError(
                "Either ysigma or ratio must be specified")
        elif xsigma is None and fwhm is None:
            raise InputParameterError(
                "Either fwhm or xsigma must be specified")
        self._amplitude = parameters.Parameter('amplitude', amplitude, self, 1)
        if xsigma is None:
            xsigma = 0.42466 * fwhm
        self._xsigma = parameters.Parameter('xsigma', xsigma, self, 1)
        if ysigma is None:
            ysigma = ratio * self._xsigma
        self._ysigma = parameters.Parameter('ysigma', ysigma, self, 1)
        self._xcen = parameters.Parameter('xcen', xcen, self, 1)
        self._ycen = parameters.Parameter('ycen', ycen, self, 1)
        self._theta = parameters.Parameter('theta', theta, self, 1)
        try:
            paramdim = len(self._amplitude)
            assert (len(self._amplitude) == len(self._xsigma) == \
                            len(self._xcen) == len(self._ycen) == \
                            len(self._theta) ), \
                            "Input parameters do not have the same dimension"
        except TypeError:
            paramdim = 1
        super(Gauss2DModel, self).__init__(self.parnames, ndim=2, outdim=1,
                                                                    paramdim=paramdim)
        self.linear = False
        if fjac:
            self.deriv = fjac
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
        x, _ = _convert_input(x, self.paramdim)
        y, fmt = _convert_input(y, self.paramdim)
        result = self.eval(x, y, self.psets)
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
    parnames = ['offsets']
    def __init__(self, offsets):
        if not isinstance(offsets, collections.Sequence):
            paramdim = 1
        else:
            paramdim = len(offsets)
        self._offsets = parameters.Parameter('offsets', offsets, self, paramdim)
        super(ShiftModel, self).__init__(self.parnames, ndim=1, outdim=1,
                                                            paramdim=paramdim)

    def __call__(self, x):
        """
        Transforms data using this model.
        """
        x, fmt = _convert_input(x, self.paramdim)
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
    parnames = ['factors']
    def __init__(self, factors):
        if not isinstance(factors, collections.Sequence):
            paramdim = 1
        else:
            paramdim = len(factors)
        self._factors = parameters.Parameter('factors', factors, self, paramdim)
        super(ScaleModel, self).__init__(self.parnames, ndim=1, outdim=1,
                                                            paramdim=paramdim)
    
    def __call__(self, x):
        """
        Transforms data using this model.
        """
        x, fmt = _convert_input(x, self.paramdim)
        result = x * self.factors
        return _convert_output(result, fmt)
