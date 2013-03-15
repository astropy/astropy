# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines base classes for all models.
The base class of all models is `Model`.
`ParametricModel` is the base class for all fittable models. Parametric 
models can be linear or nonlinear in a regression analysis sense.

All models provide a `__call__` method which performs the transformation in a
purely mathematical way, i.e. the models are unitless. In addition, when 
possible the transformation is done using multiple parameter sets, `psets`.
The number of parameter sets is stored in an attribute `paramdim`. 

Parametric models also store a flat list of all parameters as an instance of
`parameters.Parameters`. When fitting, this list-like object is modified by a
subclass of `fitting.Fitter`. When fitting nonlinear models, the values of the
parameters are used as initial guesses by the fitting class. Normally users
will not have to use the `parameters` module directly.

Input Format For Model Evaluation and Fitting

Input coordinates are passed in separate arguments, for example 2D models 
expect x and y coordinates to be passed separately as two scalars or aray-like
objects. 
The evaluation depends on the input dimensions and the number of parameter
sets but in general normal broadcasting rules apply.
For example:

-A model with one parameter set works with input in any dimensionality

-A model with N parameter sets works with 2D arrays of shape (M, N)
A parameter set is applied to each column.
 
-A model with N parameter sets works with multidimensional arrays if the
shape of the input array is (N, M, P). A parameter set is applied to each plane.
  
In all these cases the output has the same shape as the input.

-A model with N parameter sets works with 1D input arrays. The shape 
 of the output is (M, N)
 
"""
from __future__ import division, print_function
import operator
import abc
from collections import OrderedDict
import numpy as np
from . import parameters
from . import constraints
from .util import pmapdomain, InputParameterError, comb
 
__all__ = ['ChebyshevModel', 'Gauss1DModel', 'Gauss2DModel', 'ICheb2DModel', 
           'ILegend2DModel', 'LegendreModel', 'Poly1DModel', 'Poly2DModel', 
           'ScaleModel', 'ShiftModel', 'SIPModel', 
           'PCompositeModel', 'SCompositeModel', 'LabeledInput']

def _convert_input(x, pdim):
    """
    Format the input into appropriate shape
    
    Parameters
    ----------
    x : scalar, array or a sequence of numbers
        input data
    pdim : int
        number of parameter sets
        
    The meaning of the internally used format is:
    
    'N' - the format of the input was not changed
    'T' - input was transposed
    'S' - input is a scalar
    """
    x = np.asarray(x) + 0.
    fmt = 'N'
    if pdim == 1:
        if x.ndim == 0:
            fmt = 'S'
            return x, fmt
        else:
            return x, fmt
    else:
        if x.ndim < 2:
            fmt = 'N'
            return np.array([x]).T, fmt
        elif x.ndim == 2:
            assert x.shape[-1] == pdim, "Cannot broadcast with shape"\
                                            "({0}, {1})".format(x.shape[0], x.shape[1])
            return x, fmt
        elif x.ndim > 2:
            assert x.shape[0] == pdim, "Cannot broadcast with shape " \
                   "({0}, {1}, {2})".format(x.shape[0], x.shape[1], x.shape[2])
            fmt = 'T'
            return x.T, fmt
    
def _convert_output(x, fmt):
    """
    Put the output in the shpae/type of the original input
    
    Parameters
    ----------
    x : scalar, array or a sequence of numbers
        output data
    fmt : string
        original format
    """
    if fmt == 'N':
        return x
    elif fmt == 'T':
        return x.T
    elif fmt == 'S':
        return x[0]
    else:
        raise ValueError("Unrecognized output conversion format")

class _ParameterProperty(object):
    """
    Create a property for this parameter.
    """
    def __init__(self, name):
        self.aname = '_'+name
        self.name = name
        
    def __get__(self, obj, objtype):
        par = getattr(obj, self.aname)
        return par
    
    def __set__(self, obj, val):
        if self.name in obj._parcheck:
            obj._parcheck[self.name](val)
        if isinstance(obj, ParametricModel):
            if not obj._parameters._changed:
                par = parameters._Parameter(self.name, val, obj, obj.paramdim)
                oldpar = getattr(obj, self.name)
                if oldpar is not None and oldpar.parshape != par.parshape:
                    raise InputParameterError(
                        "Input parameter {0} does not "
                        "have the required shape".format(self.name))
                else:
                    setattr(obj, self.aname, par)
                obj._parameters = parameters.Parameters(obj, 
                                                        obj.parnames,
                                                         paramdim=obj.paramdim)
            else:
                setattr(obj, self.aname, val)
        else:
            par = parameters._Parameter(self.name, val, obj, obj.paramdim)
            oldpar = getattr(obj, self.name)
            if oldpar is not None and oldpar.parshape != par.parshape:
                raise InputParameterError(
                    "Input parameter {0} does not "
                    "have the required shape".format(self.name))
            else:
                setattr(obj, self.aname, par)
                
class Model(object):
    """
    Base class for all models
    
    This is an abstract class and should not be instanciated.
    
    Notes
    -----
    Models which are not meant to be fit to data should subclass this class
    
    This class sets the properties for all individual parameters and performs
    parameter validation.
    
    """
    __metaclass__ = abc.ABCMeta
    
    parnames = []

    def __init__(self, parnames, paramdim=1):
        self._paramdim = paramdim
        self.has_inverse = False
        self._parnames = parnames
        #_parcheck is a dictionary to register parameter validation funcitons
        #key: value pairs are parameter_name: parameter_validation_function_name
        #see projections.AZP for example
        self._parcheck = {}
        for par in self.parnames:
            setattr(self.__class__, par, _ParameterProperty(par))
                                                
    @property
    def paramdim(self):
        return self._paramdim
    
    @paramdim.setter
    def paramdim(self, val):
        self._paramdim = val
    
    @property
    def parnames(self):
        return self._parnames
    
    @parnames.setter
    def parnames(self, val):
        self._parnames = val
        
    def __repr__(self):
        fmt = "{0}(".format(self.__class__.__name__)
        for i in range(len(self.parnames)):
            fmt1 = """
            {0}= {1},
            """.format(self.parnames[i], getattr(self, self.parnames[i]))
            fmt += fmt1
        fmt += ")"
        
        return fmt
    
    def __str__(self):
        
        fmt = """
        Model: {0}
        Parameter sets: {1}
        Parameters: 
                   {2}
        """.format(
              self.__class__.__name__,
              self.paramdim,
              "\n                   ".join(i+':  ' + 
                str(self.__getattribute__(i)) for i in self.parnames)
                )
            
        return fmt

    @property
    def psets(self):
        """
        Return parameters as a pset
        """
        psets = np.asarray([getattr(self, attr) for attr in self.parnames])
        psets.shape = (len(self.parnames), self.paramdim)
        return psets
    
    def inverse(self):
        """
        Return a callable object which does the inverse transform
        """
        raise NotImplementedError("Subclasses should implement this")
    
    def invert(self):
        """
        Invert coordinates iteratively if possible
        """
        raise NotImplementedError("Subclasses should implement this")

    def add_model(self, newtr, mode):
        """
        Create a CompositeModel by chaining the current model with the new one
        using the specified mode.
        
        Parameters
        ----------
        newtr : an instance of a subclass of Model
        mode :  string
               'parallel', 'serial', 'p' or 's'
               a flag indicating whether to combine the models
               in series or in parallel
               
        Returns
        -------
        model : CompositeModel
            an instance of CompositeModel
        """
        if mode in ['parallel', 'p']:
            return PCompositeModel([self, newtr])
        elif mode in ['serial', 's']:
            return SCompositeModel([self, newtr])
        else:
            raise InputParameterError("Unrecognized mode {0}".format(mode))
    
    @abc.abstractmethod
    def __call__(self):
        raise NotImplementedError("Subclasses should implement this")
    
class ParametricModel(Model):
    """
    Base class for all fittable models
    
    Notes
    -----
    All models which can be fit to data and provide a `deriv` method
    should subclass this class.
    
    Sets the parameters attributes
    
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, parnames, paramdim=1, fittable=True,
                 fixed={}, tied={}, bounds={}, eqcons=[], ineqcons=[]):
        self.linear = True
        super(ParametricModel, self).__init__(parnames, paramdim=paramdim)
        self.fittable = fittable
        self._parameters = parameters.Parameters(self, self.parnames,
                                                 paramdim=paramdim)
        _fixed = {}.fromkeys(self.parnames, False)
        _tied = {}.fromkeys(self.parnames, False)
        _bounds = {}.fromkeys(self.parnames, [-1.E12, 1.E12])
        self.constraints = constraints.Constraints(self, fixed=_fixed,
                                                   tied=_tied, 
                                                   bounds=_bounds,
                                                   eqcons=eqcons, 
                                                   ineqcons=ineqcons)
        if fixed:
            for name in fixed:
                par = getattr(self, name)
                setattr(par, 'fixed', fixed[name])
        if tied:
            for name in tied:
                par = getattr(self, name)
                setattr(par, 'tied', tied[name])
        if bounds:
            for name in bounds:
                par = getattr(self, name)
                setattr(par, 'min', bounds[0])
                setattr(par, 'max', bounds[1])
        
        
    def __repr__(self):
        try:
            degree = str(self.deg)
        except AttributeError:
            degree = ""
        try:
            paramdim = str(self.paramdim)
        except AttributeError:
            paramdim = " "
            
        if degree:
            fmt = "<{0}({1},".format(self.__class__.__name__, repr(self.deg))
        else:
            fmt = "<{0}(".format(self.__class__.__name__)
        for i in range(len(self.parnames)):
            fmt1 = """
            {0}= {1},
            """.format(self.parnames[i], getattr(self, self.parnames[i]))
            fmt += fmt1.strip()
        if paramdim:
            fmt += "paramdim={0})>".format(self.paramdim)
        
        return fmt
    
    def __str__(self):
        try:
            degree = str(self.deg)
        except AttributeError:
            degree = 'N/A'
        fmt = """
        Model: {0}
        Dim:   {1}
        Degree: {2}
        Parameter sets: {3}
        Parameters: 
                   {4}
        """.format(
              self.__class__.__name__,
              self.ndim,
              degree,
              self.paramdim,
              "\n                   ".join(i+':  ' + 
                str(self.__getattribute__(i)) for i in self.parnames)
                )
            
        return fmt
    
    @property
    def parameters(self):
        return self._parameters
    
    @parameters.setter
    def parameters(self, value):
        if isinstance(value, parameters.Parameters):
            if self._parameters._is_same_length(value):
                self._parameters = value
            else:
                raise InputParameterError(
                    "Expected the list of parameters to be the same "
                    "length as the initial list.")
        elif isinstance(value, (list, np.ndarray)):
            _val = parameters._tofloat(value)
            if self._parameters._is_same_length(_val):
                self._parameters._changed = True
                self._parameters[:] = _val
            else:
                raise InputParameterError(
                    "Expected the list of parameters to be the same "
                    "length as the initial list.")
        else:
            raise TypeError("Parameters must be of type 'list' or 'Parameters'")
   
    def set_joint_parameters(self, jpars):
        """
        Used by the JointFitter class to store parameters which are 
        considered common for several models and are to be fitted together.
        """
        self.joint = jpars

class PModel(ParametricModel):
    """
    Base class for all polynomial models.
    
    Its main purpose is to determine how many coefficients are needed
    based on the polynomial order and dimension and to provide their 
    default values, names and ordering.
    
    """
    def __init__(self, degree, ndim=1, paramdim=1, **pars):
        self.deg = degree
        self.ndim = ndim
        self._order = self.get_numcoeff()
        self.parnames = self._generate_coeff_names()
        if not pars:
            self.set_coeff(pardim=paramdim)
        else:
            p = pars.get('c0', pars.get('c0_0'))
            if operator.isSequenceType(p):
                lenpars = len(p)
            else:
                lenpars = 1
            if paramdim != lenpars:
                print("Creating a model with {0} parameter sets\n".format(lenpars))
                paramdim = lenpars
            self._validate_pars(**pars)  
            self.set_coeff(pardim=paramdim, **pars)
        super(PModel, self).__init__(self.parnames, paramdim=paramdim)
    
    def _invlex(self):
        c = []
        lencoeff = self.deg + 1
        for i in range(lencoeff):
            for j in range(lencoeff):
                if i+j <= self.deg:
                    c.append((j, i))
        return c[::-1]
    
    def _generate_coeff_names(self):
        ncoeff = self.get_numcoeff()
        names = []
        if self.ndim == 1:
            for n in range(ncoeff):
                names.append('c{0}'.format(n))
        else:
            for i in range(self.deg+1):
                names.append('c{0}_{1}'.format(i, 0))
            for i in range(1, self.deg+1):
                names.append('c{0}_{1}'.format(0, i))
            for i in range(1, self.deg):
                for j in range(1, self.deg):
                    if i+j < self.deg+1:
                        names.append('c{0}_{1}'.format(i, j))
        return names
        
    def _validate_pars(self, **pars):
        numcoeff = self.get_numcoeff()
        assert(len(pars) == numcoeff)
    
    def set_coeff(self, pardim=1, **pars):
        """
        Set default values for coefficients
        """
        if not pars:
            for name in self.parnames:
                uname = '_'+name
                if pardim == 1:
                    self.__setattr__(uname, parameters._Parameter(
                                        name, 0., self, pardim))
                else:
                    self.__setattr__(uname, parameters._Parameter(
                                        name, [0.]*pardim, self, pardim))
        else:
            for name in self.parnames:
                uname = '_'+name
                self.__setattr__(uname, parameters._Parameter(
                                          name, pars[name], self, pardim))
             
    def get_numcoeff(self):
        """
        Return the number of coefficients in one parameter set
        """
        if self.deg < 1  or self.deg > 16:
            raise ValueError("Degree of polynomial must be 1< deg < 16")
        # deg+1 is used to account for the difference between iraf using 
        # degree and numpy using exact degree
        if self.ndim != 1:
            nmixed = comb(self.deg, self.ndim)
        else: 
            nmixed = 0
        numc = self.deg*self.ndim + nmixed + 1
        return numc

    def set_domain(self, x, y=None):
        """
        Map the input data into a [-1, 1] window
        """
        if self.ndim == 1:
            if not self.domain:
                self.domain = [x.min(), x.max()]
            if not self.window:
                self.window = [-1, 1]
            return pmapdomain(x, self.domain, self.window)
        if self.ndim == 2:
            assert y is not None, ("Expected 2 input coordinates")
            if not self.xdomain:
                self.xdomain = [x.min(), x.max()]
            if not self.xwindow:
                self.xwindow = [-1, 1]
            if not self.ydomain:
                self.ydomain = [y.min(), y.max()]
            if not self.ywindow:
                self.ywindow = [-1, 1]
            xnew = pmapdomain(x, self.xdomain, self.xwindow)
            ynew = pmapdomain(x, self.ydomain, self.ywindow)
            return xnew, ynew
            
class IModel(ParametricModel):
    """
    
    This is a base class for IRAF style 2D Chebyshev and Legendre models.
    
    These are polynomials which have a maximum degree in x and y.

    """
    def __init__(self, xdeg, ydeg, xdomain=None, xwindow=None, ydomain=None, 
                            ywindow=None, paramdim=1, **pars):
        """
        Parameters
        ----------
        
        xdeg : int
            degree in x
        ydeg : int
            degree in y
        xdomain : list or None
            domain of the x independent variable
        ydomain : list or None
            domain of the y independent variable
        xwindow : list or None
            range of the x independent variable
        ywindow : list or None
            range of the y independent variable
        paramdim : int
            number of parameter sets
        **pars : dict
            {keyword: value} pairs, representing {parameter_name: value}
        """
        self.ndim = 2
        self.outdim = 1
        self.xdeg = xdeg
        self.ydeg = ydeg
        self._order = self.get_numcoeff()
        self.xdomain = xdomain
        self.ydomain = ydomain
        self.xwindow = xwindow
        self.ywindow = ywindow
        self.parnames = self._generate_coeff_names()
        
        if not pars:
            self.set_coeff(pardim=paramdim)
        
        else:
            p = pars.get('c0_0')
            if operator.isSequenceType(p):
                lenpars = len(p)
            else:
                lenpars = 1
            if paramdim != lenpars:
                print("Creating a model with {0} parameter sets\n".format(lenpars))
                paramdim = lenpars
            self._validate_pars(**pars)  
            self.set_coeff(pardim=paramdim, **pars)        
        super(IModel, self).__init__(self.parnames, paramdim=paramdim)
    
    def _generate_coeff_names(self):
        names = []
        for j in range(self.ydeg+1):
            for i in range(self.xdeg+1):
                names.append('c{0}_{1}'.format(i, j))
        return names
    
    def set_coeff(self, pardim=1, **pars):
        if not pars:
            for name in self.parnames:
                uname = '_'+name
                self.__setattr__(uname, parameters._Parameter(
                                 name, [0.]*pardim, self, pardim))
        else:
            for name in self.parnames:
                uname = '_'+name
                self.__setattr__(uname, parameters._Parameter(
                    name, pars[name], self, pardim))
    
    def get_numcoeff(self):
        """
        Determine how many coefficients are needed
        
        Returns
        -------
        numc : int
            number of coefficients
        
        """
        numc = (self.xdeg+1)*(self.ydeg+1)
        return numc
    
    
    def _validate_pars(self, **pars):
        numcoeff = self.get_numcoeff()
        assert(len(pars) == numcoeff) 
 
    def _invlex(self):
        c = []
        xvar = np.arange(self.xdeg + 1)
        yvar = np.arange(self.ydeg + 1)
        for j in yvar:
            for i in xvar:
                c.append((i, j))
        return np.array(c[::-1])
    
    def invlex_coeff(self):
        coeff = []
        xvar = np.arange(self.xdeg + 1)
        yvar = np.arange(self.ydeg + 1)
        for j in yvar:
            for i in xvar:
                name = 'c'+str(i)+'_'+str(j)
                coeff.append(getattr(self, name))
        return np.array(coeff[::-1])
    
    def _alpha(self):
        invlexdeg = self._invlex()
        invlexdeg[:, 1] = invlexdeg[:, 1] + self.xdeg+1
        nx = self.xdeg + 1
        ny = self.ydeg + 1
        alpha = np.zeros((ny*nx+3, ny+nx))
        for n in range(len(invlexdeg)):
            alpha[n][invlexdeg[n]] = [1, 1]
            alpha[-2, 0] = 1
            alpha[-3, nx] = 1
        return alpha

    def imhorner(self, x, y, coeff):
        _coeff = list(coeff[:])
        _coeff.extend([0, 0, 0])
        alpha = self._alpha()
        r0 = _coeff[0]
        nalpha = len(alpha)
        
        karr = np.diff(alpha, axis=0)
        kfunc = self._fcache(x, y)
        xterms = self.xdeg+1
        yterms = self.ydeg+1
        nterms = xterms + yterms
        for n in range(1, nterms+1+3):
            setattr(self, 'r'+str(n), 0.)
        
        for n in range(1, nalpha):
            k = karr[n-1].nonzero()[0].max()+1
            rsum = 0
            for i in range(1, k+1):
                rsum = rsum + getattr(self, 'r'+str(i))
            val = kfunc[k-1] * (r0 + rsum)
            setattr(self, 'r'+str(k), val)
            r0 = _coeff[n]
            for i in range(1, k):
                setattr(self, 'r'+str(i), 0.)
        result  = r0
        for i in range(1, nterms+1+3):
            result = result + getattr(self, 'r'+str(i))
        return result
        
    def _fcache(self, x, y):
        """
        To be implemented by subclasses
        """
        raise NotImplementedError("Subclasses should implement this")

class ChebyshevModel(PModel):
    """
    
    1D Chebyshev polynomial
    
    """
    def __init__(self, degree, domain=None, window=[-1, 1], paramdim=1, **pars):
        """
        Parameters
        ----------
        degree : int
            degree of the series
        domain : list or None
        window : list or None
            If None, it is set to [-1,1]
            Fitters will remap the domain to this window
        paramdim : int
            number of parameter sets
        **pars : dict
            keyword : value pairs, representing parameter_name: value
            
        Returns
        -------
        model : ChebyshevModel
            1D Chebyshev model
            
        """
        self.domain = domain
        self.window = window
        super(ChebyshevModel, self).__init__(degree, ndim=1,
                                             paramdim=paramdim, **pars)
        self.outdim = 1
            
    def clenshaw(self, x, coeff):
        if isinstance(x, tuple) or isinstance(x, list) :
            x = np.asarray(x)
        if len(coeff) == 1 :
            c0 = coeff[0]
            c1 = 0
        elif len(coeff) == 2 :
            c0 = coeff[0]
            c1 = coeff[1]
        else :
            x2 = 2*x
            c0 = coeff[-2]
            c1 = coeff[-1]
            for i in range(3, len(coeff) + 1) :
                tmp = c0
                c0 = coeff[-i] - c1
                c1 = tmp + c1*x2
        return c0 + c1*x    
    
    def deriv(self, x):
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.deg + 1,) + x.shape, dtype=x.dtype)
        v[0] = x*0 + 1
        x2 = 2*x
        v[1] = x
        for i in range(2, self.deg + 1) :
            v[i] = v[i-1]*x2 - v[i-2]
        return np.rollaxis(v, 0, v.ndim)
 
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
        if self.domain is not None:
            x = self.set_domain(x)
        x, fmt = _convert_input(x, self.paramdim)
        result = self.clenshaw(x, self.psets)
        return _convert_output(result, fmt)
        
class LegendreModel(PModel):
    """
    
    1D Legendre polynomial

    """
    def __init__(self, degree, domain=None, window=[-1, 1], paramdim=1, **pars):
        """
        Parameters
        ----------
        degree : int
            degree of the series
        domain : list or None
        window : list or None
            If None, it is set to [-1,1]
            Fitters will remap the domain to this window
        paramdim : int
            number of parameter sets
        **pars : dict
            keyword: value pairs, representing parameter_name: value
            
        Returns
        -------
        model : LegendreModel
            1D Legendre model
            
        """
        self.domain = domain
        self.window = window
        super(LegendreModel, self).__init__(degree, ndim=1,
                                            paramdim=paramdim, **pars)
        self.outdim = 1
           
    def clenshaw(self, x, coeff):
        if isinstance(x, tuple) or isinstance(x, list) :
            x = np.asarray(x)
        if len(coeff) == 1 :
            c0 = coeff[0]
            c1 = 0
        elif len(coeff) == 2 :
            c0 = coeff[0]
            c1 = coeff[1]
        else :
            nd = len(coeff)
            c0 = coeff[-2]
            c1 = coeff[-1]
            for i in range(3, len(coeff) + 1) :
                tmp = c0
                nd = nd - 1
                c0 = coeff[-i] - (c1*(nd - 1))/nd
                c1 = tmp + (c1*x*(2*nd - 1))/nd
        return c0 + c1*x    
    
    def deriv(self, x):
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.deg + 1,) + x.shape, dtype=x.dtype)
        v[0] = x*0 + 1
        v[1] = x
        for i in range(2, self.deg + 1) :
            v[i] = (v[i-1]*x*(2*i - 1) - v[i-2]*(i - 1))/i
        return np.rollaxis(v, 0, v.ndim)
 
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
        if self.domain is not None:
            x = self.set_domain(x)
        x, fmt = _convert_input(x, self.paramdim)
        result = self.clenshaw(x, self.psets)
        return _convert_output(result, fmt)
        
class Poly1DModel(PModel):
    """
    
    1D Polynomial model
    
    """ 
    def __init__(self, degree,
                 domain=[-1, 1], window=[-1, 1],
                 paramdim=1, **pars):
        """
        Parameters
        ----------
        degree : int
            degree of the series
        domain : list or None
        window : list or None
            If None, it is set to [-1,1]
            Fitters will remap the domain to this window
        paramdim : int
            number of parameter sets
        **pars : dict
            keyword: value pairs, representing parameter_name: value
            
        Returns
        -------
        model : Poly1DModel
            1D polynomial model
        """
        self.domain = domain
        self.window = window
        super(Poly1DModel, self).__init__(degree, ndim=1,
                                          paramdim=paramdim, **pars)
        self.outdim = 1
            
    def deriv(self, x):
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.deg + 1,) + x.shape, dtype=x.dtype)
        v[0] = x*0 + 1
        v[1] = x
        for i in range(2, self.deg + 1) :
            v[i] = v[i-1]*x
        return np.rollaxis(v, 0, v.ndim)
        
    def horner(self, x, coef):
        c0 = coef[-1] + x*0
        for i in range(2, len(coef)+1):
            c0 = coef[-i] + c0 * x
        return c0

    def __call__(self, x):
        """
        Transforms data using this model.
        
        Parameters
        --------------
        x : array, of minimum dimensions 1
       
        Notes
        -----
        Rules for model evaluation are described in the module docstring
        """
        x, fmt = _convert_input(x, self.paramdim)
        result = self.horner(x, self.psets)
        return _convert_output(result, fmt)
            
class Poly2DModel(PModel):
    """
    2D Polynomial  model
    
    Represents a general polynomial of degree n:
     
    .. math:: P(x,y) = c_{0_0} + c_{1_0}x + ...+ c_{n_0}x^n + c_{0_1}y + ...+ c_{0_n}y^n + 
    c_{1_1}xy + c_{1_2}xy^2 + ... + c_{1_(n-1)}xy^{n-1}+ ... + c_{(n-1)_1}x^{n-1}y
        
    """
    def __init__(self, degree, xdomain=[-1, 1], ydomain=[-1, 1], 
                            xwindow=[-1, 1], ywindow=[-1,1], 
                            paramdim=1, **pars):
        """       
        Parameters
        ----------
        degree : int
            highest power of the polynomial, the number of terms 
            are degree+1
        xdomain : list or None
            domain of the x independent variable
        ydomain : list or None
            domain of the y independent variable
        xwindow : list or None
            range of the x independent variable
        ywindow : list or None
            range of the y independent variable
        paramdim : int
            number of parameter sets
        pars : dict
            keyword: value pairs, representing parameter_name: value
            
        Returns
        -------
        model : Poly2DModel
            2D polynomial model
        """
        self.ndim = 2
        self.outdim = 1
        super(Poly2DModel, self).__init__(degree, ndim=self.ndim,
                                          paramdim=paramdim, **pars)
        self.xdomain = xdomain
        self.ydomain = ydomain
        self.xwindow = xwindow
        self.ywindow = ywindow
    
    def mhorner(self, x, y, coeff):
        """
        Multivariate Horner's scheme
        
        Parameters
        --------------
        x, y : array
        coeff : array of coefficients in inverse lexical order
        """
        alpha = np.array(self._invlex())
        r0 = coeff[0] 
        r1 = r0 * 0.0
        r2 = r0 * 0.0
        karr = np.diff(alpha, axis=0)
        for n in range(len(karr)):
            if karr[n, 1] != 0:
                r2 = y * (r0 + r1 + r2)
                r1 = coeff[0] * 0.
            else:
                r1 = x * (r0 + r1)
            r0 = coeff[n+1] 
        return r0 + r1 + r2
        
    def deriv(self, x, y):
        """
        Derivatives with respect to parameters
        """
        if x.ndim == 2:
            x = x.flatten()
        if y.ndim == 2:
            y = y.flatten()
        if x.size != y.size:
            raise ValueError('Expected x and y to be of equal size')
        
        designx = x[:, None]**np.arange(self.deg+1)
        designy = y[:, None]**np.arange(1, self.deg+1)
        
        designmixed = []
        for i in range(1, self.deg):
            for j in range(1, self.deg):
                if i+j <= self.deg:
                    designmixed.append((x**i)*(y**j))
        designmixed = np.array(designmixed).T
        if designmixed.any():
            v = np.hstack([designx, designy, designmixed])
        else:
            v = np.hstack([designx, designy])
        return v

    def invlex_coeff(self):
        coeff = []
        lencoeff = range(self.deg + 1)
        for i in lencoeff:
            for j in lencoeff:
                if i+j <= self.deg:
                    name = 'c'+str(j)+'_'+str(i)
                    coeff.append(getattr(self, name))
        return np.array(coeff[::-1])
    
    def __call__(self, x, y):
        """
        Transforms data using this model.
        
        Parameters
        --------------
        x, y : arrays, of min dimensions 2
        
        Notes
        -----
        See the module docstring for rules for model evaluation. 
        """
        invcoeff = self.invlex_coeff()
        x, fmt = _convert_input(x, self.paramdim)
        y, fmt = _convert_input(y, self.paramdim)
        assert x.shape == y.shape, \
               "Expected input arrays to have the same shape"
        
        result = self.mhorner(x, y, invcoeff)
        return _convert_output(result, fmt)
        
class ICheb2DModel(IModel):
    """
    Chebyshev 2D polynomial:
    
    It is defined the same way as in IRAF.
    .. math:: P_{n_m}(x,y) = \sum C_{n_m}  T_n(x) T_m(y)
    """
    def __init__(self, xdeg, ydeg, xdomain=None, xwindow=[-1, 1], 
                            ydomain=None, ywindow=[-1,1], paramdim=1, **pars):
        
        """
        Parameters
        ----------
        
        xdeg : int
            degree in x
        ydeg : int
            degree in y
        xdomain : list or None
            domain of the x independent variable
        ydomain : list or None
            domain of the y independent variable
        xwindow : list or None
            range of the x independent variable
        ywindow : list or None
            range of the y independent variable
        paramdim : int
            number of parameter sets
        pars : dict
            keyword: value pairs, representing parameter_name: value
            
        Returns
        -------
        model : ICheb2DModel
            2D Chebyshev model
        """
        super(ICheb2DModel, self).__init__(xdeg, ydeg,
                                           xdomain=xdomain, ydomain=ydomain, 
                                           xwindow=xwindow, ywindow=ywindow,
                                           paramdim=paramdim, **pars)
                        
    
    def _fcache(self, x, y):
        """
        Calculate the individual Chebyshev functions once
        and store them in a dictionary to be reused.
        """
        xterms = self.xdeg+1
        yterms = self.ydeg+1
        kfunc = {}
        kfunc[0] = np.ones(x.shape)
        kfunc[1] = x.copy()
        kfunc[xterms] = np.ones(y.shape)
        kfunc[xterms+1] = y.copy()
        for n in range(2, xterms):
            kfunc[n] = 2*x*kfunc[n-1] - kfunc[n-2]
        for n in range(xterms+2, xterms+yterms):
            kfunc[n] = 2*y*kfunc[n-1] - kfunc[n-2]
        return kfunc
    
    def deriv(self, x, y):
        """
        Derivatives with respect to the coefficients.
        This is an array with Chebyshev polynomials:
        
        Tx0Ty0  Tx1Ty0...TxnTy0...TxnTym
        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")
        x = x.flatten()
        y = y.flatten()
        xderiv = self._chebderiv1d(x, self.xdeg+1).T
        yderiv = self._chebderiv1d(y, self.ydeg+1).T
        n = (self.xdeg+1)*(self.ydeg+1)
        v = np.empty((n, len(x)), dtype=x.dtype)
    
        ij = []
        for i in range(self.ydeg+1):
            for j in range(self.xdeg+1):
                ij.append(xderiv[j]*yderiv[i])
        v = np.array(ij)
        return v.T
    
    def _chebderiv1d(self, x, deg):
        """
        Derivative of 1D Chebyshev series
        """
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        d = np.empty((deg+1, len(x)), dtype=x.dtype)
        d[0] = x * 0 + 1
        if deg > 0 :
            x2 = 2*x
            d[1] = x
            for i in range(2, deg + 1) :
                d[i] = d[i-1]*x2 - d[i-2]
        return np.rollaxis(d, 0, d.ndim)
    
    def __call__(self, x, y, xdomain=None, ydomain=None):
        """
        Transforms data using this model.
        
        Parameters
        --------------
        x, y : arrays, of min dimensions 2
        xdomain, ydomain : list of two numbers
            polynomial domain for x and y variable
                    
        Notes
        -----
        See the module docstring for rules for model evaluation. 
        """
        assert x.shape == y.shape, \
               "Expected input arrays to have the same shape"
        invcoeff = self.invlex_coeff()
        x, fmt = _convert_input(x, self.paramdim)
        y, fmt = _convert_input(y, self.paramdim)
        result = self.imhorner(x, y, invcoeff)
        return _convert_output(result, fmt)
    
class ILegend2DModel(IModel):
    """
    Legendre 2D polynomial
    
    Defined as:
    .. math:: P_{nm}(x,y) = C_{n_m}  L_n(x ) L_m(y)

    """
    def __init__(self, xdeg, ydeg, xdomain=None, xwindow=[-1, 1], 
                            ydomain=None, ywindow=[-1, 1], paramdim=1, **pars):
        """
        Parameters
        ----------
        
        xdeg : int
            degree in x
        ydeg : int
            degree in y
        xdomain : list or None
            domain of the x independent variable
        ydomain : list or None
            domain of the y independent variable
        xwindow : list or None
            range of the x independent variable
        ywindow : list or None
            range of the y independent variable
        paramdim : int
            number of parameter sets
        pars : dict
            keyword: value pairs, representing parameter_name: value
            
        Returns
        -------
        model : ILegendModel
            2D Legendre model
        """
        super(ILegend2DModel, self).__init__(xdeg, ydeg,
                                             xdomain=xdomain, ydomain=ydomain, 
                                             xwindow=xwindow, ywindow=ywindow,
                                             paramdim=paramdim, **pars)

    def _fcache(self, x, y):
        """
        Calculate the individual Chebyshev functions once
        and store them in a dictionary to be reused.
        """
        xterms = self.xdeg+1
        yterms = self.ydeg+1
        kfunc = {}
        kfunc[0] = np.ones(x.shape)
        kfunc[1] = x.copy()
        kfunc[xterms] = np.ones(y.shape)
        kfunc[xterms+1] = y.copy()
        for n in range(2, xterms):
            kfunc[n] = (2*n+1)/(n+1)*kfunc[n-1] - n/(n+1)*kfunc[n-2]
        for n in range(2, yterms):
            kfunc[n+xterms] = (2*n+1)/(n+1)*kfunc[n+xterms-1] \
                                                - n/(n+1)*kfunc[n+xterms-2]
        return kfunc
    
    def deriv(self, x, y):
        """
        Derivatives with repect to the coefficients.
        So this is an array with Legender polynomials:
        
        Lx0Ly0  Lx1Ly0...LxnLy0...LxnLym
        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")
        x = x.flatten()
        y = y.flatten()
        xderiv = self._legendderiv1d(x, self.xdeg+1).T
        yderiv = self._legendderiv1d(y, self.ydeg+1).T
        n = (self.xdeg+1)*(self.ydeg+1)
        v = np.empty((n, len(x)), dtype=x.dtype)
    
        ij = []
        for i in range(self.ydeg+1):
            for j in range(self.xdeg+1):
                ij.append(xderiv[j]*yderiv[i])
                
        v = np.array(ij)
        return v.T

    def _legendderiv1d(self, x, deg):
        """
        Derivative of 1D Legendre polynomial
        """
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        d = np.empty((deg+1, len(x)), dtype=x.dtype)
        d[0] = x*0 + 1
        if deg > 0 :
            d[1] = x
            for i in range(2, deg + 1) :
                x2 = (2*i+1)*x
                d[i] = d[i-1]*x2 - d[i-2]*i/(i+1)
        return np.rollaxis(d, 0, d.ndim)
    
    def __call__(self, x, y, xdomain=None, ydomain=None):
        """
        Transforms data using this model.
        
        Parameters
        --------------
        x, y : arrays, of min dimensions 2
        xdomain, ydomain : list of two numbers
            polynomial domain for x and y variable
                    
        Notes
        -----
        See the module docstring for rules for model evaluation. 
        """
        assert x.shape == y.shape, \
               "Expected input arrays to have the same shape"
        invcoeff = self.invlex_coeff()
        x, fmt = _convert_input(x, self.paramdim)
        y, fmt = _convert_input(y, self.paramdim)
        result = self.imhorner(x, y, invcoeff)
        return _convert_output(result, fmt)
    
class Gauss1DModel(ParametricModel):
    """
    
    Implements 1D Gaussian model

    """
    parnames = ['amplitude', 'xcen', 'xsigma']
    def __init__(self, amplitude, xcen, fwhm=None, xsigma=None,
                            fjac=None, **cons):
        """
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
        fjac : callable or None
            if callable - a function to compute the Jacobian of 
            func with derivatives across the rows.
            if None - the Jacobian will be estimated
        
        Returns
        -------
        model : Gauss1DModel
            1D Gaussian
        """
        self._amplitude = parameters._Parameter('amplitude', amplitude, self, 1)
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
        self._xsigma = parameters._Parameter('xsigma', xsigmaval, self, 1)
        self._xcen = parameters._Parameter('xcen', xcen, self, 1)
        self.ndim = 1
        self.outdim = 1
        
        try:
            paramdim = len(self._amplitude)
            assert (len(amplitude) == len(xsigmaval) == len(xcen) ), \
             "Input parameters do not have the same dimension"
        except TypeError:
            paramdim = 1
        super(Gauss1DModel, self).__init__(self.parnames, paramdim=paramdim,
                                                                    **cons)
        self.linear = False
        if fjac:
            self.deriv = fjac
            
    def eval(self, x, params):
        return params[0] * np.exp((-(1/(params[2]**2)) * 
                                                (x-params[1])**2))
 
    def deriv(self, p, x, y):
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
    
    2D Gaussian
    
    """
    parnames = ['amplitude', 'xcen', 'ycen', 'xsigma', 'ysigma', 'theta']
    
    def __init__(self, amplitude, xcen, ycen, fwhm=None, xsigma=None,
                 ysigma=None, ratio=None, theta=0.0, fjac=None, **cons):
        """
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
        
        Returns
        -------
        model : Gauss2DModel
            2D Gaussian
        """
        if ysigma is None and ratio is None:
            raise InputParameterError(
                "Either ysigma or ratio must be specified")
        elif xsigma is None and fwhm is None:
            raise InputParameterError(
                "Either fwhm or xsigma must be specified")
        self._amplitude = parameters._Parameter('amplitude', amplitude, self, 1)
        if xsigma is None:
            xsigma = 0.42466 * fwhm
        self._xsigma = parameters._Parameter('xsigma', xsigma, self, 1)
        if ysigma is None:
            ysigma = ratio * self._xsigma
        self._ysigma = parameters._Parameter('ysigma', ysigma, self, 1)
        self._xcen = parameters._Parameter('xcen', xcen, self, 1)
        self._ycen = parameters._Parameter('ycen', ycen, self, 1)
        self._theta = parameters._Parameter('theta', theta, self, 1)
        
        self.ndim = 2
        self.outdim = 1
        try:
            paramdim = len(self._amplitude)
            assert (len(self._amplitude) == len(self._xsigma) == \
                            len(self._xcen) == len(self._ycen) == \
                            len(self._theta) ), \
                            "Input parameters do not have the same dimension"
        except TypeError:
            paramdim = 1
        super(Gauss2DModel, self).__init__(self.parnames, paramdim=paramdim)
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
        x, fmt = _convert_input(x, self.paramdim)
        y, fmt = _convert_input(y, self.paramdim)
        result = self.eval(x, y, self.psets)
        return _convert_output(result, fmt)

class ShiftModel(Model):
    """
    Shift a coordinate

    """
    parnames = ['offsets']
    def __init__(self, offsets): #, multiple=False):
        """
        Parameters
        ----------
        offsets : float or a list of floats
            offsets to be applied to a coordinate
            if a list - each value in the list is an offset to be applied to a
            column in the input coordinate array
            
        Returns
        -------
        model : ShiftModel
            A model representing offset
        """
        self.ndim = 1
        self.outdim = 1
        
        if not operator.isSequenceType(offsets):
            paramdim = 1
        else:
            paramdim = len(offsets)
        self._offsets = parameters._Parameter('offsets', offsets, self, paramdim)
        super(ShiftModel, self).__init__(self.parnames, paramdim=paramdim)

    def __call__(self, x):
        """
        Transforms data using this model.
        """
        x, fmt = _convert_input(x, self.paramdim)
        result = x + self.offsets
        return _convert_output(result, fmt)
 
class ScaleModel(Model):
    """
    
    Scale a coordinate

    """
    parnames = ['factors']
    def __init__(self, factors):
        """
        Parameters
        ---------------
        factors : float or a list of floats
            scale for a coordinate
        Returns
        -------
        model : ScaleModel
            Model representing scaling 
            
        """
        self.ndim = 1
        self.outdim = 1
        
        if not operator.isSequenceType(factors):
            paramdim = 1
        else:
            paramdim = len(factors)
        self._factors = parameters._Parameter('factors', factors, self, paramdim)
        super(ScaleModel, self).__init__(self.parnames, paramdim=paramdim)
    
    def __call__(self, x):
        """
        Transforms data using this model.
        """
        x, fmt = _convert_input(x, self.paramdim)
        result = x * self.factors
        return _convert_output(result, fmt)

class _SIP1D(Model):
    """
    This implements the Simple Imaging Protocol Model in 1D. 
    
    It's unlikely it will be used in 1D so this class is private
    and SIPModel should be used instead.
    
    """
    def __init__(self, order, coeffname='a', paramdim=1, **pars):
        self.order = order
        self.ndim = 2
        self.outdim = 1
        self.coeffname = coeffname.lower()
        self.parnames = self._generate_coeff_names(coeffname)
        
        if not pars:
            self.set_coeff(pardim=paramdim)
        else:
            p = pars.get('{0}02'.format(coeffname, None))
            if operator.isSequenceType(p):
                lenpars = len(p)
            else:
                lenpars = 1
            if paramdim != lenpars:
                print("Creating a model with {0} parameter sets\n".format(lenpars))
                paramdim = lenpars
            self._validate_pars(**pars)  
            self.set_coeff(pardim=paramdim, **pars)
        
        super(_SIP1D, self).__init__(self.parnames, paramdim=paramdim)
       
    def __repr__(self):
        fmt = """
        Model: {0}
        Dim:   {1}
        Order: {2}
        Parameter sets: {3}
        """.format(
              self.__class__.__name__,
              self.ndim,
              self.order,
              self.paramdim
                )
        return fmt
    
    def __str__(self):
        fmt = """
        Model: {0}
        Dim:   {1}
        Order: {2}
        Parameter sets: {3}
        Parameters: 
                   {4}
        """.format(
              self.__class__.__name__,
              self.ndim,
              self.order,
              self.paramdim,
              "\n                   ".join(i+':  ' + str(getattr(self,i)) for
                                           i in self.parnames)
                )
        return fmt
    
    def get_numcoeff(self):
        """
        Return the number of coefficients in one parset
        """
        if self.order < 2  or self.order > 9:
            raise ValueError("Degree of polynomial must be 2< deg < 9")
        nmixed = comb(self.order-1, self.ndim)
        numc = self.order * self.ndim + nmixed + 1
        return numc
    
    def _generate_coeff_names(self, coeffname):
        names = []
        for i in range(2, self.order+1):
            names.append('{0}{1}{2}'.format(coeffname, i, 0))
        for i in range(2, self.order+1):
            names.append('{0}{1}{2}'.format(coeffname, 0, i))
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i+j < self.order+1:
                    names.append('{0}{1}{2}'.format(coeffname, i, j))
        return names
    
    def set_coeff(self, pardim=1, **pars):
        if not pars:
            # default values
            for name in self.parnames:
                if pardim == 1:
                    self.__setattr__('_'+name, 
                                     parameters._Parameter(name, 0, self, 1))
                else:
                    self.__setattr__('_'+name, 
                                     parameters._Parameter(name, [0]*pardim,
                                                           self, pardim))
        else:
            for name in self.parnames:
                self.__setattr__('_'+name, 
                                 parameters._Parameter(name, pars[name],
                                                       self, pardim))
                
    def _validate_pars(self, **pars):
        numcoeff = self.get_numcoeff()
        assert(len(pars) == numcoeff)
 
    def _coef_matrix(self, coeffname):
        mat = np.zeros((self.order+1, self.order+1))
        for i in range(2, self.order+1):
            mat[i, 0] = getattr(self, '{0}{1}{2}'.format(coeffname, i, 0))[0]
        for i in range(2, self.order+1):
            mat[0, i] = getattr(self, '{0}{1}{2}'.format(coeffname, 0, i))[0]
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i+j < self.order+1:
                    mat[i, j] = getattr(self, '{0}{1}{2}'.format(coeffname, i, j))[0]
        return mat

    def _eval_sip(self, x, y, coef):
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if self.coeffname == 'a':
            result = np.zeros(x.shape)
        else:
            result = np.zeros(y.shape)
        
        for i in range(coef.shape[0]):
            for j in range(coef.shape[1]):
                if i+j > 1 and i+j < self.order+1:
                    result = result+coef[i, j]*x**i*y**j
        return result
        
    def __call__(self, x, y):
        mcoef = self._coef_matrix(self.coeffname)
        return self._eval_sip(x, y, mcoef)

class LabeledInput(dict):
    """
    Create a container with all input data arrays, assigning labels for each one.
    
    Used by CompositeModel to choose input data using labels
    
    Examples
    --------
    >>> x,y = np.mgrid[:10, :10]
    >>> l = np.arange(10)
    >>> ado = LabeledInput([x, y, l], ['x', 'y', 'pixel'])
    >>> ado.x      
    array([[0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 2],
    [3, 3, 3, 3, 3],
    [4, 4, 4, 4, 4]])    
    >>> ado['x']
    array([[0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 2],
    [3, 3, 3, 3, 3],
    [4, 4, 4, 4, 4]])
        
    """
    def __init__(self,  data, labels):
        """
        Parameters
        ----------
        data : list 
            a list of all input data
        labels : list of strings
            names matching each coordinate in data
        
        Returns
        -------
        data : LabeledData
            a dict of input data and their assigned labels
       """
        dict.__init__(self)
        assert len(labels) == len(data)
        self.labels = [l.strip() for l in labels]
        for coord, label in zip(data, labels):
            self[label] = coord
            setattr(self, '_'+label, coord)
        self._set_properties(self.labels)
    
    def _getlabel(self, name):
        par = getattr(self, '_'+name)
        return par
    
    def _setlabel(self,  name, val):
        setattr(self, '_'+name, val)
        self[name] = val
        
    def _dellabel(self,  name):
        delattr( self,  '_'+name)
        del self[name]
        
    def add(self, label=None, value=None,  **kw):
        """
        Add input data to a LabeledInput object
        
        Parameters
        --------------
        label : string
            coordinate label
        value : numerical type
            coordinate value
        kw : dictionary
            if given this is a dictionary of {label: value} pairs
            
        """
        if kw:
            if label is None or value is None:
                self.update(kw)
            else:
                kw[label] = value
                self.update(kw)
        else:
            kw = dict({label: value})
            assert(label is not None and value is not None), (
                "Expected label and value to be defined")
            self[label] = value
            
        for key in kw:
            self.__setattr__('_'+key, kw[key])
        self._set_properties(kw.keys())
            
    def _set_properties(self, attributes):
        for attr in attributes:
            setattr(self.__class__, attr , property(lambda self, attr=attr:
                                                  self._getlabel(attr),
                                                  lambda self, value, attr=attr:
                                                  self._setlabel(attr, value),
                                                  lambda self, attr=attr:
                                                  self._dellabel(attr)
                                                  )
                    )
            
    def copy(self):
        data = [self[label] for label in self.labels]
        return LabeledInput(data, self.labels)
    
class _CompositeModel(OrderedDict):
    def __init__(self, transforms, inmap=None, outmap=None):
        """
        A Base class for all composite models

        """
        OrderedDict.__init__(self)
        self.ndim = None
        self.outdim = None
        self.fittable = False
        self.has_inverse = np.array([tr.has_inverse for tr in transforms]).all()
        
    def _init_comptr(self, trans, inmap, outmap):
        # implemented by subclasses
        raise NotImplementedError("Subclasses should implement this")
    
    def __repr__(self):
        transforms = self.keys()
        fmt = """
            Model:  {0} 
            """.format(self.__class__.__name__)
        fmt1 = " %s  " * len(transforms)% tuple([repr(tr) for tr in transforms])
        fmt = fmt + fmt1
        return fmt
    
    def __str__(self):
        transforms = self.keys()
        fmt = """
            Model:  {0}
            """.format(self.__class__.__name__)
        fmt1 = " %s  " * len(transforms)% tuple([str(tr) for tr in transforms])
        fmt = fmt + fmt1
        return fmt
            
    def add_model(self, transf, inmap, outmap):
        self[transf] = [inmap, outmap]
 
    def invert(self):
        raise NotImplementedError("Subclasses should implement this")
            
    def __call__(self):
        # implemented by subclasses
        raise NotImplementedError("Subclasses should implement this")

class SCompositeModel(_CompositeModel):
    """
    
    Execute models in series
    
    Notes
    -----
    Output values of one model are used as input values of another.
    Obviously the order of the models matters.
    
    Examples
    --------
    Apply a 2D rotation followed by a shift in x and y
    
    >>> from fitting import rotations
    >>> rot = rotations.MatrixRotation2D(angle=23.5)
    >>> offx = ShiftModel(-4.23)
    >>> offy = ShiftModel(2)
    >>> linp = LabeledInput([x,y], ["x", "y"]
    >>> scomptr = SCompositeModel([rot, offx, offy], 
                                  inmap=[['x', 'y'], ['x'], ['y']],
                                  outmap=[['x', 'y'], ['x'], ['y']])
    >>> result=scomptr(linp)
        
    """
    def __init__(self, transforms, inmap=None, outmap=None):
        """
        
        Parameters
        ----------
        transforms : list
            a list of transforms in the order to be executed
        inmap : list of lists or None
            labels in an input instance of LabeledInput
            if None, the number of input coordinates is exactly what
            the transforms expect 
        outmap : list or None
            labels in an input instance of LabeledInput
            if None, the number of output coordinates is exactly what 
            the transforms expect
        
        Returns
        -------
        model : SCompositeModel
            Composite model which executes the comprising models in series
        """
        super(SCompositeModel, self).__init__(transforms, inmap, outmap)
        if transforms and inmap and outmap:
            assert len(transforms) == len(inmap) == len(outmap), \
                   "Expected sequences of transform, " \
                   "inmap and outmap to have the same length"
        if inmap is None:
            inmap = [None] * len(transforms)
        if outmap is None:
            outmap = [None]  * len(transforms)
        
        self._init_comptr(transforms, inmap, outmap)
        self.ndim = np.array([tr.ndim for tr in self]).max()
        # the output dimension is equal to the output dim of the last transform
        self.outdim = self.keys()[-1].outdim
        
    def _init_comptr(self, transforms, inmap, outmap):
        for tr, inm, outm in zip(transforms, inmap, outmap):
            self[tr] = [inm, outm]
     
    def _verify_no_mapper_input(self, *data):
        lendata = len(data)
        tr = self.keys()[0]
        
        if tr.ndim != lendata:
            
            raise ValueError("Required number of coordinates not matched for "
                             "transform # {0}: {1} required, {2} supplied ".format( 
                             self.keys().index(tr)+1, tr.ndim, lendata))

    def invert(self, inmap, outmap):
        scomptr = SCompositeModel(self[::-1], inmap=inmap, outmap=outmap)
        return scomptr
            
    def __call__(self, x, *data):
        """
        Transforms data using this model.
        """
        lendata = len(data) + 1
        if lendata == 1:
            if not isinstance(x, LabeledInput):
                data = np.asarray(x, dtype=np.float64)
                shape = data.shape
                self._verify_no_mapper_input(data)
                result = data
                for tr in self:
                    result = tr(result)
                return result
            else:
                linp = x.copy()
                # we want to return the entire labeled object because some parts
                # of it may not be used in another transform of which this 
                # one is a component
                for tr in self:
                    inmap = self[tr][0]
                    outmap = self[tr][1]
                    inlist = [getattr(linp, co) for co in inmap]
                    result = tr(*inlist)
                    if tr.outdim == 1:
                        result = [result]
                    for outcoo, res in zip(outmap, result):
                        if outcoo not in inmap:
                            linp.add(outcoo, res)
                        else:
                            linp[outcoo] = res
                            setattr(linp, outcoo, res)
                return linp
        else:
            inlist = [x]
            inlist.extend(data)
            self._verify_no_mapper_input(*inlist)
            result = self.keys()[0](*inlist)
            for tr in self.keys()[1:]:
                result = tr(result)
            return result
            
class PCompositeModel(_CompositeModel):
    """
    
    Execute models in parallel
    
    Notes
    -----
    Models are applied to input data separately and the deltas are summed.
    
    """
    def __init__(self, transforms, inmap=None, outmap=None):
        """
        Parameters
        --------------
        transforms : list
            transforms to be executed in parallel
        inmap : list or None
            labels in an input instance of LabeledInput
            if None, the number of input coordinates is exactly what the
            transforms expect 
            
        Returns
        -------
        model : PCompositeModel
            Composite model which executes the comprising models in parallel
        """
        super(PCompositeModel, self).__init__(transforms,
                                              inmap=None, outmap=None)
        self._init_comptr(transforms, inmap, outmap)
        self.ndim = self.keys()[0].ndim
        self.outdim = self.ndim
        self.inmap = inmap
        self.outmap = outmap

    def _init_comptr(self, transforms, inmap, outmap):
        for tr in transforms:
            self[tr] = [inmap, outmap]

    def _verify_no_mapper_input(self, *data):
        lendata = len(data)
        ndim = self.keys()[0].ndim
        for tr in self.keys():
            if tr.ndim != ndim:
                raise ValueError("tr.ndim ...")
    
    def invert(self, inmap, outmap):
        pcomptr = PCompositeModel(self.keys()[::-1], inmap=inmap, outmap=outmap)
        return pcomptr
    
    def __call__(self, x, *data):
        """
        Transforms data using this model.
        """
        lendata = len(data) + 1
        if lendata == 1:
            if not isinstance(x, LabeledInput):
                self._verify_no_mapper_input(x)
                result = x.copy()
                for tr in self:
                    delta = tr(x) - x
                    result = result + delta
                return result
            else:
                assert self.inmap is not None, ("Parameter 'inmap' must be "
                                    "provided when input is a labeled object")
                assert self.outmap is not None, ("Parameter 'outmap' must be "
                                    "provided when input is a labeled object")
                linp = x.copy()
                #create a list of inputs to be passed to the transforms
                inlist = [getattr(linp, co) for co in self.inmap]
                #create a list of outputs to which the deltas are applied
                result = [getattr(linp, co) for co in self.outmap]
                res = [tr(*inlist) for tr in self]
                delta = (np.asarray(res) - np.asarray(result)).sum(axis=0)
                result = np.asarray(result) + delta
                for outcoo, res in zip(self.outmap, result):
                    linp[outcoo] = res
                    setattr(linp, outcoo, res)
                #always return the entire labeled object, not just the result
                # since this may be part of another composite transform
                return linp
        else:
            self._verify_no_mapper_input(x, *data)
            inlist = [x]
            inlist.extend(list(data))
            result = inlist[:]
            for tr in self.keys():
                res = tr(*inlist)
                for i in range(len(inlist)):
                    result[i] = res[i]-inlist[i]
            return result
        
class SIPModel(SCompositeModel):
    """
    
    Simple Imaging Protocol (SIP) model [1]_
    
    References
    ----------
    .. [1] David Shupe, et al, "The SIP Convention of Representing Distortion
       In FITS Image Headers",  Astronomical Data Analysis Software And Systems, 
       ASP Conference Series, Vol. 347, 2005
       
    """
    def __init__(self, crpix, order, coeff, coeffname='a', 
                            aporder=None, apcoeff=None, paramdim=1):
        """
        Parameters
        ----------
        crpix : list or ndarray of length(2)
            CRPIX values
        order : int
            SIP polynomial order
        coeff : dict
            SIP coefficients
        coeffname : string: 'a', 'b', 'A' or 'B'
            SIP coefficient preffix
        aporder : int
            order for the inverse transformation
        apcoeff : dict
            coefficients for the inverse transform
        paramdim : int
            number of parameter sets
        multiple : boolean
            when input is 2D array, if True (default) it is to be 
            treated as multiple 1D arrays
            
        Returns
        -------
        model : SIPModel
            A model representing the Simple Imaging Protocol
        """
        self.ndim = 2
        self.outdim = 1
        self.shifta = ShiftModel(crpix[0])
        self.shiftb = ShiftModel(crpix[1])
        self.sip1d = _SIP1D(order, coeffname=coeffname,
                            paramdim=paramdim, **coeff)
        if aporder is not None and apcoeff is not None:
            self.inverse = Poly1DModel(aporder, **apcoeff)
        else:
            self.inverse = None
        super(SIPModel, self).__init__([self.shifta, self.shiftb, self.sip1d], 
                                       inmap = [['x'], ['y'], ['x', 'y']], 
                                       outmap=[['x'], ['y'], ['z']])
        
    def __repr__(self):
        models = [self.shifta, self.shiftb, self.sip1d]
        fmt = """
            Model:  {0} 
            Coeff Prefix: {1}
            """.format(self.__class__.__name__, self.sip1d.coeffname.upper())
        fmt1 = " %s  " * len(models)% tuple([repr(model) for model in models])
        fmt = fmt + fmt1
        return fmt
            
    def __str__(self):
        models = [self.shifta, self.shiftb, self.sip1d]
        fmt = """
            Model:  {0}
            Coeff Prefix: {1}
            """.format(self.__class__.__name__, self.sip1d.coeffname.upper())
        fmt1 = " %s  " * len(models)% tuple([str(model) for model in models])
        fmt = fmt + fmt1
        return fmt
    
    def __call__(self, x, y):
        """
        Transforms data using this model.
        """
        ado = LabeledInput([x, y], ['x', 'y'])
        return SCompositeModel.__call__(self, ado).z
