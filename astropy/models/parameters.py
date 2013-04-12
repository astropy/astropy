# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines two classes that deal with parameters.
It is unlikely users will need to work with these classes directly,
unless they define their own models.

"""
from __future__ import division, print_function
import numbers
import numpy as np
from ..utils import misc
from .utils import InputParameterError

__all__ = ['Parameters', 'Parameter']

def _tofloat(value):
    """
    Convert a parameter to float or float array
   
    """
    if misc.isiterable(value):
        try:
            _value = np.array(value, dtype=np.float)
        except (TypeError, ValueError): 
            #catch arrays with strings or user errors like different 
            # types of parameters in a parameter set
            raise InputParameterError(
                "Parameter of {0} could not be converted to "
                "float".format(type(value)))
    elif isinstance(value, bool):
        raise InputParameterError(
            "Expected parameter to be of numerical type, not boolean")
    elif isinstance(value, numbers.Number):
        _value = np.array(value, dtype=np.float)
    else:
        raise InputParameterError(
            "Don't know how to convert parameter of {0} to "
            "float".format(type(value)))
    return _value
    
class Parameter(list):
    """
    Wraps individual parameters.
    
    This class represents a model's parameter (in a somewhat broad
    sense). To support multiple parameter sets, a 
    parameter has a dimension (paramdim). To support some level
    of validation a parameter has also a shape (parshape).
    Parameter objects behave like numbers.

    """
    def __init__(self, name, val, mclass, paramdim, fixed=False, tied=False, 
                 minvalue=None, maxvalue=None):
        """
        Parameters
        ----------
        name : string
            parameter name
        val :  number or an iterable of numbers
        mclass : object
            an instance of a Model class
        paramdim : int
            parameter dimension      
        """
        self._paramdim = paramdim
        if isinstance(val, numbers.Number):
            if self.paramdim == 1:
                val = _tofloat(val)[()]
                super(Parameter, self).__init__([val])
                self.parshape = val.shape
            else:
                val = [_tofloat(val)[()]]
                super(Parameter, self).__init__(val)
                self.parshape = val[0].shape
        # colections.Sequence covers lists but not ndarrays
        # which are checked for in _tofloat()
        # misc.iterable allows dict which is failed in _tofloat()
        elif misc.isiterable(val):
            if paramdim == 1:
                val = [_tofloat(value)[()] for value in val]
                super(Parameter, self).__init__(val)
                self.parshape = _tofloat(val).shape
            else:
                val = [_tofloat(value)[()] for value in val]
                super(Parameter, self).__init__(val)
                self.parshape = _tofloat(val[0]).shape
        else:
            raise InputParameterError(
                "Parameter {0} is not a number".format(name))
        self._mclass = mclass
        self._name = name
        self._fixed = fixed
        self._tied = tied
        self._min = minvalue
        self._max = maxvalue
        
    @property
    def paramdim(self):
        return self._paramdim
    
    @paramdim.setter
    def model(self, val):
        self._paramdim = val
        
    @property
    def mclass(self):
        return self._mclass
    
    @mclass.setter
    def model(self, val):
        self._mclass = val
        
    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, val):
        self._name = val
        
    @property
    def fixed(self):
        return self._fixed
    
    @fixed.setter
    def fixed(self, val):
        assert isinstance(val, bool), "Fixed can be True or False"
        self._fixed = val
        self.mclass.constraints._fixed.update({self.name: val})
        self.mclass.constraints._update()
        
    @property
    def tied(self):
        return self._tied
    
    @tied.setter
    def tied(self, val):
        assert callable(val), "Tied must be a callable"
        self._tied = val
        self.mclass.constraints._tied.update({self.name:val})
        self.mclass.constraints._update()
        
    @property
    def min(self):
        return self._min
    
    @min.setter
    def min(self, val):
        assert isinstance(val, numbers.Number), "Min value must be a number"
        self._min = float(val)
        self.mclass.constraints.set_range({self.name: (val, self.max or 1E12)})
        
    @property
    def max(self):
        return self._max
    
    @max.setter
    def max(self, val):
        assert isinstance(val, numbers.Number), "Max value must be a number"
        self._max = float(val)
        self.mclass.constraints.set_range({self.name: (self.min or -1E12, val)})

    def __setslice__(self, i, j, val):
        super(Parameter, self).__setslice__(i, j, val)
        setattr(self.mclass, self.name, self)
        
    def __setitem__(self, i, val):
        super(Parameter, self).__setitem__(i, val)
        setattr(self.mclass, self.name, self)
    
    def __add__(self, val):
        return np.asarray(self) + val
    
    def __radd__(self, val):
        return np.asarray(self) + val
    
    def __sub__(self, val):
        return np.asarray(self) - val
    
    def __rsub__(self, val):
        return val - np.asarray(self)
    
    def __mul__(self, val):
        return np.asarray(self) * val
    
    def __rmul__(self, val):
        return np.asarray(self) * val
    
    def __pow__(self, val):
        return np.asarray(self)**val
    
    def __div__(self, val):
        return np.asarray(self) / val
    
    def __rdiv__(self, val):
        return val / np.asarray(self)
    
    def __truediv__(self, val):
        return np.asarray(self) / val
    
    def __rtruediv__(self, val):
        return val / np.asarray(self)
    
    def __eq__(self, val):
        return (np.asarray(self) == np.asarray(val)).all()
    
    def __ne__(self, val):
        return not (np.asarray(self) == np.asarray(val)).all()

    def __lt__(self, val):
        return (np.asarray(self) < np.asarray(val)).all()

    def __gt__(self, val):
        return (np.asarray(self) > np.asarray(val)).all()

    def __le__(self, val):
        return (np.asarray(self) <= np.asarray(val)).all()

    def __ge__(self, val):
        return (np.asarray(self) >= np.asarray(val)).all()

    def __neg__(self):
        return np.asarray(self) * (-1)
    
    def __abs__(self):
        return np.abs(np.asarray(self))
    
class Parameters(list):
    """
    Store model parameters as a flat list of floats.
    
    The `Parameters` class is a list-like object which
    stores model parameters. Only  instances of ParametricModel 
    keep an instance of this class as an attribute. The list of parameters
    can be modified by the user or by an instance of `fitting.Fitter`. 
    This list of parameters is kept in sync with single model parameter attributes.
    When more than one dimensional, a `fitting.Fitter` treats each
    set of parameters as belonging to the same model but different set of data.
    
    """
    def __init__(self, mobj, parnames, paramdim=1):
        """
        Parameters
        ----------
        mobj : object
            an instance of a subclass of `fitting.ParametricModel`
        parnames : list of strings
            parameter names
        paramdim : int
            Number of parameter sets
              
        """
        self.mobj = mobj
        self.paramdim = paramdim
        # A flag set to True by a fitter to indicate that the flat 
        # list of parameters has been changed.
        self._changed = False
        self.parinfo = {}
        parlist = [getattr(mobj, attr) for attr in parnames]
        flat = self._flatten(parnames, parlist)
        super(Parameters, self).__init__(flat)

    def __setitem__(self, ind, value):
        _val = _tofloat(value)
        super(Parameters, self).__setitem__(ind, _val[()])
        self._changed = True
        self._update_model_pars() 
        
    def __setslice__(self, istart, istop, vlist):
        super(Parameters, self).__setslice__(istart, istop, vlist)        
        self._changed = True
        self._update_model_pars()
        
    def _update_model_pars(self):
        """
        Update single parameters
        
        """
        for key in self.parinfo.keys():
            sl = self.parinfo[key][0]
            par = self[sl]
            if len(par) == 1:
                par = Parameter(key, par[0], self.mobj, self.mobj.paramdim)
            else:
                par = Parameter(key, par, self.mobj, self.mobj.paramdim)
            setattr(self.mobj, key, par)
        self._changed = False

    def _is_same_length(self, newpars):
        """
        Checks if the user supplied value of `ParametricModel.parameters`
        has the same length as the original parameters list.
        
        """
        parsize = _tofloat(newpars).size
        return parsize == self.__len__()

    def _flatten(self, parnames, parlist):
        """
        Create a list of model parameters
        """
        flatpars = []
        start = 0
        for (name, par) in zip(parnames, parlist):
            pararr = np.array(par)
            fpararr = pararr.flatten()

            stop = start+len(fpararr)  
            self.parinfo[name] = slice(start, stop, 1), pararr.shape
            start = stop
            flatpars.extend(list(fpararr))
        return flatpars
    
        
    
