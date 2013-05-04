# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines base classes for all models.
The base class of all models is `~astropy.models.Model`.
`~astropy.models.ParametricModel` is the base class for all fittable models. Parametric
models can be linear or nonlinear in a regression analysis sense.

All models provide a `__call__` method which performs the transformation in a
purely mathematical way, i.e. the models are unitless. In addition, when 
possible the transformation is done using multiple parameter sets, `param_sets`.
The number of parameter sets is stored in an attribute `param_dim`. 

Parametric models also store a flat list of all parameters as an instance of
`~astropy.models.parameters.Parameters`. When fitting, this list-like object is modified by a
subclass of `~astropy.models.fitting.Fitter`. When fitting nonlinear models, the values of the
parameters are used as initial guesses by the fitting class. Normally users
will not have to use the `~astropy.models.parameters` module directly.

Input Format For Model Evaluation and Fitting

Input coordinates are passed in separate arguments, for example 2D models 
expect x and y coordinates to be passed separately as two scalars or aray-like
objects. 
The evaluation depends on the input dimensions and the number of parameter
sets but in general normal broadcasting rules apply.
For example:

- A model with one parameter set works with input in any dimensionality

- A model with N parameter sets works with 2D arrays of shape (M, N).
  A parameter set is applied to each column.

- A model with N parameter sets works with multidimensional arrays if the
  shape of the input array is (N, M, P). A parameter set is applied to each plane.
  
In all these cases the output has the same shape as the input.

- A model with N parameter sets works with 1D input arrays. The shape 
  of the output is (M, N)
 
"""
from __future__ import division, print_function
import collections
import abc
from ..utils.compat.odict import OrderedDict
import numpy as np
from . import parameters
from . import constraints
from .utils import pmapdomain, InputParameterError, comb
 
__all__ = ['Model', 'ParametricModel', 'PCompositeModel', 'SCompositeModel',
                'LabeledInput', '_convert_input', '_convert_output']
        
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
    Create a property for a parameter.
    
    Parameters
    ----------
    name: string
        the name of the parameter
        
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
                par = parameters.Parameter(self.name, val, obj, obj.param_dim)
                oldpar = getattr(obj, self.name)
                if oldpar is not None and oldpar.parshape != par.parshape:
                    raise InputParameterError(
                        "Input parameter {0} does not "
                        "have the required shape".format(self.name))
                else:
                    setattr(obj, self.aname, par)
                obj._parameters = parameters.Parameters(obj, 
                                                        obj.param_names,
                                                         param_dim=obj.param_dim)
            else:
                setattr(obj, self.aname, val)
        else:
            par = parameters.Parameter(self.name, val, obj, obj.param_dim)
            oldpar = getattr(obj, self.name)
            if oldpar is not None and oldpar.parshape != par.parshape:
                raise InputParameterError(
                    "Input parameter {0} does not "
                    "have the required shape".format(self.name))
            else:
                setattr(obj, self.aname, par)
                
class Model(object):
    """
    Base class for all models.
    
    This is an abstract class and should not be instanciated.
    
    Notes
    -----
    Models which are not meant to be fit to data should subclass this class
    
    This class sets the properties for all individual parameters and performs
    parameter validation.
    
    """
    __metaclass__ = abc.ABCMeta
    
    param_names = []

    def __init__(self, param_names, ndim, outdim, param_dim=1):
        self._param_dim = param_dim
        self._ndim = ndim
        self._outdim = outdim
        self.has_inverse = False
        self._param_names = param_names
        #_parcheck is a dictionary to register parameter validation funcitons
        #key: value pairs are parameter_name: parameter_validation_function_name
        #see projections.AZP for example
        self._parcheck = {}
        for par in self.param_names:
            setattr(self.__class__, par, _ParameterProperty(par))
    
    @property 
    def ndim(self):
        """
        Number of input variables in model evaluation.
        """
        return self._ndim

    @property
    def outdim(self):
        """
        Number of output valiables returned when a model is evaluated.
        """
        return self._outdim
    
    @property
    def param_dim(self):
        """
        Number of parameter sets in a model.
        """
        return self._param_dim
    
    @param_dim.setter
    def param_dim(self, val):
        self._param_dim = val
    
    @property
    def param_names(self):
        """
        A list of names of the parameters defining a model.
        """
        return self._param_names
    
    @param_names.setter
    def param_names(self, val):
        self._param_names = val
        
    def __repr__(self):
        fmt = "{0}(".format(self.__class__.__name__)
        for i in range(len(self.param_names)):
            fmt1 = """
            {0}={1},
            """.format(self.param_names[i], getattr(self, self.param_names[i]))
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
              self.param_dim,
              "\n                   ".join(i+': ' + 
                str(self.__getattribute__(i)) for i in self.param_names)
                )
            
        return fmt

    @property
    def param_sets(self):
        """
        Return parameters as a pset.
        This is an array where each column represents one parameter set.
        """
        parameters = [getattr(self, attr) for attr in self.param_names]
        shapes = [par.parshape for par in parameters]
        lenshapes = np.asarray([len(p.parshape) for p in parameters])
        shapes = [p.parshape for p in parameters]
        if (lenshapes>1).any():
            if () in shapes: 
                psets = np.asarray(parameters, dtype=np.object)
            else:
                psets = np.asarray(parameters)
        else:
            psets = np.asarray(parameters)
            psets.shape = (len(self.param_names), self.param_dim)
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
    Base class for all fittable models.
    
    Notes
    -----
    All models which can be fit to data and provide a `deriv` method
    should subclass this class.
    
    Sets the parameters attributes.
    
    Parameters
    ----------
    param_names: list
        parameter names
    ndim: int
        model dimensions (number of inputs)
    outdim: int
        number of output quantities
    param_dim: int
        number of parameter sets
    fittable: boolean
        indicator if the model is fittable
    fixed: a dict
        a dictionary {parameter_name: boolean} of parameters to not be
        varied during fitting. True means the parameter is held fixed.
        Alternatively the `~astropy.models.parameters.Parameter.fixed`
        property of a parameter may be used.
    tied: dict
        a dictionary {parameter_name: callable} of parameters which are
        linked to some other parameter. The dictionary values are callables
        providing the linking relationship. 
        Alternatively the `~astropy.models.parameters.Parameter.tied`
        property of a parameter may be used.
    bounds: dict
        a dictionary {parameter_name: boolean} of lower and upper bounds of
        parameters. Keys  are parameter names. Values  are a list of length 
        2 giving the desired range for the parameter. 
        Alternatively the `~astropy.models.parameters.Parameter.min` and 
        `~astropy.models.parameters.Parameter.max` properties of a parameter 
        may be used.
    eqcons: list
        A list of functions of length n such that
        eqcons[j](x0,*args) == 0.0 in a successfully optimized
        problem.
    ineqcons : list
        A list of functions of length n such that
        ieqcons[j](x0,*args) >= 0.0 is a successfully optimized
        problem.
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, param_names, ndim, outdim, param_dim=1, fittable=True,
                 fixed={}, tied={}, bounds={}, eqcons=[], ineqcons=[]):
        self.linear = True
        super(ParametricModel, self).__init__(param_names, ndim, outdim, param_dim=param_dim)
        self.fittable = fittable
        self._parameters = parameters.Parameters(self, self.param_names,
                                                 param_dim=param_dim)
        _fixed = {}.fromkeys(self.param_names, False)
        _tied = {}.fromkeys(self.param_names, False)
        _bounds = {}.fromkeys(self.param_names, [-1.E12, 1.E12])
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
            param_dim = str(self.param_dim)
        except AttributeError:
            param_dim = " "
            
        if degree:
            fmt = "<{0}({1},".format(self.__class__.__name__, repr(self.deg))
        else:
            fmt = "<{0}(".format(self.__class__.__name__)
        for i in range(len(self.param_names)):
            fmt1 = """
            {0}={1},
            """.format(self.param_names[i], getattr(self, self.param_names[i]))
            fmt += fmt1.strip()
        if param_dim:
            fmt += "param_dim={0})>".format(self.param_dim)
        
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
              self.param_dim,
              "\n                   ".join(i+': ' + 
                str(self.__getattribute__(i)) for i in self.param_names)
                )
            
        return fmt
    
    @property
    def parameters(self):
        """
        An instance of `~astropy.models.parameters.Parameters`.
        Fittable parameters maintain this list and fitters modify it.
        """
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
            _val = parameters._tofloat(value)[0]
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

class LabeledInput(dict):
    """
    Create a container with all input data arrays, assigning labels for
    each one.
    
    Used by CompositeModel to choose input data using labels
    
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
        A Base class for all composite models.

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
    
    Execute models in series.
    
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
    
    Notes
    -----
    Output values of one model are used as input values of another.
    Obviously the order of the models matters.
    
    Examples
    --------
    Apply a 2D rotation followed by a shift in x and y
    
    >>> from astropy.models import *
    >>> rot = builtin_models.MatrixRotation2D(angle=23.5)
    >>> offx = builtin_models.ShiftModel(-4.23)
    >>> offy = builtin_models.ShiftModel(2)
    >>> linp = LabeledInput([x, y], ["x", "y"]
    >>> scomptr = SCompositeModel([rot, offx, offy], 
                                  inmap=[['x', 'y'], ['x'], ['y']],
                                  outmap=[['x', 'y'], ['x'], ['y']])
    >>> result=scomptr(linp)
        
    """
    def __init__(self, transforms, inmap=None, outmap=None):
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
    
    Execute models in parallel.
    
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
        
    Notes
    -----
    Models are applied to input data separately and the deltas are summed.
    
    """
    def __init__(self, transforms, inmap=None, outmap=None):
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
        