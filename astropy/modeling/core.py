# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines base classes for all models.
The base class of all models is `~astropy.modeling.Model`.
`~astropy.modeling.ParametricModel` is the base class for all fittable models. Parametric
models can be linear or nonlinear in a regression analysis sense.

All models provide a `__call__` method which performs the transformation in a
purely mathematical way, i.e. the models are unitless. In addition, when
possible the transformation is done using multiple parameter sets, `param_sets`.
The number of parameter sets is stored in an attribute `param_dim`.

Parametric models also store a flat list of all parameters as an instance of
`~astropy.modeling.parameters.Parameters`. When fitting, this list-like object is modified by a
subclass of `~astropy.modeling.fitting.Fitter`. When fitting nonlinear models, the values of the
parameters are used as initial guesses by the fitting class. Normally users
will not have to use the `~astropy.modeling.parameters` module directly.

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
from __future__ import division
import abc
from itertools import izip
import numpy as np
from . import parameters
from . import constraints
from .utils import InputParameterError

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
        self.aname = '_' + name
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

    def __init__(self, param_names, n_inputs, n_outputs, param_dim=1):
        self._param_dim = param_dim
        self._n_inputs = n_inputs
        self._n_outputs = n_outputs
        self.has_inverse = False
        self._param_names = param_names
        #_parcheck is a dictionary to register parameter validation funcitons
        # key: value pairs are parameter_name: parameter_validation_function_name
        # see projections.AZP for example
        self._parcheck = {}
        for par in param_names:
            setattr(self.__class__, par, _ParameterProperty(par))

    @property
    def n_inputs(self):
        """
        Number of input variables in model evaluation.
        """
        return self._n_inputs

    @property
    def n_outputs(self):
        """
        Number of output variables returned when a model is evaluated.
        """
        return self._n_outputs

    @property
    def param_dim(self):
        """
        Number of parameter sets in a model.
        """
        return self._param_dim

    @param_dim.setter
    def param_dim(self, val):
        """
        Set the number of parameter sets in a model.
        """
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
              "\n                   ".join(i + ': ' +
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
        if (lenshapes > 1).any():
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
    n_inputs: int
        number of inputs
    n_outputs: int
        number of output quantities
    param_dim: int
        number of parameter sets
    fittable: boolean
        indicator if the model is fittable
    fixed: a dict
        a dictionary {parameter_name: boolean} of parameters to not be
        varied during fitting. True means the parameter is held fixed.
        Alternatively the `~astropy.modeling.parameters.Parameter.fixed`
        property of a parameter may be used.
    tied: dict
        a dictionary {parameter_name: callable} of parameters which are
        linked to some other parameter. The dictionary values are callables
        providing the linking relationship.
        Alternatively the `~astropy.modeling.parameters.Parameter.tied`
        property of a parameter may be used.
    bounds: dict
        a dictionary {parameter_name: boolean} of lower and upper bounds of
        parameters. Keys  are parameter names. Values  are a list of length
        2 giving the desired range for the parameter.
        Alternatively the `~astropy.modeling.parameters.Parameter.min` and
        `~astropy.modeling.parameters.Parameter.max` properties of a parameter
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

    def __init__(self, param_names, n_inputs, n_outputs, param_dim=1, fittable=True,
                 fixed=None, tied=None, bounds=None, eqcons=None, ineqcons=None):
        self.linear = True
        super(ParametricModel, self).__init__(param_names, n_inputs, n_outputs,
                                              param_dim=param_dim)
        self.fittable = fittable
        self._parameters = parameters.Parameters(self, self.param_names,
                                                 param_dim=param_dim)
        # Initialize the constraints for each parameter
        _fixed = {}.fromkeys(self.param_names, False)
        _tied = {}.fromkeys(self.param_names, False)
        _bounds = {}.fromkeys(self.param_names, [-1.E12, 1.E12])
        if eqcons is None:
            eqcons = []
        if ineqcons is None:
            ineqcons = []
        self.constraints = constraints.Constraints(self, fixed=_fixed,
                                                   tied=_tied,
                                                   bounds=_bounds,
                                                   eqcons=eqcons,
                                                   ineqcons=ineqcons)
        # Set constraints
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
                setattr(par, 'min', bounds[name][0])
                setattr(par, 'max', bounds[name][1])

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
              self.n_inputs,
              degree,
              self.param_dim,
              "\n                   ".join(i + ': ' +
                                           str(self.__getattribute__(i)) for i in self.param_names)
        )

        return fmt

    @property
    def parameters(self):
        """
        An instance of `~astropy.modeling.parameters.Parameters`.
        Fittable parameters maintain this list and fitters modify it.
        """
        return self._parameters

    @parameters.setter
    def parameters(self, value):
        """
        Reset the parameters attribute as an instance of
        `~astropy.modeling.parameters.Parameters`
        """
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
    def __init__(self, data, labels):
        dict.__init__(self)
        assert len(labels) == len(data)
        self.labels = [l.strip() for l in labels]
        for coord, label in zip(data, labels):
            self[label] = coord
            setattr(self, '_' + label, coord)
        self._set_properties(self.labels)

    def _getlabel(self, name):
        par = getattr(self, '_' + name)
        return par

    def _setlabel(self, name, val):
        setattr(self, '_' + name, val)
        self[name] = val

    def _dellabel(self, name):
        delattr(self, '_' + name)
        del self[name]

    def add(self, label=None, value=None, **kw):
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
            self.__setattr__('_' + key, kw[key])
        self._set_properties(kw.keys())

    def _set_properties(self, attributes):
        for attr in attributes:
            setattr(self.__class__, attr, property(lambda self, attr=attr:
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


class _CompositeModel(Model):

    def __init__(self, transforms, n_inputs, n_outputs):
        """
        A Base class for all composite models.

        """
        self._transforms = transforms
        param_names = []
        for tr in self._transforms:
            param_names.extend(tr.param_names)
        super(_CompositeModel, self).__init__(param_names, n_inputs, n_outputs)
        self.fittable = False
        self.has_inverse = all([tr.has_inverse for tr in transforms])

    def __repr__(self):
        fmt = """
            Model:  {0}
            """.format(self.__class__.__name__)
        fmt1 = " %s  " * len(self._transforms) % tuple([repr(tr) for tr in self._transforms])
        fmt = fmt + fmt1
        return fmt

    def __str__(self):
        fmt = """
            Model:  {0}
            """.format(self.__class__.__name__)
        fmt1 = " %s  " * len(self._transforms) % tuple([str(tr) for tr in self._transforms])
        fmt = fmt + fmt1
        return fmt

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

    Notes
    -----
    Output values of one model are used as input values of another.
    Obviously the order of the models matters.

    Examples
    --------
    Apply a 2D rotation followed by a shift in x and y::

        >>> from astropy.modeling import *
        >>> rot = models.MatrixRotation2D(angle=23.5)
        >>> offx = models.ShiftModel(-4.23)
        >>> offy = models.ShiftModel(2)
        >>> linp = LabeledInput([x, y], ["x", "y"])
        >>> scomptr = SCompositeModel([rot, offx, offy],
        ...                           inmap=[['x', 'y'], ['x'], ['y']],
        ...                           outmap=[['x', 'y'], ['x'], ['y']])
        >>> result=scomptr(linp)

    """
    def __init__(self, transforms, inmap=None, outmap=None, n_inputs=None, n_outputs=None):
        if n_inputs is None:
            n_inputs = max([tr.n_inputs for tr in transforms])
            # the output dimension is equal to the output dim of the last transform
            n_outputs = transforms[-1].n_outputs
        else:
            assert n_outputs is not None, "Expected n_inputs and n_outputs"
            n_inputs = n_inputs
            n_outputs = n_outputs
        super(SCompositeModel, self).__init__(transforms, n_inputs, n_outputs)
        if transforms and inmap and outmap:
            assert len(transforms) == len(inmap) == len(outmap), \
                "Expected sequences of transform, " \
                "inmap and outmap to have the same length"
        if inmap is None:
            inmap = [None] * len(transforms)
        if outmap is None:
            outmap = [None] * len(transforms)
        self._inmap = inmap
        self._outmap = outmap

    def __call__(self, *data):
        """
        Transforms data using this model.
        """
        if len(data) == 1:
            if not isinstance(data[0], LabeledInput):
                assert self._transforms[0].n_inputs == 1, \
                    "First transform expects {0} inputs, 1 given".format(
                        self._transforms[0].n_inputs)

                result = data[0]
                for tr in self._transforms:
                    result = tr(result)
                return result
            else:
                labeled_input = data[0].copy()
                # we want to return the entire labeled object because some parts
                # of it may be used in another transform of which this
                # one is a component
                assert self._inmap is not None, ("Parameter 'inmap' must be provided when"
                                                 "input is a labeled object")
                assert self._outmap is not None, ("Parameter 'outmap' must be "
                                                  "provided when input is a labeled object")
                for i in range(len(self._transforms)):
                    inlist = [labeled_input[label] for label in self._inmap[i]]
                    result = [self._transforms[i](*inlist)]
                    output = self._outmap[i]
                    for label, res in zip(output, result):
                        labeled_input.update({label: res})
                        if label not in self._inmap[i]:
                            setattr(labeled_input, label, res)
                return labeled_input
        else:
            assert self.n_inputs == len(data), "This transform expects "
            "{0} inputs".format(self._n_inputs)

            result = self._transforms[0](*data)
            for transform in self._transforms[1:]:
                result = transform(result)
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
    outmap : list or None

    Notes
    -----
    Evaluate each model separately and add the results to the input_data.

    """
    def __init__(self, transforms, inmap=None, outmap=None):
        self._transforms = transforms
        n_inputs = self._transforms[0].n_inputs
        n_outputs = n_inputs
        for transform in self._transforms:
            assert transform.n_inputs == transform.n_outputs == n_inputs, \
                ("A PCompositeModel expects n_inputs = n_outputs for all transforms")
        super(PCompositeModel, self).__init__(transforms, n_inputs, n_outputs)

        self._inmap = inmap
        self._outmap = outmap

    def __call__(self, *data):
        """
        Transforms data using this model.
        """
        if len(data) == 1:
            if not isinstance(data[0], LabeledInput):
                result = data[0]
                x = data[0]
                deltas = sum(tr(x) for tr in self._transforms)
                return result + deltas
            else:
                assert self._inmap is not None, ("Parameter 'inmap' must be "
                                                 "provided when input is a labeled object")
                assert self._outmap is not None, ("Parameter 'outmap' must be "
                                                  "provided when input is a labeled object")
                labeled_input = data[0].copy()
                # create a list of inputs to be passed to the transforms
                inlist = [getattr(labeled_input, label) for label in self._inmap]
                deltas = [np.zeros_like(x) for x in inlist]
                for transform in self._transforms:
                    deltas = [transform(*inlist)]
                for outcoo, inp, delta in izip(self._outmap, inlist, deltas):
                    setattr(labeled_input, outcoo, inp + delta)
                # always return the entire labeled object, not just the result
                # since this may be part of another composite transform
                return labeled_input
        else:
            result = data[:]
            assert self.n_inputs == self.n_outputs
            for tr in self._transforms:
                result += tr(*data)
            return result
