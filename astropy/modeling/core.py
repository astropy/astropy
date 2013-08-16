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
from textwrap import dedent

import numpy as np

from .parameters import Parameter, Parameters, InputParameterError, _tofloat
from ..utils import indent


__all__ = ['Model', 'ParametricModel', 'ParallelCompositeModel',
           'SerialCompositeModel', 'LabeledInput', 'Parametric1DModel',
           'Parametric2DModel', 'ModelDefinitionError']


class ModelDefinitionError(Exception):
    """Used for incorrect models definitions."""


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
        try:
            return x[0]
        except IndexError:
            return x
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
        val, _ = _tofloat(val)
        oldpar = getattr(obj, self.name)
        par = Parameter(self.name, val, obj, obj.param_dim)
        if self.name in obj._parcheck:
            obj._parcheck[self.name](val)
        if isinstance(obj, ParametricModel):
            if not obj._parameters._changed:
                par = Parameter(self.name, val, obj, obj.param_dim)
                if oldpar is not None and oldpar.shape != par.shape:
                    raise InputParameterError(
                        "Input parameter {0} does not "
                        "have the required shape".format(self.name))
                else:
                    setattr(oldpar, "value", val)
                obj._parameters = Parameters(obj, obj.param_names,
                                             dim=obj.param_dim)
            else:
                setattr(oldpar, "value", val)
        else:
            if oldpar is not None and oldpar.shape != par.shape:
                raise InputParameterError(
                    "Input parameter {0} does not "
                    "have the required shape".format(self.name))
            else:
                setattr(oldpar, "value", val)


class Model(object):
    """
    Base class for all models.

    This is an abstract class and should not be instantiated.

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
        self._param_names = param_names
        self.fittable = False
        #_parcheck is a dictionary to register parameter validation funcitons
        # key: value pairs are parameter_name: parameter_validation_function_name
        # see projections.AZP for example
        self._parcheck = {}
        for par in param_names:
            setattr(self.__class__, par, _ParameterProperty(par))

    @property
    def n_inputs(self):
        """Number of input variables in model evaluation."""

        return self._n_inputs

    @property
    def n_outputs(self):
        """Number of output variables returned when a model is evaluated."""

        return self._n_outputs

    @property
    def param_dim(self):
        """Number of parameter sets in a model."""

        return self._param_dim

    @param_dim.setter
    def param_dim(self, val):
        """Set the number of parameter sets in a model."""

        self._param_dim = val

    @property
    def param_names(self):
        """A list of names of the parameters defining a model."""

        return self._param_names

    @param_names.setter
    def param_names(self, val):
        self._param_names = val

    def __repr__(self):
        fmt = "{0}(".format(self.__class__.__name__)
        for name in self.param_names:
            fmt1 = """
            {0}={1},
            """.format(name, getattr(self, name))
            fmt += fmt1
        fmt += ")"

        return fmt

    def __str__(self):
        fmt = """
        Model: {0}
        Parameter sets: {1}
        Parameters: \n{2}
        """.format(
              self.__class__.__name__,
              self.param_dim,
              indent('\n'.join('{0}: {1}'.format(n, getattr(self, n))
                               for n in self.param_names),
                     width=19))

        return dedent(fmt[1:])

    @property
    def param_sets(self):
        """
        Return parameters as a pset.

        This is an array where each column represents one parameter set.
        """

        parameters = [getattr(self, '_' + attr) for attr in self.param_names]
        values = [par.value for par in parameters]
        shapes = [par.shape for par in parameters]
        lenshapes = np.asarray([len(p.shape) for p in parameters])

        if (lenshapes > 1).any():
            if () in shapes:
                psets = np.asarray(values, dtype=np.object)
            else:
                psets = np.asarray(values)
        else:
            psets = np.asarray(values)
            psets.shape = (len(self.param_names), self.param_dim)
        return psets

    def inverse(self):
        """Returns a callable object which performs the inverse transform."""

        raise NotImplementedError("An analytical inverse transform has not "
                                  "been implemented for this model.")

    def invert(self):
        """Invert coordinates iteratively if possible."""

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
            return ParallelCompositeModel([self, newtr])
        elif mode in ['serial', 's']:
            return SerialCompositeModel([self, newtr])
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
    All models which can be fit to data should subclass this class.

    Sets the ``parameters`` attributes.

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
        parameters. Keys are parameter names. Values are a list of length 2
        giving the desired range for the parameter.  Alternatively the
        `~astropy.modeling.parameters.Parameter.min` and
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

    Examples
    --------
    >>> from astropy.modeling import models
    >>> def tie_center(model):
    ...         mean = 50 * model.stddev
    ...         return mean
    >>> tied_parameters = {'mean': tie_center}

    Specify that 'mean' is a tied parameter in one of two ways:

    >>> g1 = models.Gaussian1DModel(amplitude=10, mean=5, stddev=.3,
    ...                             tied=tied_parameters)

    or

    >>> g1 = models.Gaussian1DModel(amplitude=10, mean=5, stddev=.3)
    >>> g1.mean.tied
    False
    >>> g1.mean.tied = tie_center
    >>> g1.mean.tied
    <function tie_center at 0x...>

    Fixed parameters:

    >>> g1 = models.Gaussian1DModel(amplitude=10, mean=5, stddev=.3,
    ...                             fixed={'stddev': True})
    >>> g1.stddev.fixed
    True

    or

    >>> g1 = models.Gaussian1DModel(amplitude=10, mean=5, stddev=.3)
    >>> g1.stddev.fixed
    False
    >>> g1.stddev.fixed = True
    >>> g1.stddev.fixed
    True
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, param_names, n_inputs, n_outputs, param_dim=1, **cons):
        self.linear = False
        bounds = cons.pop('bounds', None)
        fixed = cons.pop('fixed', None)
        tied = cons.pop('tied', None)
        eqcons = cons.pop('eqcons', None)
        ineqcons = cons.pop('ineqcons', None)

        if cons:
            raise TypeError(
                "Unrecognized constraints type: {0}".format(repr(cons.keys())))

        super(ParametricModel, self).__init__(param_names, n_inputs, n_outputs,
                                              param_dim=param_dim)
        self.fittable = True
        self._parameters = Parameters(self, self.param_names, dim=param_dim)
        # Initialize the constraints for each parameter
        _bounds = {}.fromkeys(self.param_names, [None, None])
        if eqcons is None:
            self._eqcons = []
        else:
            self._eqcons = eqcons
        if ineqcons is None:
            self._ineqcons = []
        else:
            self._ineqcons = ineqcons

        # Flag that indicates if the model derivatives are given in columns
        # or rows
        self.col_deriv = 1

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
            _bounds.update(bounds)

        for name in _bounds:
            par = getattr(self, name)
            setattr(par, 'min', _bounds[name][0])
            setattr(par, 'max', _bounds[name][1])

    @property
    def fixed(self):
        """
        A dictionary mapping parameter names to their fixed constraint.
        """

        fixed = [getattr(self, name).fixed for name in self.param_names]
        return dict(zip(self.param_names, fixed))

    @property
    def tied(self):
        """
        A dictionary mapping parameter names to their tied constraint.
        """

        tied = [getattr(self, name).tied for name in self.param_names]
        return dict(zip(self.param_names, tied))

    @property
    def bounds(self):
        pars = [getattr(self, name) for name in self.param_names]
        bounds = [(par.min, par.max) for par in pars]
        return dict(zip(self.param_names, bounds))

    @property
    def eqcons(self):
        return self._eqcons

    @property
    def ineqcons(self):
        return self._ineqcons

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
        n_inputs:   {1}
        Degree: {2}
        Parameter sets: {3}
        Parameters: \n{4}
        """.format(
              self.__class__.__name__,
              self.n_inputs,
              degree,
              self.param_dim,
              indent('\n'.join('{0}: {1}'.format(n, getattr(self, n))
                     for n in self.param_names),
                     width=19))

        return dedent(fmt[1:])

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
        Resets the parameters attribute using an instance of
        `~astropy.modeling.parameters.Parameters`
        """

        if isinstance(value, Parameters):
            if self._parameters._is_same_length(value):
                self._parameters = value
            else:
                raise InputParameterError(
                    "Expected the list of parameters to be the same "
                    "length as the initial list.")
        elif isinstance(value, (list, np.ndarray)):
            _val = _tofloat(value)[0]
            if self._parameters._is_same_length(_val):
                self._parameters._changed = True
                self._parameters[:] = _val
            else:
                raise InputParameterError(
                    "Expected the list of parameters to be the same "
                    "length as the initial list.")
        else:
            raise TypeError("Parameters must be of type `list` or `Parameters`")

    def set_joint_parameters(self, jparams):
        """
        Used by the `JointFitter` class to store parameters which are
        considered common for several models and are to be fitted together.
        """
        self.joint = jparams


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
    >>> x,y = np.mgrid[:5, :5]
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

    def __repr__(self):
        fmt = """
            Model:  {0}
            """.format(self.__class__.__name__)
        fmt_args = tuple(repr(tr) for tr in self._transforms)
        fmt1 = (" %s  " * len(self._transforms)) % fmt_args
        fmt = fmt + fmt1
        return fmt

    def __str__(self):
        fmt = """
            Model:  {0}
            """.format(self.__class__.__name__)
        fmt_args = tuple(str(tr) for tr in self._transforms)
        fmt1 = (" %s  " * len(self._transforms)) % fmt_args
        fmt = fmt + fmt1
        return fmt

    def invert(self):
        raise NotImplementedError("Subclasses should implement this")

    def __call__(self):
        # implemented by subclasses
        raise NotImplementedError("Subclasses should implement this")


class SerialCompositeModel(_CompositeModel):
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

        >>> import numpy as np
        >>> from astropy.modeling import models, LabeledInput, SerialCompositeModel
        >>> x, y = np.mgrid[:5, :5]
        >>> rotation = models.MatrixRotation2D(angle=23.5)
        >>> offset_x = models.ShiftModel(-4.23)
        >>> offset_y = models.ShiftModel(2)
        >>> labeled_input = LabeledInput([x, y], ["x", "y"])
        >>> transform = SerialCompositeModel([rotation, offset_x, offset_y],
        ...                                  inmap=[['x', 'y'], ['x'], ['y']],
        ...                                  outmap=[['x', 'y'], ['x'], ['y']])
        >>> result = transform(labeled_input)
    """

    def __init__(self, transforms, inmap=None, outmap=None, n_inputs=None,
                 n_outputs=None):
        if n_inputs is None:
            n_inputs = max([tr.n_inputs for tr in transforms])
            # the output dimension is equal to the output dim of the last transform
            n_outputs = transforms[-1].n_outputs
        else:
            assert n_outputs is not None, "Expected n_inputs and n_outputs"
            n_inputs = n_inputs
            n_outputs = n_outputs

        super(SerialCompositeModel, self).__init__(transforms, n_inputs,
                                                   n_outputs)

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

    def inverse(self):
        try:
            transforms = [tr.inverse() for tr in self._transforms[::-1]]
        except NotImplementedError:
            raise
        if self._inmap is not None:
            inmap = self._inmap[::-1]
            outmap = self._outmap[::-1]
        else:
            inmap = None
            outmap = None
        return SerialCompositeModel(transforms, inmap, outmap)

    def __call__(self, *data):
        """Transforms data using this model."""

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
                # we want to return the entire labeled object because some
                # parts of it may be used in another transform of which this
                # one is a component
                assert self._inmap is not None, \
                    ("Parameter 'inmap' must be provided when "
                     "input is a labeled object.")
                assert self._outmap is not None, \
                    ("Parameter 'outmap' must be provided when input is a "
                     "labeled object")

                for transform, incoo, outcoo in izip(self._transforms,
                                                     self._inmap,
                                                     self._outmap):
                    inlist = [labeled_input[label] for label in incoo]
                    result = transform(*inlist)
                    if len(outcoo) == 1:
                        result = [result]
                    for label, res in zip(outcoo, result):

                        if label not in labeled_input.labels:
                            labeled_input[label] = res
                        setattr(labeled_input, label, res)
                return labeled_input
        else:
            assert self.n_inputs == len(data), \
                "This transform expects {0} inputs".format(self._n_inputs)

            result = self._transforms[0](*data)
            for transform in self._transforms[1:]:
                result = transform(*result)
        return result


class ParallelCompositeModel(_CompositeModel):
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
                ("A ParallelCompositeModel expects n_inputs = n_outputs for "
                 "all transforms")

        super(ParallelCompositeModel, self).__init__(transforms, n_inputs,
                                                     n_outputs)

        self._inmap = inmap
        self._outmap = outmap

    def __call__(self, *data):
        """Transforms data using this model."""

        if len(data) == 1:
            if not isinstance(data[0], LabeledInput):
                result = data[0]
                x = data[0]
                deltas = sum(tr(x) for tr in self._transforms)
                return result + deltas
            else:
                assert self._inmap is not None, \
                    ("Parameter 'inmap' must be provided when "
                     "input is a labeled object.")
                assert self._outmap is not None, \
                    ("Parameter 'outmap' must be provided when input is a "
                     "labeled object")
                labeled_input = data[0].copy()
                # create a list of inputs to be passed to the transforms
                inlist = [getattr(labeled_input, label)
                          for label in self._inmap]
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


class Parametric1DModel(ParametricModel):
    """
    Base class for one dimensional parametric models

    This class provides an easier interface to defining new models.
    Examples can be found in functional_models.py

    Parameters
    ----------
    parameter_dict : dictionary
        Dictionary of model parameters with initialisation values
        {'parameter_name': 'parameter_value'}
    """

    deriv = None

    def __init__(self, param_dict):
        # Get parameter dimension
        param_dim = np.size(param_dict[self.param_names[0]])
        cons = param_dict.pop('constraints', {})
        # Initialize model parameters. This is preliminary as long there is
        # no new parameter class. It may be more reasonable and clear to init
        # the parameters in the model constructor itself, with constraints etc.
        for param_name in self.param_names:
            setattr(self, "_" + param_name,
                    Parameter(name=param_name, val=param_dict[param_name],
                              model=self, dim=param_dim))

        super(Parametric1DModel, self).__init__(self.param_names, n_inputs=1,
                                                n_outputs=1,
                                                param_dim=param_dim, **cons)

    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        x, fmt = _convert_input(x, self.param_dim)
        result = self.eval(x, *self.param_sets)
        return _convert_output(result, fmt)


class Parametric2DModel(ParametricModel):
    """
    Base class for two dimensional parametric models

    This class provides an easier interface to defining new models.
    Examples can be found in functional_models.py

    Parameters
    ----------
    parameter_dict : dictionary
        Dictionary of model parameters with initialization values
        {'parameter_name': 'parameter_value'}
    """

    deriv = None

    def __init__(self, param_dict):
        # Get parameter dimension
        param_dim = np.size(param_dict[self.param_names[0]])
        cons = param_dict.pop('constraints', {})
        # Initialize model parameters. This is preliminary as long there is
        # no new parameter class. It may be more reasonable and clear to init
        # the parameters in the model constructor itself, with constraints etc.
        for param_name in self.param_names:
            setattr(self, "_" + param_name,
                    Parameter(name=param_name, val=param_dict[param_name],
                              model=self, dim=param_dim))

        super(Parametric2DModel, self).__init__(self.param_names, n_inputs=2,
                                                n_outputs=1,
                                                param_dim=param_dim, **cons)

    def __call__(self, x, y):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        x, _ = _convert_input(x, self.param_dim)
        y, fmt = _convert_input(y, self.param_dim)
        result = self.eval(x, y, *self.param_sets)
        return _convert_output(result, fmt)
