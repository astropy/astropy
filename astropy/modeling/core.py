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
import functools
import copy
from itertools import izip
from textwrap import dedent

import numpy as np

from .parameters import Parameter, InputParameterError
from ..utils import indent, isiterable


__all__ = ['Model', 'ParametricModel', 'SummedCompositeModel',
           'SerialCompositeModel', 'LabeledInput', 'Parametric1DModel',
           'Parametric2DModel', 'ModelDefinitionError', 'format_input']


class ModelDefinitionError(Exception):
    """Used for incorrect models definitions"""


def format_input(func):
    """
    Wraps a model's ``__call__`` method so that the input arrays are converted
    into the appropriate shape given the model's parameter dimensions.

    Wraps the result to match the shape of the last input array.
    """

    @functools.wraps(func)
    def wrapped_call(self, *args):
        converted = []

        for arg in args:
            # Reset these flags; their value only matters for the last
            # argument
            transposed = False
            scalar = False

            arg = np.asarray(arg) + 0.
            if self.param_dim == 1:
                if arg.ndim == 0:
                    scalar = True
                converted.append(arg)
                continue

            if arg.ndim < 2:
                converted.append(np.array([arg]).T)
            elif arg.ndim == 2:
                assert arg.shape[-1] == self.param_dim, \
                    ("Cannot broadcast with shape "
                     "({0}, {1})".format(arg.shape[0], arg.shape[1]))
                converted.append(arg)
            elif arg.ndim > 2:
                assert arg.shape[0] == self.param_dim, \
                    ("Cannot broadcast with shape "
                     "({0}, {1}, {2})".format(arg.shape[0], arg.shape[1],
                                              arg.shape[2]))
                transposed = True
                converted.append(arg.T)

        result = func(self, *converted)

        if transposed:
            return result.T
        elif scalar:
            try:
                return result[0]
            except IndexError:
                return result

        return result

    return wrapped_call


class _ModelMeta(abc.ABCMeta):
    """
    Metaclass for Model.

    Currently just handles auto-generating the param_names list based on
    Parameter descriptors declared at the class-level of Model subclasses.
    """

    def __new__(mcls, name, bases, members):
        param_names = members.get('param_names', [])
        parameters = dict((value.name, value) for value in members.values()
                          if isinstance(value, Parameter))

        # If no parameters were defined get out early--this is especially
        # important for PolynomialModels which take a different approach to
        # parameters, since they can have a variable number of them
        if not parameters:
            return super(_ModelMeta, mcls).__new__(mcls, name, bases, members)

        # If param_names was declared explicitly we use only the parameters
        # listed manually in param_names, but still check that all listed
        # parameters were declared
        if param_names and isiterable(param_names):
            for param_name in param_names:
                if param_name not in parameters:
                    raise RuntimeError(
                        "Parameter {0!r} listed in {1}.param_names was not "
                        "declared in the class body.".format(param_name, name))
        else:
            param_names = [param.name for param in
                           sorted(parameters.values(),
                                  key=lambda p: p._order)]
            members['param_names'] = param_names

        return super(_ModelMeta, mcls).__new__(mcls, name, bases, members)


class Model(object):
    """
    Base class for all models.

    This is an abstract class and should not be instantiated directly.

    Notes
    -----
    Models which are not meant to be fit to data should subclass this class.

    This class sets the properties for all individual parameters and performs
    parameter validation.
    """

    __metaclass__ = _ModelMeta

    param_names = []
    n_inputs = 1
    n_outputs = 1
    fittable = False
    linear = True

    def __init__(self, param_dim=1):
        self._param_dim = param_dim

    @property
    def param_dim(self):
        """Number of parameter sets in a model."""

        return self._param_dim

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

        parameters = [getattr(self, attr) for attr in self.param_names]
        values = [par.value for par in parameters]
        shapes = [par.shape for par in parameters]
        n_dims = np.asarray([len(p.shape) for p in parameters])

        if (n_dims > 1).any():
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
            return SummedCompositeModel([self, newtr])
        elif mode in ['serial', 's']:
            return SerialCompositeModel([self, newtr])
        else:
            raise InputParameterError("Unrecognized mode {0}".format(mode))

    def copy(self):
        return copy.deepcopy(self)

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
    param_dim : int
        number of parameter sets
    fittable : bool
        indicator if the model is fittable
    fixed : a dict
        a dictionary ``{parameter_name: boolean}`` of parameters to not be
        varied during fitting. True means the parameter is held fixed.
        Alternatively the `~astropy.modeling.parameters.Parameter.fixed`
        property of a parameter may be used.
    tied : dict
        a dictionary ``{parameter_name: callable}`` of parameters which are
        linked to some other parameter. The dictionary values are callables
        providing the linking relationship.
        Alternatively the `~astropy.modeling.parameters.Parameter.tied`
        property of a parameter may be used.
    bounds : dict
        a dictionary ``{parameter_name: boolean}`` of lower and upper bounds of
        parameters. Keys are parameter names. Values are a list of length 2
        giving the desired range for the parameter.  Alternatively the
        `~astropy.modeling.parameters.Parameter.min` and
        `~astropy.modeling.parameters.Parameter.max` properties of a parameter
        may be used.
    eqcons : list
        A list of functions of length n such that
        ``eqcons[j](x0, *args) == 0.0`` in a successfully optimized
        problem.
    ineqcons : list
        A list of functions of length n such that
        ``ieqcons[j](x0, *args) >= 0.0`` is a successfully optimized
        problem.

    Examples
    --------
    >>> from astropy.modeling import models
    >>> def tie_center(model):
    ...         mean = 50 * model.stddev
    ...         return mean
    >>> tied_parameters = {'mean': tie_center}

    Specify that ``'mean'`` is a tied parameter in one of two ways:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                             tied=tied_parameters)

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.mean.tied
    False
    >>> g1.mean.tied = tie_center
    >>> g1.mean.tied
    <function tie_center at 0x...>

    Fixed parameters:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                             fixed={'stddev': True})
    >>> g1.stddev.fixed
    True

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.stddev.fixed
    False
    >>> g1.stddev.fixed = True
    >>> g1.stddev.fixed
    True
    """

    linear = False
    deriv = None
    # Flag that indicates if the model derivatives are given in columns
    # or rows
    col_deriv = True
    fittable = True


    def __init__(self, **kwargs):
        # Pop any constraints off the keyword arguments
        bounds = kwargs.pop('bounds', None)
        fixed = kwargs.pop('fixed', None)
        tied = kwargs.pop('tied', None)
        eqcons = kwargs.pop('eqcons', None)
        ineqcons = kwargs.pop('ineqcons', None)

        # Pop off the param_dims
        param_dim = kwargs.pop('param_dim', None)

        # Remaining keyword args are either parameter values or invalid
        # Parameter values must be passed in as keyword arguments in order to
        # distinguish them
        params = kwargs

        # Determine the number of parameter sets: This will be based
        # on the size of any parameters whose values have been specified
        # or the default of 1 is used
        # This loop also checks that all the supplied parameter names are
        # valid
        param_names = set(self.param_names)

        max_param_dim = 1
        for name, value in params.items():
            if param_dim is None:
                # Determine the best param_dims, if not already specified,
                # based on the sizes of the input parameter values
                max_param_dim = max(max_param_dim, np.size(value))

            if name not in param_names:
                raise ValueError(
                    "Unrecognized parameter: {0}".format(name))

        if param_dim is None:
            param_dim = max_param_dim

        super(ParametricModel, self).__init__(param_dim=param_dim)

        # Initialize the constraints for each parameter
        if eqcons is None:
            self._eqcons = []
        else:
            self._eqcons = eqcons
        if ineqcons is None:
            self._ineqcons = []
        else:
            self._ineqcons = ineqcons

        # Set constraints
        self._constraints = {}

        if fixed:
            self._constraints['fixed'] = fixed
        if tied:
            self._constraints['tied'] = tied
        if bounds:
            self._constraints['bounds'] = bounds

        self._initialize_parameters(params)


    @property
    def fixed(self):
        """
        A dictionary mapping parameter names to their fixed constraint
        """

        return dict((name, getattr(self, name).fixed)
                    for name in self.param_names)

    @property
    def tied(self):
        """
        A dictionary mapping parameter names to their tied constraint
        """

        return dict((name, getattr(self, name).tied)
                    for name in self.param_names)

    @property
    def bounds(self):
        return dict((name, getattr(self, name).bounds)
                    for name in self.param_names)

    @property
    def eqcons(self):
        return self._eqcons

    @property
    def ineqcons(self):
        return self._ineqcons

    @property
    def parameters(self):
        """
        A flattened array of all parameter values in all parameter sets

        Fittable parameters maintain this list and fitters modify it.
        """

        return self._parameters

    @parameters.setter
    def parameters(self, value):
        """
        Assigning to this attribute updates the parameters array rather than
        replacing it.
        """

        try:
            value = np.array(value).reshape(self._parameters.shape)
        except ValueError as e:
            raise InputParameterError(
                "Input parameter values not compatible with the model "
                "parameters array: {0}".format(e))

        self._parameters[:] = value

    def __getattr__(self, attr):
        if len(attr) > 1 and attr[0] == '_' and attr != '_param_metrics':
            param_name = attr[1:]
            if param_name in self._param_metrics:
                param_slice, param_shape = self._param_metrics[param_name]
                value = self._parameters[param_slice]
                if self.param_dim == 1:
                    return value[0]
                else:
                    if param_shape:
                        return value.reshape(param_shape)
                    else:
                        return value[0]

        raise AttributeError(attr)

    def __setattr__(self, attr, value):
        if (len(attr) > 1 and attr[0] == '_' and
                hasattr(self, '_param_metrics')):
            param_name = attr[1:]
            if param_name in self._param_metrics:
                # TODO: Maybe handle exception on invalid input shape
                param_slice = self._param_metrics[param_name][0]
                self._parameters[param_slice] = np.array(value).ravel()
                return

        super(ParametricModel, self).__setattr__(attr, value)

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
            fmt = "<{0}({1}, ".format(self.__class__.__name__, repr(self.deg))
        else:
            fmt = "<{0}(".format(self.__class__.__name__)

        for name in self.param_names:
            fmt1 = "{0}={1}, ".format(name, getattr(self, name))
            fmt += fmt1
        if param_dim:
            fmt += "param_dim={0})>".format(self.param_dim)

        return fmt

    def __str__(self):
        try:
            degree = str(self.degree)
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

    def set_joint_parameters(self, jparams):
        """
        Used by the `JointFitter` class to store parameters which are
        considered common for several models and are to be fitted together.
        """
        self.joint = jparams

    def _model_to_fit_params(self):
        """
        Create a set of parameters to be fitted.

        These may be a subset of the model parameters, if some of them are held
        constant or tied.
        """
        fitparam_indices = range(len(self.param_names))
        if any(self.fixed.values()) or any(self.tied.values()):
            params = list(self.parameters)
            for idx, name in list(enumerate(self.param_names))[::-1]:
                if self.fixed[name] or self.tied[name]:
                    sl = self._param_metrics[name][0]
                    del params[sl]
                    del fitparam_indices[idx]
            return (np.array(params), fitparam_indices)
        else:
            return (self.parameters, fitparam_indices)

    def _fitter_to_model_params(self, fps):
        _fit_params, _fit_param_indices = self._model_to_fit_params()
        if any(self.fixed.values()) or any(self.tied.values()):
            self.parameters[_fit_param_indices] = fps
            for idx, name in enumerate(self.param_names):
                if self.tied[name] != False:
                    value = self.tied[name](self)
                    slice_ = self._param_metrics[name][0]
                    self.parameters[slice_] = value
        elif any([tuple(b) != (None, None) for b in self.bounds.values()]):
            for name, par in zip(self.param_names, _fit_params):
                if self.bounds[name] != (None, None):
                    b = self.bounds[name]
                    if b[0] is not None:
                        par = max(par, self.bounds[name][0])
                    if b[1] is not None:
                        par = min(par, self.bounds[name][1])
                    setattr(self, name, par)
        else:
            self.parameters = fps

    def _initialize_parameters(self, params):
        """
        Initialize the _parameters array that stores raw parameter values for
        all parameter sets for use with vectorized fitting algorithms; on
        ParametricModels the _param_name attributes actually just reference
        slices of this array.
        """

        # First we need to determine how much array space is needed by all the
        # parameters based on the number of parameters, the shape each input
        # parameter, and the param_dim
        self._param_metrics = {}
        total_size = 0
        for name in self.param_names:

            if params.get(name) is None:
                # parameters that were not supplied at all or that have
                # defaults of None should attempt to use the default provided
                # by their Parameter descriptor
                params[name] = getattr(self, name).default

            value = params[name]

            param_size = np.size(value)
            param_shape = np.shape(value)

            if self.param_dim == 1:
                pass
            elif self.param_dim > 1:
                if param_size == 1:
                    param_size = self.param_dim
                    # Update the value for this param to the new repeated
                    # version
                    value = params[name] = np.repeat(value, param_size)
                    param_shape = value.shape
            else:
                raise ValueError("Model param_dim must be 1 or greater.")

            if (param_size > 1 and param_size != self.param_dim and
                    len(value) != self.param_dim):
                # The len(value) == self.param_dim case is a special case
                # (see #1680) where the parameter has compound values (like [1,
                # 2]) but we're passing in two (or more) param sets, so we want
                # to make sure a value like [[1, 2], [3, 4]] is interpreted as
                # having values for 2 param sets, not 4.

                # For now all parameter values must have a number of elements
                # equal to param_dim (the number of param sets) or it is
                # invalid.  The *only* exception is, again, a scalar value is
                # allowed to be repeated across param sets
                raise InputParameterError(
                    "The input value for {0!r} has too many elements for "
                    "param_dim={1}.  Parameter values must either match the "
                    "model's param_dim in size, or may be a scalar value, in "
                    "which case the value is repeated accross parameter "
                    "sets.".format(name, self.param_dim))

            param_slice = slice(total_size, total_size + param_size)
            self._param_metrics[name] = (param_slice, param_shape)
            total_size += param_size

        self._parameters = np.empty(total_size, dtype=np.float64)

        # Now set the parameter values (this will also fill
        # self._parameters)
        for name, value in params.items():
            setattr(self, name, value)


    def _model_to_fit_params(self):
        """
        Create a set of parameters to be fitted.

        These may be a subset of the model parameters, if some of them are held
        constant or tied.

        This is an adapter of the model parameters to fitters, none of which
        currently support parameters with tied/fixed constraints.

        TODO: It may make more sense to make this adaptation on the Fitter end,
        since one could conceivably implement a fitter from scratch (for
        Astropy specifically) that understands how to use tied constraints, for
        example.
        """

        fitparam_indices = range(len(self.param_names))
        if any(self.fixed.values()) or any(self.tied.values()):
            params = list(self.parameters)
            for idx, name in list(enumerate(self.param_names))[::-1]:
                if self.fixed[name] or self.tied[name]:
                    sl = self._param_metrics[name][0]
                    del params[sl]
                    del fitparam_indices[idx]
            return (np.array(params), fitparam_indices)
        else:
            return (self.parameters, fitparam_indices)


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
        label : str
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
        """Base class for all composite models."""

        self._transforms = transforms
        param_names = []
        for tr in self._transforms:
            param_names.extend(tr.param_names)
        super(_CompositeModel, self).__init__()
        self.param_names = param_names
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs
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

    def add_model(self, transf, inmap, outmap):
        self[transf] = [inmap, outmap]

    def invert(self):
        raise NotImplementedError("Subclasses should implement this")

    def __call__(self):
        # implemented by subclasses
        raise NotImplementedError("Subclasses should implement this")


class SerialCompositeModel(_CompositeModel):
    """
    Composite model that evaluates models in series.

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
        >>> offset_x = models.Shift(-4.23)
        >>> offset_y = models.Shift(2)
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
            # the output dimension is equal to the output dim of the last
            # transform
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


class SummedCompositeModel(_CompositeModel):
    """
    Composite model that evaluates models in parallel.

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
                ("A SummedCompositeModel expects n_inputs = n_outputs for "
                 "all transforms")

        super(SummedCompositeModel, self).__init__(transforms, n_inputs,
                                                   n_outputs)

        self._inmap = inmap
        self._outmap = outmap

    def __call__(self, *data):
        """Transforms data using this model."""

        if len(data) == 1:
            if not isinstance(data[0], LabeledInput):
                x = data[0]
                deltas = sum(tr(x) for tr in self._transforms)
                return deltas
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
                sum_of_deltas = [np.zeros_like(x) for x in inlist]
                for transform in self._transforms:
                    delta = [transform(*inlist)]
                    for i in range(len(sum_of_deltas)):
                        sum_of_deltas[i] += delta[i]

                for outcoo, delta in izip(self._outmap, sum_of_deltas):
                    setattr(labeled_input, outcoo, delta)
                # always return the entire labeled object, not just the result
                # since this may be part of another composite transform
                return labeled_input
        else:
            result = self._transforms[0](*data)
            assert self.n_inputs == self.n_outputs
            for tr in self._transforms[1:]:
                result += tr(*data)
            return result


class Parametric1DModel(ParametricModel):
    """
    Base class for one dimensional parametric models.

    This class provides an easier interface to defining new models.
    Examples can be found in functional_models.py

    Parameters
    ----------
    parameters : dictionary
        Dictionary of model parameters with initialisation values
        {'parameter_name': 'parameter_value'}
    """

    @format_input
    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        return self.eval(x, *self.param_sets)


class Parametric2DModel(ParametricModel):
    """
    Base class for two dimensional parametric models.

    This class provides an easier interface to defining new models.
    Examples can be found in functional_models.py

    Parameters
    ----------
    parameter_dict : dictionary
        Dictionary of model parameters with initialization values
        {'parameter_name': 'parameter_value'}
    """

    n_inputs = 2
    n_outputs = 1

    @format_input
    def __call__(self, x, y):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array like or a number
            input
        """

        return self.eval(x, y, *self.param_sets)
