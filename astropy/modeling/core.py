# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines base classes for all models.  The base class of all
models is `~astropy.modeling.Model`. `~astropy.modeling.FittableModel` is
the base class for all fittable models. Fittable models can be linear or
nonlinear in a regression analysis sense.

All models provide a `__call__` method which performs the transformation in
a purely mathematical way, i.e. the models are unitless.  Model instances can
represent either a single model, or a "model set" representing multiple copies
of the same type of model, but with potentially different values of the
parameters in each model making up the set.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import abc
import copy
import functools
import inspect
import warnings

import numpy as np

from ..utils import indent, isiterable
from ..extern import six
from ..extern.six.moves import zip as izip
from ..extern.six.moves import range
from ..table import Table
from ..utils import deprecated
from ..utils.exceptions import AstropyDeprecationWarning
from .utils import array_repr_oneline, check_broadcast, IncompatibleShapeError

from .parameters import Parameter, InputParameterError

__all__ = ['Model', 'FittableModel', 'SummedCompositeModel',
           'SerialCompositeModel', 'LabeledInput', 'FittableModel',
           'Fittable1DModel', 'Fittable2DModel', 'ModelDefinitionError',
           'format_input']


__doctest_skip__ = ['.']


class ModelDefinitionError(Exception):
    """Used for incorrect models definitions"""


def format_input(func):
    """
    Wraps a model's ``__call__`` method so that the input arrays are converted
    into the appropriate shape given the model's parameter dimensions.

    Wraps the result to match the shape of the last input array.
    """

    argspec = inspect.getargspec(func)

    @functools.wraps(func)
    def wrapped_call(self, *inputs, **kwargs):
        model_set_axis = kwargs.pop('model_set_axis', None)

        if model_set_axis is None:
            # By default the model_set_axis for the input is assumed to be the
            # same as that for the parameters the model was defined with
            # TODO: Ensure that negative model_set_axis arguments are respected
            model_set_axis = self.model_set_axis

        n_models = len(self)

        params = [getattr(self, name) for name in self.param_names]
        inputs = [np.asanyarray(_input, dtype=float) for _input in inputs]

        scalar_params = all(not param.shape for param in params)
        scalar_inputs = all(not np.shape(_input) for _input in inputs)

        if n_models == 1 and scalar_params and scalar_inputs:
            # Simplest case is either a parameterless models (currently I don't
            # think we have any but they could exist in principle) or a single
            # model (not a model set) with all scalar paramaters and all scalar
            # inputs
            outputs = func(self, *inputs, **kwargs)

            if self.n_outputs == 1:
                return np.asscalar(outputs)
            else:
                return tuple(np.asscalar(output) for output in outputs)

        input_names = argspec.args[1:1 + len(inputs)]

        _validate_input_shapes(inputs, input_names, n_models,
                               model_set_axis, self.standard_broadcasting)

        # The input formatting required for single models versus a multiple
        # model set are different enough that they've been split into separate
        # subroutines
        if n_models == 1:
            return _format_single_model_input(func, self, params, inputs,
                                              input_names, **kwargs)
        else:
            return _format_model_set_input(func, self, params, inputs,
                                           input_names, n_models,
                                           model_set_axis, **kwargs)

    return wrapped_call


class _ModelMeta(abc.ABCMeta):
    """
    Metaclass for Model.

    Currently just handles auto-generating the param_names list based on
    Parameter descriptors declared at the class-level of Model subclasses.
    """

    def __new__(mcls, name, bases, members):
        param_names = members.get('param_names', [])
        parameters = {}
        for key, value in members.items():
            if not isinstance(value, Parameter):
                continue
            if not value.name:
                # Name not explicitly given in the constructor; add the name
                # automatically via the attribute name
                value._name = key
                value._attr = '_' + key
            if value.name != key:
                raise ModelDefinitionError(
                    "Parameters must be defined with the same name as the "
                    "class attribute they are assigned to.  Parameters may "
                    "take their name from the class attribute automatically "
                    "if the name argument is not given when initializing "
                    "them.")
            parameters[value.name] = value

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
            members['_param_orders'] = \
                    dict((name, idx) for idx, name in enumerate(param_names))

        return super(_ModelMeta, mcls).__new__(mcls, name, bases, members)


@six.add_metaclass(_ModelMeta)
class Model(object):
    """
    Base class for all models.

    This is an abstract class and should not be instantiated directly.

    This class sets the constraints and other properties for all individual
    parameters and performs parameter validation.

    Parameters
    ----------
    param_dim : int
        Number of parameter sets
    fixed : dict
        Dictionary ``{parameter_name: bool}`` setting the fixed constraint
        for one or more parameters.  `True` means the parameter is held fixed
        during fitting and is prevented from updates once an instance of the
        model has been created.

        Alternatively the `~astropy.modeling.Parameter.fixed` property of a
        parameter may be used to lock or unlock individual parameters.
    tied : dict
        Dictionary ``{parameter_name: callable}`` of parameters which are
        linked to some other parameter. The dictionary values are callables
        providing the linking relationship.

        Alternatively the `~astropy.modeling.Parameter.tied` property of a
        parameter may be used to set the ``tied`` constraint on individual
        parameters.
    bounds : dict
        Dictionary ``{parameter_name: value}`` of lower and upper bounds of
        parameters. Keys are parameter names. Values are a list of length 2
        giving the desired range for the parameter.

        Alternatively the `~astropy.modeling.Parameter.min` and
        `~astropy.modeling.Parameter.max` or
        ~astropy.modeling.Parameter.bounds` properties of a parameter may be
        used to set bounds on individual parameters.
    eqcons : list
        List of functions of length n such that ``eqcons[j](x0, *args) == 0.0``
        in a successfully optimized problem.
    ineqcons : list
        List of functions of length n such that ``ieqcons[j](x0, *args) >=
        0.0`` is a successfully optimized problem.

    Examples
    --------
    >>> from astropy.modeling import models
    >>> def tie_center(model):
    ...         mean = 50 * model.stddev
    ...         return mean
    >>> tied_parameters = {'mean': tie_center}

    Specify that ``'mean'`` is a tied parameter in one of two ways:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                        tied=tied_parameters)

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.mean.tied
    False
    >>> g1.mean.tied = tie_center
    >>> g1.mean.tied
    <function tie_center at 0x...>

    Fixed parameters:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                        fixed={'stddev': True})
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

    parameter_constraints = ['fixed', 'tied', 'bounds']
    model_constraints = ['eqcons', 'ineqcons']

    param_names = []
    """
    List of names of the parameters that describe models of this type.

    The parameters in this list are in the same order they should be passed in
    when initializing a model of a specific type.  Some types of models, such
    as polynomial models, have a different number of parameters depending on
    some other property of the model, such as the degree.
    """

    n_inputs = 1
    n_outputs = 1
    standard_broadcasting = True
    fittable = False
    linear = True

    def __init__(self, *args, **kwargs):
        super(Model, self).__init__()
        self._initialize_constraints(kwargs)
        # Remaining keyword args are either parameter values or invalid
        # Parameter values must be passed in as keyword arguments in order to
        # distinguish them
        self._initialize_parameters(args, kwargs)

    def __repr__(self):
        return self._format_repr()

    def __str__(self):
        return self._format_str()

    def __len__(self):
        return self._n_models

    @abc.abstractmethod
    def __call__(self, *args, **kwargs):
        """Evaluate the model on some input variables."""

    @property
    @deprecated('0.4', alternative='len(model)')
    def param_dim(self):
        return self._n_models

    @property
    def model_set_axis(self):
        return self._model_set_axis

    @property
    def param_sets(self):
        """
        Return parameters as a pset.

        This is an array where each column represents one parameter set.
        """

        values = [getattr(self, name).value for name in self.param_names]

        # Ensure parameter values are broadcastable
        for name, shape in six.iteritems(self._param_broadcast_shapes):
            idx = self._param_orders[name]
            values[idx] = values[idx].reshape(shape)

        shapes = [np.shape(value) for value in values]

        if len(self) == 1:
            # Add a single param set axis to the parameter's value (thus
            # converting scalars to shape (1,) array values) for consistency
            values = [np.array([value]) for value in values]

        if len(set(shapes)) != 1:
            # If the parameters are not all the same shape, converting to an
            # array is going to produce an object array
            # However the way Numpy creates object arrays is tricky in that it
            # will recurse into array objects in the list and break them up
            # into separate objects.  Doing things this way ensures a 1-D
            # object array the elements of which are the individual parameter
            # arrays.  There's not much reason to do this over returning a list
            # except for consistency
            psets = np.empty(len(values), dtype=object)
            psets[:] = values
            return psets

        return np.array(values)

    @property
    def parameters(self):
        """
        A flattened array of all parameter values in all parameter sets.

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

    @property
    def fixed(self):
        """
        A `dict` mapping parameter names to their fixed constraint.
        """

        return self._constraints['fixed']

    @property
    def tied(self):
        """
        A `dict` mapping parameter names to their tied constraint.
        """

        return self._constraints['tied']

    @property
    def bounds(self):
        """
        A `dict` mapping parameter names to their upper and lower bounds as
        ``(min, max)`` tuples.
        """

        return self._constraints['bounds']

    @property
    def eqcons(self):
        """List of parameter equality constraints."""

        return self._constraints['eqcons']

    @property
    def ineqcons(self):
        """List of parameter inequality constraints."""

        return self._constraints['ineqcons']

    def inverse(self):
        """Returns a callable object which performs the inverse transform."""

        raise NotImplementedError("An analytical inverse transform has not "
                                  "been implemented for this model.")

    def invert(self):
        """Invert coordinates iteratively if possible."""

        raise NotImplementedError("Subclasses should implement this")

    def add_model(self, model, mode):
        """
        Create a CompositeModel by chaining the current model with the new one
        using the specified mode.

        Parameters
        ----------
        model : an instance of a subclass of Model
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
            return SummedCompositeModel([self, model])
        elif mode in ['serial', 's']:
            return SerialCompositeModel([self, model])
        else:
            raise InputParameterError("Unrecognized mode {0}".format(mode))

    def copy(self):
        """
        Return a copy of this model.

        Uses a deep copy so that all model attributes, including parameter
        values, are copied as well.
        """

        return copy.deepcopy(self)

    def _initialize_constraints(self, kwargs):
        """
        Pop parameter constraint values off the keyword arguments passed to
        `Model.__init__` and store them in private instance attributes.
        """

        self._constraints = {}
        # Pop any constraints off the keyword arguments
        for constraint in self.parameter_constraints:
            values = kwargs.pop(constraint, {})
            self._constraints[constraint] = values

            # Update with default parameter constraints
            for param_name in self.param_names:
                param = getattr(self, param_name)

                # Parameters don't have all constraint types
                value = getattr(param, constraint)
                if value is not None:
                    self._constraints[constraint][param_name] = value

        for constraint in self.model_constraints:
            values = kwargs.pop(constraint, [])
            self._constraints[constraint] = values

    def _initialize_parameters(self, args, kwargs):
        """
        Initialize the _parameters array that stores raw parameter values for
        all parameter sets for use with vectorized fitting algorithms; on
        FittableModels the _param_name attributes actually just reference
        slices of this array.
        """

        n_models = None
        # Pop off param_dim and handle backwards compatibility
        if 'param_dim' in kwargs:
            n_models = kwargs.pop('param_dim')
            warnings.warn(
                'The param_dim argument to {0}.__init__ is deprecated; '
                'use n_models instead.  See also the model_set_axis argument '
                'and related discussion in the docstring for Model.'.format(
                    self.__class__.__name__), AstropyDeprecationWarning)
            if 'n_models' in kwargs:
                raise TypeError(
                    "param_dim and n_models cannot both be specified; use "
                    "n_models, as param_dim is deprecated")
        else:
            n_models = kwargs.pop('n_models', None)

        if not (n_models is None or
                    (isinstance(n_models, int) and n_models >=1)):
            raise ValueError(
                "n_models must be either None (in which case it is "
                "determined from the model_set_axis of the parameter initial "
                "values) or it must be a positive integer "
                "(got {0!r})".format(n_models))

        model_set_axis = kwargs.pop('model_set_axis', None)
        if model_set_axis is None:
            if n_models is not None and n_models > 1:
                # Default to zero
                model_set_axis = 0
            else:
                # Otherwise disable
                model_set_axis = False
        else:
            if not (model_set_axis is False or
                    (isinstance(model_set_axis, int) and
                        not isinstance(model_set_axis, bool))):
                raise ValueError(
                    "model_set_axis must be either False or an integer "
                    "specifying the parameter array axis to map to each "
                    "model in a set of models (got {0!r}).".format(
                        model_set_axis))

        # Process positional arguments by matching them up with the
        # corresponding parameters in self.param_names--if any also appear as
        # keyword arguments this presents a conflict
        params = {}
        if len(args) > len(self.param_names):
            raise TypeError(
                "{0}.__init__() takes at most {1} positional arguments ({2} "
                "given)".format(self.__class__.__name__, len(self.param_names),
                                len(args)))

        for idx, arg in enumerate(args):
            if arg is None:
                # A value of None implies using the default value, if exists
                continue
            params[self.param_names[idx]] = np.asanyarray(arg, dtype=np.float)

        # At this point the only remaining keyword arguments should be
        # parameter names; any others are in error.
        for param_name in self.param_names:
            if param_name in kwargs:
                if param_name in params:
                    raise TypeError(
                        "{0}.__init__() got multiple values for parameter "
                        "{1!r}".format(self.__class__.__name__, param_name))
                value = kwargs.pop(param_name)
                if value is None:
                    continue
                params[param_name] = np.asanyarray(value, dtype=np.float)

        if kwargs:
            # If any keyword arguments were left over at this point they are
            # invalid--the base class should only be passed the parameter
            # values, constraints, and param_dim
            for kwarg in kwargs:
                # Just raise an error on the first unrecognized argument
                raise TypeError(
                    '{0}.__init__() got an unrecognized parameter '
                    '{1!r}'.format(self.__class__.__name__, kwarg))

        # Determine the number of model sets: If the model_set_axis is
        # None then there is just one parameter set; otherwise it is determined
        # by the size of that axis on the first parameter--if the other
        # parameters don't have the right number of axes or the sizes of their
        # model_set_axis don't match an error is raised
        if model_set_axis is not False and n_models != 1 and params:
            max_ndim = 0
            if model_set_axis < 0:
                min_ndim = abs(model_set_axis)
            else:
                min_ndim = model_set_axis + 1

            for name, value in six.iteritems(params):
                param_ndim = np.ndim(value)
                if param_ndim < min_ndim:
                    raise InputParameterError(
                        "All parameter values must be arrays of dimension "
                        "at least {0} for model_set_axis={1} (the value "
                        "given for {2!r} is only {3}-dimensional)".format(
                            min_ndim, model_set_axis, name, param_ndim))

                max_ndim = max(max_ndim, param_ndim)

                if n_models is None:
                    # Use the dimensions of the first parameter to determine
                    # the number of model sets
                    n_models = value.shape[model_set_axis]
                elif value.shape[model_set_axis] != n_models:
                    raise InputParameterError(
                        "Inconsistent dimensions for parameter {0!r} for "
                        "{1} model sets.  The length of axis {2} must be the "
                        "same for all input parameter values".format(
                        name, n_models, model_set_axis))

            self._param_broadcast_shapes = self._check_param_broadcast(
                    params, max_ndim, model_set_axis)
        else:
            if n_models is None:
                n_models = 1

            self._param_broadcast_shapes = self._check_param_broadcast(
                    params, None, None)

        # First we need to determine how much array space is needed by all the
        # parameters based on the number of parameters, the shape each input
        # parameter, and the param_dim
        self._n_models = n_models
        self._model_set_axis = model_set_axis
        self._param_metrics = {}
        total_size = 0

        for name in self.param_names:
            if params.get(name) is None:
                default = getattr(self, name).default

                if default is None:
                    # No value was supplied for the parameter, and the
                    # parameter does not have a default--therefor the model is
                    # underspecified
                    raise TypeError(
                        "{0}.__init__() requires a value for parameter "
                        "{1!r}".format(self.__class__.__name__, name))

                value = params[name] = default
            else:
                value = params[name]

            param_size = np.size(value)
            param_shape = np.shape(value)

            param_slice = slice(total_size, total_size + param_size)
            self._param_metrics[name] = (param_slice, param_shape)
            total_size += param_size

        self._parameters = np.empty(total_size, dtype=np.float64)
        # Now set the parameter values (this will also fill
        # self._parameters)
        for name, value in params.items():
            setattr(self, name, value)

    def _check_param_broadcast(self, params, max_ndim, model_set_axis):
        """
        This subroutine checks that all parameter arrays can be broadcast
        against each other, and determimes the shapes parameters must have in
        order to broadcast correctly.

        If model_set_axis is None this merely checks that the parameters
        broadcast and returns an empty dict if so.  This mode is only used for
        single model sets.
        """

        broadcast_shapes = {}
        all_shapes = []
        param_names = []

        for name in self.param_names:
            # Previously this just used iteritems(params), but we loop over all
            # param_names instead just to ensure some determinism in the
            # ordering behavior
            if name not in params:
                continue

            value = params[name]
            param_names.append(name)
            # We've already checked that each parameter array is compatible in
            # the model_set_axis dimension, but now we need to check the
            # dimensions excluding that axis
            # Split the array dimensions into the axes before model_set_axis
            # and after model_set_axis
            param_shape = np.shape(value)

            param_ndim = len(param_shape)
            if max_ndim is not None and param_ndim < max_ndim:
                # All arrays have the same number of dimensions up to the
                # model_set_axis dimension, but after that they may have a
                # different number of trailing axes.  The number of trailing
                # axes must be extended for mutual compatibility.  For example
                # if max_ndim = 3 and model_set_axis = 0, an array with the
                # shape (2, 2) must be extended to (2, 1, 2).  However, an
                # array with shape (2,) is extended to (2, 1).
                new_axes = (1,) * (max_ndim - param_ndim)

                if model_set_axis < 0:
                    # Just need to prepend axes to make up the difference
                    broadcast_shape = new_axes + param_shape
                else:
                    broadcast_shape = (param_shape[:model_set_axis + 1] +
                                       new_axes +
                                       param_shape[model_set_axis + 1:])
                broadcast_shapes[name] = broadcast_shape
                all_shapes.append(broadcast_shape)
            else:
                all_shapes.append(param_shape)

        # Now check mutual broadcastability of all shapes
        try:
            check_broadcast(*all_shapes)
        except IncompatibleShapeError as exc:
            shape_a, shape_a_idx, shape_b, shape_b_idx = exc.args
            param_a = param_names[shape_a_idx]
            param_b = param_names[shape_b_idx]

            raise InputParameterError(
                "Parameter {0!r} of shape {1!r} cannot be broadcast with "
                "parameter {2!r} of shape {3!r}.  All parameter arrays "
                "must have shapes that are mutually compatible according "
                "to the broadcasting rules.".format(param_a, shape_a,
                                                    param_b, shape_b))

        return broadcast_shapes

    def _format_repr(self, args=[], kwargs={}, defaults={}):
        """
        Internal implementation of ``__repr__``.

        This is separated out for ease of use by subclasses that wish to
        override the default ``__repr__`` while keeping the same basic
        formatting.
        """

        # TODO: I think this could be reworked to preset model sets better

        parts = ['<{0}('.format(self.__class__.__name__)]

        parts.append(', '.join(repr(a) for a in args))

        if args:
            parts.append(', ')

        parts.append(', '.join(
            "{0}={1}".format(
                name, array_repr_oneline(getattr(self, name).value))
            for name in self.param_names))

        for kwarg, value in kwargs.items():
            if kwarg  in defaults and defaults[kwarg] != value:
                continue
            parts.append(', {0}={1!r}'.format(kwarg, value))

        if len(self) > 1:
            parts.append(", n_models={0}".format(len(self)))

        parts.append(')>')

        return ''.join(parts)

    def _format_str(self, keywords=[]):
        """
        Internal implementation of ``__str__``.

        This is separated out for ease of use by subclasses that wish to
        override the default ``__str__`` while keeping the same basic
        formatting.
        """

        default_keywords = [
            ('Model', self.__class__.__name__),
            ('Inputs', self.n_inputs),
            ('Outputs', self.n_outputs),
            ('Model set size', len(self))
        ]

        parts = ['{0}: {1}'.format(keyword, value)
                 for keyword, value in default_keywords + keywords]

        parts.append('Parameters:')

        if len(self) == 1:
            columns = [[getattr(self, name).value]
                       for name in self.param_names]
        else:
            columns = [getattr(self, name).value
                       for name in self.param_names]

        param_table = Table(columns, names=self.param_names)

        parts.append(indent(str(param_table), width=4))

        return '\n'.join(parts)


class FittableModel(Model):
    """
    Base class for models that can be fitted using the built-in fitting
    algorithms.
    """

    linear = False
    # derivative with respect to parameters
    fit_deriv = None
    """
    Function (similar to the model's ``eval``) to compute the derivatives of
    the model with respect to its parameters, for use by fitting algorithms.
    """
    # Flag that indicates if the model derivatives with respect to parameters
    # are given in columns or rows
    col_fit_deriv = True
    fittable = True


class LabeledInput(dict):
    """
    Used by `SerialCompositeModel` and `SummedCompositeModel` to choose input
    data using labels.

    This is a container assigning labels (names) to all input data arrays to a
    composite model.

    Parameters
    ----------
    data : list
        List of all input data
    labels : list of strings
        names matching each coordinate in data

    Examples
    --------
    >>> y, x = np.mgrid[:5, :5]
    >>> l = np.arange(10)
    >>> labeled_input = LabeledInput([x, y, l], ['x', 'y', 'pixel'])
    >>> labeled_input.x
    array([[0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4]])
    >>> labeled_input['x']
    array([[0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4]])
    """

    def __init__(self, data, labels):
        dict.__init__(self)
        if len(labels) != len(data):
            raise TypeError("Number of labels and data doesn't match")
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
            if given this is a dictionary of ``{label: value}`` pairs
        """

        if kw:
            if label is None or value is None:
                self.update(kw)
            else:
                kw[label] = value
                self.update(kw)
        else:
            kw = dict({label: value})
            if label is None or value is None:
                raise TypeError("Expected label and value to be defined")
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
        return '<{0}([\n{1}\n])>'.format(
            self.__class__.__name__,
            indent(',\n'.join(repr(tr) for tr in self._transforms),
                   width=4))

    def __str__(self):
        parts = ['Model: {0}'.format(self.__class__.__name__)]
        for tr in self._transforms:
            parts.append(indent(str(tr), width=4))
        return '\n'.join(parts)

    def add_model(self, transf, inmap, outmap):
        self[transf] = [inmap, outmap]

    def invert(self):
        raise NotImplementedError("Subclasses should implement this")

    def __call__(self):
        # implemented by subclasses
        raise NotImplementedError("Subclasses should implement this")

    @property
    def param_sets(self):
        raise NotImplementedError(
            "Composite models do not currently support multiple "
            "parameter sets.")

    @property
    def parameters(self):
        raise NotImplementedError(
            "Composite models do not currently support the .parameters "
            "array.")


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
    n_inputs : int
        dimension of input space (e.g. 2 for a spatial model)
    n_outputs : int
        dimension of output

    Notes
    -----
    Output values of one model are used as input values of another.
    Obviously the order of the models matters.

    Examples
    --------
    Apply a 2D rotation followed by a shift in x and y::

        >>> import numpy as np
        >>> from astropy.modeling import models, LabeledInput, SerialCompositeModel
        >>> y, x = np.mgrid[:5, :5]
        >>> rotation = models.Rotation2D(angle=23.5)
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
            if n_outputs is None:
                raise TypeError("Expected n_inputs and n_outputs")

        super(SerialCompositeModel, self).__init__(transforms, n_inputs,
                                                   n_outputs)

        if transforms and inmap and outmap:
            if not (len(transforms) == len(inmap) == len(outmap)):
                raise ValueError("Expected sequences of transform, "
                                 "inmap and outmap to have the same length")

        if inmap is None:
            inmap = [None] * len(transforms)

        if outmap is None:
            outmap = [None] * len(transforms)

        self._inmap = inmap
        self._outmap = outmap

    def inverse(self):
        try:
            transforms = []
            for transform in self._transforms[::-1]:
                transforms.append(transform.inverse())
        except NotImplementedError:
            raise NotImplementedError(
                "An analytical inverse has not been implemented for "
                "{0} models.".format(transform.__class__.__name__))
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
                if self._transforms[0].n_inputs != 1:
                    raise TypeError("First transform expects {0} inputs, 1 "
                                    "given".format(self._transforms[0].n_inputs))

                result = data[0]
                for tr in self._transforms:
                    result = tr(result)
                return result
            else:
                labeled_input = data[0].copy()
                # we want to return the entire labeled object because some
                # parts of it may be used in another transform of which this
                # one is a component
                if self._inmap is None:
                    raise TypeError("Parameter 'inmap' must be provided when "
                                    "input is a labeled object.")
                if self._outmap is None:
                    raise TypeError("Parameter 'outmap' must be provided when "
                                    "input is a labeled object")

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
            if self.n_inputs != len(data):
                raise TypeError("This transform expects {0} inputs".
                                format(self._n_inputs))

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
            if not (transform.n_inputs == transform.n_outputs == n_inputs):
                raise ValueError("A SummedCompositeModel expects n_inputs = "
                                 "n_outputs for all transforms")

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
                if self._inmap is None:
                    raise TypeError("Parameter 'inmap' must be provided when "
                                    "input is a labeled object.")
                if self._outmap is None:
                    raise TypeError("Parameter 'outmap' must be provided when "
                                    "input is a labeled object")
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
            if self.n_inputs != self.n_outputs:
                raise ValueError("Expected equal number of inputs and outputs")
            for tr in self._transforms[1:]:
                result += tr(*data)
            return result


class Fittable1DModel(FittableModel):
    """
    Base class for one-dimensional fittable models.

    This class provides an easier interface to defining new models.
    Examples can be found in `astropy.modeling.functional_models`.
    """

    @abc.abstractmethod
    def eval(self):
        """
        A method, `classmethod`, or `staticmethod` that implements evaluation
        of the function represented by this model.

        It must take arguments of the function's independent variables,
        followed by the function's parameters given in the same order they are
        listed by `Model.param_names`.
        """

    @format_input
    def __call__(self, x, model_set_axis=None):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array-like or numeric value
            Input coordinate values.

        model_set_axis : `int` or `False`, optional
            For `Model` instances representing a multiple-model set, this picks
            out which axis of the input array is used to map inputs to specific
            models in the set.  If `False`, this indicates that the input array
            has no such axis, and instead the same input should be broadcast to
            all models in the set.
        """

        return self.eval(x, *self.param_sets)


class Fittable2DModel(FittableModel):
    """
    Base class for one-dimensional fittable models.

    This class provides an easier interface to defining new models.
    Examples can be found in `astropy.modeling.functional_models`.
    """

    n_inputs = 2
    n_outputs = 1

    @abc.abstractmethod
    def eval(self):
        """
        A method, `classmethod`, or `staticmethod` that implements evaluation
        of the function represented by this model.

        It must take arguments of the function's independent variables,
        followed by the function's parameters given in the same order they are
        listed by `Model.param_names`.
        """

    @format_input
    def __call__(self, x, y, model_set_axis=None):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : array-like or numeric value
            First input coordinate values.

        y : array-like or numeric value
            Second input coordinate values.

        model_set_axis : `int` or `False`, optional
            For `Model` instances representing a multiple-model set, this picks
            out which axis of the input array is used to map inputs to specific
            models in the set.  If `False`, this indicates that the input array
            has no such axis, and instead the same input should be broadcast to
            all models in the set.
        """

        return self.eval(x, y, *self.param_sets)


def _format_single_model_input(func, model, params, inputs, input_names,
                               **kwargs):
    broadcasts = []

    for idx, _input in enumerate(inputs):
        _input = np.asanyarray(_input, dtype=np.float)
        input_shape = _input.shape
        max_broadcast = ()

        for param in params:
            try:
                if model.standard_broadcasting:
                    broadcast = check_broadcast(input_shape, param.shape)
                else:
                    broadcast = input_shape
            except IncompatibleShapeError:
                raise ValueError(
                    "Model input argument {0!r} of shape {1!r} cannot be "
                    "broadcast with parameter {2!r} of shape "
                    "{3!r}.".format(input_names[idx], input_shape,
                                    param.name, param.shape))

            if len(broadcast) > len(max_broadcast):
                max_broadcast = broadcast
            elif len(broadcast) == len(max_broadcast):
                max_broadcast = max(max_broadcast, broadcast)

        broadcasts.append(max_broadcast)

    outputs = func(model, *inputs, **kwargs)

    if model.n_outputs > model.n_inputs:
        if len(set(broadcasts)) > 1:
            raise ValueError(
                "For models with n_outputs > n_inputs, the combination of "
                "all inputs and parameters must broadcast to the same shape, "
                "which will be used as the shape of all outputs.  In this "
                "case some of the inputs had different shapes, so it is "
                "ambiguous how to format outputs for this model.  Try using "
                "inputs that are all the same size and shape.")
        else:
            # Extend the broadcasts list to include shapes for all outputs
            extra_outputs = model.n_outputs - model.n_inputs
            broadcasts.extend([broadcasts[0]] * extra_outputs)

    if model.n_outputs == 1:
        outputs = [outputs]
    else:
        outputs = list(outputs)

    for idx, output in enumerate(outputs):
        if broadcasts[idx]:
            outputs[idx] = output.reshape(broadcasts[idx])

    if model.n_outputs == 1:
        return outputs[0]
    else:
        return tuple(outputs)


def _format_model_set_input(func, model, params, inputs, input_names,
                            n_models, model_set_axis, **kwargs):
    reshaped = []
    pivots = []

    for idx, _input in enumerate(inputs):
        _input = np.asanyarray(_input, dtype=np.float)

        max_param_shape = ()

        if n_models > 1 and model_set_axis is not False:
            # Use the shape of the input *excluding* the model axis
            input_shape = (_input.shape[:model_set_axis] +
                           _input.shape[model_set_axis + 1:])
        else:
            input_shape = _input.shape

        for param in params:
            try:
                check_broadcast(input_shape, param.shape)
            except IncompatibleShapeError:
                raise ValueError(
                    "Model input argument {0!r} of shape {1!r} cannot be "
                    "broadcast with parameter {2!r} of shape "
                    "{3!r}.".format(input_names[idx], input_shape,
                                    param.name, param.shape))

            if len(param.shape) > len(max_param_shape):
                max_param_shape = param.shape

        # We've now determined that, excluding the model_set_axis, the
        # input can broadcast with all the parameters
        input_ndim = len(input_shape)
        if model_set_axis is False:
            if len(max_param_shape) > input_ndim:
                # Just needs to prepend new axes to the input
                n_new_axes = 1 + len(max_param_shape) - input_ndim
                new_axes = (1,) * n_new_axes
                new_shape = new_axes + _input.shape
                pivot = model.model_set_axis
            else:
                pivot = input_ndim - len(max_param_shape)
                new_shape = (_input.shape[:pivot] + (1,) +
                             _input.shape[pivot:])
            new_input = _input.reshape(new_shape)
        else:
            if len(max_param_shape) >= input_ndim:
                n_new_axes = len(max_param_shape) - input_ndim
                pivot = model.model_set_axis
                new_axes = (1,) * n_new_axes
                new_shape = (_input.shape[:pivot + 1] + new_axes +
                             _input.shape[pivot + 1:])
                new_input = _input.reshape(new_shape)
            else:
                pivot = _input.ndim - len(max_param_shape) - 1
                new_input = np.rollaxis(_input, model_set_axis,
                                        pivot + 1)

        pivots.append(pivot)
        reshaped.append(new_input)

    outputs = func(model, *reshaped, **kwargs)

    if model.n_outputs == 1:
        outputs = [outputs]
    else:
        outputs = list(outputs)

    for idx, output in enumerate(outputs):
        pivot = pivots[idx]
        if pivot < output.ndim and pivot != model.model_set_axis:
            outputs[idx] = np.rollaxis(output, pivot,
                                       model.model_set_axis)

    if model.n_outputs == 1:
        return outputs[0]
    else:
        return tuple(outputs)


def _validate_input_shapes(inputs, argnames, n_models, model_set_axis,
                           validate_broadcasting):
    """
    Perform basic validation of model inputs--that they are mutually
    broadcastable and that they have the minimum dimensions for the given
    model_set_axis.

    If validation succeeds, returns the total shape that will result from
    broadcasting the input arrays with each other.
    """

    check_model_set_axis = n_models > 1 and model_set_axis is not False

    if not (validate_broadcasting or check_model_set_axis):
        # Nothing else needed here
        return

    all_shapes = []

    for idx, _input in enumerate(inputs):
        input_shape = np.shape(_input)
        # Ensure that the input's model_set_axis matches the model's
        # n_models
        if input_shape and check_model_set_axis:
            # Note: Scalar inputs *only* get a pass on this
            if len(input_shape) < model_set_axis + 1:
                raise ValueError(
                    "For model_set_axis={0}, all inputs must be at "
                    "least {1}-dimensional.".format(
                        model_set_axis, model_set_axis + 1))
            elif input_shape[model_set_axis] != n_models:
                raise ValueError(
                    "Input argument {0!r} does not have the correct "
                    "dimensions in model_set_axis={1} for a model set with "
                    "n_models={2}.".format(argnames[idx], model_set_axis,
                                           n_models))
        all_shapes.append(input_shape)

    if not validate_broadcasting:
        return

    try:
        input_broadcast = check_broadcast(*all_shapes)
    except IncompatibleShapesError as exc:
        shape_a, shape_a_idx, shape_b, shape_b_idx = exc.args
        arg_a = argnames[shape_a_idx]
        arg_b = argnames[shape_b_idx]

        raise ValueError(
            "Model input argument {0!r} of shape {1!r} cannot "
            "be broadcast with input {2!r} of shape {3!r}".format(
                arg_a, shape_a, arg_b, shape_b))

    return input_broadcast
