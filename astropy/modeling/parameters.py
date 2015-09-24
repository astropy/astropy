# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines two classes that deal with parameters.

It is unlikely users will need to work with these classes directly, unless they
define their own models.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import functools
import numbers
import types

import numpy as np

from ..utils import isiterable, OrderedDescriptor
from ..extern import six

from .utils import get_inputs_and_params

__all__ = ['Parameter', 'InputParameterError']


class InputParameterError(ValueError):
    """Used for incorrect input parameter values and definitions."""


def _tofloat(value):
    """Convert a parameter to float or float array"""

    if isiterable(value):
        try:
            value = np.array(value, dtype=np.float)
        except (TypeError, ValueError):
            # catch arrays with strings or user errors like different
            # types of parameters in a parameter set
            raise InputParameterError(
                "Parameter of {0} could not be converted to "
                "float".format(type(value)))
    elif isinstance(value, np.ndarray):
        # A scalar/dimensionless array
        value = float(value.item())
    elif isinstance(value, (numbers.Number, np.number)):
        value = float(value)
    elif isinstance(value, bool):
        raise InputParameterError(
            "Expected parameter to be of numerical type, not boolean")
    else:
        raise InputParameterError(
            "Don't know how to convert parameter of {0} to "
            "float".format(type(value)))
    return value


class Parameter(OrderedDescriptor):
    """
    Wraps individual parameters.

    This class represents a model's parameter (in a somewhat broad sense).  It
    acts as both a descriptor that can be assigned to a class attribute to
    describe the parameters accepted by an individual model (this is called an
    "unbound parameter"), or it can act as a proxy for the parameter values on
    an individual model instance (called a "bound parameter").

    Parameter instances never store the actual value of the parameter directly.
    Rather, each instance of a model stores its own parameters parameter values
    in an array.  A *bound* Parameter simply wraps the value in a Parameter
    proxy which provides some additional information about the parameter such
    as its constraints.  In other words, this is a high-level interface to a
    model's adjustable parameter values.

    *Unbound* Parameters are not associated with any specific model instance,
    and are merely used by model classes to determine the names of their
    parameters and other information about each parameter such as their default
    values and default constraints.

    See :ref:`modeling-parameters` for more details.

    Parameters
    ----------
    name : str
        parameter name

        .. warning::

            The fact that `Parameter` accepts ``name`` as an argument is an
            implementation detail, and should not be used directly.  When
            defining a new `Model` class, parameter names are always
            automatically defined by the class attribute they're assigned to.
    description : str
        parameter description
    default : float or array
        default value to use for this parameter
    getter : callable
        a function that wraps the raw (internal) value of the parameter
        when returning the value through the parameter proxy (eg. a
        parameter may be stored internally as radians but returned to the
        user as degrees)
    setter : callable
        a function that wraps any values assigned to this parameter; should
        be the inverse of getter
    fixed : bool
        if True the parameter is not varied during fitting
    tied : callable or False
        if callable is supplied it provides a way to link the value of this
        parameter to another parameter (or some other arbitrary function)
    min : float
        the lower bound of a parameter
    max : float
        the upper bound of a parameter
    bounds : tuple
        specify min and max as a single tuple--bounds may not be specified
        simultaneously with min or max
    model : `Model` instance
        binds the the `Parameter` instance to a specific model upon
        instantiation; this should only be used internally for creating bound
        Parameters, and should not be used for `Parameter` descriptors defined
        as class attributes
    """

    constraints = ('fixed', 'tied', 'bounds')
    """
    Types of constraints a parameter can have.  Excludes 'min' and 'max'
    which are just aliases for the first and second elements of the 'bounds'
    constraint (which is represented as a 2-tuple).
    """

    # Settings for OrderedDescriptor
    _class_attribute_ = '_parameters_'
    _name_attribute_ = '_name'

    def __init__(self, name='', description='', default=None, getter=None,
                 setter=None, fixed=False, tied=False, min=None, max=None,
                 bounds=None, model=None):
        super(Parameter, self).__init__()

        self._name = name
        self.__doc__ = self._description = description.strip()
        self._default = default

        # NOTE: These are *default* constraints--on model instances constraints
        # are taken from the model if set, otherwise the defaults set here are
        # used
        if bounds is not None:
            if min is not None or max is not None:
                raise ValueError(
                    'bounds may not be specified simultaneously with min or '
                    'or max when instantiating Parameter {0}'.format(name))
        else:
            bounds = (min, max)

        self._fixed = fixed
        self._tied = tied
        self._bounds = bounds

        self._order = None
        self._model = None

        # The getter/setter functions take one or two arguments: The first
        # argument is always the value itself (either the value returned or the
        # value being set).  The second argument is optional, but if present
        # will contain a reference to the model object tied to a parameter (if
        # it exists)
        self._getter = self._create_value_wrapper(getter, None)
        self._setter = self._create_value_wrapper(setter, None)

        self._validator = None

        # Only Parameters declared as class-level descriptors require
        # and ordering ID
        if model is not None:
            self._bind(model)

    def __get__(self, obj, objtype):
        if obj is None:
            return self

        # All of the Parameter.__init__ work should already have been done for
        # the class-level descriptor; we can skip that stuff and just copy the
        # existing __dict__ and then bind to the model instance
        parameter = self.__class__.__new__(self.__class__)
        parameter.__dict__.update(self.__dict__)
        parameter._bind(obj)
        return parameter

    def __set__(self, obj, value):
        value = _tofloat(value)

        # Call the validator before the setter
        if self._validator is not None:
            self._validator(obj, value)

        if self._setter is not None:
            setter = self._create_value_wrapper(self._setter, obj)
            value = setter(value)

        self._set_model_value(obj, value)

    def __len__(self):
        if self._model is None:
            raise TypeError('Parameter definitions do not have a length.')
        return len(self._model)

    def __getitem__(self, key):
        value = self.value
        if len(self._model) == 1:
            # Wrap the value in a list so that getitem can work for sensible
            # indices like [0] and [-1]
            value = [value]
        return value[key]

    def __setitem__(self, key, value):
        # Get the existing value and check whether it even makes sense to
        # apply this index
        oldvalue = self.value
        n_models = len(self._model)

        #if n_models == 1:
        #    # Convert the single-dimension value to a list to allow some slices
        #    # that would be compatible with a length-1 array like [:] and [0:]
        #    oldvalue = [oldvalue]

        if isinstance(key, slice):
            if len(oldvalue[key]) == 0:
                raise InputParameterError(
                    "Slice assignment outside the parameter dimensions for "
                    "'{0}'".format(self.name))
            for idx, val in zip(range(*key.indices(len(self))), value):
                self.__setitem__(idx, val)
        else:
            try:
                oldvalue[key] = value
            except IndexError:
                raise InputParameterError(
                    "Input dimension {0} invalid for {1!r} parameter with "
                    "dimension {2}".format(key, self.name, n_models))

    def __repr__(self):
        args = "'{0}'".format(self._name)
        if self._model is None:
            if self._default is not None:
                args += ', default={0}'.format(self._default)
        else:
            args += ', value={0}'.format(self.value)

        for cons in self.constraints:
            val = getattr(self, cons)
            if val not in (None, False, (None, None)):
                # Maybe non-obvious, but False is the default for the fixed and
                # tied constraints
                args += ', {0}={1}'.format(cons, val)

        return "{0}({1})".format(self.__class__.__name__, args)

    @property
    def name(self):
        """Parameter name"""

        return self._name

    @property
    def default(self):
        """Parameter default value"""

        if (self._model is None or self._default is None or
                len(self._model) == 1):
            return self._default

        # Otherwise the model we are providing for has more than one parameter
        # sets, so ensure that the default is repeated the correct number of
        # times along the model_set_axis if necessary
        n_models = len(self._model)
        model_set_axis = self._model._model_set_axis
        default = self._default
        new_shape = (np.shape(default) +
                     (1,) * (model_set_axis + 1 - np.ndim(default)))
        default = np.reshape(default, new_shape)
        # Now roll the new axis into its correct position if necessary
        default = np.rollaxis(default, -1, model_set_axis)
        # Finally repeat the last newly-added axis to match n_models
        default = np.repeat(default, n_models, axis=-1)

        # NOTE: Regardless of what order the last two steps are performed in,
        # the resulting array will *look* the same, but only if the repeat is
        # performed last will it result in a *contiguous* array

        return default

    @property
    def value(self):
        """The unadorned value proxied by this parameter"""

        if self._model is None:
            raise AttributeError('Parameter definition does not have a value')

        value = self._get_model_value(self._model)

        if self._getter is None:
            return value
        else:
            return self._getter(value)

    @value.setter
    def value(self, value):
        if self._model is None:
            raise AttributeError('Cannot set a value on a parameter '
                                 'definition')

        if self._setter is not None:
            val = self._setter(value)

        self._set_model_value(self._model, value)

    @property
    def shape(self):
        """The shape of this parameter's value array."""

        if self._model is None:
            raise AttributeError('Parameter definition does not have a '
                                 'shape.')

        shape = self._model._param_metrics[self._name]['shape']

        if len(self._model) > 1:
            # If we are dealing with a model *set* the shape is the shape of
            # the parameter within a single model in the set
            model_axis = self._model._model_set_axis

            if model_axis < 0:
                model_axis = len(shape) + model_axis

            shape = shape[:model_axis] + shape[model_axis + 1:]

        return shape

    @property
    def size(self):
        """The size of this parameter's value array."""

        # TODO: Rather than using self.value this could be determined from the
        # size of the parameter in _param_metrics

        return np.size(self.value)

    @property
    def fixed(self):
        """
        Boolean indicating if the parameter is kept fixed during fitting.
        """

        if self._model is not None:
            fixed = self._model._constraints['fixed']
            return fixed.get(self._name, self._fixed)
        else:
            return self._fixed

    @fixed.setter
    def fixed(self, value):
        """Fix a parameter"""
        if self._model is not None:
            if not isinstance(value, bool):
                raise TypeError("Fixed can be True or False")
            self._model._constraints['fixed'][self._name] = value
        else:
            raise AttributeError("can't set attribute 'fixed' on Parameter "
                                 "definition")

    @property
    def tied(self):
        """
        Indicates that this parameter is linked to another one.

        A callable which provides the relationship of the two parameters.
        """

        if self._model is not None:
            tied = self._model._constraints['tied']
            return tied.get(self._name, self._tied)
        else:
            return self._tied

    @tied.setter
    def tied(self, value):
        """Tie a parameter"""

        if self._model is not None:
            if not six.callable(value) and value not in (False, None):
                    raise TypeError("Tied must be a callable")
            self._model._constraints['tied'][self._name] = value
        else:
            raise AttributeError("can't set attribute 'tied' on Parameter "
                                 "definition")

    @property
    def bounds(self):
        """The minimum and maximum values of a parameter as a tuple"""

        if self._model is not None:
            bounds = self._model._constraints['bounds']
            return bounds.get(self._name, self._bounds)
        else:
            return self._bounds

    @bounds.setter
    def bounds(self, value):
        """Set the minimum and maximum values of a parameter from a tuple"""

        if self._model is not None:
            _min, _max = value
            if _min is not None:
                if not isinstance(_min, numbers.Number):
                        raise TypeError("Min value must be a number")
                _min = float(_min)

            if _max is not None:
                if not isinstance(_max, numbers.Number):
                        raise TypeError("Max value must be a number")
                _max = float(_max)

            bounds = self._model._constraints.setdefault('bounds', {})
            self._model._constraints['bounds'][self._name] = (_min, _max)
        else:
            raise AttributeError("can't set attribute 'bounds' on Parameter "
                                 "definition")

    @property
    def min(self):
        """A value used as a lower bound when fitting a parameter"""

        return self.bounds[0]

    @min.setter
    def min(self, value):
        """Set a minimum value of a parameter"""

        if self._model is not None:
            self.bounds = (value, self.max)
        else:
            raise AttributeError("can't set attribute 'min' on Parameter "
                                 "definition")

    @property
    def max(self):
        """A value used as an upper bound when fitting a parameter"""

        return self.bounds[1]

    @max.setter
    def max(self, value):
        """Set a maximum value of a parameter."""

        if self._model is not None:
            self.bounds = (self.min, value)
        else:
            raise AttributeError("can't set attribute 'max' on Parameter "
                                 "definition")

    @property
    def validator(self):
        """
        Used as a decorator to set the validator method for a `Parameter`.
        The validator method validates any value set for that parameter.
        It takes two arguments--``self``, which refers to the `Model`
        instance (remember, this is a method defined on a `Model`), and
        the value being set for this parameter.  The validator method's
        return value is ignored, but it may raise an exception if the value
        set on the parameter is invalid (typically an `InputParameterError`
        should be raised, though this is not currently a requirement).

        The decorator *returns* the `Parameter` instance that the validator
        is set on, so the underlying validator method should have the same
        name as the `Parameter` itself (think of this as analogous to
        ``property.setter``).  For example::

            >>> from astropy.modeling import Fittable1DModel
            >>> class TestModel(Fittable1DModel):
            ...     a = Parameter()
            ...     b = Parameter()
            ...
            ...     @a.validator
            ...     def a(self, value):
            ...         # Remember, the value can be an array
            ...         if np.any(value < self.b):
            ...             raise InputParameterError(
            ...                 "parameter 'a' must be greater than or equal "
            ...                 "to parameter 'b'")
            ...
            ...     @staticmethod
            ...     def evaluate(x, a, b):
            ...         return a * x + b
            ...
            >>> m = TestModel(a=1, b=2)  # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
            ...
            InputParameterError: parameter 'a' must be greater than or equal
            to parameter 'b'
            >>> m = TestModel(a=2, b=2)
            >>> m.a = 0  # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
            ...
            InputParameterError: parameter 'a' must be greater than or equal
            to parameter 'b'

        On bound parameters this property returns the validator method itself,
        as a bound method on the `Parameter`.  This is not often as useful, but
        it allows validating a parameter value without setting that parameter::

            >>> m.a.validator(42)  # Passes
            >>> m.a.validator(-42)  # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
            ...
            InputParameterError: parameter 'a' must be greater than or equal
            to parameter 'b'
        """

        if self._model is None:
            # For unbound parameters return the validator setter
            def validator(func, self=self):
                self._validator = func
                return self

            return validator
        else:
            # Return the validator method, bound to the Parameter instance with
            # the name "validator"
            def validator(self, value):
                if self._validator is not None:
                    return self._validator(self._model, value)

            if six.PY2:
                return types.MethodType(validator, self, type(self))
            else:
                return types.MethodType(validator, self)

    def copy(self, name=None, description=None, default=None, getter=None,
             setter=None, fixed=False, tied=False, min=None, max=None,
             bounds=None):
        """
        Make a copy of this `Parameter`, overriding any of its core attributes
        in the process (or an exact copy).

        The arguments to this method are the same as those for the `Parameter`
        initializer.  This simply returns a new `Parameter` instance with any
        or all of the attributes overridden, and so returns the equivalent of:

        .. code:: python

            Parameter(self.name, self.description, ...)

        """

        kwargs = locals().copy()
        del kwargs['self']

        for key, value in six.iteritems(kwargs):
            if value is None:
                # Annoying special cases for min/max where are just aliases for
                # the components of bounds
                if key in ('min', 'max'):
                    continue
                else:
                    if hasattr(self, key):
                        value = getattr(self, key)
                    elif hasattr(self, '_' + key):
                        value = getattr(self, '_' + key)
                kwargs[key] = value

        return self.__class__(**kwargs)

    @property
    def _raw_value(self):
        """
        Currently for internal use only.

        Like Parameter.value but does not pass the result through
        Parameter.getter.  By design this should only be used from bound
        parameters.

        This will probably be removed are retweaked at some point in the
        process of rethinking how parameter values are stored/updated.
        """

        return self._get_model_value(self._model)

    def _bind(self, model):
        """
        Bind the `Parameter` to a specific `Model` instance; don't use this
        directly on *unbound* parameters, i.e. `Parameter` descriptors that
        are defined in class bodies.
        """

        self._model = model
        self._getter = self._create_value_wrapper(self._getter, model)
        self._setter = self._create_value_wrapper(self._setter, model)

    def _get_model_value(self, model):
        """
        This method implements how to retrieve the value of this parameter from
        the model instance.  See also `Parameter._set_model_value`.

        These methods take an explicit model argument rather than using
        self._model so that they can be used from unbound `Parameter`
        instances.
        """

        if not hasattr(model, '_parameters'):
            # The _parameters array hasn't been initialized yet; just translate
            # this to an AttributeError
            raise AttributeError(self._name)

        # Use the _param_metrics to extract the parameter value from the
        # _parameters array
        param_slice = model._param_metrics[self._name]['slice']
        param_shape = model._param_metrics[self._name]['shape']
        value = model._parameters[param_slice]
        if param_shape:
            value = value.reshape(param_shape)
        else:
            value = value[0]
        return value

    def _set_model_value(self, model, value):
        """
        This method implements how to store the value of a parameter on the
        model instance.

        Currently there is only one storage mechanism (via the ._parameters
        array) but other mechanisms may be desireable, in which case really the
        model class itself should dictate this and *not* `Parameter` itself.
        """

        # TODO: Maybe handle exception on invalid input shape
        param_slice = model._param_metrics[self._name]['slice']
        param_shape = model._param_metrics[self._name]['shape']
        param_size = np.prod(param_shape)

        if np.size(value) != param_size:
            raise InputParameterError(
                "Input value for parameter {0!r} does not have {1} elements "
                "as the current value does".format(self._name, param_size))

        model._parameters[param_slice] = np.array(value).ravel()

    @staticmethod
    def _create_value_wrapper(wrapper, model):
        """Wraps a getter/setter function to support optionally passing in
        a reference to the model object as the second argument.

        If a model is tied to this parameter and its getter/setter supports
        a second argument then this creates a partial function using the model
        instance as the second argument.
        """

        if isinstance(wrapper, np.ufunc):
            if wrapper.nin != 1:
                raise TypeError("A numpy.ufunc used for Parameter "
                                "getter/setter may only take one input "
                                "argument")
        elif wrapper is None:
            # Just allow non-wrappers to fall through silently, for convenience
            return None
        else:
            inputs, params = get_inputs_and_params(wrapper)
            nargs = len(inputs)

            if nargs == 1:
                pass
            elif nargs == 2:
                if model is not None:
                    # Don't make a partial function unless we're tied to a
                    # specific model instance
                    model_arg = inputs.args[1].name
                    wrapper = functools.partial(wrapper, **{model_arg: model})
            else:
                raise TypeError("Parameter getter/setter must be a function "
                                "of either one or two arguments")

        return wrapper

    def __array__(self, dtype=None):
        # Make np.asarray(self) work a little more straightforwardly
        return np.asarray(self.value, dtype=dtype)

    def __nonzero__(self):
        if self._model is None:
            return True
        else:
            return bool(self.value)

    __bool__ = __nonzero__

    def __add__(self, val):
        if self._model is None:
            # If we don't do this, __add__ will raise an AttributeError instead
            # (from self.value) which is strange and unexpected
            return NotImplemented
        return self.value + val

    def __radd__(self, val):
        if self._model is None:
            return NotImplemented
        return val + self.value

    def __sub__(self, val):
        if self._model is None:
            return NotImplemented
        return self.value - val

    def __rsub__(self, val):
        if self._model is None:
            return NotImplemented
        return val - self.value

    def __mul__(self, val):
        if self._model is None:
            return NotImplemented
        return self.value * val

    def __rmul__(self, val):
        if self._model is None:
            return NotImplemented
        return val * self.value

    def __pow__(self, val):
        if self._model is None:
            return NotImplemented
        return self.value ** val

    def __rpow__(self, val):
        if self._model is None:
            return NotImplemented
        return val ** self.value

    def __div__(self, val):
        if self._model is None:
            return NotImplemented
        return self.value / val

    def __rdiv__(self, val):
        if self._model is None:
            return NotImplemented
        return val / self.value

    def __truediv__(self, val):
        if self._model is None:
            return NotImplemented
        return self.value / val

    def __rtruediv__(self, val):
        if self._model is None:
            return NotImplemented
        return val / self.value

    def __eq__(self, val):
        if self._model is None:
            return NotImplemented

        return self.__array__() == val

    def __ne__(self, val):
        if self._model is None:
            return NotImplemented

        return self.__array__() != val

    def __lt__(self, val):
        # Because OrderedDescriptor uses __lt__ to work, we need to call the
        # super method, but only when not bound to an instance anyways
        if self._model is None:
            return super(Parameter, self).__lt__(val)

        return self.__array__() < val

    def __gt__(self, val):
        if self._model is None:
            return NotImplemented

        return self.__array__() > val

    def __le__(self, val):
        if self._model is None:
            return NotImplemented

        return self.__array__() <= val

    def __ge__(self, val):
        if self._model is None:
            return NotImplemented

        return self.__array__() >= val

    def __neg__(self):
        if self._model is None:
            return NotImplemented

        return -self.value

    def __abs__(self):
        if self._model is None:
            return NotImplemented

        return np.abs(self.value)
