# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines two classes that deal with parameters.

It is unlikely users will need to work with these classes directly, unless they
define their own models.
"""

from __future__ import division

import inspect
import functools
import numbers

import numpy as np

from ..utils import isiterable


__all__ = ['Parameter', 'InputParameterError']


class InputParameterError(ValueError):
    """Used for incorrect input parameter values and definitions."""


def _tofloat(value):
    """Convert a parameter to float or float array"""

    if isiterable(value):
        try:
            value = np.array(value, dtype=np.float)
            shape = value.shape
        except (TypeError, ValueError):
            # catch arrays with strings or user errors like different
            # types of parameters in a parameter set
            raise InputParameterError(
                "Parameter of {0} could not be converted to "
                "float".format(type(value)))
    elif isinstance(value, bool):
        raise InputParameterError(
            "Expected parameter to be of numerical type, not boolean")
    elif isinstance(value, (numbers.Number, np.number)):
        value = float(value)
        shape = ()
    else:
        raise InputParameterError(
            "Don't know how to convert parameter of {0} to "
            "float".format(type(value)))
    return value, shape


class Parameter(object):
    """
    Wraps individual parameters.

    This class represents a model's parameter (in a somewhat broad sense).  It
    acts as both a descriptor that can be assigned to a class attribute to
    describe the parameters accepted by an individual model (this is called an
    "unbound parameter"), or it can act as a proxy for the parameter values on
    an individual model instance (called a "bound parameter").

    Parameter instances never store the actual value of the parameter
    directly.  Rather, each instance of a model stores its own parameters
    as either hidden attributes or (in the case of
    `~astropy.modeling.core.ParametricModel`) in an array.  A *bound*
    Parameter simply wraps the value in a Parameter proxy which provides some
    additional information about the parameter such as its constraints.

    *Unbound* Parameters are not associated with any specific model instance,
    and are merely used by model classes to determine the names of their
    parameters and other information about each parameter such as their default
    values and default constraints.

    Parameters
    ----------
    name : str
        parameter name
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
    model : object
        an instance of a Model class; this should only be used internally for
        creating bound Parameters
    """

    # See the _nextid classmethod
    _nextid = 1

    def __init__(self, name, description='', default=None, getter=None,
                 setter=None, fixed=False, tied=False, min=None, max=None,
                 model=None):
        super(Parameter, self).__init__()
        self._name = name
        self.__doc__ = description.strip()
        self._default = default
        self._attr = '_' + name

        self._default_fixed = fixed
        self._default_tied = tied
        self._default_min = min
        self._default_max = max

        self._order = None
        self._shape = None
        self._model = model

        # The getter/setter functions take one or two arguments: The first
        # argument is always the value itself (either the value returned or the
        # value being set).  The second argument is optional, but if present
        # will contain a reference to the model object tied to a parameter (if
        # it exists)
        if getter is not None:
            self._getter = self._create_value_wrapper(getter, model)
        else:
            self._getter = None
        if setter is not None:
            self._setter = self._create_value_wrapper(setter, model)
        else:
            self._setter = None

        if model is not None:
            try:
                _, self._shape = self._validate_value(model, self.value)
            except AttributeError:
                # This can happen if the paramter's value has not been set yet
                pass
        else:
            # Only Parameters declared as class-level descriptors require
            # and ordering ID
            self._order = self._get_nextid()


    def __get__(self, obj, objtype):
        if obj is None:
            return self

        return self.__class__(self._name, default=self._default,
                              getter=self._getter,
                              setter=self._setter, model=obj)

    def __set__(self, obj, value):
        value, shape = self._validate_value(obj, value)
        # Compare the shape against the previous value's shape, if it exists
        if hasattr(obj, self._attr):
            current_shape = getattr(obj, self.name).shape
            if shape != current_shape:
                raise InputParameterError(
                    "Input value for parameter {0!r} does not have the "
                    "required shape {1}".format(self.name, current_shape))

        if self._setter is not None:
            setter = self._create_value_wrapper(self._setter, obj)
            value = setter(value)

        setattr(obj, self._attr, value)

    def __len__(self):
        if self._model is None:
            raise TypeError('Parameter definitions do not have a length.')
        return self._model.param_dim

    def __getitem__(self, key):
        value = self.value
        if self._model.param_dim == 1:
            # Wrap the value in a list so that getitem can work for sensible
            # indcies like [0] and [-1]
            value = [value]
        return value[key]

    def __setitem__(self, key, value):
        # Get the existing value and check whether it even makes sense to
        # apply this index
        oldvalue = self.value
        param_dim = self._model.param_dim

        if param_dim == 1:
            # Convert the single-dimension value to a list to allow some slices
            # that would be compatible with a length-1 array like [:] and [0:]
            oldvalue = [oldvalue]

        if isinstance(key, slice):
            if len(oldvalue[key]) == 0:
                raise InputParameterError(
                    "Slice assignment outside the parameter dimensions for "
                    "{0!r}".format(self.name))
            for idx, val in zip(range(*key.indices(len(self))), value):
                self.__setitem__(idx, val)
        else:
            try:
                oldvalue[key] = value
                if param_dim == 1:
                    self.value = value
            except IndexError:
                raise InputParameterError(
                    "Input dimension {0} invalid for {1!r} parameter with "
                    "dimension {2}".format(key, self.name, param_dim))

    def __repr__(self):
        if self._model is None:
            return 'Parameter({0!r})'.format(self._name)
        else:
            return 'Parameter({0!r}, value={1!r})'.format(
                self._name, self.value)

    @property
    def name(self):
        """Parameter name"""

        return self._name

    @property
    def default(self):
        """Parameter default value"""

        if self._model is None:
            return self._default

        if self._model.param_dim == 1:
            return self._default
        else:
            return np.repeat(self._default, self._model.param_dim)

    @property
    def value(self):
        """The unadorned value proxied by this parameter"""

        if self._model is not None:
            if not hasattr(self._model, self._attr):
                if self._default is not None:
                    value = self.default
                else:
                    raise AttributeError(
                        'Parameter value for {0!r} not set'.format(self._name))
            else:
                value = getattr(self._model, self._attr)
            if self._getter is None:
                return value
            else:
                return self._getter(value)
        raise AttributeError('Parameter definition does not have a value')

    @value.setter
    def value(self, val):
        if self._model is not None:
            if self._setter is not None:
                val = self._setter(val)
            setattr(self._model, self._attr, val)
        raise AttributeError('Cannot set a value on a parameter definition')

    @property
    def shape(self):
        """The shape of this parameter's value array."""

        return self._shape

    @property
    def size(self):
        """The size of this parameter's value array."""

        return np.size(self.value)

    @property
    def fixed(self):
        """
        Boolean indicating if the parameter is kept fixed during fitting.
        """

        if self._model is not None:
            fixed = self._model._constraints.setdefault('fixed', {})
            return fixed.setdefault(self._name, self._default_fixed)
        else:
            return self._default_fixed

    @fixed.setter
    def fixed(self, value):
        """Fix a parameter"""
        if self._model is not None:
            assert isinstance(value, bool), "Fixed can be True or False"
            fixed = self._model._constraints.setdefault('fixed', {})
            fixed[self._name] = value
            self._model._fit_parameters, self._model._fit_param_indices, \
                = self._model._model_to_fit_params()
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
            tied = self._model._constraints.setdefault('tied', {})
            return tied.setdefault(self._name, self._default_tied)
        else:
            return self._default_tied

    @tied.setter
    def tied(self, value):
        """Tie a parameter"""

        if self._model is not None:
            assert callable(value) or value in (False, None), \
                    "Tied must be a callable"
            tied = self._model._constraints.setdefault('tied', {})
            tied[self._name] = value
            self._model._fit_parameters, self._model._fit_param_indices, \
                = self._model._model_to_fit_params()
        else:
            raise AttributeError("can't set attribute 'tied' on Parameter "
                                 "definition")

    @property
    def bounds(self):
        """The minimum and maximum values of a parameter as a tuple"""

        if self._model is not None:
            bounds = self._model._constraints.setdefault('bounds', {})
            return bounds.setdefault(self._name,
                                     (self._default_min, self._default_max))
        else:
            return (self._default_min, self._default_max)

    @bounds.setter
    def bounds(self, value):
        """Set the minimum and maximum values of a parameter from a tuple"""

        if self._model is not None:
            _min, _max = value
            if _min is not None:
                assert isinstance(_min, numbers.Number), \
                        "Min value must be a number"
                _min = float(_min)

            if _max is not None:
                assert isinstance(_max, numbers.Number), \
                        "Max value must be a number"
                _max = float(_max)

            bounds = self._model._constraints.setdefault('bounds', {})
            bounds[self._name] = (_min, _max)
            self._model._model_to_fit_params()
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

    @classmethod
    def _get_nextid(cls):
        """Returns a monotonically increasing ID used to order Parameter
        descriptors declared at the class-level of Model subclasses.

        This allows the desired parameter order to be determined without
        having to list it manually in the param_names class attribute.
        """

        nextid = cls._nextid
        cls._nextid += 1
        return nextid

    def _validate_value(self, model, value):
        if model is None:
            return

        param_dim = model.param_dim
        if param_dim == 1:
            # Just validate the value with _tofloat
            return _tofloat(value)
        else:
            try:
                # Validate each value
                value, shape = _tofloat(value)
            except (TypeError, IndexError):
                raise InputParameterError(
                    "Expected a multivalued input of dimension {0} "
                    "for parameter {1!r}".format(param_dim, self.name))

            return value, shape

    def _create_value_wrapper(self, wrapper, model):
        """Wrappers a getter/setter function to support optionally passing in
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
        else:
            wrapper_args = inspect.getargspec(wrapper)
            nargs = len(wrapper_args.args)

            if nargs == 1:
                pass
            elif nargs == 2:
                if model is not None:
                    # Don't make a partial function unless we're tied to a
                    # specific model instance
                    model_arg = wrapper_args.args[1]
                    wrapper = functools.partial(wrapper, **{model_arg: model})
            else:
                raise TypeError("Parameter getter/setter must be a function "
                                "of either one or two arguments")

        return wrapper

    def __add__(self, val):
        return self.value + val

    def __radd__(self, val):
        return self.value + val

    def __sub__(self, val):
        return self.value - val

    def __rsub__(self, val):
        return val - self.value

    def __mul__(self, val):
        return self.value * val

    def __rmul__(self, val):
        return self.value * val

    def __pow__(self, val):
        return self.value ** val

    def __rpow__(self, val):
        return val ** self.value

    def __div__(self, val):
        return self.value / val

    def __rdiv__(self, val):
        return val / self.value

    def __truediv__(self, val):
        return self.value / val

    def __rtruediv__(self, val):
        return val / self.value

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
        return -self.value

    def __abs__(self):
        return np.abs(self.value)
