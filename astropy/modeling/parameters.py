# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines two classes that deal with parameters.

It is unlikely users will need to work with these classes directly, unless they
define their own models.
"""

from __future__ import division

import bisect
import inspect
import functools
import numbers

import numpy as np

from ..utils import isiterable, lazyproperty


__all__ = ['Parameters', 'Parameter', 'InputParameterError']


class InputParameterError(Exception):
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


# TODO: Is there really any need for this special exception over just, say, a
# ValueError? Or maybe this should be a subclass of ValueError?
class InputParameterError(Exception):
    """
    Called when there's a problem with input parameters.
    """


class Parameter(object):
    """
    Wraps individual parameters.

    This class represents a model's parameter (in a somewhat broad sense). To
    support multiple parameter sets, a parameter has a dimension (dim).
    Parameter objects behave like numbers.

    Parameters
    ----------
    name : string
        parameter name
    val :  number or an iterable of numbers
    model : object
        an instance of a Model class
    dim : int
        parameter dimension
    fixed: boolean
        if True the parameter is not varied during fitting
    tied: callable or False
        if callable is suplplied it provides a way to link to another parameter
    min: float
        the lower bound of a parameter
    max: float
        the upper bound of a parameter
    """

    # See the _nextid classmethod
    _nextid = 1

    def __init__(self, name, getter=None, setter=None, fixed=False,
                 tied=False, min=None, max=None, model=None):
        super(Parameter, self).__init__()
        self._name = name
        self._attr = '_' + name

        self._default_fixed = fixed
        self._default_tied = tied
        self._default_min = min
        self._default_max = max

        self._order = None

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
            _, self._shape = self._validate_value(model, self.value)
        else:
            # Only Parameters declared as class-level descriptors require
            # and ordering ID
            self._order = self._get_nextid()


    def __get__(self, obj, objtype):
        return self.__class__(self._name, getter=self._getter,
                              setter=self._setter, model=obj)

    def __set__(self, obj, value):
        value, shape = self._validate_value(obj, value)
        # Compare the shape against the previous value's shape, if it exists
        if hasattr(obj, self._attr):
            oldvalue = getattr(obj, self._attr)
            if self._getter is not None:
                getter = self._create_value_wrapper(self._getter, obj)
                oldvalue = getter(oldvalue)
            oldvalue, oldshape = self._validate_value(obj, oldvalue)
            if shape != oldshape:
                raise InputParameterError(
                    "Input value for parameter {0!r} does not have the "
                    "required shape {1}".format(self.name, oldshape))

        if self._setter is not None:
            setter = self._create_value_wrapper(self._setter, obj)
            value = setter(value)

        if (hasattr(obj, '_parameters') and
                isinstance(obj._parameters, Parameters)):
            # Model instance has a Parameters list (in general this means it's
            # a ParametricModel)
            setattr(obj, self._attr, value)
            # Rebuild the Parameters list
            # TODO: Is this really necessary? Could we get away with just
            # updating this parameter's value in the Parameters list?
            obj._parameters = Parameters(obj)
        else:
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
    def value(self):
        if self._model is not None:
            if not hasattr(self._model, self._attr):
                raise AttributeError(
                    'Parameter value for {0!r} not set'.format(self._name))
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

        if isinstance(self.value, np.ndarray):
            return self.value.size
        else:
            # A scalar value
            return 1

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
        descriptors delcared at the class-level of Model subclasses.

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
            # If there are more parameter dimensions the value should
            # be a sequence with at least one item
            try:
                if len(value) != param_dim:
                    raise InputParameterError(
                        "Expected parameter {0!r} to be of dimension "
                        "{1}".format(self.name, param_dim))
                # Validate each value
                values = [_tofloat(v) for v in value]
            except (TypeError, IndexError):
                raise InputParameterError(
                    "Expected a multivalued input of dimension {0} "
                    "for parameter {1!r}".format(param_dim, self.name))

            # Check that the value for each dimension has the same shape
            shapes = set(v[1] for v in values)
            if len(shapes) != 1:
                raise InputParameterError(
                    "The value for parameter {0!r} does not have the same "
                    "shape for every dimension".format(self.name))

            # Return the value for each dimension as a list, along with the
            # shape
            return np.array([v[0] for v in values]), shapes.pop()

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
        return np.asarray(self) ** val

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


class Parameters(object):
    """
    View model parameters as a flat list of floats.

    This is a sequence object which provides a view of model parameters. Only
    instances of `~astropy.modeling.core.ParametricModel` keep an instance of
    this class as an attribute. The list of parameters can be modified by the
    user or by an instance of `~astropy.modeling.fitting.Fitter`.

    This is a list-like object which stores model parameters. Only  instances
    of `~astropy.modeling.core.ParametricModel` keep an instance of this class
    as an attribute. The list of parameters can be modified by the user or by
    an instance of `~astropy.modeling.fitting.Fitter`.  This list of parameters
    is kept in sync with single model parameter attributes.  When more than one
    dimensional, a `~astropy.modeling.fitting.Fitter` treats each set of
    parameters as belonging to the same model but different set of data.

    Parameters
    ----------
    model : object
        an instance of a subclass of `~astropy.modeling.core.ParametricModel`
    """

    def __init__(self, model):
        self._model = model

    def __repr__(self):
        return repr(self._flatten())

    def __len__(self):
        return self.size

    def __getitem__(self, key):
        if isinstance(key, slice):
            # TODO: For now slices are handled by generating the full flattened
            # list and slicing it, but this can be made much more efficient
            return self._flatten()[key]
        else:
            if key >= self.size:
                raise IndexError
            starts_idx = bisect.bisect_left(self.starts, (key, ''))
            if starts_idx == len(self.starts):
                starts_idx -= 1
            param_start, param_name = self.starts[starts_idx]

            if param_start > key:
                # We want the previous parameter bin
                param_start, param_name = self.starts[starts_idx - 1]

            param = getattr(self._model, param_name)
            if param.size == 1:
                # Quick out for scalar parameter values
                return param.value

            # key is the index of the virtual flattened array; here we want to
            # convert it into the correct index for this parameter's array
            param_idx = key - param_start
            return getattr(self._model, param_name).value.flat[param_idx]

    def __setitem__(self, key, value):
        if isinstance(key, slice):
            for idx, val in zip(range(*key.indices(len(self))), value):
                self.__setitem__(idx, val)
        else:
            if key >= self.size:
                raise IndexError
            starts_idx = bisect.bisect_left(self.starts, (key, ''))
            if starts_idx == len(self.starts):
                starts_idx -= 1
            param_start, param_name = self.starts[starts_idx]
            if param_start > key:
                # We want the previous parameter bin
                param_start, param_name = self.starts[starts_idx - 1]

            param = getattr(self._model, param_name)

            if param.size == 1:
                setattr(self._model, param_name, value)
            else:
                param_idx = key - param_start
                param.value.flat[param_idx] = value

    def __eq__(self, other):
        """Implement list equality."""

        return self[:] == other

    def __ne__(self, other):
        """Implement list inequality."""

        return self[:] != other

    @lazyproperty
    def size(self):
        """The number of items in the flattened array."""

        # This is just the end index of the last parameter
        last = self._model.param_names[-1]
        return self.slices[last].stop

    @lazyproperty
    def slices(self):
        """The start and stop indces in the flattened array for each parameter.

        Each parameter can have multiple items in the flattened array depending
        on the number of parameter dimensions and the shape of each parameter.
        This number is already fixed once a model has been created, so this
        mapping only needs to be generated once.
        """

        slices = {}
        index = 0

        for name in self._model.param_names:
            parameter = getattr(self._model, name)
            size = parameter.size
            slices[parameter.name] = slice(index, index + size)
            index += size

        return slices

    @lazyproperty
    def starts(self):
        """A list of tuples pairing the start index of a parameter with
        the parameter name.

        This is most useful for mapping some index in the "flattened" parameter
        value list to a specific parameter (see `__getitem__` for example).
        """

        return sorted((slc.start, name) for name, slc in self.slices.items())

    def _is_same_length(self, newparams):
        """
        Checks if the user supplied value of
        `~astropy.modeling.core.ParametricModel.parameters` has the same length
        as the original parameters list.
        """

        # TODO: Would len(newpars) == len(parsize) really not be sufficient?
        parsize = _tofloat(newparams)[0].size
        return parsize == len(self)

    def _flatten(self):
        """
        Create a flat list of model parameters.
        """

        params = []
        for name in self._model.param_names:
            value = getattr(self._model, name).value
            if isinstance(value, np.ndarray):
                params.extend(value.flatten())
            else:
                params.append(value)

        return params
