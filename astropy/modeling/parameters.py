# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name

"""
This module defines classes that deal with parameters.

It is unlikely users will need to work with these classes directly,
unless they define their own models.
"""


import functools
import numbers
import operator

import numpy as np

from astropy.units import MagUnit, Quantity
from astropy.utils import isiterable
from astropy.utils.compat import COPY_IF_NEEDED

from .utils import array_repr_oneline, get_inputs_and_params

__all__ = ["Parameter", "InputParameterError", "ParameterError"]


class ParameterError(Exception):
    """Generic exception class for all exceptions pertaining to Parameters."""


class InputParameterError(ValueError, ParameterError):
    """Used for incorrect input parameter values and definitions."""


class ParameterDefinitionError(ParameterError):
    """Exception in declaration of class-level Parameters."""


def _tofloat(value):
    """Convert a parameter to float or float array."""
    if isiterable(value):
        try:
            value = np.asanyarray(value, dtype=float)
        except (TypeError, ValueError):
            # catch arrays with strings or user errors like different
            # types of parameters in a parameter set
            raise InputParameterError(
                f"Parameter of {type(value)} could not be converted to float"
            )
    elif isinstance(value, Quantity):
        # Quantities are fine as is
        pass
    elif isinstance(value, np.ndarray):
        # A scalar/dimensionless array
        value = float(value.item())
    elif isinstance(value, (numbers.Number, np.number)) and not isinstance(value, bool):
        value = float(value)
    elif isinstance(value, bool):
        raise InputParameterError(
            "Expected parameter to be of numerical type, not boolean"
        )
    else:
        raise InputParameterError(
            f"Don't know how to convert parameter of {type(value)} to float"
        )
    return value


# Helpers for implementing operator overloading on Parameter


def _binary_arithmetic_operation(op, reflected=False):
    @functools.wraps(op)
    def wrapper(self, val):
        if self.unit is not None:
            self_value = Quantity(self.value, self.unit)
        else:
            self_value = self.value

        if reflected:
            return op(val, self_value)
        else:
            return op(self_value, val)

    return wrapper


def _binary_comparison_operation(op):
    @functools.wraps(op)
    def wrapper(self, val):
        if self.unit is not None:
            self_value = Quantity(self.value, self.unit)
        else:
            self_value = self.value

        return op(self_value, val)

    return wrapper


def _unary_arithmetic_operation(op):
    @functools.wraps(op)
    def wrapper(self):
        if self.unit is not None:
            self_value = Quantity(self.value, self.unit)
        else:
            self_value = self.value

        return op(self_value)

    return wrapper


class Parameter:
    """
    Wraps individual parameters.

    Since 4.0 Parameters are no longer descriptors and are based on a new
    implementation of the Parameter class. Parameters now  (as of 4.0) store
    values locally (as instead previously in the associated model)

    This class represents a model's parameter (in a somewhat broad sense). It
    serves a number of purposes:

    1) A type to be recognized by models and treated specially at class
    initialization (i.e., if it is found that there is a class definition
    of a Parameter, the model initializer makes a copy at the instance level).

    2) Managing the handling of allowable parameter values and once defined,
    ensuring updates are consistent with the Parameter definition. This
    includes the optional use of units and quantities as well as transforming
    values to an internally consistent representation (e.g., from degrees to
    radians through the use of getters and setters).

    3) Holding attributes of parameters relevant to fitting, such as whether
    the parameter may be varied in fitting, or whether there are constraints
    that must be satisfied.



    See :ref:`astropy:modeling-parameters` for more details.

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
    unit : `~astropy.units.Unit`
        if specified, the parameter will be in these units, and when the
        parameter is updated in future, it should be set to a
        :class:`~astropy.units.Quantity` that has equivalent units.
    getter : callable or `None`, optional
        A function that wraps the raw (internal) value of the parameter
        when returning the value through the parameter proxy (e.g., a
        parameter may be stored internally as radians but returned to
        the user as degrees). The internal value is what is used for
        computations while the proxy value is what users will interact
        with (passing and viewing). If ``getter`` is not `None`, then a
        ``setter`` must also be input.
    setter : callable or `None`, optional
        A function that wraps any values assigned to this parameter; should
        be the inverse of ``getter``.  If ``setter`` is not `None`, then a
        ``getter`` must also be input.
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
    mag : bool
        Specify if the unit of the parameter can be a Magnitude unit or not
    """

    constraints = ("fixed", "tied", "bounds")
    """
    Types of constraints a parameter can have.  Excludes 'min' and 'max'
    which are just aliases for the first and second elements of the 'bounds'
    constraint (which is represented as a 2-tuple). 'prior' and 'posterior'
    are available for use by user fitters but are not used by any built-in
    fitters as of this writing.
    """

    def __init__(
        self,
        name="",
        description="",
        default=None,
        unit=None,
        getter=None,
        setter=None,
        fixed=False,
        tied=False,
        min=None,
        max=None,
        bounds=None,
        prior=None,
        posterior=None,
        mag=False,
    ):
        super().__init__()

        self._model = None
        self._model_required = False

        if (setter is not None and getter is None) or (
            getter is not None and setter is None
        ):
            raise ValueError("setter and getter must both be input")
        self._setter = self._create_value_wrapper(setter, None)
        self._getter = self._create_value_wrapper(getter, None)
        self._name = name
        self.__doc__ = self._description = description.strip()

        # We only need to perform this check on unbound parameters
        if isinstance(default, Quantity):
            if unit is not None and not unit.is_equivalent(default.unit):
                raise ParameterDefinitionError(
                    f"parameter default {default} does not have units equivalent to "
                    f"the required unit {unit}"
                )
            unit = default.unit
            default = default.value

        self._default = default

        self._mag = mag
        self._set_unit(unit, force=True)
        # Internal units correspond to raw_units held by the model in the
        # previous implementation. The private _getter and _setter methods
        # use this to convert to and from the public unit defined for the
        # parameter.
        self._internal_unit = None
        if not self._model_required:
            if self._default is not None:
                self.value = self._default
            else:
                self._value = None

        # NOTE: These are *default* constraints--on model instances constraints
        # are taken from the model if set, otherwise the defaults set here are
        # used
        if bounds is not None:
            if min is not None or max is not None:
                raise ValueError(
                    "bounds may not be specified simultaneously with min or "
                    f"max when instantiating Parameter {name}"
                )
        else:
            bounds = (min, max)

        self._fixed = fixed
        self._tied = tied
        self._bounds = bounds
        self._order = None

        self._validator = None
        self._prior = prior
        self._posterior = posterior

        self._std = None

    def __set_name__(self, owner, name):
        self._name = name

    def __len__(self):
        val = self.value
        if val.shape == ():
            return 1
        else:
            return val.shape[0]

    def __getitem__(self, key):
        value = self.value
        if len(value.shape) == 0:
            # Wrap the value in a list so that getitem can work for sensible
            # indices like [0] and [-1]
            value = [value]
        return value[key]

    def __setitem__(self, key, value):
        # Get the existing value and check whether it even makes sense to
        # apply this index
        oldvalue = self.value
        if isinstance(key, slice):
            if len(oldvalue[key]) == 0:
                raise InputParameterError(
                    "Slice assignment outside the parameter dimensions for "
                    f"'{self.name}'"
                )
            for idx, val in zip(range(*key.indices(len(self))), value):
                self.__setitem__(idx, val)
        else:
            try:
                oldvalue[key] = value
            except IndexError:
                raise InputParameterError(
                    f"Input dimension {key} invalid for {self.name!r} parameter with "
                    f"dimension {value.shape[0]}"
                )  # likely wrong

    def __repr__(self):
        args = f"'{self._name}'"
        args += f", value={self.value}"

        if self.unit is not None:
            args += f", unit={self.unit}"

        for cons in self.constraints:
            val = getattr(self, cons)
            if val not in (None, False, (None, None)):
                # Maybe non-obvious, but False is the default for the fixed and
                # tied constraints
                args += f", {cons}={val}"

        return f"{self.__class__.__name__}({args})"

    @property
    def name(self):
        """Parameter name."""
        return self._name

    @property
    def default(self):
        """Parameter default value."""
        return self._default

    @property
    def value(self):
        """The unadorned value proxied by this parameter."""
        if self._getter is None and self._setter is None:
            value = self._value
        else:
            # This new implementation uses the names of internal_unit
            # in place of raw_unit used previously. The contrast between
            # internal values and units is that between the public
            # units that the parameter advertises to what it actually
            # uses internally.
            if self.internal_unit:
                value = self._getter(
                    self._internal_value, self.internal_unit, self.unit
                ).value
            else:
                value = self._getter(self._internal_value)

        if value.size == 1:
            # return scalar number as np.float64 object
            return np.float64(value.item())

        return np.float64(value)

    @value.setter
    def value(self, value):
        if isinstance(value, Quantity):
            raise TypeError(
                "The .value property on parameters should be set"
                " to unitless values, not Quantity objects. To set"
                "a parameter to a quantity simply set the "
                "parameter directly without using .value"
            )
        if self._setter is None:
            self._value = np.array(value, dtype=np.float64)
        else:
            self._internal_value = np.array(self._setter(value), dtype=np.float64)

    @property
    def unit(self):
        """
        The unit attached to this parameter, if any.

        On unbound parameters (i.e. parameters accessed through the
        model class, rather than a model instance) this is the required/
        default unit for the parameter.
        """
        return self._unit

    @unit.setter
    def unit(self, unit):
        if self.unit is None:
            raise ValueError(
                "Cannot attach units to parameters that were "
                "not initially specified with units"
            )
        else:
            raise ValueError(
                "Cannot change the unit attribute directly, "
                "instead change the parameter to a new quantity"
            )

    def _set_unit(self, unit, force=False):
        if force:
            if isinstance(unit, MagUnit) and not self._mag:
                raise ValueError(
                    "This parameter does not support the magnitude units such as"
                    f" {unit}"
                )
            self._unit = unit
        else:
            self.unit = unit

    @property
    def internal_unit(self):
        """
        Return the internal unit the parameter uses for the internal value stored.
        """
        return self._internal_unit

    @internal_unit.setter
    def internal_unit(self, internal_unit):
        """
        Set the unit the parameter will convert the supplied value to the
        representation used internally.
        """
        self._internal_unit = internal_unit

    @property
    def input_unit(self):
        """Unit for the input value."""
        if self.internal_unit is not None:
            return self.internal_unit
        elif self.unit is not None:
            return self.unit
        else:
            return None

    @property
    def quantity(self):
        """
        This parameter, as a :class:`~astropy.units.Quantity` instance.
        """
        if self.unit is None:
            return None
        return self.value * self.unit

    @quantity.setter
    def quantity(self, quantity):
        if not isinstance(quantity, Quantity):
            raise TypeError(
                "The .quantity attribute should be set to a Quantity object"
            )
        self.value = quantity.value
        self._set_unit(quantity.unit, force=True)

    @property
    def shape(self):
        """The shape of this parameter's value array."""
        if self._setter is None:
            return self._value.shape
        return self._internal_value.shape

    @shape.setter
    def shape(self, value):
        if isinstance(self.value, np.generic):
            if value not in ((), (1,)):
                raise ValueError("Cannot assign this shape to a scalar quantity")
        else:
            self.value.shape = value

    @property
    def size(self):
        """The size of this parameter's value array."""
        return np.size(self.value)

    @property
    def std(self):
        """Standard deviation, if available from fit."""
        return self._std

    @std.setter
    def std(self, value):
        self._std = value

    @property
    def prior(self):
        return self._prior

    @prior.setter
    def prior(self, val):
        self._prior = val

    @property
    def posterior(self):
        return self._posterior

    @posterior.setter
    def posterior(self, val):
        self._posterior = val

    @property
    def fixed(self):
        """
        Boolean indicating if the parameter is kept fixed during fitting.
        """
        return self._fixed

    @fixed.setter
    def fixed(self, value):
        """Fix a parameter."""
        if not isinstance(value, bool):
            raise ValueError("Value must be boolean")
        self._fixed = value

    @property
    def tied(self):
        """
        Indicates that this parameter is linked to another one.

        A callable which provides the relationship of the two parameters.
        """
        return self._tied

    @tied.setter
    def tied(self, value):
        """Tie a parameter."""
        if not callable(value) and value not in (False, None):
            raise TypeError("Tied must be a callable or set to False or None")
        self._tied = value

    @property
    def bounds(self):
        """The minimum and maximum values of a parameter as a tuple."""
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        """Set the minimum and maximum values of a parameter from a tuple."""
        _min, _max = value
        if _min is not None:
            if not isinstance(_min, (numbers.Number, Quantity)):
                raise TypeError("Min value must be a number or a Quantity")
            if isinstance(_min, Quantity):
                _min = float(_min.value)
            else:
                _min = float(_min)

        if _max is not None:
            if not isinstance(_max, (numbers.Number, Quantity)):
                raise TypeError("Max value must be a number or a Quantity")
            if isinstance(_max, Quantity):
                _max = float(_max.value)
            else:
                _max = float(_max)

        self._bounds = (_min, _max)

    @property
    def min(self):
        """A value used as a lower bound when fitting a parameter."""
        return self.bounds[0]

    @min.setter
    def min(self, value):
        """Set a minimum value of a parameter."""
        self.bounds = (value, self.max)

    @property
    def max(self):
        """A value used as an upper bound when fitting a parameter."""
        return self.bounds[1]

    @max.setter
    def max(self, value):
        """Set a maximum value of a parameter."""
        self.bounds = (self.min, value)

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

        Note: Using this method as a decorator will cause problems with
        pickling the model. An alternative is to assign the actual validator
        function to ``Parameter._validator`` (see examples in modeling).

        """

        def validator(func, self=self):
            if callable(func):
                self._validator = func
                return self
            else:
                raise ValueError(
                    "This decorator method expects a callable.\n"
                    "The use of this method as a direct validator is\n"
                    "deprecated; use the new validate method instead\n"
                )

        return validator

    def validate(self, value):
        """Run the validator on this parameter."""
        if self._validator is not None and self._model is not None:
            self._validator(self._model, value)

    def copy(
        self,
        name=None,
        description=None,
        default=None,
        unit=None,
        getter=None,
        setter=None,
        fixed=False,
        tied=False,
        min=None,
        max=None,
        bounds=None,
        prior=None,
        posterior=None,
    ):
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
        del kwargs["self"]

        for key, value in kwargs.items():
            if value is None:
                # Annoying special cases for min/max where are just aliases for
                # the components of bounds
                if key in ("min", "max"):
                    continue
                else:
                    if hasattr(self, key):
                        value = getattr(self, key)
                    elif hasattr(self, "_" + key):
                        value = getattr(self, "_" + key)
                kwargs[key] = value

        return self.__class__(**kwargs)

    @property
    def model(self):
        """Return the model this  parameter is associated with."""
        return self._model

    @model.setter
    def model(self, value):
        self._model = value
        self._setter = self._create_value_wrapper(self._setter, value)
        self._getter = self._create_value_wrapper(self._getter, value)
        if self._model_required:
            if self._default is not None:
                self.value = self._default
            else:
                self._value = None

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
        if self._setter:
            return self._internal_value
        return self.value

    def _create_value_wrapper(self, wrapper, model):
        """Wraps a getter/setter function to support optionally passing in
        a reference to the model object as the second argument.
        If a model is tied to this parameter and its getter/setter supports
        a second argument then this creates a partial function using the model
        instance as the second argument.
        """
        if isinstance(wrapper, np.ufunc):
            if wrapper.nin != 1:
                raise TypeError(
                    "A numpy.ufunc used for Parameter "
                    "getter/setter may only take one input "
                    "argument"
                )

            return _wrap_ufunc(wrapper)
        elif wrapper is None:
            # Just allow non-wrappers to fall through silently, for convenience
            return None
        else:
            inputs, _ = get_inputs_and_params(wrapper)
            nargs = len(inputs)

            if nargs == 1:
                pass
            elif nargs == 2:
                self._model_required = True
                if model is not None:
                    # Don't make a partial function unless we're tied to a
                    # specific model instance
                    model_arg = inputs[1].name
                    wrapper = functools.partial(wrapper, **{model_arg: model})
            else:
                raise TypeError(
                    "Parameter getter/setter must be a function "
                    "of either one or two arguments"
                )

        return wrapper

    def __array__(self, dtype=None, copy=COPY_IF_NEEDED):
        # Make np.asarray(self) work a little more straightforwardly
        arr = np.asarray(self.value, dtype=dtype)

        if self.unit is not None:
            arr = Quantity(arr, self.unit, copy=copy, subok=True)

        return arr

    def __bool__(self):
        return bool(np.all(self.value))

    __add__ = _binary_arithmetic_operation(operator.add)
    __radd__ = _binary_arithmetic_operation(operator.add, reflected=True)
    __sub__ = _binary_arithmetic_operation(operator.sub)
    __rsub__ = _binary_arithmetic_operation(operator.sub, reflected=True)
    __mul__ = _binary_arithmetic_operation(operator.mul)
    __rmul__ = _binary_arithmetic_operation(operator.mul, reflected=True)
    __pow__ = _binary_arithmetic_operation(operator.pow)
    __rpow__ = _binary_arithmetic_operation(operator.pow, reflected=True)
    __truediv__ = _binary_arithmetic_operation(operator.truediv)
    __rtruediv__ = _binary_arithmetic_operation(operator.truediv, reflected=True)
    __eq__ = _binary_comparison_operation(operator.eq)
    __ne__ = _binary_comparison_operation(operator.ne)
    __lt__ = _binary_comparison_operation(operator.lt)
    __gt__ = _binary_comparison_operation(operator.gt)
    __le__ = _binary_comparison_operation(operator.le)
    __ge__ = _binary_comparison_operation(operator.ge)
    __neg__ = _unary_arithmetic_operation(operator.neg)
    __abs__ = _unary_arithmetic_operation(operator.abs)


def param_repr_oneline(param):
    """
    Like array_repr_oneline but works on `Parameter` objects and supports
    rendering parameters with units like quantities.
    """
    out = array_repr_oneline(param.value)
    if param.unit is not None:
        out = f"{out} {param.unit!s}"
    return out


def _wrap_ufunc(ufunc):
    def _wrapper(value, raw_unit=None, orig_unit=None):
        """
        Wrap ufuncs to support passing in units
            raw_unit is the unit of the value
            orig_unit is the value after the ufunc has been applied
            it is assumed ufunc(raw_unit) == orig_unit
        """
        if orig_unit is not None:
            return ufunc(value) * orig_unit
        elif raw_unit is not None:
            return ufunc(value * raw_unit)

        return ufunc(value)

    return _wrapper
