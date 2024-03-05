# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines the `Quantity` object, which represents a number with some
associated units. `Quantity` objects support operations like ordinary numbers,
but will deal with unit conversions internally.
"""

# STDLIB
import numbers
import operator
import re
import warnings
from fractions import Fraction

# THIRD PARTY
import numpy as np

# LOCAL
from astropy import config as _config
from astropy.utils.compat.numpycompat import COPY_IF_NEEDED, NUMPY_LT_2_0
from astropy.utils.data_info import ParentDtypeInfo
from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.misc import isiterable

from .core import (
    Unit,
    UnitBase,
    UnitConversionError,
    UnitsError,
    UnitTypeError,
    dimensionless_unscaled,
    get_current_unit_registry,
)
from .format import Base, Latex
from .quantity_helper import can_have_arbitrary_unit, check_output, converters_and_unit
from .quantity_helper.function_helpers import (
    DISPATCHED_FUNCTIONS,
    FUNCTION_HELPERS,
    SUBCLASS_SAFE_FUNCTIONS,
    UNSUPPORTED_FUNCTIONS,
)
from .structured import StructuredUnit, _structured_unit_like_dtype
from .utils import is_effectively_unity

__all__ = [
    "Quantity",
    "SpecificTypeQuantity",
    "QuantityInfoBase",
    "QuantityInfo",
    "allclose",
    "isclose",
]


# We don't want to run doctests in the docstrings we inherit from Numpy
__doctest_skip__ = ["Quantity.*"]

_UNIT_NOT_INITIALISED = "(Unit not initialised)"
_UFUNCS_FILTER_WARNINGS = {np.arcsin, np.arccos, np.arccosh, np.arctanh}


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for Quantity.
    """

    latex_array_threshold = _config.ConfigItem(
        100,
        "The maximum size an array Quantity can be before its LaTeX "
        'representation for IPython gets "summarized" (meaning only the first '
        'and last few elements are shown with "..." between). Setting this to a '
        "negative number means that the value will instead be whatever numpy "
        "gets from get_printoptions.",
    )


conf = Conf()


class QuantityIterator:
    """
    Flat iterator object to iterate over Quantities.

    A `QuantityIterator` iterator is returned by ``q.flat`` for any Quantity
    ``q``.  It allows iterating over the array as if it were a 1-D array,
    either in a for-loop or by calling its `next` method.

    Iteration is done in C-contiguous style, with the last index varying the
    fastest. The iterator can also be indexed using basic slicing or
    advanced indexing.

    See Also
    --------
    Quantity.flatten : Returns a flattened copy of an array.

    Notes
    -----
    `QuantityIterator` is inspired by `~numpy.ma.core.MaskedIterator`.  It
    is not exported by the `~astropy.units` module.  Instead of
    instantiating a `QuantityIterator` directly, use `Quantity.flat`.
    """

    def __init__(self, q):
        self._quantity = q
        self._dataiter = q.view(np.ndarray).flat

    def __iter__(self):
        return self

    def __getitem__(self, indx):
        out = self._dataiter.__getitem__(indx)
        # For single elements, ndarray.flat.__getitem__ returns scalars; these
        # need a new view as a Quantity.
        if isinstance(out, type(self._quantity)):
            return out
        else:
            return self._quantity._new_view(out)

    def __setitem__(self, index, value):
        self._dataiter[index] = self._quantity._to_own_unit(value)

    def __next__(self):
        """
        Return the next value, or raise StopIteration.
        """
        out = next(self._dataiter)
        # ndarray.flat._dataiter returns scalars, so need a view as a Quantity.
        return self._quantity._new_view(out)

    next = __next__

    def __len__(self):
        return len(self._dataiter)

    # Properties and methods to match `numpy.ndarray.flatiter`

    @property
    def base(self):
        """A reference to the array that is iterated over."""
        return self._quantity

    @property
    def coords(self):
        """An N-dimensional tuple of current coordinates."""
        return self._dataiter.coords

    @property
    def index(self):
        """Current flat index into the array."""
        return self._dataiter.index

    def copy(self):
        """Get a copy of the iterator as a 1-D array."""
        return self._quantity.flatten()


class QuantityInfoBase(ParentDtypeInfo):
    # This is on a base class rather than QuantityInfo directly, so that
    # it can be used for EarthLocationInfo yet make clear that that class
    # should not be considered a typical Quantity subclass by Table.
    attrs_from_parent = {"dtype", "unit"}  # dtype and unit taken from parent
    _supports_indexing = True

    @staticmethod
    def default_format(val):
        return f"{val.value}"

    @staticmethod
    def possible_string_format_functions(format_):
        """Iterate through possible string-derived format functions.

        A string can either be a format specifier for the format built-in,
        a new-style format string, or an old-style format string.

        This method is overridden in order to suppress printing the unit
        in each row since it is already at the top in the column header.
        """
        yield lambda format_, val: format(val.value, format_)
        yield lambda format_, val: format_.format(val.value)
        yield lambda format_, val: format_ % val.value


class QuantityInfo(QuantityInfoBase):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """

    _represent_as_dict_attrs = ("value", "unit")
    _construct_from_dict_args = ["value"]
    _represent_as_dict_primary_data = "value"

    def new_like(self, cols, length, metadata_conflicts="warn", name=None):
        """
        Return a new Quantity instance which is consistent with the
        input ``cols`` and has ``length`` rows.

        This is intended for creating an empty column object whose elements can
        be set in-place for table operations like join or vstack.

        Parameters
        ----------
        cols : list
            List of input columns
        length : int
            Length of the output column object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : `~astropy.units.Quantity` (or subclass)
            Empty instance of this class consistent with ``cols``

        """
        # Get merged info attributes like shape, dtype, format, description, etc.
        attrs = self.merge_cols_attributes(
            cols, metadata_conflicts, name, ("meta", "format", "description")
        )

        # Make an empty quantity using the unit of the last one.
        shape = (length,) + attrs.pop("shape")
        dtype = attrs.pop("dtype")
        # Use zeros so we do not get problems for Quantity subclasses such
        # as Longitude and Latitude, which cannot take arbitrary values.
        data = np.zeros(shape=shape, dtype=dtype)
        # Get arguments needed to reconstruct class
        map = {
            key: (data if key == "value" else getattr(cols[-1], key))
            for key in self._represent_as_dict_attrs
        }
        map["copy"] = False
        out = self._construct_from_dict(map)

        # Set remaining info attributes
        for attr, value in attrs.items():
            setattr(out.info, attr, value)

        return out

    def get_sortable_arrays(self):
        """
        Return a list of arrays which can be lexically sorted to represent
        the order of the parent column.

        For Quantity this is just the quantity itself.


        Returns
        -------
        arrays : list of ndarray
        """
        return [self._parent]


class Quantity(np.ndarray):
    """A `~astropy.units.Quantity` represents a number with some associated unit.

    See also: https://docs.astropy.org/en/stable/units/quantity.html

    Parameters
    ----------
    value : number, `~numpy.ndarray`, `~astropy.units.Quantity` (sequence), or str
        The numerical value of this quantity in the units given by unit.  If a
        `Quantity` or sequence of them (or any other valid object with a
        ``unit`` attribute), creates a new `Quantity` object, converting to
        `unit` units as needed.  If a string, it is converted to a number or
        `Quantity`, depending on whether a unit is present.

    unit : unit-like
        An object that represents the unit associated with the input value.
        Must be an `~astropy.units.UnitBase` object or a string parseable by
        the :mod:`~astropy.units` package.

    dtype : ~numpy.dtype, optional
        The dtype of the resulting Numpy array or scalar that will
        hold the value.  If not provided, it is determined from the input,
        except that any integer and (non-Quantity) object inputs are converted
        to float by default.
        If `None`, the normal `numpy.dtype` introspection is used, e.g.
        preventing upcasting of integers.

    copy : bool, optional
        If `True` (default), then the value is copied.  Otherwise, a copy will
        only be made if ``__array__`` returns a copy, if value is a nested
        sequence, or if a copy is needed to satisfy an explicitly given
        ``dtype``.  (The `False` option is intended mostly for internal use,
        to speed up initialization where a copy is known to have been made.
        Use with care.)

    order : {'C', 'F', 'A'}, optional
        Specify the order of the array.  As in `~numpy.array`.  This parameter
        is ignored if the input is a `Quantity` and ``copy=False``.

    subok : bool, optional
        If `False` (default), the returned array will be forced to be a
        `Quantity`.  Otherwise, `Quantity` subclasses will be passed through,
        or a subclass appropriate for the unit will be used (such as
        `~astropy.units.Dex` for ``u.dex(u.AA)``).

    ndmin : int, optional
        Specifies the minimum number of dimensions that the resulting array
        should have.  Ones will be prepended to the shape as needed to meet
        this requirement.  This parameter is ignored if the input is a
        `Quantity` and ``copy=False``.

    Raises
    ------
    TypeError
        If the value provided is not a Python numeric type.
    TypeError
        If the unit provided is not either a :class:`~astropy.units.Unit`
        object or a parseable string unit.

    Notes
    -----
    Quantities can also be created by multiplying a number or array with a
    :class:`~astropy.units.Unit`. See https://docs.astropy.org/en/latest/units/

    Unless the ``dtype`` argument is explicitly specified, integer
    or (non-Quantity) object inputs are converted to `float` by default.
    """

    # Need to set a class-level default for _equivalencies, or
    # Constants can not initialize properly
    _equivalencies = []

    # Default unit for initialization; can be overridden by subclasses,
    # possibly to `None` to indicate there is no default unit.
    _default_unit = dimensionless_unscaled

    # Ensures views have an undefined unit.
    _unit = None

    __array_priority__ = 10000

    def __class_getitem__(cls, unit_shape_dtype):
        """Quantity Type Hints.

        Unit-aware type hints are ``Annotated`` objects that encode the class,
        the unit, and possibly shape and dtype information, depending on the
        python and :mod:`numpy` versions.

        Schematically, ``Annotated[cls[shape, dtype], unit]``

        As a classmethod, the type is the class, ie ``Quantity``
        produces an ``Annotated[Quantity, ...]`` while a subclass
        like :class:`~astropy.coordinates.Angle` returns
        ``Annotated[Angle, ...]``.

        Parameters
        ----------
        unit_shape_dtype : :class:`~astropy.units.UnitBase`, str, `~astropy.units.PhysicalType`, or tuple
            Unit specification, can be the physical type (ie str or class).
            If tuple, then the first element is the unit specification
            and all other elements are for `numpy.ndarray` type annotations.
            Whether they are included depends on the python and :mod:`numpy`
            versions.

        Returns
        -------
        `typing.Annotated`, `astropy.units.Unit`, or `astropy.units.PhysicalType`
            Return type in this preference order:
            * `typing.Annotated`
            * `astropy.units.Unit` or `astropy.units.PhysicalType`

        Raises
        ------
        TypeError
            If the unit/physical_type annotation is not Unit-like or
            PhysicalType-like.

        Examples
        --------
        Create a unit-aware Quantity type annotation

            >>> Quantity[Unit("s")]
            Annotated[Quantity, Unit("s")]

        See Also
        --------
        `~astropy.units.quantity_input`
            Use annotations for unit checks on function arguments and results.

        Notes
        -----
        With Python 3.9+ or :mod:`typing_extensions`, |Quantity| types are also
        static-type compatible.
        """
        from typing import Annotated

        # process whether [unit] or [unit, shape, ptype]
        if isinstance(unit_shape_dtype, tuple):  # unit, shape, dtype
            target = unit_shape_dtype[0]
            shape_dtype = unit_shape_dtype[1:]
        else:  # just unit
            target = unit_shape_dtype
            shape_dtype = ()

        # Allowed unit/physical types. Errors if neither.
        try:
            unit = Unit(target)
        except (TypeError, ValueError):
            from astropy.units.physical import get_physical_type

            try:
                unit = get_physical_type(target)
            except (TypeError, ValueError, KeyError):  # KeyError for Enum
                raise TypeError(
                    "unit annotation is not a Unit or PhysicalType"
                ) from None

        # Quantity does not (yet) properly extend the NumPy generics types,
        # introduced in numpy v1.22+, instead just including the unit info as
        # metadata using Annotated.
        # TODO: ensure we do interact with NDArray.__class_getitem__.
        return Annotated[cls, unit]

    def __new__(
        cls,
        value,
        unit=None,
        dtype=np.inexact,
        copy=True,
        order=None,
        subok=False,
        ndmin=0,
    ):
        if unit is not None:
            # convert unit first, to avoid multiple string->unit conversions
            unit = Unit(unit)

        # inexact -> upcast to float dtype
        float_default = dtype is np.inexact
        if float_default:
            dtype = None

        if not NUMPY_LT_2_0 and copy is False:
            # strict backward compatibility
            # see https://github.com/astropy/astropy/pull/16142
            copy = None

        # optimize speed for Quantity with no dtype given, copy=COPY_IF_NEEDED
        if isinstance(value, Quantity):
            if unit is not None and unit is not value.unit:
                value = value.to(unit)
                # the above already makes a copy (with float dtype)
                copy = COPY_IF_NEEDED

            if type(value) is not cls and not (subok and isinstance(value, cls)):
                value = value.view(cls)

            if float_default and value.dtype.kind in "iu":
                dtype = float

            return np.array(
                value, dtype=dtype, copy=copy, order=order, subok=True, ndmin=ndmin
            )

        # Maybe str, or list/tuple of Quantity? If so, this may set value_unit.
        # To ensure array remains fast, we short-circuit it.
        value_unit = None
        if not isinstance(value, np.ndarray):
            if isinstance(value, str):
                # The first part of the regex string matches any integer/float;
                # the second parts adds possible trailing .+-, which will break
                # the float function below and ensure things like 1.2.3deg
                # will not work.
                pattern = (
                    r"\s*[+-]?"
                    r"((\d+\.?\d*)|(\.\d+)|([nN][aA][nN])|"
                    r"([iI][nN][fF]([iI][nN][iI][tT][yY]){0,1}))"
                    r"([eE][+-]?\d+)?"
                    r"[.+-]?"
                )

                v = re.match(pattern, value)
                unit_string = None
                try:
                    value = float(v.group())

                except Exception:
                    raise TypeError(
                        f'Cannot parse "{value}" as a {cls.__name__}. It does not '
                        "start with a number."
                    )

                unit_string = v.string[v.end() :].strip()
                if unit_string:
                    value_unit = Unit(unit_string)
                    if unit is None:
                        unit = value_unit  # signal no conversion needed below.

            elif isiterable(value) and len(value) > 0:
                # Iterables like lists and tuples.
                if all(isinstance(v, Quantity) for v in value):
                    # If a list/tuple containing only quantities, convert all
                    # to the same unit.
                    if unit is None:
                        unit = value[0].unit
                    value = [q.to_value(unit) for q in value]
                    value_unit = unit  # signal below that conversion has been done
                elif (
                    dtype is None
                    and not hasattr(value, "dtype")
                    and isinstance(unit, StructuredUnit)
                ):
                    # Special case for list/tuple of values and a structured unit:
                    # ``np.array(value, dtype=None)`` would treat tuples as lower
                    # levels of the array, rather than as elements of a structured
                    # array, so we use the structure of the unit to help infer the
                    # structured dtype of the value.
                    dtype = unit._recursively_get_dtype(value)

        using_default_unit = False
        if value_unit is None:
            # If the value has a `unit` attribute and if not None
            # (for Columns with uninitialized unit), treat it like a quantity.
            value_unit = getattr(value, "unit", None)
            if value_unit is None:
                # Default to dimensionless for no (initialized) unit attribute.
                if unit is None:
                    using_default_unit = True
                    unit = cls._default_unit
                value_unit = unit  # signal below that no conversion is needed
            else:
                try:
                    value_unit = Unit(value_unit)
                except Exception as exc:
                    raise TypeError(
                        f"The unit attribute {value.unit!r} of the input could "
                        "not be parsed as an astropy Unit."
                    ) from exc

                if unit is None:
                    unit = value_unit
                elif unit is not value_unit:
                    copy = COPY_IF_NEEDED  # copy will be made in conversion at end

        value = np.array(
            value, dtype=dtype, copy=copy, order=order, subok=True, ndmin=ndmin
        )

        # For no-user-input unit, make sure the constructed unit matches the
        # structure of the data.
        if using_default_unit and value.dtype.names is not None:
            unit = value_unit = _structured_unit_like_dtype(value_unit, value.dtype)

        # check that array contains numbers or long int objects
        if value.dtype.kind in "OSU" and not (
            value.dtype.kind == "O" and isinstance(value.item(0), numbers.Number)
        ):
            raise TypeError("The value must be a valid Python or Numpy numeric type.")

        # by default, cast any integer, boolean, etc., to float
        if float_default and value.dtype.kind in "iuO":
            value = value.astype(float)

        # if we allow subclasses, allow a class from the unit.
        if subok:
            qcls = getattr(unit, "_quantity_class", cls)
            if issubclass(qcls, cls):
                cls = qcls

        value = value.view(cls)
        value._set_unit(value_unit)
        if unit is value_unit:
            return value
        else:
            # here we had non-Quantity input that had a "unit" attribute
            # with a unit different from the desired one.  So, convert.
            return value.to(unit)

    def __array_finalize__(self, obj):
        # Check whether super().__array_finalize should be called
        # (sadly, ndarray.__array_finalize__ is None; we cannot be sure
        # what is above us).
        super_array_finalize = super().__array_finalize__
        if super_array_finalize is not None:
            super_array_finalize(obj)

        # If we're a new object or viewing an ndarray, nothing has to be done.
        if obj is None or obj.__class__ is np.ndarray:
            return

        # Copy over the unit and possibly info.  Note that the only way the
        # unit can already be set is if one enters via _new_view(), where the
        # unit is often different from that of self, and where propagation of
        # info is not always desirable.
        if self._unit is None:
            unit = getattr(obj, "_unit", None)
            if unit is not None:
                self._set_unit(unit)

            # Copy info if the original had `info` defined.  Because of the way the
            # DataInfo works, `'info' in obj.__dict__` is False until the
            # `info` attribute is accessed or set.
            if "info" in obj.__dict__:
                self.info = obj.info

    def __array_wrap__(self, obj, context=None, return_scalar=False):
        if context is None:
            # Methods like .squeeze() created a new `ndarray` and then call
            # __array_wrap__ to turn the array into self's subclass.
            return self._new_view(obj)

        raise NotImplementedError(
            "__array_wrap__ should not be used with a context any more since all "
            "use should go through array_function. Please raise an issue on "
            "https://github.com/astropy/astropy"
        )

    def __array_ufunc__(self, function, method, *inputs, **kwargs):
        """Wrap numpy ufuncs, taking care of units.

        Parameters
        ----------
        function : callable
            ufunc to wrap.
        method : str
            Ufunc method: ``__call__``, ``at``, ``reduce``, etc.
        inputs : tuple
            Input arrays.
        kwargs : keyword arguments
            As passed on, with ``out`` containing possible quantity output.

        Returns
        -------
        result : `~astropy.units.Quantity` or `NotImplemented`
            Results of the ufunc, with the unit set properly.
        """
        # Determine required conversion functions -- to bring the unit of the
        # input to that expected (e.g., radian for np.sin), or to get
        # consistent units between two inputs (e.g., in np.add) --
        # and the unit of the result (or tuple of units for nout > 1).
        try:
            converters, unit = converters_and_unit(function, method, *inputs)

            out = kwargs.get("out", None)
            # Avoid loop back by turning any Quantity output into array views.
            if out is not None:
                # If pre-allocated output is used, check it is suitable.
                # This also returns array view, to ensure we don't loop back.
                if function.nout == 1:
                    out = out[0]
                out_array = check_output(out, unit, inputs, function=function)
                # Ensure output argument remains a tuple.
                kwargs["out"] = (out_array,) if function.nout == 1 else out_array

            if method == "reduce" and "initial" in kwargs and unit is not None:
                # Special-case for initial argument for reductions like
                # np.add.reduce.  This should be converted to the output unit as
                # well, which is typically the same as the input unit (but can
                # in principle be different: unitless for np.equal, radian
                # for np.arctan2, though those are not necessarily useful!)
                kwargs["initial"] = self._to_own_unit(
                    kwargs["initial"], check_precision=False, unit=unit
                )

            # Same for inputs, but here also convert if necessary.
            arrays = []
            for input_, converter in zip(inputs, converters):
                input_ = getattr(input_, "value", input_)
                arrays.append(converter(input_) if converter else input_)

            # Call our superclass's __array_ufunc__
            result = super().__array_ufunc__(function, method, *arrays, **kwargs)
            # If unit is None, a plain array is expected (e.g., comparisons), which
            # means we're done.
            # We're also done if the result was None (for method 'at') or
            # NotImplemented, which can happen if other inputs/outputs override
            # __array_ufunc__; hopefully, they can then deal with us.
            if unit is None or result is None or result is NotImplemented:
                return result

            return self._result_as_quantity(result, unit, out)

        except (TypeError, ValueError, AttributeError) as e:
            out_normalized = kwargs.get("out", tuple())
            inputs_and_outputs = inputs + out_normalized
            ignored_ufunc = (
                None,
                np.ndarray.__array_ufunc__,
                type(self).__array_ufunc__,
            )
            if not all(
                getattr(type(io), "__array_ufunc__", None) in ignored_ufunc
                for io in inputs_and_outputs
            ):
                return NotImplemented
            else:
                raise e

    def _result_as_quantity(self, result, unit, out):
        """Turn result into a quantity with the given unit.

        If no output is given, it will take a view of the array as a quantity,
        and set the unit.  If output is given, those should be quantity views
        of the result arrays, and the function will just set the unit.

        Parameters
        ----------
        result : ndarray or tuple thereof
            Array(s) which need to be turned into quantity.
        unit : `~astropy.units.Unit`
            Unit for the quantities to be returned (or `None` if the result
            should not be a quantity).  Should be tuple if result is a tuple.
        out : `~astropy.units.Quantity` or None
            Possible output quantity. Should be `None` or a tuple if result
            is a tuple.

        Returns
        -------
        out : `~astropy.units.Quantity`
           With units set.
        """
        if isinstance(result, (tuple, list)):
            if out is None:
                out = (None,) * len(result)
            # Some np.linalg functions return namedtuple, which is handy to access
            # elements by name, but cannot be directly initialized with an iterator.
            result_cls = getattr(result, "_make", result.__class__)
            return result_cls(
                self._result_as_quantity(result_, unit_, out_)
                for (result_, unit_, out_) in zip(result, unit, out)
            )

        if out is None:
            # View the result array as a Quantity with the proper unit.
            return (
                result
                if unit is None
                else self._new_view(result, unit, propagate_info=False)
            )

        elif isinstance(out, Quantity):
            # For given Quantity output, just set the unit. We know the unit
            # is not None and the output is of the correct Quantity subclass,
            # as it was passed through check_output.
            # (We cannot do this unconditionally, though, since it is possible
            # for out to be ndarray and the unit to be dimensionless.)
            out._set_unit(unit)

        return out

    def __quantity_subclass__(self, unit):
        """
        Overridden by subclasses to change what kind of view is
        created based on the output unit of an operation.

        Parameters
        ----------
        unit : UnitBase
            The unit for which the appropriate class should be returned

        Returns
        -------
        tuple :
            - `~astropy.units.Quantity` subclass
            - bool: True if subclasses of the given class are ok
        """
        return Quantity, True

    def _new_view(self, obj=None, unit=None, propagate_info=True):
        """Create a Quantity view of some array-like input, and set the unit.

        By default, return a view of ``obj`` of the same class as ``self`` and
        with the same unit.  Subclasses can override the type of class for a
        given unit using ``__quantity_subclass__``, and can ensure properties
        other than the unit are copied using ``__array_finalize__``.

        If the given unit defines a ``_quantity_class`` of which ``self``
        is not an instance, a view using this class is taken.

        Parameters
        ----------
        obj : ndarray or scalar, optional
            The array to create a view of.  If obj is a numpy or python scalar,
            it will be converted to an array scalar.  By default, ``self``
            is converted.

        unit : unit-like, optional
            The unit of the resulting object.  It is used to select a
            subclass, and explicitly assigned to the view if given.
            If not given, the subclass and unit will be that of ``self``.

        propagate_info : bool, optional
            Whether to transfer ``info`` if present.  Default: `True`, as
            appropriate for, e.g., unit conversions or slicing, where the
            nature of the object does not change.

        Returns
        -------
        view : `~astropy.units.Quantity` subclass

        """
        # Determine the unit and quantity subclass that we need for the view.
        if unit is None:
            unit = self.unit
            quantity_subclass = self.__class__
        elif unit is self.unit and self.__class__ is Quantity:
            # The second part is because we should not presume what other
            # classes want to do for the same unit.  E.g., Constant will
            # always want to fall back to Quantity, and relies on going
            # through `__quantity_subclass__`.
            quantity_subclass = Quantity
        else:
            unit = Unit(unit)
            quantity_subclass = getattr(unit, "_quantity_class", Quantity)
            if isinstance(self, quantity_subclass):
                quantity_subclass, subok = self.__quantity_subclass__(unit)
                if subok:
                    quantity_subclass = self.__class__

        # We only want to propagate information from ``self`` to our new view,
        # so obj should be a regular array.  By using ``np.array``, we also
        # convert python and numpy scalars, which cannot be viewed as arrays
        # and thus not as Quantity either, to zero-dimensional arrays.
        # (These are turned back into scalar in `.value`)
        # Note that for an ndarray input, the np.array call takes only double
        # ``obj.__class is np.ndarray``. So, not worth special-casing.
        if obj is None:
            obj = self.view(np.ndarray)
        else:
            obj = np.array(obj, copy=COPY_IF_NEEDED, subok=True)

        # Take the view, set the unit, and update possible other properties
        # such as ``info``, ``wrap_angle`` in `Longitude`, etc.
        view = obj.view(quantity_subclass)
        view._set_unit(unit)
        view.__array_finalize__(self)
        if propagate_info and "info" in self.__dict__:
            view.info = self.info
        return view

    def _set_unit(self, unit):
        """Set the unit.

        This is used anywhere the unit is set or modified, i.e., in the
        initializer, in ``__imul__`` and ``__itruediv__`` for in-place
        multiplication and division by another unit, as well as in
        ``__array_finalize__`` for wrapping up views.  For Quantity, it just
        sets the unit, but subclasses can override it to check that, e.g.,
        a unit is consistent.
        """
        if not isinstance(unit, UnitBase):
            if isinstance(self._unit, StructuredUnit) or isinstance(
                unit, StructuredUnit
            ):
                unit = StructuredUnit(unit, self.dtype)
            else:
                # Trying to go through a string ensures that, e.g., Magnitudes with
                # dimensionless physical unit become Quantity with units of mag.
                unit = Unit(str(unit), parse_strict="silent")
                if not isinstance(unit, (UnitBase, StructuredUnit)):
                    raise UnitTypeError(
                        f"{self.__class__.__name__} instances require normal units, "
                        f"not {unit.__class__} instances."
                    )

        self._unit = unit

    def __deepcopy__(self, memo):
        # If we don't define this, ``copy.deepcopy(quantity)`` will
        # return a bare Numpy array.
        return self.copy()

    def __reduce__(self):
        # patch to pickle Quantity objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        object_state = list(super().__reduce__())
        object_state[2] = (object_state[2], self.__dict__)
        return tuple(object_state)

    def __setstate__(self, state):
        # patch to unpickle Quantity objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        nd_state, own_state = state
        super().__setstate__(nd_state)
        self.__dict__.update(own_state)

    info = QuantityInfo()

    def _to_value(self, unit, equivalencies=[]):
        """Helper method for to and to_value."""
        if equivalencies == []:
            equivalencies = self._equivalencies
        if not self.dtype.names or isinstance(self.unit, StructuredUnit):
            # Standard path, let unit to do work.
            return self.unit.to(
                unit, self.view(np.ndarray), equivalencies=equivalencies
            )

        else:
            # The .to() method of a simple unit cannot convert a structured
            # dtype, so we work around it, by recursing.
            # TODO: deprecate this?
            # Convert simple to Structured on initialization?
            result = np.empty_like(self.view(np.ndarray))
            for name in self.dtype.names:
                result[name] = self[name]._to_value(unit, equivalencies)
            return result

    def to(self, unit, equivalencies=[], copy=True):
        """
        Return a new `~astropy.units.Quantity` object with the specified unit.

        Parameters
        ----------
        unit : unit-like
            An object that represents the unit to convert to. Must be
            an `~astropy.units.UnitBase` object or a string parseable
            by the `~astropy.units` package.

        equivalencies : list of tuple
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`astropy:unit_equivalencies`.
            If not provided or ``[]``, class default equivalencies will be used
            (none for `~astropy.units.Quantity`, but may be set for subclasses)
            If `None`, no equivalencies will be applied at all, not even any
            set globally or within a context.

        copy : bool, optional
            If `True` (default), then the value is copied.  Otherwise, a copy
            will only be made if necessary.

        See Also
        --------
        to_value : get the numerical value in a given unit.
        """
        # We don't use `to_value` below since we always want to make a copy
        # and don't want to slow down this method (esp. the scalar case).
        unit = Unit(unit)
        if copy:
            # Avoid using to_value to ensure that we make a copy. We also
            # don't want to slow down this method (esp. the scalar case).
            value = self._to_value(unit, equivalencies)
        else:
            # to_value only copies if necessary
            value = self.to_value(unit, equivalencies)
        return self._new_view(value, unit)

    def to_value(self, unit=None, equivalencies=[]):
        """
        The numerical value, possibly in a different unit.

        Parameters
        ----------
        unit : unit-like, optional
            The unit in which the value should be given. If not given or `None`,
            use the current unit.

        equivalencies : list of tuple, optional
            A list of equivalence pairs to try if the units are not directly
            convertible (see :ref:`astropy:unit_equivalencies`). If not provided
            or ``[]``, class default equivalencies will be used (none for
            `~astropy.units.Quantity`, but may be set for subclasses).
            If `None`, no equivalencies will be applied at all, not even any
            set globally or within a context.

        Returns
        -------
        value : ndarray or scalar
            The value in the units specified. For arrays, this will be a view
            of the data if no unit conversion was necessary.

        See Also
        --------
        to : Get a new instance in a different unit.
        """
        if unit is None or unit is self.unit:
            value = self.view(np.ndarray)
        elif not self.dtype.names:
            # For non-structured, we attempt a short-cut, where we just get
            # the scale.  If that is 1, we do not have to do anything.
            unit = Unit(unit)
            # We want a view if the unit does not change.  One could check
            # with "==", but that calculates the scale that we need anyway.
            # TODO: would be better for `unit.to` to have an in-place flag.
            try:
                scale = self.unit._to(unit)
            except Exception:
                # Short-cut failed; try default (maybe equivalencies help).
                value = self._to_value(unit, equivalencies)
            else:
                value = self.view(np.ndarray)
                if not is_effectively_unity(scale):
                    # not in-place!
                    value = value * scale
        else:
            # For structured arrays, we go the default route.
            value = self._to_value(unit, equivalencies)

        # Index with empty tuple to decay array scalars in to numpy scalars.
        return value if value.shape else value[()]

    value = property(
        to_value,
        doc="""The numerical value of this instance.

    See also
    --------
    to_value : Get the numerical value in a given unit.
    """,
    )

    @property
    def unit(self):
        """
        A `~astropy.units.UnitBase` object representing the unit of this
        quantity.
        """
        return self._unit

    @property
    def equivalencies(self):
        """
        A list of equivalencies that will be applied by default during
        unit conversions.
        """
        return self._equivalencies

    def _recursively_apply(self, func):
        """Apply function recursively to every field.

        Returns a copy with the result.
        """
        result = np.empty_like(self)
        result_value = result.view(np.ndarray)
        result_unit = ()
        for name in self.dtype.names:
            part = func(self[name])
            result_value[name] = part.value
            result_unit += (part.unit,)

        result._set_unit(result_unit)
        return result

    @property
    def si(self):
        """
        Returns a copy of the current `Quantity` instance with SI units. The
        value of the resulting object will be scaled.
        """
        if self.dtype.names:
            return self._recursively_apply(operator.attrgetter("si"))
        si_unit = self.unit.si
        return self._new_view(self.value * si_unit.scale, si_unit / si_unit.scale)

    @property
    def cgs(self):
        """
        Returns a copy of the current `Quantity` instance with CGS units. The
        value of the resulting object will be scaled.
        """
        if self.dtype.names:
            return self._recursively_apply(operator.attrgetter("cgs"))
        cgs_unit = self.unit.cgs
        return self._new_view(self.value * cgs_unit.scale, cgs_unit / cgs_unit.scale)

    @property
    def isscalar(self):
        """
        True if the `value` of this quantity is a scalar, or False if it
        is an array-like object.

        .. note::
            This is subtly different from `numpy.isscalar` in that
            `numpy.isscalar` returns False for a zero-dimensional array
            (e.g. ``np.array(1)``), while this is True for quantities,
            since quantities cannot represent true numpy scalars.
        """
        return not self.shape

    # This flag controls whether convenience conversion members, such
    # as `q.m` equivalent to `q.to_value(u.m)` are available.  This is
    # not turned on on Quantity itself, but is on some subclasses of
    # Quantity, such as `astropy.coordinates.Angle`.
    _include_easy_conversion_members = False

    def __dir__(self):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.  This function is implemented in
        order to make autocompletion still work correctly in IPython.
        """
        if not self._include_easy_conversion_members:
            return super().__dir__()

        dir_values = set(super().__dir__())
        equivalencies = Unit._normalize_equivalencies(self.equivalencies)
        for equivalent in self.unit._get_units_with_same_physical_type(equivalencies):
            dir_values.update(equivalent.names)
        return sorted(dir_values)

    def __getattr__(self, attr):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.
        """
        if not self._include_easy_conversion_members:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no '{attr}' member"
            )

        def get_virtual_unit_attribute():
            registry = get_current_unit_registry().registry
            to_unit = registry.get(attr, None)
            if to_unit is None:
                return None

            try:
                return self.unit.to(
                    to_unit, self.value, equivalencies=self.equivalencies
                )
            except UnitsError:
                return None

        value = get_virtual_unit_attribute()

        if value is None:
            raise AttributeError(
                f"{self.__class__.__name__} instance has no attribute '{attr}'"
            )
        else:
            return value

    # Equality needs to be handled explicitly as ndarray.__eq__ gives
    # DeprecationWarnings on any error, which is distracting, and does not
    # deal well with structured arrays (nor does the ufunc).
    def __eq__(self, other):
        try:
            other_value = self._to_own_unit(other)
        except UnitsError:
            return False
        except Exception:
            return NotImplemented
        return self.value.__eq__(other_value)

    def __ne__(self, other):
        try:
            other_value = self._to_own_unit(other)
        except UnitsError:
            return True
        except Exception:
            return NotImplemented
        return self.value.__ne__(other_value)

    # Unit conversion operator (<<).
    def __lshift__(self, other):
        try:
            other = Unit(other, parse_strict="silent")
        except UnitTypeError:
            return NotImplemented

        return self.__class__(self, other, copy=False, subok=True)

    def __ilshift__(self, other):
        try:
            other = Unit(other, parse_strict="silent")
        except UnitTypeError:
            return NotImplemented  # try other.__rlshift__(self)

        try:
            factor = self.unit._to(other)
        except UnitConversionError:  # incompatible, or requires an Equivalency
            return NotImplemented
        except AttributeError:  # StructuredUnit does not have `_to`
            # In principle, in-place might be possible.
            return NotImplemented

        view = self.view(np.ndarray)
        try:
            view *= factor  # operates on view
        except TypeError:
            # The error is `numpy.core._exceptions._UFuncOutputCastingError`,
            # which inherits from `TypeError`.
            return NotImplemented

        self._set_unit(other)
        return self

    def __rlshift__(self, other):
        if not self.isscalar:
            return NotImplemented
        return Unit(self).__rlshift__(other)

    # Give warning for other >> self, since probably other << self was meant.
    def __rrshift__(self, other):
        warnings.warn(
            ">> is not implemented. Did you mean to convert "
            "something to this quantity as a unit using '<<'?",
            AstropyWarning,
        )
        return NotImplemented

    # Also define __rshift__ and __irshift__ so we override default ndarray
    # behaviour, but instead of emitting a warning here, let it be done by
    # other (which likely is a unit if this was a mistake).
    def __rshift__(self, other):
        return NotImplemented

    def __irshift__(self, other):
        return NotImplemented

    # Arithmetic operations
    def __mul__(self, other):
        """Multiplication between `Quantity` objects and other objects."""
        if isinstance(other, (UnitBase, str)):
            try:
                return self._new_view(
                    self.value.copy(), other * self.unit, propagate_info=False
                )
            except UnitsError:  # let other try to deal with it
                return NotImplemented

        return super().__mul__(other)

    def __imul__(self, other):
        """In-place multiplication between `Quantity` objects and others."""
        if isinstance(other, (UnitBase, str)):
            self._set_unit(other * self.unit)
            return self

        return super().__imul__(other)

    def __rmul__(self, other):
        """
        Right Multiplication between `Quantity` objects and other objects.
        """
        return self.__mul__(other)

    def __truediv__(self, other):
        """Division between `Quantity` objects and other objects."""
        if isinstance(other, (UnitBase, str)):
            try:
                return self._new_view(
                    self.value.copy(), self.unit / other, propagate_info=False
                )
            except UnitsError:  # let other try to deal with it
                return NotImplemented

        return super().__truediv__(other)

    def __itruediv__(self, other):
        """Inplace division between `Quantity` objects and other objects."""
        if isinstance(other, (UnitBase, str)):
            self._set_unit(self.unit / other)
            return self

        return super().__itruediv__(other)

    def __rtruediv__(self, other):
        """Right Division between `Quantity` objects and other objects."""
        if isinstance(other, (UnitBase, str)):
            return self._new_view(
                1.0 / self.value, other / self.unit, propagate_info=False
            )

        return super().__rtruediv__(other)

    def __pow__(self, other):
        if isinstance(other, Fraction):
            # Avoid getting object arrays by raising the value to a Fraction.
            return self._new_view(
                self.value ** float(other), self.unit**other, propagate_info=False
            )

        return super().__pow__(other)

    # other overrides of special functions
    def __hash__(self):
        return hash(self.value) ^ hash(self.unit)

    def __iter__(self):
        if self.isscalar:
            raise TypeError(
                f"'{self.__class__.__name__}' object with a scalar value is not"
                " iterable"
            )

        # Otherwise return a generator
        def quantity_iter():
            for val in self.value:
                yield self._new_view(val)

        return quantity_iter()

    def __getitem__(self, key):
        if isinstance(key, str) and isinstance(self.unit, StructuredUnit):
            return self._new_view(
                self.view(np.ndarray)[key], self.unit[key], propagate_info=False
            )

        try:
            out = super().__getitem__(key)
        except IndexError:
            # We want zero-dimensional Quantity objects to behave like scalars,
            # so they should raise a TypeError rather than an IndexError.
            if self.isscalar:
                raise TypeError(
                    f"'{self.__class__.__name__}' object with a scalar value "
                    "does not support indexing"
                )
            else:
                raise
        # For single elements, ndarray.__getitem__ returns scalars; these
        # need a new view as a Quantity.
        if not isinstance(out, np.ndarray):
            out = self._new_view(out)
        return out

    def __setitem__(self, i, value):
        if isinstance(i, str):
            # Indexing will cause a different unit, so by doing this in
            # two steps we effectively try with the right unit.
            self[i][...] = value
            return

        # update indices in info if the info property has been accessed
        # (in which case 'info' in self.__dict__ is True; this is guaranteed
        # to be the case if we're part of a table).
        if not self.isscalar and "info" in self.__dict__:
            self.info.adjust_indices(i, value, len(self))
        self.view(np.ndarray).__setitem__(i, self._to_own_unit(value))

    # __contains__ is OK

    def __bool__(self):
        """This method raises ValueError, since truthiness of quantities is ambiguous,
        especially for logarithmic units and temperatures. Use explicit comparisons.
        """
        raise ValueError(
            f"{type(self).__name__} truthiness is ambiguous, especially for logarithmic units"
            " and temperatures. Use explicit comparisons."
        )

    def __len__(self):
        if self.isscalar:
            raise TypeError(
                f"'{self.__class__.__name__}' object with a scalar value has no len()"
            )
        else:
            return len(self.value)

    # Numerical types
    def __float__(self):
        try:
            return float(self.to_value(dimensionless_unscaled))
        except (UnitsError, TypeError):
            raise TypeError(
                "only dimensionless scalar quantities can be "
                "converted to Python scalars"
            )

    def __int__(self):
        try:
            return int(self.to_value(dimensionless_unscaled))
        except (UnitsError, TypeError):
            raise TypeError(
                "only dimensionless scalar quantities can be "
                "converted to Python scalars"
            )

    def __index__(self):
        # for indices, we do not want to mess around with scaling at all,
        # so unlike for float, int, we insist here on unscaled dimensionless
        try:
            assert self.unit.is_unity()
            return self.value.__index__()
        except Exception:
            raise TypeError(
                "only integer dimensionless scalar quantities "
                "can be converted to a Python index"
            )

    # TODO: we may want to add a hook for dimensionless quantities?
    @property
    def _unitstr(self):
        if self.unit is None:
            unitstr = _UNIT_NOT_INITIALISED
        else:
            unitstr = str(self.unit)

        if unitstr:
            unitstr = " " + unitstr

        return unitstr

    def to_string(self, unit=None, precision=None, format=None, subfmt=None):
        """
        Generate a string representation of the quantity and its unit.

        The behavior of this function can be altered via the
        `numpy.set_printoptions` function and its various keywords.  The
        exception to this is the ``threshold`` keyword, which is controlled via
        the ``[units.quantity]`` configuration item ``latex_array_threshold``.
        This is treated separately because the numpy default of 1000 is too big
        for most browsers to handle.

        Parameters
        ----------
        unit : unit-like, optional
            Specifies the unit.  If not provided,
            the unit used to initialize the quantity will be used.

        precision : number, optional
            The level of decimal precision. If `None`, or not provided,
            it will be determined from NumPy print options.

        format : str, optional
            The format of the result. If not provided, an unadorned
            string is returned. Supported values are:

            - 'latex': Return a LaTeX-formatted string

            - 'latex_inline': Return a LaTeX-formatted string that uses
              negative exponents instead of fractions

        subfmt : str, optional
            Subformat of the result. For the moment, only used for
            ``format='latex'`` and ``format='latex_inline'``. Supported
            values are:

            - 'inline': Use ``$ ... $`` as delimiters.

            - 'display': Use ``$\\displaystyle ... $`` as delimiters.

        Returns
        -------
        str
            A string with the contents of this Quantity
        """
        if unit is not None and unit != self.unit:
            return self.to(unit).to_string(
                unit=None, precision=precision, format=format, subfmt=subfmt
            )

        formats = {
            None: None,
            "latex": {
                None: ("$", "$"),
                "inline": ("$", "$"),
                "display": (r"$\displaystyle ", r"$"),
            },
        }
        formats["latex_inline"] = formats["latex"]

        if format not in formats:
            raise ValueError(f"Unknown format '{format}'")
        elif format is None:
            if precision is None:
                # Use default formatting settings
                return f"{self.value}{self._unitstr:s}"
            else:
                # np.array2string properly formats arrays as well as scalars
                return (
                    np.array2string(self.value, precision=precision, floatmode="fixed")
                    + self._unitstr
                )

        # else, for the moment we assume format="latex" or "latex_inline".

        # Set the precision if set, otherwise use numpy default
        pops = np.get_printoptions()
        format_spec = f".{precision if precision is not None else pops['precision']}g"

        def float_formatter(value):
            return Latex.format_exponential_notation(value, format_spec=format_spec)

        def complex_formatter(value):
            return "({}{}i)".format(
                Latex.format_exponential_notation(value.real, format_spec=format_spec),
                Latex.format_exponential_notation(
                    value.imag, format_spec="+" + format_spec
                ),
            )

        # The view is needed for the scalar case - self.value might be float.
        latex_value = np.array2string(
            self.view(np.ndarray),
            threshold=(
                conf.latex_array_threshold
                if conf.latex_array_threshold > -1
                else pops["threshold"]
            ),
            formatter={
                "float_kind": float_formatter,
                "complex_kind": complex_formatter,
            },
            max_line_width=np.inf,
            separator=",~",
        )

        latex_value = latex_value.replace("...", r"\dots")

        # Format unit
        # [1:-1] strips the '$' on either side needed for math mode
        if self.unit is None:
            latex_unit = _UNIT_NOT_INITIALISED
        elif format == "latex":
            latex_unit = self.unit._repr_latex_()[1:-1]  # note this is unicode
        elif format == "latex_inline":
            latex_unit = self.unit.to_string(format="latex_inline")[1:-1]

        delimiter_left, delimiter_right = formats[format][subfmt]

        return rf"{delimiter_left}{latex_value} \; {latex_unit}{delimiter_right}"

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        prefixstr = "<" + self.__class__.__name__ + " "
        arrstr = np.array2string(
            self.view(np.ndarray), separator=", ", prefix=prefixstr
        )
        return f"{prefixstr}{arrstr}{self._unitstr:s}>"

    def _repr_latex_(self):
        """
        Generate a latex representation of the quantity and its unit.

        Returns
        -------
        lstr
            A LaTeX string with the contents of this Quantity
        """
        # NOTE: This should change to display format in a future release
        return self.to_string(format="latex", subfmt="inline")

    def __format__(self, format_spec):
        try:
            return self.to_string(format=format_spec)
        except ValueError:
            # We might have a unit format not implemented in `to_string()`.
            if format_spec in Base.registry:
                if self.unit is dimensionless_unscaled:
                    return f"{self.value}"
                else:
                    return f"{self.value} {format(self.unit, format_spec)}"
            # Can the value be formatted on its own?
            try:
                return f"{format(self.value, format_spec)}{self._unitstr:s}"
            except ValueError:
                # Format the whole thing as a single string.
                return format(f"{self.value}{self._unitstr:s}", format_spec)

    def decompose(self, bases=[]):
        """
        Generates a new `Quantity` with the units
        decomposed. Decomposed units have only irreducible units in
        them (see `astropy.units.UnitBase.decompose`).

        Parameters
        ----------
        bases : sequence of `~astropy.units.UnitBase`, optional
            The bases to decompose into.  When not provided,
            decomposes down to any irreducible units.  When provided,
            the decomposed result will only contain the given units.
            This will raises a `~astropy.units.UnitsError` if it's not possible
            to do so.

        Returns
        -------
        newq : `~astropy.units.Quantity`
            A new object equal to this quantity with units decomposed.
        """
        return self._decompose(False, bases=bases)

    def _decompose(self, allowscaledunits=False, bases=[]):
        """
        Generates a new `Quantity` with the units decomposed. Decomposed
        units have only irreducible units in them (see
        `astropy.units.UnitBase.decompose`).

        Parameters
        ----------
        allowscaledunits : bool
            If True, the resulting `Quantity` may have a scale factor
            associated with it.  If False, any scaling in the unit will
            be subsumed into the value of the resulting `Quantity`

        bases : sequence of UnitBase, optional
            The bases to decompose into.  When not provided,
            decomposes down to any irreducible units.  When provided,
            the decomposed result will only contain the given units.
            This will raises a `~astropy.units.UnitsError` if it's not possible
            to do so.

        Returns
        -------
        newq : `~astropy.units.Quantity`
            A new object equal to this quantity with units decomposed.

        """
        new_unit = self.unit.decompose(bases=bases)

        # Be careful here because self.value usually is a view of self;
        # be sure that the original value is not being modified.
        if not allowscaledunits and hasattr(new_unit, "scale"):
            new_value = self.value * new_unit.scale
            new_unit = new_unit / new_unit.scale
            return self._new_view(new_value, new_unit)
        else:
            return self._new_view(self.copy(), new_unit)

    # These functions need to be overridden to take into account the units
    # Array conversion
    # https://numpy.org/doc/stable/reference/arrays.ndarray.html#array-conversion

    def item(self, *args):
        """Copy an element of an array to a scalar Quantity and return it.

        Like :meth:`~numpy.ndarray.item` except that it always
        returns a `Quantity`, not a Python scalar.

        """
        return self._new_view(super().item(*args))

    def tolist(self):
        raise NotImplementedError(
            "cannot make a list of Quantities. Get list of values with"
            " q.value.tolist()."
        )

    def _to_own_unit(self, value, check_precision=True, *, unit=None):
        """Convert value to one's own unit (or that given).

        Here, non-quantities are treated as dimensionless, and care is taken
        for values of 0, infinity or nan, which are allowed to have any unit.

        Parameters
        ----------
        value : anything convertible to `~astropy.units.Quantity`
            The value to be converted to the requested unit.
        check_precision : bool
            Whether to forbid conversion of float to integer if that changes
            the input number.  Default: `True`.
        unit : `~astropy.units.Unit` or None
            The unit to convert to.  By default, the unit of ``self``.

        Returns
        -------
        value : number or `~numpy.ndarray`
            In the requested units.

        """
        if unit is None:
            unit = self.unit
        try:
            _value = value.to_value(unit)
        except AttributeError:
            # We're not a Quantity.
            # First remove two special cases (with a fast test):
            # 1) Maybe masked printing? MaskedArray with quantities does not
            # work very well, but no reason to break even repr and str.
            # 2) np.ma.masked? useful if we're a MaskedQuantity.
            if value is np.ma.masked or (
                value is np.ma.masked_print_option and self.dtype.kind == "O"
            ):
                return value
            # Now, let's try a more general conversion.
            # Plain arrays will be converted to dimensionless in the process,
            # but anything with a unit attribute will use that.
            try:
                as_quantity = Quantity(value)
                _value = as_quantity.to_value(unit)
            except UnitsError:
                # last chance: if this was not something with a unit
                # and is all 0, inf, or nan, we treat it as arbitrary unit.
                if not hasattr(value, "unit") and can_have_arbitrary_unit(
                    as_quantity.value
                ):
                    _value = as_quantity.value
                else:
                    raise

        if self.dtype.kind == "i" and check_precision:
            # If, e.g., we are casting float to int, we want to fail if
            # precision is lost, but let things pass if it works.
            _value = np.array(_value, copy=COPY_IF_NEEDED, subok=True)
            if not np.can_cast(_value.dtype, self.dtype):
                self_dtype_array = np.array(_value, self.dtype, subok=True)
                if not np.all((self_dtype_array == _value) | np.isnan(_value)):
                    raise TypeError(
                        "cannot convert value type to array type without precision loss"
                    )

        # Setting names to ensure things like equality work (note that
        # above will have failed already if units did not match).
        # TODO: is this the best place to do this?
        if _value.dtype.names is not None:
            _value = _value.astype(self.dtype, copy=False)
        return _value

    if NUMPY_LT_2_0:

        def itemset(self, *args):
            if len(args) == 0:
                raise ValueError("itemset must have at least one argument")

            self.view(np.ndarray).itemset(*(args[:-1] + (self._to_own_unit(args[-1]),)))

    def tostring(self, order="C"):
        """Not implemented, use ``.value.tostring()`` instead."""
        raise NotImplementedError(
            "cannot write Quantities to string.  Write array with"
            " q.value.tostring(...)."
        )

    def tobytes(self, order="C"):
        """Not implemented, use ``.value.tobytes()`` instead."""
        raise NotImplementedError(
            "cannot write Quantities to bytes.  Write array with q.value.tobytes(...)."
        )

    def tofile(self, fid, sep="", format="%s"):
        """Not implemented, use ``.value.tofile()`` instead."""
        raise NotImplementedError(
            "cannot write Quantities to file.  Write array with q.value.tofile(...)"
        )

    def dump(self, file):
        """Not implemented, use ``.value.dump()`` instead."""
        raise NotImplementedError(
            "cannot dump Quantities to file.  Write array with q.value.dump()"
        )

    def dumps(self):
        """Not implemented, use ``.value.dumps()`` instead."""
        raise NotImplementedError(
            "cannot dump Quantities to string.  Write array with q.value.dumps()"
        )

    # astype, byteswap, copy, view, getfield, setflags OK as is

    def fill(self, value):
        self.view(np.ndarray).fill(self._to_own_unit(value))

    # Shape manipulation: resize cannot be done (does not own data), but
    # shape, transpose, swapaxes, flatten, ravel, squeeze all OK.  Only
    # the flat iterator needs to be overwritten, otherwise single items are
    # returned as numbers.
    @property
    def flat(self):
        """A 1-D iterator over the Quantity array.

        This returns a ``QuantityIterator`` instance, which behaves the same
        as the `~numpy.flatiter` instance returned by `~numpy.ndarray.flat`,
        and is similar to, but not a subclass of, Python's built-in iterator
        object.
        """
        return QuantityIterator(self)

    @flat.setter
    def flat(self, value):
        y = self.ravel()
        y[:] = value

    # Item selection and manipulation
    # repeat, sort, compress, diagonal OK
    def take(self, indices, axis=None, out=None, mode="raise"):
        out = super().take(indices, axis=axis, out=out, mode=mode)
        # For single elements, ndarray.take returns scalars; these
        # need a new view as a Quantity.
        if type(out) is not type(self):
            out = self._new_view(out)
        return out

    def put(self, indices, values, mode="raise"):
        self.view(np.ndarray).put(indices, self._to_own_unit(values), mode)

    def choose(self, choices, out=None, mode="raise"):
        raise NotImplementedError(
            "cannot choose based on quantity.  Choose using array with"
            " q.value.choose(...)"
        )

    # ensure we do not return indices as quantities
    if NUMPY_LT_2_0:

        def argsort(self, axis=-1, kind=None, order=None):
            return self.view(np.ndarray).argsort(axis=axis, kind=kind, order=order)

    else:

        def argsort(self, axis=-1, kind=None, order=None, *, stable=None):
            return self.view(np.ndarray).argsort(
                axis=axis, kind=kind, order=order, stable=stable
            )

    def searchsorted(self, v, *args, **kwargs):
        return np.searchsorted(
            np.array(self), self._to_own_unit(v, check_precision=False), *args, **kwargs
        )  # avoid numpy 1.6 problem

    def argmax(self, axis=None, out=None, *, keepdims=False):
        return self.view(np.ndarray).argmax(axis=axis, out=out, keepdims=keepdims)

    def argmin(self, axis=None, out=None, *, keepdims=False):
        return self.view(np.ndarray).argmin(axis=axis, out=out, keepdims=keepdims)

    def __array_function__(self, function, types, args, kwargs):
        """Wrap numpy functions, taking care of units.

        Parameters
        ----------
        function : callable
            Numpy function to wrap
        types : iterable of classes
            Classes that provide an ``__array_function__`` override. Can
            in principle be used to interact with other classes. Below,
            mostly passed on to `~numpy.ndarray`, which can only interact
            with subclasses.
        args : tuple
            Positional arguments provided in the function call.
        kwargs : dict
            Keyword arguments provided in the function call.

        Returns
        -------
        result: `~astropy.units.Quantity`, `~numpy.ndarray`
            As appropriate for the function.  If the function is not
            supported, `NotImplemented` is returned, which will lead to
            a `TypeError` unless another argument overrode the function.

        Raises
        ------
        ~astropy.units.UnitsError
            If operands have incompatible units.
        """
        # A function should be in one of the following sets or dicts:
        # 1. SUBCLASS_SAFE_FUNCTIONS (set), if the numpy implementation
        #    supports Quantity; we pass on to ndarray.__array_function__.
        # 2. FUNCTION_HELPERS (dict), if the numpy implementation is usable
        #    after converting quantities to arrays with suitable units,
        #    and possibly setting units on the result.
        # 3. DISPATCHED_FUNCTIONS (dict), if the function makes sense but
        #    requires a Quantity-specific implementation.
        # 4. UNSUPPORTED_FUNCTIONS (set), if the function does not make sense.
        # For now, since we may not yet have complete coverage, if a
        # function is in none of the above, we simply call the numpy
        # implementation.
        if function in SUBCLASS_SAFE_FUNCTIONS:
            return super().__array_function__(function, types, args, kwargs)

        elif function in FUNCTION_HELPERS:
            function_helper = FUNCTION_HELPERS[function]
            try:
                args, kwargs, unit, out = function_helper(*args, **kwargs)
            except NotImplementedError:
                return self._not_implemented_or_raise(function, types)

            result = super().__array_function__(function, types, args, kwargs)
            # Fall through to return section

        elif function in DISPATCHED_FUNCTIONS:
            dispatched_function = DISPATCHED_FUNCTIONS[function]
            try:
                result, unit, out = dispatched_function(*args, **kwargs)
            except NotImplementedError:
                return self._not_implemented_or_raise(function, types)

            # Fall through to return section

        elif function in UNSUPPORTED_FUNCTIONS:
            return NotImplemented

        else:
            warnings.warn(
                f"function '{function.__name__}' is not known to astropy's Quantity."
                " Will run it anyway, hoping it will treat ndarray subclasses"
                " correctly. Please raise an issue at"
                " https://github.com/astropy/astropy/issues.",
                AstropyWarning,
            )
            return super().__array_function__(function, types, args, kwargs)

        # If unit is None, a plain array is expected (e.g., boolean), which
        # means we're done.
        # We're also done if the result was NotImplemented, which can happen
        # if other inputs/outputs override __array_function__;
        # hopefully, they can then deal with us.
        if unit is None or result is NotImplemented:
            return result

        return self._result_as_quantity(result, unit, out=out)

    def _not_implemented_or_raise(self, function, types):
        # Our function helper or dispatcher found that the function does not
        # work with Quantity.  In principle, there may be another class that
        # knows what to do with us, for which we should return NotImplemented.
        # But if there is ndarray (or a non-Quantity subclass of it) around,
        # it quite likely coerces, so we should just break.
        if any(
            issubclass(t, np.ndarray) and not issubclass(t, Quantity) for t in types
        ):
            raise TypeError(
                f"the Quantity implementation cannot handle {function} "
                "with the given arguments."
            ) from None
        else:
            return NotImplemented

    # Calculation -- override ndarray methods to take into account units.
    # We use the corresponding numpy functions to evaluate the results, since
    # the methods do not always allow calling with keyword arguments.
    # For instance, np.array([0.,2.]).clip(a_min=0., a_max=1.) gives
    # TypeError: 'a_max' is an invalid keyword argument for this function.
    def _wrap_function(self, function, *args, unit=None, out=None, **kwargs):
        """Wrap a numpy function that processes self, returning a Quantity.

        Parameters
        ----------
        function : callable
            Numpy function to wrap.
        args : positional arguments
            Any positional arguments to the function beyond the first argument
            (which will be set to ``self``).
        kwargs : keyword arguments
            Keyword arguments to the function.

        If present, the following arguments are treated specially:

        unit : `~astropy.units.Unit`
            Unit of the output result.  If not given, the unit of ``self``.
        out : `~astropy.units.Quantity`
            A Quantity instance in which to store the output.

        Notes
        -----
        Output should always be assigned via a keyword argument, otherwise
        no proper account of the unit is taken.

        Returns
        -------
        out : `~astropy.units.Quantity`
            Result of the function call, with the unit set properly.
        """
        if unit is None:
            unit = self.unit
        # Ensure we don't loop back by turning any Quantity into array views.
        args = (self.value,) + tuple(
            (arg.value if isinstance(arg, Quantity) else arg) for arg in args
        )
        if out is not None:
            # If pre-allocated output is used, check it is suitable.
            # This also returns array view, to ensure we don't loop back.
            arrays = tuple(arg for arg in args if isinstance(arg, np.ndarray))
            kwargs["out"] = check_output(out, unit, arrays, function=function)
        # Apply the function and turn it back into a Quantity.
        result = function(*args, **kwargs)
        return self._result_as_quantity(result, unit, out)

    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        return self._wrap_function(np.trace, offset, axis1, axis2, dtype, out=out)

    def var(
        self, axis=None, dtype=None, out=None, ddof=0, keepdims=False, *, where=True
    ):
        return self._wrap_function(
            np.var,
            axis,
            dtype,
            out=out,
            ddof=ddof,
            keepdims=keepdims,
            where=where,
            unit=self.unit**2,
        )

    def std(
        self, axis=None, dtype=None, out=None, ddof=0, keepdims=False, *, where=True
    ):
        return self._wrap_function(
            np.std, axis, dtype, out=out, ddof=ddof, keepdims=keepdims, where=where
        )

    def mean(self, axis=None, dtype=None, out=None, keepdims=False, *, where=True):
        return self._wrap_function(
            np.mean, axis, dtype, out=out, keepdims=keepdims, where=where
        )

    def round(self, decimals=0, out=None):
        return self._wrap_function(np.round, decimals, out=out)

    def dot(self, b, out=None):
        result_unit = self.unit * getattr(b, "unit", dimensionless_unscaled)
        return self._wrap_function(np.dot, b, out=out, unit=result_unit)

    # Calculation: override methods that do not make sense.

    def all(self, axis=None, out=None):
        raise TypeError(
            "cannot evaluate truth value of quantities. "
            "Evaluate array with q.value.all(...)"
        )

    def any(self, axis=None, out=None):
        raise TypeError(
            "cannot evaluate truth value of quantities. "
            "Evaluate array with q.value.any(...)"
        )

    # Calculation: numpy functions that can be overridden with methods.

    def diff(self, n=1, axis=-1):
        return self._wrap_function(np.diff, n, axis)

    def ediff1d(self, to_end=None, to_begin=None):
        return self._wrap_function(np.ediff1d, to_end, to_begin)

    @deprecated("5.3", alternative="np.nansum", obj_type="method")
    def nansum(self, axis=None, out=None, keepdims=False, *, initial=None, where=True):
        if initial is not None:
            initial = self._to_own_unit(initial)
        return self._wrap_function(
            np.nansum,
            axis,
            out=out,
            keepdims=keepdims,
            initial=initial,
            where=where,
        )

    def insert(self, obj, values, axis=None):
        """
        Insert values along the given axis before the given indices and return
        a new `~astropy.units.Quantity` object.

        This is a thin wrapper around the `numpy.insert` function.

        Parameters
        ----------
        obj : int, slice or sequence of int
            Object that defines the index or indices before which ``values`` is
            inserted.
        values : array-like
            Values to insert.  If the type of ``values`` is different
            from that of quantity, ``values`` is converted to the matching type.
            ``values`` should be shaped so that it can be broadcast appropriately
            The unit of ``values`` must be consistent with this quantity.
        axis : int, optional
            Axis along which to insert ``values``.  If ``axis`` is None then
            the quantity array is flattened before insertion.

        Returns
        -------
        out : `~astropy.units.Quantity`
            A copy of quantity with ``values`` inserted.  Note that the
            insertion does not occur in-place: a new quantity array is returned.

        Examples
        --------
        >>> import astropy.units as u
        >>> q = [1, 2] * u.m
        >>> q.insert(0, 50 * u.cm)
        <Quantity [ 0.5,  1.,  2.] m>

        >>> q = [[1, 2], [3, 4]] * u.m
        >>> q.insert(1, [10, 20] * u.m, axis=0)
        <Quantity [[  1.,  2.],
                   [ 10., 20.],
                   [  3.,  4.]] m>

        >>> q.insert(1, 10 * u.m, axis=1)
        <Quantity [[  1., 10.,  2.],
                   [  3., 10.,  4.]] m>

        """
        out_array = np.insert(self.value, obj, self._to_own_unit(values), axis)
        return self._new_view(out_array)


class SpecificTypeQuantity(Quantity):
    """Superclass for Quantities of specific physical type.

    Subclasses of these work just like :class:`~astropy.units.Quantity`, except
    that they are for specific physical types (and may have methods that are
    only appropriate for that type).  Astropy examples are
    :class:`~astropy.coordinates.Angle` and
    :class:`~astropy.coordinates.Distance`

    At a minimum, subclasses should set ``_equivalent_unit`` to the unit
    associated with the physical type.
    """

    # The unit for the specific physical type.  Instances can only be created
    # with units that are equivalent to this.
    _equivalent_unit = None

    # The default unit used for views.  Even with `None`, views of arrays
    # without units are possible, but will have an uninitialized unit.
    _unit = None

    # Default unit for initialization through the constructor.
    _default_unit = None

    # ensure that we get precedence over our superclass.
    __array_priority__ = Quantity.__array_priority__ + 10

    def __quantity_subclass__(self, unit):
        if unit.is_equivalent(self._equivalent_unit):
            return type(self), True
        else:
            return super().__quantity_subclass__(unit)[0], False

    def _set_unit(self, unit):
        if unit is None or not unit.is_equivalent(self._equivalent_unit):
            raise UnitTypeError(
                "{} instances require units equivalent to '{}'".format(
                    type(self).__name__, self._equivalent_unit
                )
                + (
                    ", but no unit was given."
                    if unit is None
                    else f", so cannot set it to '{unit}'."
                )
            )

        super()._set_unit(unit)


def isclose(a, b, rtol=1.0e-5, atol=None, equal_nan=False):
    """
    Return a boolean array where two arrays are element-wise equal
    within a tolerance.

    Parameters
    ----------
    a, b : array-like or `~astropy.units.Quantity`
        Input values or arrays to compare
    rtol : array-like or `~astropy.units.Quantity`
        The relative tolerance for the comparison, which defaults to
        ``1e-5``.  If ``rtol`` is a :class:`~astropy.units.Quantity`,
        then it must be dimensionless.
    atol : number or `~astropy.units.Quantity`
        The absolute tolerance for the comparison.  The units (or lack
        thereof) of ``a``, ``b``, and ``atol`` must be consistent with
        each other.  If `None`, ``atol`` defaults to zero in the
        appropriate units.
    equal_nan : `bool`
        Whether to compare NaN’s as equal. If `True`, NaNs in ``a`` will
        be considered equal to NaN’s in ``b``.

    Notes
    -----
    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.isclose`. However, this differs from the `numpy` function in
    that the default for the absolute tolerance here is zero instead of
    ``atol=1e-8`` in `numpy`, as there is no natural way to set a default
    *absolute* tolerance given two inputs that may have differently scaled
    units.

    Raises
    ------
    `~astropy.units.UnitsError`
        If the dimensions of ``a``, ``b``, or ``atol`` are incompatible,
        or if ``rtol`` is not dimensionless.

    See Also
    --------
    allclose
    """
    return np.isclose(*_unquantify_allclose_arguments(a, b, rtol, atol), equal_nan)


def allclose(a, b, rtol=1.0e-5, atol=None, equal_nan=False) -> bool:
    """
    Whether two arrays are element-wise equal within a tolerance.

    Parameters
    ----------
    a, b : array-like or `~astropy.units.Quantity`
        Input values or arrays to compare
    rtol : array-like or `~astropy.units.Quantity`
        The relative tolerance for the comparison, which defaults to
        ``1e-5``.  If ``rtol`` is a :class:`~astropy.units.Quantity`,
        then it must be dimensionless.
    atol : number or `~astropy.units.Quantity`
        The absolute tolerance for the comparison.  The units (or lack
        thereof) of ``a``, ``b``, and ``atol`` must be consistent with
        each other.  If `None`, ``atol`` defaults to zero in the
        appropriate units.
    equal_nan : `bool`
        Whether to compare NaN’s as equal. If `True`, NaNs in ``a`` will
        be considered equal to NaN’s in ``b``.

    Notes
    -----
    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.allclose`. However, this differs from the `numpy` function in
    that the default for the absolute tolerance here is zero instead of
    ``atol=1e-8`` in `numpy`, as there is no natural way to set a default
    *absolute* tolerance given two inputs that may have differently scaled
    units.

    Raises
    ------
    `~astropy.units.UnitsError`
        If the dimensions of ``a``, ``b``, or ``atol`` are incompatible,
        or if ``rtol`` is not dimensionless.

    See Also
    --------
    isclose
    """
    return np.allclose(*_unquantify_allclose_arguments(a, b, rtol, atol), equal_nan)


def _unquantify_allclose_arguments(actual, desired, rtol, atol):
    actual = Quantity(actual, subok=True, copy=False)

    desired = Quantity(desired, subok=True, copy=False)
    try:
        desired = desired.to(actual.unit)
    except UnitsError:
        raise UnitsError(
            f"Units for 'desired' ({desired.unit}) and 'actual' "
            f"({actual.unit}) are not convertible"
        )

    if atol is None:
        # By default, we assume an absolute tolerance of zero in the
        # appropriate units.  The default value of None for atol is
        # needed because the units of atol must be consistent with the
        # units for a and b.
        atol = Quantity(0)
    else:
        atol = Quantity(atol, subok=True, copy=False)
        try:
            atol = atol.to(actual.unit)
        except UnitsError:
            raise UnitsError(
                f"Units for 'atol' ({atol.unit}) and 'actual' "
                f"({actual.unit}) are not convertible"
            )

    rtol = Quantity(rtol, subok=True, copy=False)
    try:
        rtol = rtol.to(dimensionless_unscaled)
    except Exception:
        raise UnitsError("'rtol' should be dimensionless")

    return actual.value, desired.value, rtol.value, atol.value
