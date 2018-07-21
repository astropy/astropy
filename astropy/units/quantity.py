# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines the `Quantity` object, which represents a number with some
associated units. `Quantity` objects support operations like ordinary numbers,
but will deal with unit conversions internally.
"""


# Standard library
import re
import numbers
from fractions import Fraction
import warnings

import numpy as np

# AstroPy
from .core import (Unit, dimensionless_unscaled, get_current_unit_registry,
                   UnitBase, UnitsError, UnitConversionError, UnitTypeError)
from .utils import is_effectively_unity
from .format.latex import Latex
from ..utils.compat import NUMPY_LT_1_14
from ..utils.compat.misc import override__dir__
from ..utils.exceptions import AstropyDeprecationWarning, AstropyWarning
from ..utils.misc import isiterable, InheritDocstrings
from ..utils.data_info import ParentDtypeInfo
from .. import config as _config
from .quantity_helper import (converters_and_unit, can_have_arbitrary_unit,
                              check_output)

__all__ = ["Quantity", "SpecificTypeQuantity",
           "QuantityInfoBase", "QuantityInfo", "allclose", "isclose"]


# We don't want to run doctests in the docstrings we inherit from Numpy
__doctest_skip__ = ['Quantity.*']

_UNIT_NOT_INITIALISED = "(Unit not initialised)"
_UFUNCS_FILTER_WARNINGS = {np.arcsin, np.arccos, np.arccosh, np.arctanh}


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for Quantity
    """
    latex_array_threshold = _config.ConfigItem(100,
        'The maximum size an array Quantity can be before its LaTeX '
        'representation for IPython gets "summarized" (meaning only the first '
        'and last few elements are shown with "..." between). Setting this to a '
        'negative number means that the value will instead be whatever numpy '
        'gets from get_printoptions.')


conf = Conf()


class QuantityIterator:
    """
    Flat iterator object to iterate over Quantities

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


class QuantityInfoBase(ParentDtypeInfo):
    # This is on a base class rather than QuantityInfo directly, so that
    # it can be used for EarthLocationInfo yet make clear that that class
    # should not be considered a typical Quantity subclass by Table.
    attrs_from_parent = {'dtype', 'unit'}  # dtype and unit taken from parent
    _supports_indexing = True

    @staticmethod
    def default_format(val):
        return '{0.value:}'.format(val)

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
    _represent_as_dict_attrs = ('value', 'unit')
    _construct_from_dict_args = ['value']
    _represent_as_dict_primary_data = 'value'

    def new_like(self, cols, length, metadata_conflicts='warn', name=None):
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
        col : Quantity (or subclass)
            Empty instance of this class consistent with ``cols``

        """

        # Get merged info attributes like shape, dtype, format, description, etc.
        attrs = self.merge_cols_attributes(cols, metadata_conflicts, name,
                                           ('meta', 'format', 'description'))

        # Make an empty quantity using the unit of the last one.
        shape = (length,) + attrs.pop('shape')
        dtype = attrs.pop('dtype')
        # Use zeros so we do not get problems for Quantity subclasses such
        # as Longitude and Latitude, which cannot take arbitrary values.
        data = np.zeros(shape=shape, dtype=dtype)
        # Get arguments needed to reconstruct class
        map = {key: (data if key == 'value' else getattr(cols[-1], key))
               for key in self._represent_as_dict_attrs}
        map['copy'] = False
        out = self._construct_from_dict(map)

        # Set remaining info attributes
        for attr, value in attrs.items():
            setattr(out.info, attr, value)

        return out


class Quantity(np.ndarray, metaclass=InheritDocstrings):
    """A `~astropy.units.Quantity` represents a number with some associated unit.

    Parameters
    ----------
    value : number, `~numpy.ndarray`, `Quantity` object (sequence), str
        The numerical value of this quantity in the units given by unit.  If a
        `Quantity` or sequence of them (or any other valid object with a
        ``unit`` attribute), creates a new `Quantity` object, converting to
        `unit` units as needed.  If a string, it is converted to a number or
        `Quantity`, depending on whether a unit is present.

    unit : `~astropy.units.UnitBase` instance, str
        An object that represents the unit associated with the input value.
        Must be an `~astropy.units.UnitBase` object or a string parseable by
        the :mod:`~astropy.units` package.

    dtype : ~numpy.dtype, optional
        The dtype of the resulting Numpy array or scalar that will
        hold the value.  If not provided, it is determined from the input,
        except that any input that cannot represent float (integer and bool)
        is converted to float.

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
        should have.  Ones will be pre-pended to the shape as needed to meet
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
    :class:`~astropy.units.Unit`. See http://docs.astropy.org/en/latest/units/

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

    def __new__(cls, value, unit=None, dtype=None, copy=True, order=None,
                subok=False, ndmin=0):

        if unit is not None:
            # convert unit first, to avoid multiple string->unit conversions
            unit = Unit(unit)
            # if we allow subclasses, allow a class from the unit.
            if subok:
                qcls = getattr(unit, '_quantity_class', cls)
                if issubclass(qcls, cls):
                    cls = qcls

        # optimize speed for Quantity with no dtype given, copy=False
        if isinstance(value, Quantity):
            if unit is not None and unit is not value.unit:
                value = value.to(unit)
                # the above already makes a copy (with float dtype)
                copy = False

            if type(value) is not cls and not (subok and
                                               isinstance(value, cls)):
                value = value.view(cls)

            if dtype is None:
                if not copy:
                    return value

                if not np.can_cast(np.float32, value.dtype):
                    dtype = float

            return np.array(value, dtype=dtype, copy=copy, order=order,
                            subok=True, ndmin=ndmin)

        # Maybe str, or list/tuple of Quantity? If so, this may set value_unit.
        # To ensure array remains fast, we short-circuit it.
        value_unit = None
        if not isinstance(value, np.ndarray):
            if isinstance(value, str):
                # The first part of the regex string matches any integer/float;
                # the second parts adds possible trailing .+-, which will break
                # the float function below and ensure things like 1.2.3deg
                # will not work.
                pattern = (r'\s*[+-]?'
                           r'((\d+\.?\d*)|(\.\d+)|([nN][aA][nN])|'
                           r'([iI][nN][fF]([iI][nN][iI][tT][yY]){0,1}))'
                           r'([eE][+-]?\d+)?'
                           r'[.+-]?')

                v = re.match(pattern, value)
                unit_string = None
                try:
                    value = float(v.group())

                except Exception:
                    raise TypeError('Cannot parse "{0}" as a {1}. It does not '
                                    'start with a number.'
                                    .format(value, cls.__name__))

                unit_string = v.string[v.end():].strip()
                if unit_string:
                    value_unit = Unit(unit_string)
                    if unit is None:
                        unit = value_unit  # signal no conversion needed below.

            elif (isiterable(value) and len(value) > 0 and
                  all(isinstance(v, Quantity) for v in value)):
                # Convert all quantities to the same unit.
                if unit is None:
                    unit = value[0].unit
                value = [q.to_value(unit) for q in value]
                value_unit = unit  # signal below that conversion has been done

        if value_unit is None:
            # If the value has a `unit` attribute and if not None
            # (for Columns with uninitialized unit), treat it like a quantity.
            value_unit = getattr(value, 'unit', None)
            if value_unit is None:
                # Default to dimensionless for no (initialized) unit attribute.
                if unit is None:
                    unit = cls._default_unit
                value_unit = unit  # signal below that no conversion is needed
            else:
                try:
                    value_unit = Unit(value_unit)
                except Exception as exc:
                    raise TypeError("The unit attribute {0!r} of the input could "
                                    "not be parsed as an astropy Unit, raising "
                                    "the following exception:\n{1}"
                                    .format(value.unit, exc))

                if unit is None:
                    unit = value_unit
                elif unit is not value_unit:
                    copy = False  # copy will be made in conversion at end

        value = np.array(value, dtype=dtype, copy=copy, order=order,
                         subok=False, ndmin=ndmin)

        # check that array contains numbers or long int objects
        if (value.dtype.kind in 'OSU' and
            not (value.dtype.kind == 'O' and
                 isinstance(value.item(() if value.ndim == 0 else 0),
                            numbers.Number))):
            raise TypeError("The value must be a valid Python or "
                            "Numpy numeric type.")

        # by default, cast any integer, boolean, etc., to float
        if dtype is None and (not np.can_cast(np.float32, value.dtype)
                              or value.dtype.kind == 'O'):
            value = value.astype(float)

        value = value.view(cls)
        value._set_unit(value_unit)
        if unit is value_unit:
            return value
        else:
            # here we had non-Quantity input that had a "unit" attribute
            # with a unit different from the desired one.  So, convert.
            return value.to(unit)

    def __array_finalize__(self, obj):
        # If we're a new object or viewing an ndarray, nothing has to be done.
        if obj is None or obj.__class__ is np.ndarray:
            return

        # If our unit is not set and obj has a valid one, use it.
        if self._unit is None:
            unit = getattr(obj, '_unit', None)
            if unit is not None:
                self._set_unit(unit)

        # Copy info if the original had `info` defined.  Because of the way the
        # DataInfo works, `'info' in obj.__dict__` is False until the
        # `info` attribute is accessed or set.
        if 'info' in obj.__dict__:
            self.info = obj.info

    def __array_wrap__(self, obj, context=None):

        if context is None:
            # Methods like .squeeze() created a new `ndarray` and then call
            # __array_wrap__ to turn the array into self's subclass.
            return self._new_view(obj)

        raise NotImplementedError('__array_wrap__ should not be used '
                                  'with a context any more, since we require '
                                  'numpy >=1.13.  Please raise an issue on '
                                  'https://github.com/astropy/astropy')

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
        result : `~astropy.units.Quantity`
            Results of the ufunc, with the unit set properly.
        """
        # Determine required conversion functions -- to bring the unit of the
        # input to that expected (e.g., radian for np.sin), or to get
        # consistent units between two inputs (e.g., in np.add) --
        # and the unit of the result (or tuple of units for nout > 1).
        converters, unit = converters_and_unit(function, method, *inputs)

        out = kwargs.get('out', None)
        # Avoid loop back by turning any Quantity output into array views.
        if out is not None:
            # If pre-allocated output is used, check it is suitable.
            # This also returns array view, to ensure we don't loop back.
            if function.nout == 1:
                out = out[0]
            out_array = check_output(out, unit, inputs, function=function)
            # Ensure output argument remains a tuple.
            kwargs['out'] = (out_array,) if function.nout == 1 else out_array

        # Same for inputs, but here also convert if necessary.
        arrays = [(converter(input_.value) if converter else
                   getattr(input_, 'value', input_))
                  for input_, converter in zip(inputs, converters)]

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

    def _result_as_quantity(self, result, unit, out):
        """Turn result into a quantity with the given unit.

        If no output is given, it will take a view of the array as a quantity,
        and set the unit.  If output is given, those should be quantity views
        of the result arrays, and the function will just set the unit.

        Parameters
        ----------
        result : `~numpy.ndarray` or tuple of `~numpy.ndarray`
            Array(s) which need to be turned into quantity.
        unit : `~astropy.units.Unit` or None
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
        if isinstance(result, tuple):
            if out is None:
                out = (None,) * len(result)
            return tuple(self._result_as_quantity(result_, unit_, out_)
                         for (result_, unit_, out_) in
                         zip(result, unit, out))

        if out is None:
            # View the result array as a Quantity with the proper unit.
            return result if unit is None else self._new_view(result, unit)

        # For given output, just set the unit. We know the unit is not None and
        # the output is of the correct Quantity subclass, as it was passed
        # through check_output.
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
            - `Quantity` subclass
            - bool: True if subclasses of the given class are ok
        """
        return Quantity, True

    def _new_view(self, obj=None, unit=None):
        """
        Create a Quantity view of some array-like input, and set the unit

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

        unit : `UnitBase`, or anything convertible to a :class:`~astropy.units.Unit`, optional
            The unit of the resulting object.  It is used to select a
            subclass, and explicitly assigned to the view if given.
            If not given, the subclass and unit will be that of ``self``.

        Returns
        -------
        view : Quantity subclass
        """
        # Determine the unit and quantity subclass that we need for the view.
        if unit is None:
            unit = self.unit
            quantity_subclass = self.__class__
        else:
            # In principle, could gain time by testing unit is self.unit
            # as well, and then quantity_subclass = self.__class__, but
            # Constant relies on going through `__quantity_subclass__`.
            unit = Unit(unit)
            quantity_subclass = getattr(unit, '_quantity_class', Quantity)
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
            obj = np.array(obj, copy=False)

        # Take the view, set the unit, and update possible other properties
        # such as ``info``, ``wrap_angle`` in `Longitude`, etc.
        view = obj.view(quantity_subclass)
        view._set_unit(unit)
        view.__array_finalize__(self)
        return view

    def _set_unit(self, unit):
        """Set the unit.

        This is used anywhere the unit is set or modified, i.e., in the
        initilizer, in ``__imul__`` and ``__itruediv__`` for in-place
        multiplication and division by another unit, as well as in
        ``__array_finalize__`` for wrapping up views.  For Quantity, it just
        sets the unit, but subclasses can override it to check that, e.g.,
        a unit is consistent.
        """
        if not isinstance(unit, UnitBase):
            # Trying to go through a string ensures that, e.g., Magnitudes with
            # dimensionless physical unit become Quantity with units of mag.
            unit = Unit(str(unit), parse_strict='silent')
            if not isinstance(unit, UnitBase):
                raise UnitTypeError(
                    "{0} instances require {1} units, not {2} instances."
                    .format(type(self).__name__, UnitBase, type(unit)))

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
        return self.unit.to(unit, self.view(np.ndarray),
                            equivalencies=equivalencies)

    def to(self, unit, equivalencies=[]):
        """
        Return a new `~astropy.units.Quantity` object with the specified unit.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to convert to. Must be
            an `~astropy.units.UnitBase` object or a string parseable
            by the `~astropy.units` package.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            If not provided or ``[]``, class default equivalencies will be used
            (none for `~astropy.units.Quantity`, but may be set for subclasses)
            If `None`, no equivalencies will be applied at all, not even any
            set globally or within a context.

        See also
        --------
        to_value : get the numerical value in a given unit.
        """
        # We don't use `to_value` below since we always want to make a copy
        # and don't want to slow down this method (esp. the scalar case).
        unit = Unit(unit)
        return self._new_view(self._to_value(unit, equivalencies), unit)

    def to_value(self, unit=None, equivalencies=[]):
        """
        The numerical value, possibly in a different unit.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance or str, optional
            The unit in which the value should be given. If not given or `None`,
            use the current unit.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not directly
            convertible (see :ref:`unit_equivalencies`). If not provided or
            ``[]``, class default equivalencies will be used (none for
            `~astropy.units.Quantity`, but may be set for subclasses).
            If `None`, no equivalencies will be applied at all, not even any
            set globally or within a context.

        Returns
        -------
        value : `~numpy.ndarray` or scalar
            The value in the units specified. For arrays, this will be a view
            of the data if no unit conversion was necessary.

        See also
        --------
        to : Get a new instance in a different unit.
        """
        if unit is None or unit is self.unit:
            value = self.view(np.ndarray)
        else:
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

        return value if self.shape else (value[()] if self.dtype.fields
                                         else value.item())

    value = property(to_value,
                     doc="""The numerical value of this instance.

    See also
    --------
    to_value : Get the numerical value in a given unit.
    """)

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

    @property
    def si(self):
        """
        Returns a copy of the current `Quantity` instance with SI units. The
        value of the resulting object will be scaled.
        """
        si_unit = self.unit.si
        return self._new_view(self.value * si_unit.scale,
                              si_unit / si_unit.scale)

    @property
    def cgs(self):
        """
        Returns a copy of the current `Quantity` instance with CGS units. The
        value of the resulting object will be scaled.
        """
        cgs_unit = self.unit.cgs
        return self._new_view(self.value * cgs_unit.scale,
                              cgs_unit / cgs_unit.scale)

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

    @override__dir__
    def __dir__(self):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.  This function is implemented in
        order to make autocompletion still work correctly in IPython.
        """
        if not self._include_easy_conversion_members:
            return []
        extra_members = set()
        equivalencies = Unit._normalize_equivalencies(self.equivalencies)
        for equivalent in self.unit._get_units_with_same_physical_type(
                equivalencies):
            extra_members.update(equivalent.names)
        return extra_members

    def __getattr__(self, attr):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.
        """
        if not self._include_easy_conversion_members:
            raise AttributeError(
                "'{0}' object has no '{1}' member".format(
                    self.__class__.__name__,
                    attr))

        def get_virtual_unit_attribute():
            registry = get_current_unit_registry().registry
            to_unit = registry.get(attr, None)
            if to_unit is None:
                return None

            try:
                return self.unit.to(
                    to_unit, self.value, equivalencies=self.equivalencies)
            except UnitsError:
                return None

        value = get_virtual_unit_attribute()

        if value is None:
            raise AttributeError(
                "{0} instance has no attribute '{1}'".format(
                    self.__class__.__name__, attr))
        else:
            return value

    # Equality (return False if units do not match) needs to be handled
    # explicitly for numpy >=1.9, since it no longer traps errors.
    def __eq__(self, other):
        try:
            try:
                return super().__eq__(other)
            except DeprecationWarning:
                # We treat the DeprecationWarning separately, since it may
                # mask another Exception.  But we do not want to just use
                # np.equal, since super's __eq__ treats recarrays correctly.
                return np.equal(self, other)
        except UnitsError:
            return False
        except TypeError:
            return NotImplemented

    def __ne__(self, other):
        try:
            try:
                return super().__ne__(other)
            except DeprecationWarning:
                return np.not_equal(self, other)
        except UnitsError:
            return True
        except TypeError:
            return NotImplemented

    # Unit conversion operator (<<).
    def __lshift__(self, other):
        try:
            other = Unit(other, parse_strict='silent')
        except UnitTypeError:
            return NotImplemented

        return self.__class__(self, other, copy=False, subok=True)

    def __ilshift__(self, other):
        try:
            other = Unit(other, parse_strict='silent')
        except UnitTypeError:
            return NotImplemented

        try:
            factor = self.unit._to(other)
        except UnitConversionError:
            # Maybe via equivalencies?  Now we do make a temporary copy.
            try:
                value = self._to_value(other)
            except UnitConversionError:
                return NotImplemented

            self.view(np.ndarray)[...] = value

        else:
            self.view(np.ndarray)[...] *= factor

        self._set_unit(other)
        return self

    def __rlshift__(self, other):
        if not self.isscalar:
            return NotImplemented
        return Unit(self).__rlshift__(other)

    # Give warning for other >> self, since probably other << self was meant.
    def __rrshift__(self, other):
        warnings.warn(">> is not implemented. Did you mean to convert "
                      "something to this quantity as a unit using '<<'?",
                      AstropyWarning)
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
        """ Multiplication between `Quantity` objects and other objects."""

        if isinstance(other, (UnitBase, str)):
            try:
                return self._new_view(self.copy(), other * self.unit)
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
        """ Right Multiplication between `Quantity` objects and other
        objects.
        """

        return self.__mul__(other)

    def __truediv__(self, other):
        """ Division between `Quantity` objects and other objects."""

        if isinstance(other, (UnitBase, str)):
            try:
                return self._new_view(self.copy(), self.unit / other)
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
        """ Right Division between `Quantity` objects and other objects."""

        if isinstance(other, (UnitBase, str)):
            return self._new_view(1. / self.value, other / self.unit)

        return super().__rtruediv__(other)

    def __div__(self, other):
        """ Division between `Quantity` objects. """
        return self.__truediv__(other)

    def __idiv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__itruediv__(other)

    def __rdiv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__rtruediv__(other)

    def __pow__(self, other):
        if isinstance(other, Fraction):
            # Avoid getting object arrays by raising the value to a Fraction.
            return self._new_view(self.value ** float(other),
                                  self.unit ** other)

        return super().__pow__(other)

    # For Py>=3.5
    def __matmul__(self, other, reverse=False):
        result_unit = self.unit * getattr(other, 'unit', dimensionless_unscaled)
        result_array = np.matmul(self.value, getattr(other, 'value', other))
        return self._new_view(result_array, result_unit)

    def __rmatmul__(self, other):
        result_unit = self.unit * getattr(other, 'unit', dimensionless_unscaled)
        result_array = np.matmul(getattr(other, 'value', other), self.value)
        return self._new_view(result_array, result_unit)

    # In numpy 1.13, 1.14, a np.positive ufunc exists, but ndarray.__pos__
    # does not go through it, so we define it, to allow subclasses to override
    # it inside __array_ufunc__. This can be removed if a solution to
    # https://github.com/numpy/numpy/issues/9081 is merged.
    def __pos__(self):
        """Plus the quantity."""
        return np.positive(self)

    # other overrides of special functions
    def __hash__(self):
        return hash(self.value) ^ hash(self.unit)

    def __iter__(self):
        if self.isscalar:
            raise TypeError(
                "'{cls}' object with a scalar value is not iterable"
                .format(cls=self.__class__.__name__))

        # Otherwise return a generator
        def quantity_iter():
            for val in self.value:
                yield self._new_view(val)

        return quantity_iter()

    def __getitem__(self, key):
        try:
            out = super().__getitem__(key)
        except IndexError:
            # We want zero-dimensional Quantity objects to behave like scalars,
            # so they should raise a TypeError rather than an IndexError.
            if self.isscalar:
                raise TypeError(
                    "'{cls}' object with a scalar value does not support "
                    "indexing".format(cls=self.__class__.__name__))
            else:
                raise
        # For single elements, ndarray.__getitem__ returns scalars; these
        # need a new view as a Quantity.
        if type(out) is not type(self):
            out = self._new_view(out)
        return out

    def __setitem__(self, i, value):
        # update indices in info if the info property has been accessed
        # (in which case 'info' in self.__dict__ is True; this is guaranteed
        # to be the case if we're part of a table).
        if not self.isscalar and 'info' in self.__dict__:
            self.info.adjust_indices(i, value, len(self))
        self.view(np.ndarray).__setitem__(i, self._to_own_unit(value))

    # __contains__ is OK

    def __bool__(self):
        """Quantities should always be treated as non-False; there is too much
        potential for ambiguity otherwise.
        """
        warnings.warn('The truth value of a Quantity is ambiguous. '
                      'In the future this will raise a ValueError.',
                      AstropyDeprecationWarning)
        return True

    def __len__(self):
        if self.isscalar:
            raise TypeError("'{cls}' object with a scalar value has no "
                            "len()".format(cls=self.__class__.__name__))
        else:
            return len(self.value)

    # Numerical types
    def __float__(self):
        try:
            return float(self.to_value(dimensionless_unscaled))
        except (UnitsError, TypeError):
            raise TypeError('only dimensionless scalar quantities can be '
                            'converted to Python scalars')

    def __int__(self):
        try:
            return int(self.to_value(dimensionless_unscaled))
        except (UnitsError, TypeError):
            raise TypeError('only dimensionless scalar quantities can be '
                            'converted to Python scalars')

    def __index__(self):
        # for indices, we do not want to mess around with scaling at all,
        # so unlike for float, int, we insist here on unscaled dimensionless
        try:
            assert self.unit.is_unity()
            return self.value.__index__()
        except Exception:
            raise TypeError('only integer dimensionless scalar quantities '
                            'can be converted to a Python index')

    @property
    def _unitstr(self):
        if self.unit is None:
            unitstr = _UNIT_NOT_INITIALISED
        else:
            unitstr = str(self.unit)

        if unitstr:
            unitstr = ' ' + unitstr

        return unitstr

    # Display
    # TODO: we may want to add a hook for dimensionless quantities?
    def __str__(self):
        return '{0}{1:s}'.format(self.value, self._unitstr)

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        sep = ',' if NUMPY_LT_1_14 else ', '
        arrstr = np.array2string(self.view(np.ndarray), separator=sep,
                                 prefix=prefixstr)
        return '{0}{1}{2:s}>'.format(prefixstr, arrstr, self._unitstr)

    def _repr_latex_(self):
        """
        Generate a latex representation of the quantity and its unit.

        The behavior of this function can be altered via the
        `numpy.set_printoptions` function and its various keywords.  The
        exception to this is the ``threshold`` keyword, which is controlled via
        the ``[units.quantity]`` configuration item ``latex_array_threshold``.
        This is treated separately because the numpy default of 1000 is too big
        for most browsers to handle.

        Returns
        -------
        lstr
            A LaTeX string with the contents of this Quantity
        """
        # need to do try/finally because "threshold" cannot be overridden
        # with array2string
        pops = np.get_printoptions()

        format_spec = '.{}g'.format(pops['precision'])

        def float_formatter(value):
            return Latex.format_exponential_notation(value,
                                                     format_spec=format_spec)

        def complex_formatter(value):
            return '({0}{1}i)'.format(
                Latex.format_exponential_notation(value.real,
                                                  format_spec=format_spec),
                Latex.format_exponential_notation(value.imag,
                                                  format_spec='+' + format_spec))

        try:
            formatter = {'float_kind': float_formatter,
                         'complex_kind': complex_formatter}
            if conf.latex_array_threshold > -1:
                np.set_printoptions(threshold=conf.latex_array_threshold,
                                    formatter=formatter)

            # the view is needed for the scalar case - value might be float
            if NUMPY_LT_1_14:   # style deprecated in 1.14
                latex_value = np.array2string(
                    self.view(np.ndarray),
                    style=(float_formatter if self.dtype.kind == 'f'
                           else complex_formatter if self.dtype.kind == 'c'
                           else repr),
                    max_line_width=np.inf, separator=',~')
            else:
                latex_value = np.array2string(
                    self.view(np.ndarray),
                    max_line_width=np.inf, separator=',~')

            latex_value = latex_value.replace('...', r'\dots')
        finally:
            np.set_printoptions(**pops)

        # Format unit
        # [1:-1] strips the '$' on either side needed for math mode
        latex_unit = (self.unit._repr_latex_()[1:-1]  # note this is unicode
                      if self.unit is not None
                      else _UNIT_NOT_INITIALISED)

        return r'${0} \; {1}$'.format(latex_value, latex_unit)

    def __format__(self, format_spec):
        """
        Format quantities using the new-style python formatting codes
        as specifiers for the number.

        If the format specifier correctly applies itself to the value,
        then it is used to format only the value. If it cannot be
        applied to the value, then it is applied to the whole string.

        """
        try:
            value = format(self.value, format_spec)
            full_format_spec = "s"
        except ValueError:
            value = self.value
            full_format_spec = format_spec

        return format("{0}{1:s}".format(value, self._unitstr),
                      full_format_spec)

    def decompose(self, bases=[]):
        """
        Generates a new `Quantity` with the units
        decomposed. Decomposed units have only irreducible units in
        them (see `astropy.units.UnitBase.decompose`).

        Parameters
        ----------
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
        if not allowscaledunits and hasattr(new_unit, 'scale'):
            new_value = self.value * new_unit.scale
            new_unit = new_unit / new_unit.scale
            return self._new_view(new_value, new_unit)
        else:
            return self._new_view(self.copy(), new_unit)

    # These functions need to be overridden to take into account the units
    # Array conversion
    # http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html#array-conversion

    def item(self, *args):
        return self._new_view(super().item(*args))

    def tolist(self):
        raise NotImplementedError("cannot make a list of Quantities.  Get "
                                  "list of values with q.value.list()")

    def _to_own_unit(self, value, check_precision=True):
        try:
            _value = value.to_value(self.unit)
        except AttributeError:
            # We're not a Quantity, so let's try a more general conversion.
            # Plain arrays will be converted to dimensionless in the process,
            # but anything with a unit attribute will use that.
            try:
                _value = Quantity(value).to_value(self.unit)
            except UnitsError as exc:
                # last chance: if this was not something with a unit
                # and is all 0, inf, or nan, we treat it as arbitrary unit.
                if (not hasattr(value, 'unit') and
                        can_have_arbitrary_unit(value)):
                    _value = value
                else:
                    raise exc

        if check_precision:
            value_dtype = getattr(value, 'dtype', None)
            if self.dtype != value_dtype:
                self_dtype_array = np.array(_value, self.dtype)
                value_dtype_array = np.array(_value, dtype=value_dtype,
                                             copy=False)
                if not np.all(np.logical_or(self_dtype_array ==
                                            value_dtype_array,
                                            np.isnan(value_dtype_array))):
                    raise TypeError("cannot convert value type to array type "
                                    "without precision loss")
        return _value

    def itemset(self, *args):
        if len(args) == 0:
            raise ValueError("itemset must have at least one argument")

        self.view(np.ndarray).itemset(*(args[:-1] +
                                        (self._to_own_unit(args[-1]),)))

    def tostring(self, order='C'):
        raise NotImplementedError("cannot write Quantities to string.  Write "
                                  "array with q.value.tostring(...).")

    def tofile(self, fid, sep="", format="%s"):
        raise NotImplementedError("cannot write Quantities to file.  Write "
                                  "array with q.value.tofile(...)")

    def dump(self, file):
        raise NotImplementedError("cannot dump Quantities to file.  Write "
                                  "array with q.value.dump()")

    def dumps(self):
        raise NotImplementedError("cannot dump Quantities to string.  Write "
                                  "array with q.value.dumps()")

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
    # take, repeat, sort, compress, diagonal OK
    def put(self, indices, values, mode='raise'):
        self.view(np.ndarray).put(indices, self._to_own_unit(values), mode)

    def choose(self, choices, out=None, mode='raise'):
        raise NotImplementedError("cannot choose based on quantity.  Choose "
                                  "using array with q.value.choose(...)")

    # ensure we do not return indices as quantities
    def argsort(self, axis=-1, kind='quicksort', order=None):
        return self.view(np.ndarray).argsort(axis=axis, kind=kind, order=order)

    def searchsorted(self, v, *args, **kwargs):
        return np.searchsorted(np.array(self),
                               self._to_own_unit(v, check_precision=False),
                               *args, **kwargs)  # avoid numpy 1.6 problem

    def argmax(self, axis=None, out=None):
        return self.view(np.ndarray).argmax(axis, out=out)

    def argmin(self, axis=None, out=None):
        return self.view(np.ndarray).argmin(axis, out=out)

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
        args = (self.value,) + tuple((arg.value if isinstance(arg, Quantity)
                                      else arg) for arg in args)
        if out is not None:
            # If pre-allocated output is used, check it is suitable.
            # This also returns array view, to ensure we don't loop back.
            arrays = tuple(arg for arg in args if isinstance(arg, np.ndarray))
            kwargs['out'] = check_output(out, unit, arrays, function=function)
        # Apply the function and turn it back into a Quantity.
        result = function(*args, **kwargs)
        return self._result_as_quantity(result, unit, out)

    def clip(self, a_min, a_max, out=None):
        return self._wrap_function(np.clip, self._to_own_unit(a_min),
                                   self._to_own_unit(a_max), out=out)

    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        return self._wrap_function(np.trace, offset, axis1, axis2, dtype,
                                   out=out)

    def var(self, axis=None, dtype=None, out=None, ddof=0):
        return self._wrap_function(np.var, axis, dtype,
                                   out=out, ddof=ddof, unit=self.unit**2)

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        return self._wrap_function(np.std, axis, dtype, out=out, ddof=ddof)

    def mean(self, axis=None, dtype=None, out=None):
        return self._wrap_function(np.mean, axis, dtype, out=out)

    def ptp(self, axis=None, out=None):
        return self._wrap_function(np.ptp, axis, out=out)

    def round(self, decimals=0, out=None):
        return self._wrap_function(np.round, decimals, out=out)

    def max(self, axis=None, out=None, keepdims=False):
        return self._wrap_function(np.max, axis, out=out, keepdims=keepdims)

    def min(self, axis=None, out=None, keepdims=False):
        return self._wrap_function(np.min, axis, out=out, keepdims=keepdims)

    def sum(self, axis=None, dtype=None, out=None, keepdims=False):
        return self._wrap_function(np.sum, axis, dtype, out=out,
                                   keepdims=keepdims)

    def prod(self, axis=None, dtype=None, out=None, keepdims=False):
        if not self.unit.is_unity():
            raise ValueError("cannot use prod on scaled or "
                             "non-dimensionless Quantity arrays")
        return self._wrap_function(np.prod, axis, dtype, out=out,
                                   keepdims=keepdims)

    def dot(self, b, out=None):
        result_unit = self.unit * getattr(b, 'unit', dimensionless_unscaled)
        return self._wrap_function(np.dot, b, out=out, unit=result_unit)

    def cumsum(self, axis=None, dtype=None, out=None):
        return self._wrap_function(np.cumsum, axis, dtype, out=out)

    def cumprod(self, axis=None, dtype=None, out=None):
        if not self.unit.is_unity():
            raise ValueError("cannot use cumprod on scaled or "
                             "non-dimensionless Quantity arrays")
        return self._wrap_function(np.cumprod, axis, dtype, out=out)

    # Calculation: override methods that do not make sense.

    def all(self, axis=None, out=None):
        raise NotImplementedError("cannot evaluate truth value of quantities. "
                                  "Evaluate array with q.value.all(...)")

    def any(self, axis=None, out=None):
        raise NotImplementedError("cannot evaluate truth value of quantities. "
                                  "Evaluate array with q.value.any(...)")

    # Calculation: numpy functions that can be overridden with methods.

    def diff(self, n=1, axis=-1):
        return self._wrap_function(np.diff, n, axis)

    def ediff1d(self, to_end=None, to_begin=None):
        return self._wrap_function(np.ediff1d, to_end, to_begin)

    def nansum(self, axis=None, out=None, keepdims=False):
        return self._wrap_function(np.nansum, axis,
                                   out=out, keepdims=keepdims)

    def insert(self, obj, values, axis=None):
        """
        Insert values along the given axis before the given indices and return
        a new `~astropy.units.Quantity` object.

        This is a thin wrapper around the `numpy.insert` function.

        Parameters
        ----------
        obj : int, slice or sequence of ints
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
    # without units are possible, but will have an uninitalized unit.
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
                "{0} instances require units equivalent to '{1}'"
                .format(type(self).__name__, self._equivalent_unit) +
                (", but no unit was given." if unit is None else
                 ", so cannot set it to '{0}'.".format(unit)))

        super()._set_unit(unit)


def isclose(a, b, rtol=1.e-5, atol=None, **kwargs):
    """
    Notes
    -----
    Returns True if two arrays are element-wise equal within a tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.isclose`.
    """
    return np.isclose(*_unquantify_allclose_arguments(a, b, rtol, atol),
                      **kwargs)


def allclose(a, b, rtol=1.e-5, atol=None, **kwargs):
    """
    Notes
    -----
    Returns True if two arrays are element-wise equal within a tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.allclose`.
    """
    return np.allclose(*_unquantify_allclose_arguments(a, b, rtol, atol),
                       **kwargs)


def _unquantify_allclose_arguments(actual, desired, rtol, atol):
    actual = Quantity(actual, subok=True, copy=False)

    desired = Quantity(desired, subok=True, copy=False)
    try:
        desired = desired.to(actual.unit)
    except UnitsError:
        raise UnitsError("Units for 'desired' ({0}) and 'actual' ({1}) "
                         "are not convertible"
                         .format(desired.unit, actual.unit))

    if atol is None:
        # by default, we assume an absolute tolerance of 0
        atol = Quantity(0)
    else:
        atol = Quantity(atol, subok=True, copy=False)
        try:
            atol = atol.to(actual.unit)
        except UnitsError:
            raise UnitsError("Units for 'atol' ({0}) and 'actual' ({1}) "
                             "are not convertible"
                             .format(atol.unit, actual.unit))

    rtol = Quantity(rtol, subok=True, copy=False)
    try:
        rtol = rtol.to(dimensionless_unscaled)
    except Exception:
        raise UnitsError("`rtol` should be dimensionless")

    return actual.value, desired.value, rtol.value, atol.value
