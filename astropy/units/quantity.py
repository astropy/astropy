# Coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines the `Quantity` object, which represents a number with some
associated units. `Quantity` objects support operations like ordinary numbers,
but will deal with unit conversions internally.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import numbers

import numpy as np

# AstroPy
from .core import Unit, UnitBase, UnitsException
from ..utils import lazyproperty
from ..utils.compat.misc import override__dir__


__all__ = ["Quantity"]


from .quantity_helper import _is_unity

UNIT_NOT_INITIALISED = "(Unit not initialised)"


def _validate_value(value):
    """ Make sure that the input is a Python or Numpy numeric type.

    Parameters
    ----------
    value : number
        An object that will be checked whether it is a numeric type or not.

    Returns
    -------
    newval
        The new value either as an array or a scalar
    """

    from ..utils.misc import isiterable

    if (isinstance(value, (numbers.Number, np.number, np.ndarray)) or
            isiterable(value)):
        value_obj = np.array(value, copy=True)
    else:
        raise TypeError("The value must be a valid Python or Numpy numeric "
                        "type.")

    return value_obj


class Quantity(np.ndarray):
    """ A `Quantity` represents a number with some associated unit.

    Parameters
    ----------
    value : number, `Quantity` object, or sequence of `Quantity` objects.
        The numerical value of this quantity in the units given by
        unit.  If a `Quantity` or sequence of them, creates a new
        `Quantity` object, converting to `unit` units as needed.

    unit : `~astropy.units.UnitBase` instance, str
        An object that represents the unit associated with the input value.
        Must be an `~astropy.units.UnitBase` object or a string parseable by
        the `units` package.

    equivalencies : list of equivalence pairs, optional
        A list of equivalence pairs. See :ref:`unit_equivalencies`.

    Raises
    ------
    TypeError
        If the value provided is not a Python numeric type.
    TypeError
        If the unit provided is not either a `Unit` object or a parseable
        string unit.
    """
    # Need to set a class-level default for _equivalencies, or
    # Constants can not initialize properly
    _equivalencies = []

    __array_priority__ = 10000

    def __new__(cls, value, unit=None, dtype=None, equivalencies=[]):

        from ..utils.misc import isiterable

        if isinstance(value, Quantity):
            _value = _validate_value(value.to(unit).value)
        elif isiterable(value) and all(isinstance(v, Quantity) for v in value):
            _value = _validate_value([q.to(unit).value for q in value])
        else:
            _value = _validate_value(value)

        if dtype is not None and dtype != _value.dtype:
            _value = _value.astype(dtype)
        else:
            dtype = _value.dtype

        self = super(Quantity, cls).__new__(cls, _value.shape, dtype=dtype,
                                            buffer=_value.data)
        if unit is None:
            self._unit = Unit(1)
        else:
            self._unit = Unit(unit)
        self._equivalencies = Unit._normalize_equivalencies(equivalencies)

        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return
        elif isinstance(obj, Quantity):
            self._unit = obj._unit

    def __array_prepare__(self, obj, context=None):

        # This method gets called by Numpy whenever a ufunc is called on the
        # array. The object passed in ``obj`` is an empty version of the
        # output array which we can e.g. change to an array sub-class, add
        # attributes to, etc. After this is called, then the ufunc is called
        # and the values in this empty array are set.

        # If no context is set, just return the input
        if context is None:
            return obj

        # Find out which ufunc is being used
        function = context[0]

        from .quantity_helper import UNSUPPORTED_UFUNCS, UFUNC_HELPERS

        # Check whether we even support this ufunc
        if function in UNSUPPORTED_UFUNCS:
            raise TypeError("Cannot use function '{0}' with quantities"
                            .format(function.__name__))

        # Now find out what arguments were passed to the ufunc, usually, this
        # will include at least the present object, and another, which could
        # be a Quantity, or a Numpy array, etc. when using two-argument ufuncs.
        args = context[1][:function.nin]
        units = [getattr(arg, 'unit', None) for arg in args]

        # If the ufunc is supported, then we call a helper function (defined
        # in quantity_helper.py) which returns the scale by which the inputs
        # should be multiplied before being passed to the ufunc, as well as
        # the unit the output from the ufunc will have.
        if function in UFUNC_HELPERS:
            scales, result_unit = UFUNC_HELPERS[function](function, *units)
        else:
            raise TypeError("Unknown ufunc {0}.  Please raise issue on "
                            "https://github.com/astropy/astropy"
                            .format(function.__name__))

        # In the case of np.power, the unit itself needs to be modified by an
        # amount that depends on one of the input values, so we need to treat
        # this as a special case.
        # TODO: find a better way to deal with this case
        if function is np.power and result_unit is not None:
            if units[1] is None:
                result_unit = result_unit ** args[1]
            else:
                from . import dimensionless_unscaled
                result_unit = result_unit ** args[1].to(dimensionless_unscaled)

        # We now prepare the output object

        if self is obj:  # happens if the output object is self, which happens
                         # for in-place operations such as q1 += q2

            # In some cases, the result of a ufunc should be a plain Numpy
            # array, which we can't do if we are doing an in-place operation.
            if result_unit is None:
                raise TypeError("Cannot store non-quantity output from {0} "
                                "function in Quantity object"
                                .format(function.__name__))

            # If the Quantity has an integer dtype, in-place operations are
            # dangerous because in some cases the quantity will be e.g.
            # decomposed, which involves being scaled by a float, but since
            # the array is an integer the output then gets converted to an int
            # and truncated.
            if(any(not np.can_cast(arg, obj.dtype) for arg in args) or
               np.any(np.array(scales, dtype=obj.dtype) != np.array(scales))):
                raise TypeError("Arguments cannot be cast safely to inplace "
                                "output with dtype={0}".format(self.dtype))

            result = self  # no need for a view since we are returning the object itself

            # in principle, if self is also an argument, it could be rescaled
            # here, since it won't be needed anymore.  But maybe not change
            # inputs before the calculation even if they will get destroyed

        else:  # normal case: set up output as a Quantity

            # We should use Quantity here otherwise we might re-cast as e.g.
            # an EMConstant, which is not what is wanted.
            result = obj.view(Quantity)

        # We now need to treat the case where the inputs have to be scaled -
        # the issue is that we can't actually scale the inputs since that
        # would be changing the objects passed to the ufunc, which would not
        # be expected by the user.
        if any(scale != 1. for scale in scales):

            # If self is both output and input (which happens for in-place
            # operations), input will get overwritten with junk. To avoid
            # that, hide it in a new object
            if self is obj and any(self is arg for arg in args):
                # but with two outputs it would become unhidden too soon
                # [ie., np.modf(q1, q1, other)].  Bail.
                if context[2] < function.nout-1:
                    raise TypeError("Cannot apply multi-output {0} function "
                                    "to quantities with in-place replacement "
                                    "of an input by any but the last output."
                                    .format(function.__name__))
                # cannot use self.copy(), as it can change dtype
                result = self.view(np.ndarray).copy().view(Quantity)
                result._result = self

            # ensure we remember the scales we need
            result._scales = scales

        # unit output will get (setting _unit could prematurely change input)
        result._result_unit = result_unit
        return result

    def __array_wrap__(self, obj, context=None):

        if context is not None:

            if hasattr(obj, '_result_unit'):
                result_unit = obj._result_unit
                del obj._result_unit
            else:
                result_unit = None

            # We now need to re-calculate quantities for which the input
            # needed to be scaled.
            if hasattr(obj, '_scales'):

                scales = obj._scales
                del obj._scales

                # since values were not initially scaled, the output is junk,
                # so we can use it as storage for scaled input; (array view,
                # to ensure pure ndarray recalculation)

                if hasattr(obj, '_result'):
                    # real output was also one of the inputs (happend when
                    # using an in-place ufunc), so was hidden. junked copy
                    # still useful as storage for scaled inputs and possible
                    # junk output
                    junk = obj.view(np.ndarray)
                    storage = [junk]
                    # Retrieve real object
                    obj = obj._result
                else:
                    junk = None
                    storage = []

                # take array view to which output can be written without
                # getting back here
                obj_array = obj.view(np.ndarray)
                storage.append(obj_array)

                # Find out which ufunc was called and with which inputs
                function = context[0]
                args = context[1][:function.nin]

                # set the inputs, rescaling as necessary
                inputs = []
                for arg, scale in zip(args, scales):
                    if scale != 1.:
                        # if possible, use storage for scaled input Quantity
                        if(storage != [] and
                           storage[0].shape == arg.shape and
                           storage[0].dtype == arg.dtype and
                           np.array(scale, storage[0].dtype) == scale):
                            inputs.append(np.multiply(arg.view(np.ndarray),
                                                      scale,
                                                      out=storage.pop(0)))
                        else:
                            inputs.append(arg.value * scale)
                    else:  # for scale==1, input is not necessarily a Quantity
                        inputs.append(getattr(arg, 'value', arg))

                # scaling may change whether or not output can be written to
                if(result_unit is not None and
                   any(not np.can_cast(scaled_arg, obj_array.dtype)
                       for scaled_arg in inputs)):
                    obj_array = None

                if function.nin == 1:
                    if function.nout == 1:
                        out = function(inputs[0], obj_array)
                    else:  # 2-output function (np.modf, np.frexp); 1 input
                        if context[2] == 0:
                            out, _ = function(inputs[0], obj_array, junk)
                        else:
                            _, out = function(inputs[0], junk, obj_array)
                else:
                    out = function(inputs[0], inputs[1], obj_array)

                if obj_array is None:
                    if type(out) != np.ndarray:  # array scalar; cannot view
                        return Quantity(out, result_unit)
                    else:
                        obj = out.view(Quantity)

            if result_unit is None:  # we're returning a plain array
                obj = obj.view(np.ndarray)
            else:
                obj._unit = result_unit

        return obj

    def to(self, unit, equivalencies=None):
        """ Returns a new `Quantity` object with the specified units.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to convert to. Must be
            an `~astropy.units.UnitBase` object or a string parseable
            by the `units` package.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.  If
            not provided, the equivalencies that were provided in the
            constructor will be used.
        """
        if equivalencies is None:
            equivalencies = self._equivalencies
        new_val = self.unit.to(unit, self.value, equivalencies=equivalencies)
        new_unit = Unit(unit)
        return Quantity(new_val, new_unit)

    @property
    def value(self):
        """ The numerical value of this quantity. """
        value = self.view(np.ndarray)
        if self.shape:
            return value
        else:
            return value.item()

    @property
    def unit(self):
        """
        A `~astropy.units.UnitBase` object representing the unit of this
        quantity.
        """

        return self._unit

    # this ensures that if we do a view, __repr__ and __str__ do not balk
    _unit = None

    @property
    def equivalencies(self):
        """
        A list of equivalencies that will be applied implicitly during
        unit conversions.
        """

        return self._equivalencies

    @property
    def si(self):
        """
        Returns a copy of the current `Quantity` instance with SI units. The
        value of the resulting object will be scaled.
        """

        from . import si
        si_unit = self.unit.to_system(si)[0]
        return Quantity(self.value * si_unit.scale, si_unit / si_unit.scale)

    @property
    def cgs(self):
        """
        Returns a copy of the current `Quantity` instance with CGS units. The
        value of the resulting object will be scaled.
        """

        from . import cgs
        cgs_unit = self.unit.to_system(cgs)[0]
        return Quantity(self.value * cgs_unit.scale, cgs_unit / cgs_unit.scale)

    @lazyproperty
    def isscalar(self):
        """
        True if the `value` of this quantity is a scalar, or False if it
        is an array-like object.

        .. note::
            This is subtly different from `numpy.isscalar` in that
            `numpy.isscalar` returns False for a zero-dimensional array
            (e.g. ``np.array(1)``), while this is True in that case.
        """

        from ..utils.misc import isiterable

        return not isiterable(self.value)

    def copy(self):
        """ Return a copy of this `Quantity` instance """

        return self.__class__(self, unit=self.unit)

    @override__dir__
    def __dir__(self):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.  This function is implemented in
        order to make autocompletion still work correctly in IPython.
        """
        extra_members = set()
        for equivalent in self.unit._get_units_with_same_physical_type(
                self._equivalencies):
            if len(equivalent.aliases):
                name = equivalent.aliases[0]
            else:
                name = equivalent.name
            extra_members.add(name)
        return extra_members

    def __getattr__(self, attr):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.
        """
        def get_virtual_unit_attribute():
            try:
                to_unit = Unit(attr)
            except ValueError:
                return None

            if len(to_unit.aliases):
                if to_unit.aliases[0] != attr:
                    return None
            else:
                if to_unit.name != attr:
                    return None

            try:
                return self.unit.to(
                    to_unit, self.value, equivalencies=self.equivalencies)
            except UnitsException:
                return None

        value = get_virtual_unit_attribute()

        if value is None:
            raise AttributeError(
                "{0} instance has no attribute '{1}'".format(
                    self.__class__.__name__, attr))
        else:
            return value

    # Arithmetic operations
    def __mul__(self, other):
        """ Multiplication between `Quantity` objects and other objects."""

        if isinstance(other, basestring):
            return Quantity(self.value, unit=Unit(other) * self.unit)
        elif isinstance(other, UnitBase):
            return Quantity(self.value, unit=other * self.unit)
        else:
            return np.multiply(self, other)

    def __imul__(self, other):
        """In-place multiplication between `Quantity` objects and others."""
        if isinstance(other, basestring):
            self._unit = Unit(other) * self.unit
        elif isinstance(other, UnitBase):
            self._unit = other * self.unit
        else:
            return np.multiply(self, other, out=self)

        return self

    def __rmul__(self, other):
        """ Right Multiplication between `Quantity` objects and other
        objects.
        """

        return self.__mul__(other)

    def __div__(self, other):
        """ Division between `Quantity` objects and other objects."""

        if isinstance(other, basestring):
            return Quantity(self.value, unit=self.unit / Unit(other))
        elif isinstance(other, UnitBase):
            return Quantity(self.value, unit=self.unit / other)
        else:
            return np.true_divide(self, other)

    def __idiv__(self, other):
        """Inplace division between `Quantity` objects and other objects."""

        if isinstance(other, basestring):
            self._unit = self.unit / Unit(other)
        elif isinstance(other, UnitBase):
            self._unit = self.unit / other
        else:
            return np.true_divide(self, other, out=self)

        return self

    def __rdiv__(self, other):
        """ Right Division between `Quantity` objects and other objects."""

        if isinstance(other, basestring):
            return Quantity(self.value, unit=Unit(other) / self.unit)
        elif isinstance(other, UnitBase):
            return Quantity(1. / self.value, unit=other / self.unit)
        else:
            return np.divide(other, self)

    def __truediv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__div__(other)

    def __itruediv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__idiv__(other)

    def __rtruediv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__rdiv__(other)

    def __divmod__(self, other):
        from . import dimensionless_unscaled
        other_value = self._to_own_unit(other)
        result_tuple = super(Quantity, self.__class__).__divmod__(
            self.view(np.ndarray), other_value)

        return (Quantity(result_tuple[0], dimensionless_unscaled),
                Quantity(result_tuple[1], self.unit))

    def __pos__(self):
        """
        Plus the quantity. This is implemented in case users use +q where q is
        a quantity.  (Required for scalar case.)
        """

        return Quantity(self.value, unit=self.unit)

    # Comparison operations
    def __eq__(self, other):
        try:
            return np.equal(self, other)
        except Exception as exc:
            if isinstance(other, Quantity):
                raise exc
            return False

    def __ne__(self, other):
        try:
            return np.not_equal(self, other)
        except Exception as exc:
            if isinstance(other, Quantity):
                raise exc
            return True

    #other overrides of special functions
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
                yield Quantity(val, unit=self.unit)

        return quantity_iter()

    def __getitem__(self, key):
        if self.isscalar:
            raise TypeError(
                "'{cls}' object with a scalar value does not support "
                "indexing".format(cls=self.__class__.__name__))
        else:
            return Quantity(self.value[key], unit=self.unit)

    def __setitem__(self, i, value):
        self.view(np.ndarray).__setitem__(i, self._to_own_unit(value))

    def __setslice__(self, i, j, value):
        self.view(np.ndarray).__setslice__(i, j, self._to_own_unit(value))

    # __contains__ is OK

    def __nonzero__(self):
        """Quantities should always be treated as non-False; there is too much
        potential for ambiguity otherwise.
        """

        return True

    def __len__(self):
        if self.isscalar:
            raise TypeError("'{cls}' object with a scalar value has no "
                            "len()".format(cls=self.__class__.__name__))
        else:
            return len(self.value)

    # Numerical types
    def __float__(self):
        if not self.isscalar or not _is_unity(self.unit):
            raise TypeError('Only dimensionless scalar quantities can be '
                            'converted to Python scalars')
        else:
            return float(self.value)

    def __int__(self):
        if not self.isscalar or not _is_unity(self.unit):
            raise TypeError('Only dimensionless scalar quantities can be '
                            'converted to Python scalars')
        else:
            return int(self.value)

    def __long__(self):
        if not self.isscalar or not _is_unity(self.unit):
            raise TypeError('Only dimensionless scalar quantities can be '
                            'converted to Python scalars')
        else:
            return long(self.value)

    # Display
    # TODO: we may want to add a hook for dimensionless quantities?
    def __str__(self):
        return "{0} {1:s}".format(self.value,
                                  self.unit.to_string() if
                                  self.unit is not None
                                  else UNIT_NOT_INITIALISED)

    def __repr__(self):
        return "<Quantity {0} {1:s}>".format(self.value,
                                             self.unit.to_string() if
                                             self.unit is not None
                                             else UNIT_NOT_INITIALISED)

    def _repr_latex_(self):
        """
        Generate latex representation of unit name.  This is used by
        the IPython notebook to show it all latexified.

        Returns
        -------
        lstr
            LaTeX string
        """

        # Format value
        latex_value = "{0:g}".format(self.value)
        if "e" in latex_value:
            latex_value = latex_value.replace('e', '\\times 10^{') + '}'

        # Format unit
        # [1:-1] strips the '$' on either side needed for math mode
        latex_unit = (self.unit._repr_latex_()[1:-1]  # note this is unicode
                      if self.unit is not None
                      else UNIT_NOT_INITIALISED)

        return u'${0} \; {1}$'.format(latex_value, latex_unit)

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
            This will raises a `UnitsException` if it's not possible
            to do so.

        Returns
        -------
        newq : `~astropy.units.quantity.Quantity`
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
            This will raises a `UnitsException` if it's not possible
            to do so.

        Returns
        -------
        newq : `~astropy.units.quantity.Quantity`
            A new object equal to this quantity with units decomposed.

        """

        new_unit = self.unit.decompose(bases=bases)

        if not allowscaledunits and hasattr(new_unit, 'scale'):
            # Be careful here because self.value might be an array, so if the
            # following is changed, always be sure that the original value is
            # not being modified.
            new_value = self.value * new_unit.scale
            new_unit = new_unit / Unit(new_unit.scale)
        else:
            new_value = self.value

        return Quantity(new_value, new_unit)

    # These functions need to be overridden to take into account the units
    # Array conversion
    # http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html#array-conversion

    def item(self, *args):
        return Quantity(self.view(np.ndarray).item(*args), self.unit)

    def list(self):
        raise NotImplementedError("cannot make a list of Quantities.  Get "
                                  "list of values with q.value.list()")

    def _to_own_unit(self, value, check_precision=True):
        from . import dimensionless_unscaled
        try:
            value = value.to(self.unit).value
        except AttributeError:
            value = dimensionless_unscaled.to(self.unit, value)

        if(check_precision and
           np.any(np.array(value, self.dtype) != np.array(value))):
            raise TypeError("cannot convert value type to array type without "
                            "precision loss")
        return value

    def itemset(self, *args):
        if len(args) == 0:
            raise ValueError("itemset must have at least one argument")

        self.view(np.ndarray).itemset(*(args[:-1] +
                                        (self._to_own_unit(args[1]),)))

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

    # astype, byteswap OK as is
    # copy done above
    # view, getfield, setflags OK as is

    def fill(self, value):
        self.view(np.ndarray).fill(self._to_own_unit(value))

    # Shape manipulation: resize cannot be done (does not own data), but
    # shape, transpose, swapaxes, flatten, ravel, squeeze all OK.

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

    # Calculation

    # ensure we do not return indices as quantities
    # conj OK
    def argmax(self, axis=None, out=None):
        return self.view(np.ndarray).argmax(axis=axis, out=out)

    def argmin(self, axis=None, out=None):
        return self.view(np.ndarray).argmin(axis=axis, out=out)

    def clip(self, a_min, a_max, out=None):
        if out is not None:
            out = out.view(Quantity)
            out._unit = self.unit
        return super(Quantity, self.__class__).clip(self,
                                                    self._to_own_unit(a_min),
                                                    self._to_own_unit(a_max),
                                                    out=out)

    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        if out is not None:
            out = out.view(Quantity)
            out._unit = self.unit
        return Quantity(np.trace(self.value, offset=offset, axis1=axis1,
                                 axis2=axis2, dtype=None, out=out),
                        self.unit)

    def var(self, axis=None, dtype=None, out=None, ddof=0):
        result_unit = self.unit ** 2
        if out is not None:
            out = out.view(Quantity)
            out._unit = result_unit
        return Quantity(np.var(self.value, axis=axis, dtype=dtype, ddof=ddof),
                        result_unit)

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        out = out and out.view(Quantity)
        value = np.std(self.value, axis=axis, dtype=dtype, out=out, ddof=ddof)
        return Quantity(value, self.unit)

    def mean(self, axis=None, dtype=None, out=None):
        out = out and out.view(Quantity)
        value = np.mean(self.value, axis=axis, dtype=dtype, out=out)
        return Quantity(value, self.unit)

    def ptp(self, axis=None, out=None):
        out = out and out.view(Quantity)
        value = np.ptp(self.value, axis=axis, out=out)
        return Quantity(value, self.unit)

    def max(self, axis=None, out=None, keepdims=False):
        out = out and out.view(Quantity)
        try:
            value = np.max(self.value, axis=axis, out=out, keepdims=keepdims)
        except:  # numpy < 1.7
            value = np.max(self.value, axis=axis, out=out)
        return Quantity(value, self.unit)

    def min(self, axis=None, out=None, keepdims=False):
        out = out and out.view(Quantity)
        try:
            value = np.min(self.value, axis=axis, out=out, keepdims=keepdims)
        except:  # numpy < 1.7
            value = np.min(self.value, axis=axis, out=out)
        return Quantity(value, self.unit)

    def dot(self, b, out=None):
        result_unit = self.unit * getattr(b, 'unit', 1.)
        if out is not None:
            out = out.view(Quantity)
            out._unit = result_unit
        value = np.ndarray.dot(self, b, out=out)
        return Quantity(value, result_unit)

    def diff(self, n=1, axis=-1):
        value = np.diff(self.value, n=n, axis=axis)
        return Quantity(value, self.unit)

    def ediff1d(self, to_end=None, to_begin=None):
        value = np.ediff1d(self.value, to_end=to_end, to_begin=to_begin)
        return Quantity(value, self.unit)

    def nansum(self, axis=None):
        value = np.nansum(self.value, axis=axis)
        return Quantity(value, self.unit)

    def sum(self, axis=None, dtype=None, out=None, keepdims=False):
        out = out and out.view(Quantity)
        try:
            value = np.sum(self.value, axis=axis, dtype=dtype,
                           out=out, keepdims=keepdims)
        except:  # numpy < 1.7
            value = np.sum(self.value, axis=axis, dtype=dtype,
                           out=out)
        return Quantity(value, self.unit)

    def cumsum(self, axis=None, dtype=None, out=None):
        out = out and out.view(Quantity)
        value = np.cumsum(self.value, axis=axis, dtype=dtype, out=out)
        return Quantity(value, self.unit)

    def prod(self, axis=None, dtype=None, out=None, keepdims=False):
        if _is_unity(self.unit):
            out = out and out.view(Quantity)
            try:
                value = np.prod(self.value, axis=axis, dtype=dtype,
                                out=out, keepdims=keepdims)
            except:  # numpy < 1.7
                value = np.prod(self.value, axis=axis, dtype=dtype,
                                out=out)
            return Quantity(value, self.unit)
        else:
            raise ValueError("cannot use prod on scaled or "
                             "non-dimensionless Quantity arrays")

    def cumprod(self, axis=None, dtype=None, out=None):
        if _is_unity(self.unit):
            out = out and out.view(Quantity)
            value = np.cumprod(self.value, axis=axis, dtype=dtype, out=out)
            return Quantity(value, self.unit)
        else:
            raise ValueError("cannot use cumprod on scaled or "
                             "non-dimensionless Quantity arrays")

    def all(self, axis=None, out=None):
        raise NotImplementedError("cannot evaluate truth value of quantities. "
                                  "Evaluate array with q.value.all(...)")

    def any(self, axis=None, out=None):
        raise NotImplementedError("cannot evaluate truth value of quantities. "
                                  "Evaluate array with {0}.value.any(...)")
