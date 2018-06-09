# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines structured units and quantities.
"""

# Standard library
import operator

import numpy as np

from .core import Unit, UnitBase
from .quantity import Quantity


__all__ = ['StructuredUnit', 'StructuredQuantity']


def _replace_dtype_fields_recursive(dtype, primitive_dtype):
    "Private function allowing recursion in replace_dtype_fields."
    # Copied from numpy.ma.core, except that we don't do subarrays.
    _recurse = _replace_dtype_fields_recursive

    # Do we have some name fields ?
    if dtype.names:
        descr = []
        for name in dtype.names:
            field = dtype.fields[name]
            if len(field) == 3:
                # Prepend the title to the name
                name = (field[-1], name)
            descr.append((name, _recurse(field[0], primitive_dtype)))
        new_dtype = np.dtype(descr)

    # this is a sub-array or primitive type, so do a direct replacement
    else:
        new_dtype = primitive_dtype

    # preserve identity of dtypes
    if new_dtype == dtype:
        new_dtype = dtype

    return new_dtype


def _replace_dtype_fields(dtype, primitive_dtype):
    """Construct a dtype description list from a given dtype.

    Returns a new dtype object, with all fields and subtypes in the
    given type recursively replaced with ``primitive_dtype``, but with
    sub-arrays removed.

    Arguments are coerced to dtypes first.
    """
    # Copied from numpy.ma.core.
    dtype = np.dtype(dtype)
    primitive_dtype = np.dtype(primitive_dtype)
    return _replace_dtype_fields_recursive(dtype, primitive_dtype)


class StructuredUnit(np.void):
    """Container for units for a Structured Quantity.

    Parameters
    ----------
    units : nested tuple or `~numpy.void`
        Containing items that can initialize regular `~astropy.units.Unit`.
    dtype : `~numpy.dtype`, optional
        A dtype of the `~astropy.units.StructuredQuantity` the units
        are associated with.  Each element will be replaced with a unit.
        Optional only if ``units`` is alread a structured `~numpy.void`.
    """

    def __new__(cls, units, dtype=None):
        if dtype is None:
            dtype = getattr(units, 'dtype', None)
            if dtype is None:
                # Assume single field is wanted.
                dtype = [('f0', 'O')]

        if not isinstance(dtype, cls):
            dtype = _replace_dtype_fields(dtype, 'O')
            dtype = np.dtype((cls, dtype))
        self = np.array(units, dtype)[()]
        self._recursively_check()
        return self

    def _recursively_check(self):
        for field in self.dtype.names:
            # For nested items, __getitem__ ensures this is a StructuredUnit.
            unit = self[field]
            if isinstance(unit, type(self)):
                unit._recursively_check()
            else:
                self[field] = Unit(unit)

    @property
    def _quantity_class(self):
        return StructuredQuantity

    def __getitem__(self, item):
        val = super().__getitem__(item)
        if isinstance(val, np.void) and val.dtype.fields:
            return val.view((self.__class__, val.dtype.fields))
        else:
            return val

    def _recursively_apply(self, func, cls=None):
        r = self.copy()
        if cls is not None:
            r = r.view((cls, self.dtype))
        for field in self.dtype.fields:
            r[field] = func(self[field])

        return r

    @property
    def si(self):
        """Returns a copy of the current `Unit` instance in SI units."""
        return self._recursively_apply(operator.attrgetter('si'))

    @property
    def cgs(self):
        """Returns a copy of the current `Unit` instance in SI units."""
        return self._recursively_apply(operator.attrgetter('cgs'))

    # Needed to pass through Unit initializer, so might as well use it.
    def _get_physical_type_id(self):
        return self._recursively_apply(
            operator.methodcaller('_get_physical_type_id'), np.void)

    @property
    def physical_type(self):
        """Physical types of all the fields."""
        return self._recursively_apply(
            operator.attrgetter('physical_type'), np.void)

    def decompose(self, bases=set()):
        """Return a unit object composed of only irreducible units.

        Parameters
        ----------
        bases : sequence of UnitBase, optional
            The bases to decompose into.  When not provided,
            decomposes down to any irreducible units.  When provided,
            the decomposed result will only contain the given units.
            This will raises a `UnitsError` if it's not possible
            to do so.

        Returns
        -------
        unit : New `~astropy.units.StructuredUnit` instance
            With the unit for each field containing only irreducible units.
        """
        return self._recursively_apply(
            operator.methodcaller('decompose', bases=bases))

    def is_equivalent(self, other, equivalencies=[]):
        """`True` if all fields are equivalent to the other's fields.

        Parameters
        ----------
        other : `~astropy.units.StructuredUnit` or (nested) tuple
            The unit to compare with.
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            The list will be applied to all fields.

        Returns
        -------
        bool
        """
        if not isinstance(other, type(self)):
            try:
                other = self.__class__(other, self.dtype)
            except Exception:
                return False

        is_equivalent = self.dtype.fields == other.dtype.fields
        for field in self.dtype.fields:
            is_equivalent = is_equivalent and self[field].is_equivalent(
                other[field], equivalencies=equivalencies)
        return is_equivalent

    def _get_converter(self, unit, equivalencies=[]):
        if not isinstance(unit, type(self)):
            unit = self.__class__(unit, self.dtype)

        def converter(value):
            if not hasattr(value, 'dtype'):
                dtype = _replace_dtype_fields(self.dtype, 'f8')
                value = np.array(value, dtype)
            result = np.empty(value.shape, value.dtype)
            for name, r_name, u_name in zip(self.dtype.names,
                                            result.dtype.names,
                                            unit.dtype.names):
                result[r_name] = self[name].to(unit[u_name], value[r_name],
                                               equivalencies=equivalencies)
            return result
        return converter

    def to(self, other, value, equivalencies=[]):
        """Return values converted to the specified unit.

        Parameters
        ----------
        other : `~astropy.units.StructuredUnit` or (nested) tuple of str
            The unit to convert to.  If necessary, will be converted to
            a Structured Unit using the dtype of ``value``.
        value : scalar void or array, or sequence convertible to array.
            Value(s) in the current unit to be converted to the
            specified unit.
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            This list is in addition to possible global defaults set by, e.g.,
            `set_enabled_equivalencies`.
            Use `None` to turn off all equivalencies.

        Returns
        -------
        values : scalar or array
            Converted value(s).

        Raises
        ------
        UnitsError
            If units are inconsistent
        """
        return self._get_converter(other, equivalencies=equivalencies)(value)


class StructuredQuantity(Quantity):
    """Structured array with associated units.

    Parameters
    ----------
    value : structured `~numpy.ndarray` or (nested) tuple
        Numerical values of the various elements of the structure.
    unit : `~astropy.units.StructuredUnit` or (nested) tuple
        Associated units.
    dtype : `~numpy.dtype`, optional
        If not given, inferred from ``value``
    *args, **kwargs
        All other arguments are passed on to the `~numpy.array` constructor.
    """
    def __new__(cls, value, unit, dtype=None, *args, **kwargs):
        if dtype is None and not hasattr(value, 'dtype'):
            try:
                dtype = _replace_dtype_fields(unit.dtype, 'f8')
            except AttributeError:
                pass
        value = np.array(value, dtype, *args, **kwargs)
        if not value.dtype.fields:
            raise TypeError("Values should be structured.")
        if not isinstance(value, cls):
            value = value.view(cls)
        value._set_unit(unit)
        return value

    def __getitem__(self, item):
        out = super().__getitem__(item)
        if out.dtype is not self.dtype:
            # item caused dtype change -> indexed with string-like.
            out_unit = self.unit[item]
            if not out.dtype.fields:
                out = out.view(getattr(out_unit, '_quantity_cls',
                                       Quantity))
            out._unit = out_unit

        return out

    def _set_unit(self, unit):
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit, self.dtype)
        self._unit = unit

    def to(self, unit, equivalencies=[]):
        """Convert to the specific structured unit.

        Parameters
        ----------
        unit : `~astropy.units.StructuredUnit` or (nested) tuple of units
            Will be converted if necessary.
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            These will apply to all elements in the structured quantity.
        """
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit, self.dtype)
        return self._new_view(self._to_value(unit, equivalencies), unit)

    def to_value(self, unit=None, equivalencies=[]):
        """The numerical value, possibly in a different unit.

        Parameters
        ----------
        unit : `~astropy.units.StructuredUnit` or (nested) tuple of units
            Will be converted if necessary.
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            These will apply to all elements in the structured quantity.
        """
        if unit is None or unit is self.unit:
            result = self.view(np.ndarray)
        else:
            result = self._to_value(unit, equivalencies)
        # [()] to return a void rather than a tuple.
        return result if result.shape else result[()]

    value = property(to_value,
                     doc="""The numerical value of this instance.

    See also
    --------
    to_value : Get the numerical value in a given unit.
    """)
