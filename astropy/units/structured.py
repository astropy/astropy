# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines structured units and quantities.
"""

# Standard library
import operator

import numpy as np

from .core import Unit, UnitConversionError
from .quantity import Quantity


__all__ = ['StructuredUnit', 'StructuredQuantity']


def _names_from_dtype(dtype):
    names = []
    for field, (subdtype, offset) in dtype.fields.items():
        if subdtype.fields:
            names.append([field, _names_from_dtype(subdtype)])
        else:
            names.append(field)
    return tuple(names)


class StructuredUnit(np.void):
    """Container for units for a Structured Quantity.

    Parameters
    ----------
    units : nested tuple or `~numpy.void`
        Containing items that can initialize regular `~astropy.units.Unit`.
    names : nested tuple of names or `~numpy.dtype`, optional
        For nested tuples, by default the name of the upper entry will just
        be the concatenation of the names of the lower levels.  One can
        pass in a list with the upper-level name and a tuple of lower-level
        names to avoid this.
    """
    def __new__(cls, units=(), names=None):
        if isinstance(units, str):
            units = Unit(units)

        if isinstance(units, StructuredUnit):
            if names is None:
                return units
            else:
                units = units.item()

        if not isinstance(units, tuple):
            units = (units,)
        elif units == ():
            return

        if names is None:
            names = tuple('f{}'.format(i) for i in range(len(units)))
        elif isinstance(names, str):
            names = (names,)
        elif isinstance(names, np.dtype):
            names = _names_from_dtype(names)

        if len(units) != len(names):
            raise ValueError("lengths of units and names must match.")

        dtype = []
        converted = []
        for unit, name in zip(units, names):
            if isinstance(name, list):
                assert len(name) == 2
                unit = cls(unit, name[1])
                name = name[0]
            elif isinstance(name, tuple):
                unit = cls(unit, name)
                name = ''.join(name)
            elif isinstance(unit, str):
                unit = Unit(unit)
            elif isinstance(unit, tuple):
                unit = cls(unit)

            converted.append(unit)
            dtype.append((name, 'O'))

        dtype = np.dtype((cls, dtype))
        self = np.array(tuple(converted), dtype)[()]
        return self

    @property
    def _quantity_class(self):
        return StructuredQuantity

    def parts(self):
        for name in self.dtype.names:
            yield self[name]

    def _recursively_apply(self, func, cls=None):
        if cls is None:
            dtype = self.dtype
        else:
            dtype = np.dtype((cls, self.dtype))

        results = tuple([func(part) for part in self.parts()])

        return np.array(results, dtype)[()]

    def _recursively_get_dtype(self, value):
        if not isinstance(value, tuple) or len(self.dtype.names) != len(value):
            raise ValueError("cannot interpret value for unit {}"
                             .format(self))
        new_dtype = []
        for name, value_part in zip(self.dtype.names, value):
            part = self[name]
            if isinstance(part, type(self)):
                new_dtype.append((name, part._recursively_get_dtype(value_part)))
            else:
                value_part = np.array(value_part)
                new_dtype.append((name, 'f8', value_part.shape))
        return np.dtype(new_dtype)

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

        is_equivalent = self.dtype == other.dtype
        for self_part, other_part in zip(self.parts(), other.parts()):
            is_equivalent = is_equivalent and self_part.is_equivalent(
                other_part, equivalencies=equivalencies)
        return is_equivalent

    def _get_converter(self, other, equivalencies=[]):
        if not isinstance(other, type(self)):
            other = self.__class__(other, self.dtype)

        converters = [self_part._get_converter(other_part,
                                               equivalencies=equivalencies)
                      for (self_part, other_part) in zip(self.parts(),
                                                         other.parts())]

        def converter(value):
            if not hasattr(value, 'dtype'):
                # Get inner tuple.
                tmp = value
                while isinstance(tmp, list):
                    tmp = tmp[0]
                # Get correct dtype using first element.
                dtype = self._recursively_get_dtype(tmp)
                # Use that for second conversion.
                value = np.array(value, dtype)
            result = np.empty(value.shape, value.dtype)
            for name, converter_ in zip(result.dtype.names, converters):
                result[name] = converter_(value[name])
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

    def to_string(self, format='generic'):
        """Output the unit in the given format as a string.

        Units are separated by commas.

        Parameters
        ----------
        format : `astropy.units.format.Base` instance or str
            The name of a format or a formatter object.  If not
            provided, defaults to the generic format.
        """
        if format not in ('generic', 'unscaled', 'latex'):
            raise ValueError("Structured units cannot be written in {0} format. "
                             "Only 'generic', 'unscaled' and 'latex' are "
                             "supported.".format(format))
        items = self.item()
        out = '({})'.format(', '.join([item.to_string(format)
                                       for item in items]))
        return out if len(items) > 1 else out[:-1] + ',)'

    def _repr_latex_(self):
        return self.to_string('latex')

    def __str__(self):
        return self.to_string('generic')

    def __repr__(self):
        return 'Unit("{}")'.format(self.to_string())


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
                dtype = unit.dtype('f8')
            except AttributeError:
                pass
        value = np.array(value, dtype, *args, **kwargs)
        if not value.dtype.fields:
            raise TypeError("Values should be structured.")
        if not isinstance(value, cls):
            value = value.view(cls)
        value._set_unit(unit)
        return value

    def parts(self):
        for name in self.dtype.names:
            yield self[name]

    def __quantity_subclass__(self, unit):
        if isinstance(unit, StructuredUnit):
            return type(self), True
        else:
            return super().__quantity_subclass__(unit)[0], False

    def __getitem__(self, item):
        out = super().__getitem__(item)
        if out.dtype is not self.dtype:
            # item caused dtype change -> indexed with string-like.
            out_unit = self.unit[item]
            if not out.dtype.fields:
                out = self._new_view(out, out_unit)
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
