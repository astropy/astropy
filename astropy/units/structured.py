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


def _names_from_dtype(dtype):
    names = []
    for field, (subdtype, offset) in dtype.fields.items():
        if subdtype.fields:
            names.append([field, _names_from_dtype(subdtype)])
        else:
            names.append(field)
    return tuple(names)


def _normalize_tuples(names):
    return tuple(
        ([''.join(name), _normalize_tuples(name)]
         if isinstance(name, tuple) else name) for name in names)


class StructuredUnit:
    """Container for units for a Structured Quantity.

    Parameters
    ----------
    units : nested tuple of units, or str
        Tuple elements should contain items that can initialize regular units.
    names : nested tuple of str or `~numpy.dtype`, optional
        For nested tuples, by default the name of the upper entry will just
        be the concatenation of the names of the lower levels.  One can
        pass in a list with the upper-level name and a tuple of lower-level
        names to avoid this.  Not all levels have to be given, except when
        a `~numpy.dtype` is passed in.
    """

    def __new__(cls, units, names=None):
        allow_default_names = True
        if names is not None:
            if isinstance(names, np.dtype):
                if not names.names:
                    raise ValueError('dtype should be structured, with names.')
                names = _names_from_dtype(names)
                allow_default_names = False
            else:
                if not isinstance(names, tuple):
                    names = (names,)
                names = _normalize_tuples(names)

        if not isinstance(units, tuple):
            units = Unit(units)
            if isinstance(units, StructuredUnit):
                # Avoid constructing a new StructuredUnit if no names
                # are given, or if all names are the same already anyway.
                if names is None or units.names == names:
                    return units

                # Otherwise, turn (the upper level) into a tuple.
                units = units.values()
            else:
                # Single regular unit: make a tuple for iteration below.
                units = (units,)

        if names is None:
            names = tuple('f{}'.format(i) for i in range(len(units)))

        elif len(units) != len(names):
            raise ValueError("lengths of units and names must match.")

        dtype = []
        converted = []
        for unit, name in zip(units, names):
            if isinstance(name, list):
                # For list, the first item is the name of our level,
                # and the second another tuple, i.e., we recurse.
                assert len(name) == 2
                unit = cls(unit, name[1])
                name = name[0]
            else:
                # We are at the lowest level.  Check unit.
                if isinstance(unit, str):
                    unit = Unit(unit)
                elif isinstance(unit, tuple):
                    unit = cls(unit)

                if not allow_default_names and isinstance(unit,
                                                          StructuredUnit):
                    raise ValueError('units and names from dtype do not '
                                     'match in depth.')

            converted.append(unit)
            dtype.append((name, 'O'))

        self = super().__new__(cls)
        dtype = np.dtype(dtype)
        self._units = np.array(tuple(converted), dtype)
        return self

    @property
    def _quantity_class(self):
        return StructuredQuantity

    @property
    def names(self):
        """Possibly nested tuple of the names of the parts."""
        return tuple(([name, part.names]
                      if isinstance(part, StructuredUnit) else name)
                     for name, part in self.items())

    # Allow StructuredUnit to be treated as an (ordered) mapping.
    def __len__(self):
        return len(self._units.dtype.names)

    def __getitem__(self, item):
        return self._units[item].item()

    def values(self):
        return self._units.item()

    def keys(self):
        return self._units.dtype.names

    def items(self):
        return zip(self._units.dtype.names, self._units.item())

    def __iter__(self):
        yield from self._units.dtype.names

    # Helpers for methods below.
    def _recursively_apply(self, func, as_void=False):
        """Apply func recursively, storing the results in an instance of cls"""
        results = np.array(tuple([func(part) for part in self.values()]),
                           self._units.dtype)
        if as_void:
            return results[()]

        # Short-cut; no need to interpret names, etc.
        result = super().__new__(self.__class__)
        result._units = results
        return result

    def _recursively_get_dtype(self, value):
        """Get structured dtype according to value, using our names."""
        # Helper for _create_array below.
        if not isinstance(value, tuple) or len(self) != len(value):
            raise ValueError("cannot interpret value for unit {}"
                             .format(self))
        new_dtype = []
        for (name, unit), part in zip(self.items(), value):
            if isinstance(unit, type(self)):
                new_dtype.append((name, unit._recursively_get_dtype(part)))
            else:
                value_part = np.array(part)
                new_dtype.append((name, value_part.dtype, value_part.shape))
        return np.dtype(new_dtype)

    def _create_array(self, value):
        """Turn a (nested) list of properly nested tuples into an array.

        This routine is needed since ``np.array(value)`` would treat the
        tuples as lower levels of the array, rather than as elements
        of a structured array.  The routine does presume that the type
        of the first tuple is representative of the rest.
        """
        # Get inner tuple.
        tmp = value
        while isinstance(tmp, list):
            tmp = tmp[0]
        # Get correct dtype using first element.
        dtype = self._recursively_get_dtype(tmp)
        # Use that for second conversion.
        return np.array(value, dtype)

    @property
    def si(self):
        """The `StructuredUnit` instance in SI units."""
        return self._recursively_apply(operator.attrgetter('si'))

    @property
    def cgs(self):
        """The `StructuredUnit` instance in cgs units."""
        return self._recursively_apply(operator.attrgetter('cgs'))

    # Needed to pass through Unit initializer, so might as well use it.
    def _get_physical_type_id(self):
        return self._recursively_apply(
            operator.methodcaller('_get_physical_type_id'), as_void=True)

    @property
    def physical_type(self):
        """Physical types of all the fields."""
        return self._recursively_apply(
            operator.attrgetter('physical_type'), as_void=True)

    def decompose(self, bases=set()):
        """Return a StructuredUnit object composed of only irreducible units.

        Parameters
        ----------
        bases : sequence of `~astropy.units.UnitBase`, optional
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
        other : `~astropy.units.StructuredUnit`, or what can initialize one
            The structured unit to compare with.
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
                other = self.__class__(other, self.names)
            except Exception:
                return False

        if len(self) != len(other):
            return False

        for self_part, other_part in zip(self.values(), other.values()):
            if not self_part.is_equivalent(other_part,
                                           equivalencies=equivalencies):
                return False

        return True

    def _get_converter(self, other, equivalencies=[]):
        if not isinstance(other, type(self)):
            other = self.__class__(other, self.names)

        converters = [self_part._get_converter(other_part,
                                               equivalencies=equivalencies)
                      for (self_part, other_part) in zip(self.values(),
                                                         other.values())]

        def converter(value):
            if not hasattr(value, 'dtype'):
                value = self._create_array(value)
            result = np.empty(value.shape, value.dtype)
            for name, converter_ in zip(result.dtype.names, converters):
                result[name] = converter_(value[name])
            return result

        return converter

    def to(self, other, value, equivalencies=[]):
        """Return values converted to the specified unit.

        Parameters
        ----------
        other : `~astropy.units.StructuredUnit`, or what can initialize one
            The unit to convert to.  If necessary, will be converted to
            a Structured Unit using the dtype of ``value``.
        value : scalar void or array, or sequence convertible to array.
            Value(s) in the current unit to be converted to the
            specified unit.  If a sequence, the first element must have
            entries of the correct type to represent all elements (i.e.,
            not have, e.g., an ``int`` where other elements have ``float``).
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
            raise ValueError("Structured units cannot be written in {0} "
                             "format. Only 'generic', 'unscaled' and "
                             "'latex' are supported.".format(format))
        out = '({})'.format(', '.join([part.to_string(format)
                                       for part in self.values()]))
        return out if len(self) > 1 else out[:-1] + ',)'

    def _repr_latex_(self):
        return self.to_string('latex')

    __array_ufunc__ = None

    def __mul__(self, other):
        if isinstance(other, str):
            try:
                other = Unit(other, parse_strict='silent')
            except Exception:
                return NotImplemented
        if isinstance(other, UnitBase):
            new_units = tuple(part * other for part in self.values())
            return self.__class__(new_units, self.names)
        if isinstance(other, StructuredUnit):
            return NotImplemented

        # Anything not like a unit, try initialising as a structured quantity.
        try:
            return self._quantity_class(other, unit=self)
        except Exception:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, str):
            try:
                other = Unit(other, parse_strict='silent')
            except Exception:
                return NotImplemented

        if isinstance(other, UnitBase):
            new_units = tuple(part / other for part in self.values())
            return self.__class__(new_units, self.names)
        return NotImplemented

    def __rlshift__(self, m):
        try:
            return StructuredQuantity(m, self, copy=False, subok=True)
        except Exception:
            return NotImplemented

    def __str__(self):
        return self.to_string('generic')

    def __repr__(self):
        return 'Unit("{}")'.format(self.to_string())

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            try:
                other = self.__class__(other, self.names)
            except Exception:
                return NotImplemented

        return self.values() == other.values()

    def __ne__(self, other):
        eq = self.__eq__(other)
        return eq if eq is NotImplemented else not eq


class StructuredQuantity(Quantity):
    """Structured array with associated units.

    Parameters
    ----------
    value : structured `~numpy.ndarray` or (nested) tuple
        Numerical values of the various elements of the structure.
    unit : `~astropy.units.StructuredUnit` or (nested) tuple of units or str
        Associated units.
    dtype : `~numpy.dtype`, optional
        If not given, inferred from ``value``
    *args, **kwargs
        All other arguments are passed on to the `~numpy.array` constructor.
    """
    def __new__(cls, value, unit, dtype=None, *args, **kwargs):
        # TODO: special-case StructuredQuantity input.
        # Tuples of tuples get interpreted as an array rather than
        # as a structured void, so avoid that.
        if (dtype is None and not hasattr(value, 'dtype') and
                isinstance(unit, StructuredUnit)):
            value = unit._create_array(value)
        else:
            value = np.array(value, dtype, *args, **kwargs)
        if not value.dtype.fields:
            raise TypeError("Values should be structured.")
        if not isinstance(value, cls):
            value = value.view(cls)
        value._set_unit(unit)
        return value

    def __quantity_subclass__(self, unit):
        if isinstance(unit, StructuredUnit):
            return type(self), True
        else:
            return super().__quantity_subclass__(unit)[0], False

    def __getitem__(self, item):
        out_value = self.value[item]
        if out_value.dtype is self.dtype:
            out_unit = self.unit
        else:
            # item caused dtype change -> indexed with string-like.
            # Index unit as well.
            out_unit = self.unit[item]

        return self._new_view(out_value, out_unit)

    def __setitem__(self, item, value):
        out_item = self[item]
        if out_item.dtype is self.dtype:
            super(StructuredQuantity, out_item).__setitem__(Ellipsis, value)
        else:
            # item caused dtype change -> indexed with string-like, so we
            # have a different unit (and may well be a Quantity); try again.
            out_item[...] = value

    def _set_unit(self, unit):
        unit = StructuredUnit(unit, self.dtype)
        self._unit = unit

    def to(self, unit, equivalencies=[]):
        """Convert to the specific structured unit.

        Parameters
        ----------
        unit : `~astropy.units.StructuredUnit`, or what can initialize one
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

    def _to_own_unit(self, other, check_precision=False):
        other_value = super()._to_own_unit(other, check_precision)
        # Setting names to ensure things like equality work (note that
        # above will have failed already if units did not match).
        other_value.dtype.names = self.dtype.names
        return other_value

    def _recursively_apply(self, func):
        """Apply function recursively to every field, a copy with the result."""
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
        return self._recursively_apply(operator.attrgetter('si'))

    @property
    def cgs(self):
        return self._recursively_apply(operator.attrgetter('cgs'))
