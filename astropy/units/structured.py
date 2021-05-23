# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines structured units and quantities.
"""

# Standard library
import operator

import numpy as np

from .core import Unit, UnitBase


__all__ = ['StructuredUnit']


def _names_from_dtype(dtype):
    """Recursively extract names from a dtype."""
    names = []
    for field, (subdtype, offset) in dtype.fields.items():
        if subdtype.fields:
            names.append([field, _names_from_dtype(subdtype)])
        else:
            names.append(field)
    return tuple(names)


def _normalize_tuples(names):
    """Infer upper level names from unadorned tuples.

    Generally, we want the names to be organized like dtypes, as in
    ``(['pv', ('p', 'v')], 't')``.  But we automatically infer upper
    names if the list is absent from items like ``(('p', 'v'), 't')``,
    by concatenating the names inside the tuple.
    """
    result = []
    for name in names:
        if isinstance(name, str) and len(name) > 0:
            result.append(name)
        elif (isinstance(name, list)
              and len(name) == 2
              and isinstance(name[0], str) and len(name[0]) > 0
              and isinstance(name[1], tuple) and len(name[1]) > 0):
            result.append([name[0], _normalize_tuples(name[1])])
        elif isinstance(name, tuple) and len(name) > 0:
            new_tuple = _normalize_tuples(name)
            result.append([''.join([(i[0] if isinstance(i, list) else i)
                                    for i in new_tuple]), new_tuple])
        else:
            raise ValueError(f'invalid entry {name!r}. Should be a name, '
                             'tuple of names, or 2-element list of the '
                             'form [name, tuple of names].')

    return tuple(result)


class StructuredUnit:
    """Container for units for a Structured Quantity.

    Parameters
    ----------
    units : tuple of unit-like
        Tuple elements should contain items that can initialize regular units.
        They can be nested tuples.
    names : tuple of str, or `~numpy.dtype`, optional
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
            names = tuple(f'f{i}' for i in range(len(units)))

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
                unit = Unit(unit)
                if not allow_default_names and isinstance(unit,
                                                          StructuredUnit):
                    raise ValueError('units and names from dtype do not '
                                     'match in depth.')

            converted.append(unit)
            dtype.append((name, 'O'))

        self = super().__new__(cls)
        dtype = np.dtype(dtype)
        # Decay array to void so we can access by field name and number.
        self._units = np.array(tuple(converted), dtype)[()]
        return self

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
        return self._units[item]

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
        """Apply func recursively.

        The result is stored in an instance of cls or a void.
        """
        results = np.array(tuple([func(part) for part in self.values()]),
                           self._units.dtype)[()]
        if as_void:
            return results

        # Short-cut; no need to interpret names, etc.
        result = super().__new__(self.__class__)
        result._units = results
        return result

    def _recursively_get_dtype(self, value, enter_lists=True):
        """Get structured dtype according to value, using our names.

        This is useful since ``np.array(value)`` would treat tuples as lower
        levels of the array, rather than as elements of a structured array.
        The routine does presume that the type of the first tuple is
        representative of the rest.

        """
        # Used in _get_converter below.
        if enter_lists:
            while isinstance(value, list):
                value = value[0]
        if not isinstance(value, tuple) or len(self) != len(value):
            raise ValueError("cannot interpret value for unit {}"
                             .format(self))
        new_dtype = []
        for (name, unit), part in zip(self.items(), value):
            if isinstance(unit, type(self)):
                new_dtype.append(
                    (name, unit._recursively_get_dtype(part, enter_lists=False)))
            else:
                part = np.array(part)
                part_dtype = part.dtype
                if part_dtype.kind in 'iu':
                    part_dtype = np.dtype(float)
                new_dtype.append((name, part_dtype, part.shape))
        return np.dtype(new_dtype)

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
        `~astropy.units.StructuredUnit`
            With the unit for each field containing only irreducible units.
        """
        return self._recursively_apply(
            operator.methodcaller('decompose', bases=bases))

    def is_equivalent(self, other, equivalencies=[]):
        """`True` if all fields are equivalent to the other's fields.

        Parameters
        ----------
        other : `~astropy.units.StructuredUnit`
            The structured unit to compare with, or what can initialize one.
        equivalencies : list of tuple, optional
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
                value = np.array(value, self._recursively_get_dtype(value))
            result = np.empty_like(value)
            for name, converter_ in zip(result.dtype.names, converters):
                result[name] = converter_(value[name])
            return result

        return converter

    def to(self, other, value, equivalencies=[]):
        """Return values converted to the specified unit.

        Parameters
        ----------
        other : `~astropy.units.StructuredUnit`
            The unit to convert to.  If necessary, will be converted to
            a Structured Unit using the dtype of ``value``.
        value : array-like
            Value(s) in the current unit to be converted to the
            specified unit.  If a sequence, the first element must have
            entries of the correct type to represent all elements (i.e.,
            not have, e.g., a ``float`` where other elements have ``complex``).
        equivalencies : list of tuple, optional
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
            from .quantity import Quantity
            return Quantity(other, unit=self)
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
            from .quantity import Quantity
            return Quantity(m, self, copy=False, subok=True)
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
