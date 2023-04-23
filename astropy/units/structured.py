# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines structured units and quantities.
"""

from __future__ import annotations  # For python < 3.10

# Standard library
import operator

import numpy as np

from .core import UNITY, Unit, UnitBase

__all__ = ["StructuredUnit"]


DTYPE_OBJECT = np.dtype("O")


def _names_from_dtype(dtype):
    """Recursively extract field names from a dtype."""
    names = []
    for name in dtype.names:
        subdtype = dtype.fields[name][0].base
        if subdtype.names:
            names.append([name, _names_from_dtype(subdtype)])
        else:
            names.append(name)
    return tuple(names)


def _normalize_names(names):
    """Recursively normalize, inferring upper level names for unadorned tuples.

    Generally, we want the field names to be organized like dtypes, as in
    ``(['pv', ('p', 'v')], 't')``.  But we automatically infer upper
    field names if the list is absent from items like ``(('p', 'v'), 't')``,
    by concatenating the names inside the tuple.
    """
    result = []
    for name in names:
        if isinstance(name, str) and len(name) > 0:
            result.append(name)
        elif (
            isinstance(name, list)
            and len(name) == 2
            and isinstance(name[0], str)
            and len(name[0]) > 0
            and isinstance(name[1], tuple)
            and len(name[1]) > 0
        ):
            result.append([name[0], _normalize_names(name[1])])
        elif isinstance(name, tuple) and len(name) > 0:
            new_tuple = _normalize_names(name)
            name = "".join([(i[0] if isinstance(i, list) else i) for i in new_tuple])
            result.append([name, new_tuple])
        else:
            raise ValueError(
                f"invalid entry {name!r}. Should be a name, "
                "tuple of names, or 2-element list of the "
                "form [name, tuple of names]."
            )

    return tuple(result)


class StructuredUnit:
    """Container for units for a structured Quantity.

    Parameters
    ----------
    units : unit-like, tuple of unit-like, or `~astropy.units.StructuredUnit`
        Tuples can be nested.  If a `~astropy.units.StructuredUnit` is passed
        in, it will be returned unchanged unless different names are requested.
    names : tuple of str, tuple or list; `~numpy.dtype`; or `~astropy.units.StructuredUnit`, optional
        Field names for the units, possibly nested. Can be inferred from a
        structured `~numpy.dtype` or another `~astropy.units.StructuredUnit`.
        For nested tuples, by default the name of the upper entry will be the
        concatenation of the names of the lower levels.  One can pass in a
        list with the upper-level name and a tuple of lower-level names to
        avoid this.  For tuples, not all levels have to be given; for any level
        not passed in, default field names of 'f0', 'f1', etc., will be used.

    Notes
    -----
    It is recommended to initialize the class indirectly, using
    `~astropy.units.Unit`.  E.g., ``u.Unit('AU,AU/day')``.

    When combined with a structured array to produce a structured
    `~astropy.units.Quantity`, array field names will take precedence.
    Generally, passing in ``names`` is needed only if the unit is used
    unattached to a `~astropy.units.Quantity` and one needs to access its
    fields.

    Examples
    --------
    Various ways to initialize a `~astropy.units.StructuredUnit`::

        >>> import astropy.units as u
        >>> su = u.Unit('(AU,AU/day),yr')
        >>> su
        Unit("((AU, AU / d), yr)")
        >>> su.field_names
        (['f0', ('f0', 'f1')], 'f1')
        >>> su['f1']
        Unit("yr")
        >>> su2 = u.StructuredUnit(((u.AU, u.AU/u.day), u.yr), names=(('p', 'v'), 't'))
        >>> su2 == su
        True
        >>> su2.field_names
        (['pv', ('p', 'v')], 't')
        >>> su3 = u.StructuredUnit((su2['pv'], u.day), names=(['p_v', ('p', 'v')], 't'))
        >>> su3.field_names
        (['p_v', ('p', 'v')], 't')
        >>> su3.keys()
        ('p_v', 't')
        >>> su3.values()
        (Unit("(AU, AU / d)"), Unit("d"))

    Structured units share most methods with regular units::

        >>> su.physical_type
        ((PhysicalType('length'), PhysicalType({'speed', 'velocity'})), PhysicalType('time'))
        >>> su.si
        Unit("((1.49598e+11 m, 1.73146e+06 m / s), 3.15576e+07 s)")

    """

    def __new__(cls, units, names=None):
        dtype = None
        if names is not None:
            if isinstance(names, StructuredUnit):
                dtype = names._units.dtype
                names = names.field_names
            elif isinstance(names, np.dtype):
                if not names.fields:
                    raise ValueError("dtype should be structured, with fields.")
                dtype = np.dtype([(name, DTYPE_OBJECT) for name in names.names])
                names = _names_from_dtype(names)
            else:
                if not isinstance(names, tuple):
                    names = (names,)
                names = _normalize_names(names)

        if not isinstance(units, tuple):
            units = Unit(units)
            if isinstance(units, StructuredUnit):
                # Avoid constructing a new StructuredUnit if no field names
                # are given, or if all field names are the same already anyway.
                if names is None or units.field_names == names:
                    return units

                # Otherwise, turn (the upper level) into a tuple, for renaming.
                units = units.values()
            else:
                # Single regular unit: make a tuple for iteration below.
                units = (units,)

        if names is None:
            names = tuple(f"f{i}" for i in range(len(units)))

        elif len(units) != len(names):
            raise ValueError("lengths of units and field names must match.")

        converted = []
        for unit, name in zip(units, names):
            if isinstance(name, list):
                # For list, the first item is the name of our level,
                # and the second another tuple of names, i.e., we recurse.
                unit = cls(unit, name[1])
                name = name[0]
            else:
                # We are at the lowest level.  Check unit.
                unit = Unit(unit)
                if dtype is not None and isinstance(unit, StructuredUnit):
                    raise ValueError(
                        "units do not match in depth with field "
                        "names from dtype or structured unit."
                    )

            converted.append(unit)

        self = super().__new__(cls)
        if dtype is None:
            dtype = np.dtype(
                [
                    ((name[0] if isinstance(name, list) else name), DTYPE_OBJECT)
                    for name in names
                ]
            )
        # Decay array to void so we can access by field name and number.
        self._units = np.array(tuple(converted), dtype)[()]
        return self

    def __getnewargs__(self):
        """When de-serializing, e.g. pickle, start with a blank structure."""
        return (), None

    @property
    def field_names(self):
        """Possibly nested tuple of the field names of the parts."""
        return tuple(
            ([name, unit.field_names] if isinstance(unit, StructuredUnit) else name)
            for name, unit in self.items()
        )

    # Allow StructuredUnit to be treated as an (ordered) mapping.
    def __len__(self):
        return len(self._units.dtype.names)

    def __getitem__(self, item):
        # Since we are based on np.void, indexing by field number works too.
        return self._units[item]

    def values(self):
        return self._units.item()

    def keys(self):
        return self._units.dtype.names

    def items(self):
        return tuple(zip(self._units.dtype.names, self._units.item()))

    def __iter__(self):
        yield from self._units.dtype.names

    # Helpers for methods below.
    def _recursively_apply(self, func, cls=None):
        """Apply func recursively.

        Parameters
        ----------
        func : callable
            Function to apply to all parts of the structured unit,
            recursing as needed.
        cls : type, optional
            If given, should be a subclass of `~numpy.void`. By default,
            will return a new `~astropy.units.StructuredUnit` instance.
        """
        applied = tuple(func(part) for part in self.values())
        # Once not NUMPY_LT_1_23: results = np.void(applied, self._units.dtype).
        results = np.array(applied, self._units.dtype)[()]
        if cls is not None:
            return results.view((cls, results.dtype))

        # Short-cut; no need to interpret field names, etc.
        result = super().__new__(self.__class__)
        result._units = results
        return result

    def _recursively_get_dtype(self, value, enter_lists=True):
        """Get structured dtype according to value, using our field names.

        This is useful since ``np.array(value)`` would treat tuples as lower
        levels of the array, rather than as elements of a structured array.
        The routine does presume that the type of the first tuple is
        representative of the rest.  Used in ``_get_converter``.

        For the special value of ``UNITY``, all fields are assumed to be 1.0,
        and hence this will return an all-float dtype.

        """
        if enter_lists:
            while isinstance(value, list):
                value = value[0]
        if value is UNITY:
            value = (UNITY,) * len(self)
        elif not isinstance(value, tuple) or len(self) != len(value):
            raise ValueError(f"cannot interpret value {value} for unit {self}.")
        descr = []
        for (name, unit), part in zip(self.items(), value):
            if isinstance(unit, StructuredUnit):
                descr.append(
                    (name, unit._recursively_get_dtype(part, enter_lists=False))
                )
            else:
                # Got a part associated with a regular unit. Gets its dtype.
                # Like for Quantity, we cast integers to float.
                part = np.array(part)
                part_dtype = part.dtype
                if part_dtype.kind in "iu":
                    part_dtype = np.dtype(float)
                descr.append((name, part_dtype, part.shape))
        return np.dtype(descr)

    @property
    def si(self):
        """The `StructuredUnit` instance in SI units."""
        return self._recursively_apply(operator.attrgetter("si"))

    @property
    def cgs(self):
        """The `StructuredUnit` instance in cgs units."""
        return self._recursively_apply(operator.attrgetter("cgs"))

    # Needed to pass through Unit initializer, so might as well use it.
    def _get_physical_type_id(self):
        return self._recursively_apply(
            operator.methodcaller("_get_physical_type_id"), cls=Structure
        )

    @property
    def physical_type(self):
        """Physical types of all the fields."""
        return self._recursively_apply(
            operator.attrgetter("physical_type"), cls=Structure
        )

    def decompose(self, bases=set()):
        """The `StructuredUnit` composed of only irreducible units.

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
        return self._recursively_apply(operator.methodcaller("decompose", bases=bases))

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
        try:
            other = StructuredUnit(other)
        except Exception:
            return False

        if len(self) != len(other):
            return False

        for self_part, other_part in zip(self.values(), other.values()):
            if not self_part.is_equivalent(other_part, equivalencies=equivalencies):
                return False

        return True

    def _get_converter(self, other, equivalencies=[]):
        if not isinstance(other, type(self)):
            other = self.__class__(other, names=self)

        converters = [
            self_part._get_converter(other_part, equivalencies=equivalencies)
            for (self_part, other_part) in zip(self.values(), other.values())
        ]

        def converter(value):
            if not hasattr(value, "dtype"):
                value = np.array(value, self._recursively_get_dtype(value))
            result = np.empty_like(value)
            for name, converter_ in zip(result.dtype.names, converters):
                result[name] = converter_(value[name])
            # Index with empty tuple to decay array scalars to numpy void.
            return result if result.shape else result[()]

        return converter

    def to(self, other, value=np._NoValue, equivalencies=[]):
        """Return values converted to the specified unit.

        Parameters
        ----------
        other : `~astropy.units.StructuredUnit`
            The unit to convert to.  If necessary, will be converted to
            a `~astropy.units.StructuredUnit` using the dtype of ``value``.
        value : array-like, optional
            Value(s) in the current unit to be converted to the
            specified unit.  If a sequence, the first element must have
            entries of the correct type to represent all elements (i.e.,
            not have, e.g., a ``float`` where other elements have ``complex``).
            If not given, assumed to have 1. in all fields.
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
        if value is np._NoValue:
            # We do not have UNITY as a default, since then the docstring
            # would list 1.0 as default, yet one could not pass that in.
            value = UNITY
        return self._get_converter(other, equivalencies=equivalencies)(value)

    def to_string(self, format="generic"):
        """Output the unit in the given format as a string.

        Units are separated by commas.

        Parameters
        ----------
        format : `astropy.units.format.Base` instance or str
            The name of a format or a formatter object.  If not
            provided, defaults to the generic format.

        Notes
        -----
        Structured units can be written to all formats, but can be
        re-read only with 'generic'.

        """
        parts = [part.to_string(format) for part in self.values()]
        out_fmt = "({})" if len(self) > 1 else "({},)"
        if format.startswith("latex"):
            # Strip $ from parts and add them on the outside.
            parts = [part[1:-1] for part in parts]
            out_fmt = "$" + out_fmt + "$"
        return out_fmt.format(", ".join(parts))

    def _repr_latex_(self):
        return self.to_string("latex")

    __array_ufunc__ = None

    def __mul__(self, other):
        if isinstance(other, str):
            try:
                other = Unit(other, parse_strict="silent")
            except Exception:
                return NotImplemented
        if isinstance(other, UnitBase):
            new_units = tuple(part * other for part in self.values())
            return self.__class__(new_units, names=self)
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
                other = Unit(other, parse_strict="silent")
            except Exception:
                return NotImplemented

        if isinstance(other, UnitBase):
            new_units = tuple(part / other for part in self.values())
            return self.__class__(new_units, names=self)
        return NotImplemented

    def __rlshift__(self, m):
        try:
            from .quantity import Quantity

            return Quantity(m, self, copy=False, subok=True)
        except Exception:
            return NotImplemented

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return f'Unit("{self.to_string()}")'

    def __eq__(self, other):
        try:
            other = StructuredUnit(other)
        except Exception:
            return NotImplemented

        return self.values() == other.values()

    def __ne__(self, other):
        if not isinstance(other, type(self)):
            try:
                other = StructuredUnit(other)
            except Exception:
                return NotImplemented

        return self.values() != other.values()


class Structure(np.void):
    """Single element structure for physical type IDs, etc.

    Behaves like a `~numpy.void` and thus mostly like a tuple which can also
    be indexed with field names, but overrides ``__eq__`` and ``__ne__`` to
    compare only the contents, not the field names.  Furthermore, this way no
    `FutureWarning` about comparisons is given.

    """

    # Note that it is important for physical type IDs to not be stored in a
    # tuple, since then the physical types would be treated as alternatives in
    # :meth:`~astropy.units.UnitBase.is_equivalent`.  (Of course, in that
    # case, they could also not be indexed by name.)

    def __eq__(self, other):
        if isinstance(other, np.void):
            other = other.item()

        return self.item() == other

    def __ne__(self, other):
        if isinstance(other, np.void):
            other = other.item()

        return self.item() != other


def _structured_unit_like_dtype(
    unit: UnitBase | StructuredUnit, dtype: np.dtype
) -> StructuredUnit:
    """Make a `StructuredUnit` of one unit, with the structure of a `numpy.dtype`.

    Parameters
    ----------
    unit : UnitBase
        The unit that will be filled into the structure.
    dtype : `numpy.dtype`
        The structure for the StructuredUnit.

    Returns
    -------
    StructuredUnit
    """
    if isinstance(unit, StructuredUnit):
        # If unit is structured, it should match the dtype. This function is
        # only used in Quantity, which performs this check, so it's fine to
        # return as is.
        return unit

    # Make a structured unit
    units = []
    for name in dtype.names:
        subdtype = dtype.fields[name][0]
        if subdtype.names is not None:
            units.append(_structured_unit_like_dtype(unit, subdtype))
        else:
            units.append(unit)
    return StructuredUnit(tuple(units), names=dtype.names)
