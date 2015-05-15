# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Core units classes and functions
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six
from ..extern.six.moves import zip
if six.PY2:
    import cmath

import inspect
import textwrap
import warnings
import numpy as np

from ..utils.decorators import lazyproperty, deprecated
from ..utils.exceptions import AstropyWarning
from ..utils.misc import isiterable, InheritDocstrings
from .utils import (is_effectively_unity, sanitize_scale, validate_power,
                    add_powers)
from . import format as unit_format

# TODO: Support functional units, e.g. log(x), ln(x)

__all__ = [
    'UnitsError', 'UnitsWarning', 'UnitConversionError', 'UnitBase',
    'NamedUnit', 'IrreducibleUnit', 'Unit', 'def_unit', 'CompositeUnit',
    'PrefixUnit', 'UnrecognizedUnit', 'get_current_unit_registry',
    'set_enabled_units', 'add_enabled_units',
    'set_enabled_equivalencies', 'add_enabled_equivalencies',
    'dimensionless_unscaled', 'one']


def _flatten_units_collection(items):
    """
    Given a list of sequences, modules or dictionaries of units, or
    single units, return a flat set of all the units found.
    """
    if not isinstance(items, list):
        items = [items]

    result = set()
    for item in items:
        if isinstance(item, UnitBase):
            result.add(item)
        else:
            if isinstance(item, dict):
                units = item.values()
            elif inspect.ismodule(item):
                units = vars(item).values()
            elif isiterable(item):
                units = item
            else:
                continue

            for unit in units:
                if isinstance(unit, UnitBase):
                    result.add(unit)

    return result


def _normalize_equivalencies(equivalencies):
    """
    Normalizes equivalencies, ensuring each is a 4-tuple of the form::

    (from_unit, to_unit, forward_func, backward_func)

    Parameters
    ----------
    equivalencies : list of equivalency pairs

    Raises
    ------
    ValueError if an equivalency cannot be interpreted
    """
    if equivalencies is None:
        return []

    normalized = []

    for i, equiv in enumerate(equivalencies):
        if len(equiv) == 2:
            funit, tunit = equiv
            a = b = lambda x: x
        elif len(equiv) == 3:
            funit, tunit, a = equiv
            b = a
        elif len(equiv) == 4:
            funit, tunit, a, b = equiv
        else:
            raise ValueError(
                "Invalid equivalence entry {0}: {1!r}".format(i, equiv))
        if not (isinstance(funit, UnitBase) and
                (isinstance(tunit, UnitBase) or tunit is None) and
                six.callable(a) and
                six.callable(b)):
            raise ValueError(
                "Invalid equivalence entry {0}: {1!r}".format(i, equiv))
        normalized.append((funit, tunit, a, b))

    return normalized


class _UnitRegistry(object):
    """
    Manages a registry of the enabled units.
    """
    def __init__(self, init=[], equivalencies=[]):
        if isinstance(init, _UnitRegistry):
            equivalencies = init.equivalencies
            init = init.all_units

        self._reset_units()
        self._reset_equivalencies()
        self.add_enabled_units(init)
        self.add_enabled_equivalencies(equivalencies)

    def _reset_units(self):
        self._all_units = set()
        self._non_prefix_units = set()
        self._registry = {}
        self._by_physical_type = {}

    def _reset_equivalencies(self):
        self._equivalencies = set()

    @property
    def registry(self):
        return self._registry

    @property
    def all_units(self):
        return self._all_units

    @property
    def non_prefix_units(self):
        return self._non_prefix_units

    def set_enabled_units(self, units):
        """
        Sets the units enabled in the unit registry.

        These units are searched when using
        `UnitBase.find_equivalent_units`, for example.

        Parameters
        ----------
        units : list of sequences, dicts, or modules containing units, or units
            This is a list of things in which units may be found
            (sequences, dicts or modules), or units themselves.  The
            entire set will be "enabled" for searching through by
            methods like `UnitBase.find_equivalent_units` and
            `UnitBase.compose`.
        """
        self._reset_units()
        return self.add_enabled_units(units)

    def add_enabled_units(self, units):
        """
        Adds to the set of units enabled in the unit registry.

        These units are searched when using
        `UnitBase.find_equivalent_units`, for example.

        Parameters
        ----------
        units : list of sequences, dicts, or modules containing units, or units
            This is a list of things in which units may be found
            (sequences, dicts or modules), or units themselves.  The
            entire set will be added to the "enabled" set for
            searching through by methods like
            `UnitBase.find_equivalent_units` and `UnitBase.compose`.
        """
        units = _flatten_units_collection(units)

        for unit in units:
            # Loop through all of the names first, to ensure all of them
            # are new, then add them all as a single "transaction" below.
            for st in unit._names:
                if (st in self._registry and unit != self._registry[st]):
                    raise ValueError(
                        "Object with name {0!r} already exists in namespace. "
                        "Filter the set of units to avoid name classes before "
                        "enabling them.".format(st))

            for st in unit._names:
                self._registry[st] = unit

            self._all_units.add(unit)
            if not isinstance(unit, PrefixUnit):
                self._non_prefix_units.add(unit)

            hash = unit._get_physical_type_id()
            self._by_physical_type.setdefault(hash, set()).add(unit)

    def get_units_with_physical_type(self, unit):
        """
        Get all units in the registry with the same physical type as
        the given unit.

        Parameters
        ----------
        unit : UnitBase instance
        """
        return self._by_physical_type.get(
            unit._get_physical_type_id(),
            set())

    @property
    def equivalencies(self):
        return list(self._equivalencies)

    def set_enabled_equivalencies(self, equivalencies):
        """
        Sets the equivalencies enabled in the unit registry.

        These equivalencies are used if no explicit equivalencies are given,
        both in unit conversion and in finding equivalent units.

        This is meant in particular for allowing angles to be dimensionless.
        Use with care.

        Parameters
        ----------
        equivalencies : list of equivalent pairs
            E.g., as returned by
            `~astropy.units.equivalencies.dimensionless_angles`.
        """
        self._reset_equivalencies()
        return self.add_enabled_equivalencies(equivalencies)

    def add_enabled_equivalencies(self, equivalencies):
        """
        Adds to the set of equivalencies enabled in the unit registry.

        These equivalencies are used if no explicit equivalencies are given,
        both in unit conversion and in finding equivalent units.

        This is meant in particular for allowing angles to be dimensionless.
        Use with care.

        Parameters
        ----------
        equivalencies : list of equivalent pairs
            E.g., as returned by
            `~astropy.units.equivalencies.dimensionless_angles`.
        """
        # pre-normalize list to help catch mistakes
        equivalencies = _normalize_equivalencies(equivalencies)
        self._equivalencies |= set(equivalencies)


class _UnitContext(object):
    def __init__(self, init=[], equivalencies=[]):
        _unit_registries.append(
            _UnitRegistry(
                init=init,
                equivalencies=equivalencies))

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        _unit_registries.pop()


_unit_registries = [_UnitRegistry()]


def get_current_unit_registry():
    return _unit_registries[-1]


def set_enabled_units(units):
    """
    Sets the units enabled in the unit registry.

    These units are searched when using
    `UnitBase.find_equivalent_units`, for example.

    This may be used either permanently, or as a context manager using
    the ``with`` statement (see example below).

    Parameters
    ----------
    units : list of sequences, dicts, or modules containing units, or units
        This is a list of things in which units may be found
        (sequences, dicts or modules), or units themselves.  The
        entire set will be "enabled" for searching through by methods
        like `UnitBase.find_equivalent_units` and `UnitBase.compose`.

    Examples
    --------

    >>> from astropy import units as u
    >>> with u.set_enabled_units([u.pc]):
    ...     u.m.find_equivalent_units()
    ...
      Primary name | Unit definition | Aliases
    [
      pc           | 3.08568e+16 m   | parsec  ,
    ]
    >>> u.m.find_equivalent_units()
      Primary name | Unit definition | Aliases
    [
      AU           | 1.49598e+11 m   | au, astronomical_unit ,
      Angstrom     | 1e-10 m         | AA, angstrom          ,
      cm           | 0.01 m          | centimeter            ,
      lyr          | 9.46073e+15 m   | lightyear             ,
      m            | irreducible     | meter                 ,
      micron       | 1e-06 m         |                       ,
      pc           | 3.08568e+16 m   | parsec                ,
      solRad       | 6.95508e+08 m   | R_sun, Rsun           ,
    ]
    """
    # get a context with a new registry, using equivalencies of the current one
    context = _UnitContext(
        equivalencies=get_current_unit_registry().equivalencies)
    # in this new current registry, enable the units requested
    get_current_unit_registry().set_enabled_units(units)
    return context


def add_enabled_units(units):
    """
    Adds to the set of units enabled in the unit registry.

    These units are searched when using
    `UnitBase.find_equivalent_units`, for example.

    This may be used either permanently, or as a context manager using
    the ``with`` statement (see example below).

    Parameters
    ----------
    units : list of sequences, dicts, or modules containing units, or units
        This is a list of things in which units may be found
        (sequences, dicts or modules), or units themselves.  The
        entire set will be added to the "enabled" set for searching
        through by methods like `UnitBase.find_equivalent_units` and
        `UnitBase.compose`.

    Examples
    --------

    >>> from astropy import units as u
    >>> from astropy.units import imperial
    >>> with u.add_enabled_units(imperial):
    ...     u.m.find_equivalent_units()
    ...
      Primary name | Unit definition | Aliases
    [
      AU           | 1.49598e+11 m   | au, astronomical_unit ,
      Angstrom     | 1e-10 m         | AA, angstrom          ,
      cm           | 0.01 m          | centimeter            ,
      ft           | 0.3048 m        | foot                  ,
      fur          | 201.168 m       | furlong               ,
      inch         | 0.0254 m        |                       ,
      lyr          | 9.46073e+15 m   | lightyear             ,
      m            | irreducible     | meter                 ,
      mi           | 1609.34 m       | mile                  ,
      micron       | 1e-06 m         |                       ,
      mil          | 2.54e-05 m      | thou                  ,
      nmi          | 1852 m          | nauticalmile, NM      ,
      pc           | 3.08568e+16 m   | parsec                ,
      solRad       | 6.95508e+08 m   | R_sun, Rsun           ,
      yd           | 0.9144 m        | yard                  ,
    ]
    """
    # get a context with a new registry, which is a copy of the current one
    context = _UnitContext(get_current_unit_registry())
    # in this new current registry, enable the further units requested
    get_current_unit_registry().add_enabled_units(units)
    return context


def set_enabled_equivalencies(equivalencies):
    """
    Sets the equivalencies enabled in the unit registry.

    These equivalencies are used if no explicit equivalencies are given,
    both in unit conversion and in finding equivalent units.

    This is meant in particular for allowing angles to be dimensionless.
    Use with care.

    Parameters
    ----------
    equivalencies : list of equivalent pairs
        E.g., as returned by
        `~astropy.units.equivalencies.dimensionless_angles`.

    Examples
    --------
    Exponentiation normally requires dimensionless quantities.  To avoid
    problems with complex phases::

        >>> from astropy import units as u
        >>> with u.set_enabled_equivalencies(u.dimensionless_angles()):
        ...     phase = 0.5 * u.cycle
        ...     np.exp(1j*phase)  # doctest: +FLOAT_CMP
        <Quantity (-1+1.2246063538223773e-16j)>
    """
    # get a context with a new registry, using all units of the current one
    context = _UnitContext(get_current_unit_registry().all_units)
    # in this new current registry, enable the equivalencies requested
    get_current_unit_registry().set_enabled_equivalencies(equivalencies)
    return context


def add_enabled_equivalencies(equivalencies):
    """
    Adds to the equivalencies enabled in the unit registry.

    These equivalencies are used if no explicit equivalencies are given,
    both in unit conversion and in finding equivalent units.

    This is meant in particular for allowing angles to be dimensionless.
    Since no equivalencies are enabled by default, generally it is recommended
    to use `set_enabled_equivalencies`.

    Parameters
    ----------
    equivalencies : list of equivalent pairs
        E.g., as returned by
        `~astropy.units.equivalencies.dimensionless_angles`.
    """
    # get a context with a new registry, which is a copy of the current one
    context = _UnitContext(get_current_unit_registry())
    # in this new current registry, enable the further equivalencies requested
    get_current_unit_registry().add_enabled_equivalencies(equivalencies)
    return context


class UnitsError(Exception):
    """
    The base class for unit-specific exceptions.
    """


class UnitScaleError(UnitsError, ValueError):
    """
    Used to catch the errors involving scaled units,
    which are not recognized by FITS format.
    """
    pass


class UnitConversionError(UnitsError, ValueError):
    """
    Used specifically for errors related to converting between units or
    interpreting units in terms of other units.
    """


# Maintain error in old location for backward compatibility
from .format import fits as _fits
_fits.UnitScaleError = UnitScaleError


class UnitsWarning(AstropyWarning):
    """
    The base class for unit-specific exceptions.
    """


@six.add_metaclass(InheritDocstrings)
class UnitBase(object):
    """
    Abstract base class for units.

    Most of the arithmetic operations on units are defined in this
    base class.

    Should not be instantiated by users directly.
    """
    # Make sure that __rmul__ of units gets called over the __mul__ of Numpy
    # arrays to avoid element-wise multiplication.
    __array_priority__ = 1000

    def __deepcopy__(self, memo):
        # This may look odd, but the units conversion will be very
        # broken after deep-copying if we don't guarantee that a given
        # physical unit corresponds to only one instance
        return self

    def _repr_latex_(self):
        """
        Generate latex representation of unit name.  This is used by
        the IPython notebook to print a unit with a nice layout.

        Returns
        -------
        Latex string
        """
        return unit_format.Latex().to_string(self)

    def __bytes__(self):
        """Return string representation for unit"""
        return unit_format.Generic().to_string(self).encode('ascii')
    if six.PY2:
        __str__ = __bytes__

    def __unicode__(self):
        """Return string representation for unit"""
        return unit_format.Generic().to_string(self)
    if six.PY3:
        __str__ = __unicode__

    def __repr__(self):
        string = unit_format.Generic().to_string(self)
        if six.PY2:
            string = string.encode('unicode_escape')

        return 'Unit("{0}")'.format(string)

    def _get_physical_type_id(self):
        """
        Returns an identifier that uniquely identifies the physical
        type of this unit.  It is comprised of the bases and powers of
        this unit, without the scale.  Since it is hashable, it is
        useful as a dictionary key.
        """
        unit = self.decompose()
        r = zip([x.name for x in unit.bases], unit.powers)
        # bases and powers are already sorted in a unique way
        # r.sort()
        r = tuple(r)
        return r

    @property
    def names(self):
        """
        Returns all of the names associated with this unit.
        """
        raise AttributeError(
            "Can not get names from unnamed units. "
            "Perhaps you meant to_string()?")
        return self._names

    @property
    def name(self):
        """
        Returns the canonical (short) name associated with this unit.
        """
        raise AttributeError(
            "Can not get names from unnamed units. "
            "Perhaps you meant to_string()?")

    @property
    def aliases(self):
        """
        Returns the alias (long) names for this unit.
        """
        raise AttributeError(
            "Can not get aliases from unnamed units. "
            "Perhaps you meant to_string()?")

    @property
    def scale(self):
        """
        Return the scale of the unit.
        """
        return 1.0

    @property
    def bases(self):
        """
        Return the bases of the unit.
        """
        return [self]

    @property
    def powers(self):
        """
        Return the powers of the unit.
        """
        return [1]

    def to_string(self, format=unit_format.Generic):
        """
        Output the unit in the given format as a string.

        Parameters
        ----------
        format : `astropy.units.format.Base` instance or str
            The name of a format or a formatter object.  If not
            provided, defaults to the generic format.
        """
        f = unit_format.get_format(format)
        return f.to_string(self)

    def __format__(self, format_spec):
        """Try to format units using a formatter."""
        try:
            return self.to_string(format=format_spec)
        except ValueError:
            return format(six.text_type(self), format_spec)

    @staticmethod
    def _normalize_equivalencies(equivalencies):
        """
        Normalizes equivalencies, ensuring each is a 4-tuple of the form::

        (from_unit, to_unit, forward_func, backward_func)

        Parameters
        ----------
        equivalencies : list of equivalency pairs, or `None`

        Returns
        -------
        A normalized list, including possible global defaults set by, e.g.,
        `set_enabled_equivalencies`, except when `equivalencies`=`None`,
        in which case the returned list is always empty.

        Raises
        ------
        ValueError if an equivalency cannot be interpreted
        """
        normalized = _normalize_equivalencies(equivalencies)
        if equivalencies is not None:
            normalized += get_current_unit_registry().equivalencies

        return normalized

    def __pow__(self, p):
        return CompositeUnit(1, [self], [p])

    def __div__(self, m):
        if isinstance(m, (bytes, six.text_type)):
            m = Unit(m)

        if isinstance(m, UnitBase):
            if m.is_unity():
                return self
            return CompositeUnit(1, [self, m], [1, -1], _error_check=False)

        # Cannot handle this as Unit, re-try as Quantity
        from .quantity import Quantity
        return Quantity(1, self) / m

    def __rdiv__(self, m):
        if isinstance(m, (bytes, six.text_type)):
            return Unit(m) / self

        # Cannot handle this as Unit.  Here, m cannot be a Quantity,
        # so we make it into one, fasttracking when it does not have a unit,
        # for the common case of <array> / <unit>.
        from .quantity import Quantity
        if hasattr(m, 'unit'):
            result = Quantity(m)
            result /= self
            return result
        else:
            return Quantity(m, self**(-1))

    __truediv__ = __div__

    __rtruediv__ = __rdiv__

    def __mul__(self, m):
        if isinstance(m, (bytes, six.text_type)):
            m = Unit(m)

        if isinstance(m, UnitBase):
            if m.is_unity():
                return self
            elif self.is_unity():
                return m
            return CompositeUnit(1, [self, m], [1, 1], _error_check=False)

        # Cannot handle this as Unit, re-try as Quantity.
        from .quantity import Quantity
        return Quantity(1, self) * m

    def __rmul__(self, m):
        if isinstance(m, (bytes, six.text_type)):
            return Unit(m) * self

        # Cannot handle this as Unit.  Here, m cannot be a Quantity,
        # so we make it into one, fasttracking when it does not have a unit
        # for the common case of <array> * <unit>.
        from .quantity import Quantity
        if hasattr(m, 'unit'):
            result = Quantity(m)
            result *= self
            return result
        else:
            return Quantity(m, self)

    def __hash__(self):
        # This must match the hash used in CompositeUnit for a unit
        # with only one base and no scale or power.
        return hash((str(self.scale), self.name, str('1')))

    def __eq__(self, other):
        if self is other:
            return True

        try:
            other = Unit(other, parse_strict='silent')
        except (ValueError, UnitsError, TypeError):
            return False
        try:
            return is_effectively_unity(self._to(other))
        except UnitsError:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __le__(self, other):
        scale = self._to(Unit(other))
        return scale <= 1. or is_effectively_unity(scale)

    def __ge__(self, other):
        scale = self._to(Unit(other))
        return scale >= 1. or is_effectively_unity(scale)

    def __lt__(self, other):
        return not (self >= other)

    def __gt__(self, other):
        return not (self <= other)

    def __neg__(self):
        return self * -1.

    def is_equivalent(self, other, equivalencies=[]):
        """
        Returns `True` if this unit is equivalent to ``other``.

        Parameters
        ----------
        other : unit object or string or tuple
            The unit to convert to. If a tuple of units is specified, this
            method returns true if the unit matches any of those in the tuple.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            This list is in addition to possible global defaults set by, e.g.,
            `set_enabled_equivalencies`.
            Use `None` to turn off all equivalencies.

        Returns
        -------
        bool
        """
        equivalencies = self._normalize_equivalencies(equivalencies)

        if isinstance(other, tuple):
            return any(self.is_equivalent(u, equivalencies=equivalencies)
                       for u in other)
        else:
            other = Unit(other, parse_strict='silent')

        return self._is_equivalent(other, equivalencies)

    def _is_equivalent(self, other, equivalencies=[]):
        """Returns `True` if this unit is equivalent to `other`.
        See `is_equivalent`, except that a proper Unit object should be
        given (i.e., no string) and that the equivalency list should be
        normalized using `_normalize_equivalencies`.
        """
        if isinstance(other, UnrecognizedUnit):
            return False

        if (self._get_physical_type_id() ==
                other._get_physical_type_id()):
            return True
        elif len(equivalencies):
            unit = self.decompose()
            other = other.decompose()
            for a, b, forward, backward in equivalencies:
                if b is None:
                    # after canceling, is what's left convertable
                    # to dimensionless (according to the equivalency)?
                    try:
                        (other/unit).decompose([a])
                        return True
                    except:
                        pass
                else:
                    if(a._is_equivalent(unit) and b._is_equivalent(other) or
                       b._is_equivalent(unit) and a._is_equivalent(other)):
                        return True

        return False

    def _apply_equivalencies(self, unit, other, equivalencies):
        """
        Internal function (used from `_get_converter`) to apply
        equivalence pairs.
        """
        def make_converter(scale1, func, scale2):
            def convert(v):
                return func(_condition_arg(v) / scale1) * scale2
            return convert

        orig_unit = unit
        orig_other = other

        unit = self.decompose()
        other = other.decompose()

        for funit, tunit, a, b in equivalencies:
            if tunit is None:
                try:
                    ratio_in_funit = (other/unit).decompose([funit])
                    return make_converter(ratio_in_funit.scale, a, 1.)
                except UnitsError:
                    pass
            else:
                try:
                    scale1 = funit._to(unit)
                    scale2 = tunit._to(other)
                    return make_converter(scale1, a, scale2)
                except UnitsError:
                    pass
                try:
                    scale1 = tunit._to(unit)
                    scale2 = funit._to(other)
                    return make_converter(scale1, b, scale2)
                except UnitsError:
                    pass

        def get_err_str(unit):
            unit_str = unit.to_string('unscaled')
            physical_type = unit.physical_type
            if physical_type != 'unknown':
                unit_str = "'{0}' ({1})".format(
                    unit_str, physical_type)
            else:
                unit_str = "'{0}'".format(unit_str)
            return unit_str

        unit_str = get_err_str(orig_unit)
        other_str = get_err_str(orig_other)

        raise UnitConversionError(
            "{0} and {1} are not convertible".format(
                unit_str, other_str))

    def _get_converter(self, other, equivalencies=[]):
        other = Unit(other)

        try:
            scale = self._to(other)
        except UnitsError:
            return self._apply_equivalencies(
                self, other, self._normalize_equivalencies(equivalencies))
        return lambda val: scale * _condition_arg(val)

    @deprecated('1.0')
    def get_converter(self, other, equivalencies=[]):
        """
        Return the conversion function to convert values from ``self``
        to the specified unit.

        Parameters
        ----------
        other : unit object or string
            The unit to convert to.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            This list is in addition to possible global defaults set by, e.g.,
            `set_enabled_equivalencies`.
            Use `None` to turn off all equivalencies.

        Returns
        -------
        func : callable
            A callable that normally expects a single argument that is
            a scalar value or an array of values (or anything that may
            be converted to an array).

        Raises
        ------
        UnitsError
            If units are inconsistent
        """
        return self._get_converter(other, equivalencies=equivalencies)

    def _to(self, other):
        """
        Returns the scale to the specified unit.

        See `to`, except that a Unit object should be given (i.e., no
        string), and that all defaults are used, i.e., no
        equivalencies and value=1.
        """
        # There are many cases where we just want to ensure a Quantity is
        # of a particular unit, without checking whether it's already in
        # a particular unit.  If we're being asked to convert from a unit
        # to itself, we can short-circuit all of this.
        if self is other:
            return 1.0

        self_decomposed = self.decompose()
        other_decomposed = other.decompose()

        # Check quickly whether equivalent.  This is faster than
        # `is_equivalent`, because it doesn't generate the entire
        # physical type list of both units.  In other words it "fails
        # fast".
        if(self_decomposed.powers == other_decomposed.powers and
           all(self_base is other_base for (self_base, other_base)
               in zip(self_decomposed.bases, other_decomposed.bases))):
            return self_decomposed.scale / other_decomposed.scale

        raise UnitConversionError(
            "'{0!r}' is not a scaled version of '{1!r}'".format(self, other))

    def to(self, other, value=1.0, equivalencies=[]):
        """
        Return the converted values in the specified unit.

        Parameters
        ----------
        other : unit object or string
            The unit to convert to.

        value : scalar int or float, or sequence convertable to array, optional
            Value(s) in the current unit to be converted to the
            specified unit.  If not provided, defaults to 1.0

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            This list is in addition to possible global defaults set by, e.g.,
            `set_enabled_equivalencies`.
            Use `None` to turn off all equivalencies.

        Returns
        -------
        values : scalar or array
            Converted value(s). Input value sequences are returned as
            numpy arrays.

        Raises
        ------
        UnitsError
            If units are inconsistent
        """
        return self._get_converter(other, equivalencies=equivalencies)(value)

    def in_units(self, other, value=1.0, equivalencies=[]):
        """
        Alias for `to` for backward compatibility with pynbody.
        """
        return self.to(
            other, value=value, equivalencies=equivalencies)

    def decompose(self, bases=set()):
        """
        Return a unit object composed of only irreducible units.

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
        unit : CompositeUnit object
            New object containing only irreducible unit objects.
        """
        raise NotImplementedError()

    def _compose(self, equivalencies=[], namespace=[], max_depth=2, depth=0,
                 cached_results=None):
        def is_final_result(unit):
            # Returns True if this result contains only the expected
            # units
            for base in unit.bases:
                if base not in namespace:
                    return False
            return True

        unit = self.decompose()
        key = hash(unit)

        cached = cached_results.get(key)
        if cached is not None:
            if isinstance(cached, Exception):
                raise cached
            return cached

        # Prevent too many levels of recursion
        # And special case for dimensionless unit
        if depth >= max_depth:
            cached_results[key] = [unit]
            return [unit]

        # Make a list including all of the equivalent units
        units = [unit]
        for funit, tunit, a, b in equivalencies:
            if tunit is not None:
                if self._is_equivalent(funit):
                    scale = funit.decompose().scale / unit.scale
                    units.append(Unit(a(1.0 / scale) * tunit).decompose())
                elif self._is_equivalent(tunit):
                    scale = tunit.decompose().scale / unit.scale
                    units.append(Unit(b(1.0 / scale) * funit).decompose())
            else:
                if self._is_equivalent(funit):
                    units.append(Unit(unit.scale))

        # Store partial results
        partial_results = []
        # Store final results that reduce to a single unit or pair of
        # units
        if len(unit.bases) == 0:
            final_results = [set([unit]), set()]
        else:
            final_results = [set(), set()]

        for tunit in namespace:
            tunit_decomposed = tunit.decompose()
            for u in units:
                # If the unit is a base unit, look for an exact match
                # to one of the bases of the target unit.  If found,
                # factor by the same power as the target unit's base.
                # This allows us to factor out fractional powers
                # without needing to do an exhaustive search.
                if len(tunit_decomposed.bases) == 1:
                    for base, power in zip(u.bases, u.powers):
                        if tunit_decomposed._is_equivalent(base):
                            tunit = tunit ** power
                            tunit_decomposed = tunit_decomposed ** power
                            break

                composed = (u / tunit_decomposed).decompose()
                factored = composed * tunit
                len_bases = len(composed.bases)
                if is_final_result(factored) and len_bases <= 1:
                    final_results[len_bases].add(factored)
                else:
                    partial_results.append(
                        (len_bases, composed, tunit))

        # Do we have any minimal results?
        for final_result in final_results:
            if len(final_result):
                results = final_results[0].union(final_results[1])
                cached_results[key] = results
                return results

        partial_results.sort(key=lambda x: x[0])

        # ...we have to recurse and try to further compose
        results = []
        for len_bases, composed, tunit in partial_results:
            try:
                composed_list = composed._compose(
                    equivalencies=equivalencies,
                    namespace=namespace,
                    max_depth=max_depth, depth=depth + 1,
                    cached_results=cached_results)
            except UnitsError:
                composed_list = []
            for subcomposed in composed_list:
                results.append(
                    (len(subcomposed.bases), subcomposed, tunit))

        if len(results):
            results.sort(key=lambda x: x[0])

            min_length = results[0][0]
            subresults = set()
            for len_bases, composed, tunit in results:
                if len_bases > min_length:
                    break
                else:
                    factored = composed * tunit
                    if is_final_result(factored):
                        subresults.add(factored)

            if len(subresults):
                cached_results[key] = subresults
                return subresults

        if not is_final_result(self):
            result = UnitsError(
                "Cannot represent unit {0} in terms of the given "
                "units".format(self))
            cached_results[key] = result
            raise result

        cached_results[key] = [self]
        return [self]

    def compose(self, equivalencies=[], units=None, max_depth=2,
                include_prefix_units=False):
        """
        Return the simplest possible composite unit(s) that represent
        the given unit.  Since there may be multiple equally simple
        compositions of the unit, a list of units is always returned.

        Parameters
        ----------
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to also list.  See
            :ref:`unit_equivalencies`.
            This list is in addition to possible global defaults set by, e.g.,
            `set_enabled_equivalencies`.
            Use `None` to turn off all equivalencies.

        units : set of units to compose to, optional
            If not provided, any known units may be used to compose
            into.  Otherwise, ``units`` is a dict, module or sequence
            containing the units to compose into.

        max_depth : int, optional
            The maximum recursion depth to use when composing into
            composite units.

        include_prefix_units : bool, optional
            When `True`, include prefixed units in the result.
            Default is `False`.

        Returns
        -------
        units : list of `CompositeUnit`
            A list of candidate compositions.  These will all be
            equally simple, but it may not be possible to
            automatically determine which of the candidates are
            better.
        """
        # Pre-normalize the equivalencies list
        equivalencies = self._normalize_equivalencies(equivalencies)

        # The namespace of units to compose into should be filtered to
        # only include units with bases in common with self, otherwise
        # they can't possibly provide useful results.  Having too many
        # destination units greatly increases the search space.

        def has_bases_in_common(a, b):
            if len(a.bases) == 0 and len(b.bases) == 0:
                return True
            for ab in a.bases:
                for bb in b.bases:
                    if ab == bb:
                        return True
            return False

        def has_bases_in_common_with_equiv(unit, other):
            if has_bases_in_common(unit, other):
                return True
            for funit, tunit, a, b in equivalencies:
                if tunit is not None:
                    if unit._is_equivalent(funit):
                        if has_bases_in_common(tunit.decompose(), other):
                            return True
                    elif unit._is_equivalent(tunit):
                        if has_bases_in_common(funit.decompose(), other):
                            return True
                else:
                    if unit._is_equivalent(funit):
                        if has_bases_in_common(dimensionless_unscaled, other):
                            return True
            return False

        def filter_units(units):
            filtered_namespace = set()
            for tunit in units:
                if (isinstance(tunit, UnitBase) and
                    (include_prefix_units or
                     not isinstance(tunit, PrefixUnit)) and
                    has_bases_in_common_with_equiv(
                        decomposed, tunit.decompose())):
                    filtered_namespace.add(tunit)
            return filtered_namespace

        decomposed = self.decompose()

        if units is None:
            units = filter_units(self._get_units_with_same_physical_type(
                equivalencies=equivalencies))
            if len(units) == 0:
                units = get_current_unit_registry().non_prefix_units
        elif isinstance(units, dict):
            units = set(filter_units(six.itervalues(units)))
        elif inspect.ismodule(units):
            units = filter_units(six.itervalues(vars(units)))
        else:
            units = filter_units(_flatten_units_collection(units))

        def sort_results(results):
            if not len(results):
                return []

            # Sort the results so the simplest ones appear first.
            # Simplest is defined as "the minimum sum of absolute
            # powers" (i.e. the fewest bases), and preference should
            # be given to results where the sum of powers is positive
            # and the scale is exactly equal to 1.0
            results = list(results)
            results.sort(key=lambda x: np.abs(x.scale))
            results.sort(key=lambda x: np.sum(np.abs(x.powers)))
            results.sort(key=lambda x: np.sum(x.powers) < 0.0)
            results.sort(key=lambda x: not is_effectively_unity(x.scale))

            last_result = results[0]
            filtered = [last_result]
            for result in results[1:]:
                if str(result) != str(last_result):
                    filtered.append(result)
                last_result = result

            return filtered

        return sort_results(self._compose(
            equivalencies=equivalencies, namespace=units,
            max_depth=max_depth, depth=0, cached_results={}))

    def to_system(self, system):
        """
        Converts this unit into ones belonging to the given system.
        Since more than one result may be possible, a list is always
        returned.

        Parameters
        ----------
        system : module
            The module that defines the unit system.  Commonly used
            ones include `astropy.units.si` and `astropy.units.cgs`.

            To use your own module it must contain unit objects and a
            sequence member named ``bases`` containing the base units of
            the system.

        Returns
        -------
        units : list of `CompositeUnit`
            The list is ranked so that units containing only the base
            units of that system will appear first.
        """
        bases = set(system.bases)

        def score(compose):
            # In case that compose._bases has no elements we return
            # 'np.inf' as 'score value'.  It does not really matter which
            # number we would return. This case occurs for instance for
            # dimensionless quantities:
            compose_bases = compose.bases
            if len(compose_bases) == 0:
                return np.inf
            else:
                sum = 0
                for base in compose_bases:
                    if base in bases:
                        sum += 1

                return sum / float(len(compose_bases))

        x = self.decompose(bases=bases)
        composed = x.compose(units=system)
        composed = sorted(composed, key=score, reverse=True)
        return composed

    @lazyproperty
    def si(self):
        """
        Returns a copy of the current `Unit` instance in SI units.
        """

        from . import si
        return self.to_system(si)[0]

    @lazyproperty
    def cgs(self):
        """
        Returns a copy of the current `Unit` instance with CGS units.
        """
        from . import cgs
        return self.to_system(cgs)[0]

    @property
    def physical_type(self):
        """
        Return the physical type on the unit.

        Examples
        --------
        >>> from astropy import units as u
        >>> print(u.m.physical_type)
        length

        """
        from . import physical
        return physical.get_physical_type(self)

    def _get_units_with_same_physical_type(self, equivalencies=[]):
        """
        Return a list of registered units with the same physical type
        as this unit.

        This function is used by Quantity to add its built-in
        conversions to equivalent units.

        This is a private method, since end users should be encouraged
        to use the more powerful `compose` and `find_equivalent_units`
        methods (which use this under the hood).

        Parameters
        ----------
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to also pull options from.
            See :ref:`unit_equivalencies`.  It must already be
            normalized using `_normalize_equivalencies`.
        """
        unit_registry = get_current_unit_registry()
        units = set(unit_registry.get_units_with_physical_type(self))
        for funit, tunit, a, b in equivalencies:
            if tunit is not None:
                if self.is_equivalent(funit) and tunit not in units:
                    units.update(
                        unit_registry.get_units_with_physical_type(tunit))
                if self._is_equivalent(tunit) and funit not in units:
                    units.update(
                        unit_registry.get_units_with_physical_type(funit))
            else:
                if self.is_equivalent(funit):
                    units.add(dimensionless_unscaled)
        return units

    class EquivalentUnitsList(list):
        """
        A class to handle pretty-printing the result of
        `find_equivalent_units`.
        """
        def __repr__(self):
            if len(self) == 0:
                return "[]"
            else:
                lines = []
                for u in self:
                    irred = u.decompose().to_string()
                    if irred == u.name:
                        irred = "irreducible"
                    lines.append((u.name, irred, ', '.join(u.aliases)))

                lines.sort()
                lines.insert(0, ('Primary name', 'Unit definition', 'Aliases'))
                widths = [0, 0, 0]
                for line in lines:
                    for i, col in enumerate(line):
                        widths[i] = max(widths[i], len(col))

                f = "  {{0:<{0}s}} | {{1:<{1}s}} | {{2:<{2}s}}".format(*widths)
                lines = [f.format(*line) for line in lines]
                lines = (lines[0:1] +
                         ['['] +
                         ['{0} ,'.format(x) for x in lines[1:]] +
                         [']'])
                return '\n'.join(lines)

    def find_equivalent_units(self, equivalencies=[], units=None,
                              include_prefix_units=False):
        """
        Return a list of all the units that are the same type as ``self``.

        Parameters
        ----------
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to also list.  See
            :ref:`unit_equivalencies`.
            Any list given, including an empty one, supercedes global defaults
            that may be in effect (as set by `set_enabled_equivalencies`)

        units : set of units to search in, optional
            If not provided, all defined units will be searched for
            equivalencies.  Otherwise, may be a dict, module or
            sequence containing the units to search for equivalencies.

        include_prefix_units : bool, optional
            When `True`, include prefixed units in the result.
            Default is `False`.

        Returns
        -------
        units : list of `UnitBase`
            A list of unit objects that match ``u``.  A subclass of
            `list` (``EquivalentUnitsList``) is returned that
            pretty-prints the list of units when output.
        """
        results = self.compose(
            equivalencies=equivalencies, units=units, max_depth=1,
            include_prefix_units=include_prefix_units)
        results = set(
            x.bases[0] for x in results if len(x.bases) == 1)
        return self.EquivalentUnitsList(results)

    def is_unity(self):
        """
        Returns `True` if the unit is unscaled and dimensionless.
        """
        return False


class NamedUnit(UnitBase):
    """
    The base class of units that have a name.

    Parameters
    ----------
    st : str, list of str, 2-tuple
        The name of the unit.  If a list of strings, the first element
        is the canonical (short) name, and the rest of the elements
        are aliases.  If a tuple of lists, the first element is a list
        of short names, and the second element is a list of long
        names; all but the first short name are considered "aliases".
        Each name *should* be a valid Python identifier to make it
        easy to access, but this is not required.

    namespace : dict, optional
        When provided, inject the unit, and all of its aliases, in the
        given namespace dictionary.  If a unit by the same name is
        already in the namespace, a ValueError is raised.

    doc : str, optional
        A docstring describing the unit.

    format : dict, optional
        A mapping to format-specific representations of this unit.
        For example, for the ``Ohm`` unit, it might be nice to have it
        displayed as ``\\Omega`` by the ``latex`` formatter.  In that
        case, `format` argument should be set to::

            {'latex': r'\\Omega'}

    Raises
    ------
    ValueError
        If any of the given unit names are already in the registry.

    ValueError
        If any of the given unit names are not valid Python tokens.
    """
    def __init__(self, st, register=None, doc=None, format=None,
                 namespace=None):
        if register is not None:
            warnings.warn(
                "The registry kwarg was removed in astropy 0.3. "
                "Use the namespace kwarg to inject the unit into "
                "a namespace and add_enabled_units() to enable it "
                "in the global unit registry.",
                DeprecationWarning)

        UnitBase.__init__(self)

        if isinstance(st, (bytes, six.text_type)):
            self._names = [st]
            self._short_names = [st]
            self._long_names = []
        elif isinstance(st, tuple):
            if not len(st) == 2:
                raise ValueError("st must be string, list or 2-tuple")
            self._names = st[0] + st[1]
            if not len(self._names):
                raise ValueError("must provide at least one name")
            self._short_names = st[0][:]
            self._long_names = st[1][:]
        else:
            if len(st) == 0:
                raise ValueError(
                    "st list must have at least one entry")
            self._names = st[:]
            self._short_names = [st[0]]
            self._long_names = st[1:]

        if format is None:
            format = {}
        self._format = format

        if doc is None:
            doc = self._generate_doc()
        else:
            doc = textwrap.dedent(doc)
            doc = textwrap.fill(doc)

        self.__doc__ = doc

        self._inject(namespace)

    def _generate_doc(self):
        """
        Generate a docstring for the unit if the user didn't supply
        one.  This is only used from the constructor and may be
        overridden in subclasses.
        """
        names = self.names
        if len(self.names) > 1:
            return "{1} ({0})".format(*names[:2])
        else:
            return names[0]

    def get_format_name(self, format):
        """
        Get a name for this unit that is specific to a particular
        format.

        Uses the dictionary passed into the `format` kwarg in the
        constructor.

        Parameters
        ----------
        format : str
            The name of the format

        Returns
        -------
        name : str
            The name of the unit for the given format.
        """
        return self._format.get(format, self.name)

    @property
    def names(self):
        """
        Returns all of the names associated with this unit.
        """
        return self._names

    @property
    def name(self):
        """
        Returns the canonical (short) name associated with this unit.
        """
        return self._names[0]

    @property
    def aliases(self):
        """
        Returns the alias (long) names for this unit.
        """
        return self._names[1:]

    @property
    def short_names(self):
        """
        Returns all of the short names associated with this unit.
        """
        return self._short_names

    @property
    def long_names(self):
        """
        Returns all of the long names associated with this unit.
        """
        return self._long_names

    def register(self, add_to_namespace=False):
        raise NotImplementedError(
            "The register method has been removed in astropy 0.3. "
            "Use add_enabled_units/set_enabled_units instead.")

    def deregister(self, remove_from_namespace=False):
        raise NotImplementedError(
            "The deregister method has been removed in astropy 0.3. "
            "Use add_enabled_units/set_enabled_units instead.")

    def _inject(self, namespace=None):
        """
        Injects the unit, and all of its aliases, in the given
        namespace dictionary.
        """
        if namespace is None:
            return

        # Loop through all of the names first, to ensure all of them
        # are new, then add them all as a single "transaction" below.
        for name in self._names:
            if name in namespace and self != namespace[name]:
                raise ValueError(
                    "Object with name {0!r} already exists in "
                    "given namespace ({1!r}).".format(
                        name, namespace[name]))

        for name in self._names:
            namespace[name] = self


def _recreate_irreducible_unit(cls, names, registered):
    """
    This is used to reconstruct units when passed around by
    multiprocessing.
    """
    registry = get_current_unit_registry().registry
    if names[0] in registry:
        return registry[names[0]]
    else:
        unit = cls(names)
        if registered:
            get_current_unit_registry().add_enabled_units([unit])


class IrreducibleUnit(NamedUnit):
    """
    Irreducible units are the units that all other units are defined
    in terms of.

    Examples are meters, seconds, kilograms, amperes, etc.  There is
    only once instance of such a unit per type.
    """
    def __reduce__(self):
        # When IrreducibleUnit objects are passed to other processes
        # over multiprocessing, they need to be recreated to be the
        # ones already in the subprocesses' namespace, not new
        # objects, or they will be considered "unconvertible".
        # Therefore, we have a custom pickler/unpickler that
        # understands how to recreate the Unit on the other side.
        registry = get_current_unit_registry().registry
        return (_recreate_irreducible_unit,
                (self.__class__, list(self.names), self.name in registry),
                self.__dict__)

    def decompose(self, bases=set()):
        if len(bases) and not self in bases:
            for base in bases:
                try:
                    scale = self._to(base)
                except UnitsError:
                    pass
                else:
                    if is_effectively_unity(scale):
                        return base
                    else:
                        return CompositeUnit(scale, [base], [1],
                                             _error_check=False)

            raise UnitConversionError(
                "Unit {0} can not be decomposed into the requested "
                "bases".format(self))

        return self


class UnrecognizedUnit(IrreducibleUnit):
    """
    A unit that did not parse correctly.  This allows for
    roundtripping it as a string, but no unit operations actually work
    on it.

    Parameters
    ----------
    st : str
        The name of the unit.
    """
    # For UnrecognizedUnits, we want to use "standard" Python
    # pickling, not the special case that is used for
    # IrreducibleUnits.
    __reduce__ = object.__reduce__

    def __repr__(self):
        return "UnrecognizedUnit({0})".format(str(self))

    def __bytes__(self):
        return self.name.encode('ascii', 'replace')
    if six.PY2:
        __str__ = __bytes__

    def __unicode__(self):
        return self.name
    if six.PY3:
        __str__ = __unicode__

    def to_string(self, format=unit_format.Generic):
        return self.name

    def _unrecognized_operator(self, *args, **kwargs):
        raise ValueError(
            "The unit {0!r} is unrecognized, so all arithmetic operations "
            "with it are invalid.".format(self.name))

    __pow__ = __div__ = __rdiv__ = __truediv__ = __rtruediv__ = __mul__ = \
        __rmul__ = __lt__ = __gt__ = __le__ = __ge__ = __neg__ = \
        _unrecognized_operator

    def __eq__(self, other):
        other = Unit(other, parse_strict='silent')
        return isinstance(other, UnrecognizedUnit) and self.name == other.name

    def __ne__(self, other):
        return not (self == other)

    def is_equivalent(self, other, equivalencies=None):
        self._normalize_equivalencies(equivalencies)
        return self == other

    def _get_converter(self, other, equivalencies=None):
        self._normalize_equivalencies(equivalencies)
        raise ValueError(
            "The unit {0!r} is unrecognized.  It can not be converted "
            "to other units.".format(self.name))

    def get_format_name(self, format):
        return self.name

    def is_unity(self):
        return False


class _UnitMetaClass(InheritDocstrings):
    """
    This metaclass exists because the Unit constructor should
    sometimes return instances that already exist.  This "overrides"
    the constructor before the new instance is actually created, so we
    can return an existing one.
    """

    def __call__(self, s, represents=None, format=None, namespace=None,
                 doc=None, parse_strict='raise'):

        # Short-circuit if we're already a unit
        if isinstance(s, UnitBase):
            return s

        # turn possible Quantity input for s or represents into a Unit
        from .quantity import Quantity

        if isinstance(represents, Quantity):
            if is_effectively_unity(represents.value):
                represents = represents.unit
            else:
                # cannot use _error_check=False: scale may be effectively unity
                represents = CompositeUnit(represents.value *
                                           represents.unit.scale,
                                           bases=represents.unit.bases,
                                           powers=represents.unit.powers)

        if isinstance(s, Quantity):
            if is_effectively_unity(s.value):
                s = s.unit
            else:
                s = CompositeUnit(s.value * s.unit.scale,
                                  bases=s.unit.bases,
                                  powers=s.unit.powers)

        # now decide what we really need to do; define derived Unit?
        if isinstance(represents, UnitBase):
            # This has the effect of calling the real __new__ and
            # __init__ on the Unit class.
            return super(_UnitMetaClass, self).__call__(
                s, represents, format=format, namespace=namespace, doc=doc)

        # or interpret a Quantity (now became unit), string or number?
        if isinstance(s, UnitBase):
            return s

        elif isinstance(s, (bytes, six.text_type)):
            if len(s.strip()) == 0:
                # Return the NULL unit
                return dimensionless_unscaled

            if format is None:
                format = unit_format.Generic

            f = unit_format.get_format(format)
            if six.PY3 and isinstance(s, bytes):
                s = s.decode('ascii')

            try:
                return f.parse(s)
            except Exception as e:
                if parse_strict == 'silent':
                    pass
                else:
                    # Deliberately not isinstance here. Subclasses
                    # should use their name.
                    if type(f) is not unit_format.Generic:
                        format_clause = f.name + ' '
                    else:
                        format_clause = ''
                    msg = ("'{0}' did not parse as {1}unit: {2}"
                           .format(s, format_clause, six.text_type(e)))
                    if parse_strict == 'raise':
                        raise ValueError(msg)
                    elif parse_strict == 'warn':
                        warnings.warn(msg, UnitsWarning)
                    else:
                        raise ValueError("'parse_strict' must be 'warn', "
                                         "'raise' or 'silent'")
                return UnrecognizedUnit(s)

        elif isinstance(s, (int, float, np.floating, np.integer)):
            return CompositeUnit(s, [], [])

        elif s is None:
            raise TypeError("None is not a valid Unit")

        else:
            raise TypeError("{0} can not be converted to a Unit".format(s))


@six.add_metaclass(_UnitMetaClass)
class Unit(NamedUnit):
    """
    The main unit class.

    There are a number of different ways to construct a Unit, but
    always returns a `UnitBase` instance.  If the arguments refer to
    an already-existing unit, that existing unit instance is returned,
    rather than a new one.

    - From a string::

        Unit(s, format=None, parse_strict='silent')

      Construct from a string representing a (possibly compound) unit.

      The optional `format` keyword argument specifies the format the
      string is in, by default ``"generic"``.  For a description of
      the available formats, see `astropy.units.format`.

      The optional ``parse_strict`` keyword controls what happens when an
      unrecognized unit string is passed in.  It may be one of the following:

         - ``'raise'``: (default) raise a ValueError exception.

         - ``'warn'``: emit a Warning, and return an
           `UnrecognizedUnit` instance.

         - ``'silent'``: return an `UnrecognizedUnit` instance.

    - From a number::

        Unit(number)

      Creates a dimensionless unit.

    - From a `UnitBase` instance::

        Unit(unit)

      Returns the given unit unchanged.

    - From `None`::

        Unit()

      Returns the null unit.

    - The last form, which creates a new `Unit` is described in detail
      below.

    Parameters
    ----------
    st : str or list of str
        The name of the unit.  If a list, the first element is the
        canonical (short) name, and the rest of the elements are
        aliases.

    represents : UnitBase instance
        The unit that this named unit represents.

    doc : str, optional
        A docstring describing the unit.

    format : dict, optional
        A mapping to format-specific representations of this unit.
        For example, for the ``Ohm`` unit, it might be nice to have it
        displayed as ``\\Omega`` by the ``latex`` formatter.  In that
        case, `format` argument should be set to::

            {'latex': r'\\Omega'}

    namespace : dictionary, optional
        When provided, inject the unit (and all of its aliases) into
        the given namespace.

    Raises
    ------
    ValueError
        If any of the given unit names are already in the registry.

    ValueError
        If any of the given unit names are not valid Python tokens.
    """

    def __init__(self, st, represents=None, register=None, doc=None,
                 format=None, namespace=None):

        if register is not None:
            warnings.warn(
                "The registry kwarg was removed in astropy 0.3. "
                "Use the namespace kwarg to inject the unit into "
                "a namespace and add_enabled_units() to enable it "
                "in the global unit registry.",
                DeprecationWarning)

        represents = Unit(represents)
        self._represents = represents

        NamedUnit.__init__(self, st, namespace=namespace, doc=doc,
                           format=format)

    def decompose(self, bases=set()):
        return self._represents.decompose(bases=bases)

    def is_unity(self):
        return self._represents.is_unity()

    def __hash__(self):
        return hash(self.name) + hash(self._represents)


class PrefixUnit(Unit):
    """
    A unit that is simply a SI-prefixed version of another unit.

    For example, ``mm`` is a `PrefixUnit` of ``.001 * m``.

    The constructor is the same as for `Unit`.
    """


class CompositeUnit(UnitBase):
    """
    Create a composite unit using expressions of previously defined
    units.

    Direct use of this class is not recommended. Instead use the
    factory function `Unit` and arithmetic operators to compose
    units.

    Parameters
    ----------
    scale : number
        A scaling factor for the unit.

    bases : sequence of `UnitBase`
        A sequence of units this unit is composed of.

    powers : sequence of numbers
        A sequence of powers (in parallel with ``bases``) for each
        of the base units.
    """
    def __init__(self, scale, bases, powers, decompose=False,
                 decompose_bases=set(), _error_check=True):
        # There are many cases internal to astropy.units where we
        # already know that all the bases are Unit objects, and the
        # powers have been validated.  In those cases, we can skip the
        # error checking for performance reasons.  When the private
        # kwarg `_error_check` is False, the error checking is turned
        # off.
        if _error_check:
            scale = sanitize_scale(scale)
            for base in bases:
                if not isinstance(base, UnitBase):
                    raise TypeError(
                        "bases must be sequence of UnitBase instances")
            powers = [validate_power(p, support_tuples=True) for p in powers]

        self._scale = scale
        self._bases = bases
        self._powers = powers
        self._decomposed_cache = None
        self._expand_and_gather(decompose=decompose, bases=decompose_bases)
        self._hash = None

    def __repr__(self):
        if len(self._bases):
            return super(CompositeUnit, self).__repr__()
        else:
            if self._scale != 1.0:
                return 'Unit(dimensionless with a scale of {0})'.format(
                    self._scale)
            else:
                return 'Unit(dimensionless)'

    def __hash__(self):
        if self._hash is None:
            parts = ([str(self._scale)] +
                     [x.name for x in self._bases] +
                     [str(x) for x in self._powers])
            self._hash = hash(tuple(parts))
        return self._hash

    @property
    def scale(self):
        """
        Return the scale of the composite unit.
        """
        return self._scale

    @property
    def bases(self):
        """
        Return the bases of the composite unit.
        """
        return self._bases

    @property
    def powers(self):
        """
        Return the powers of the composite unit.
        """
        return self._powers

    def _expand_and_gather(self, decompose=False, bases=set()):
        def add_unit(unit, power, scale):
            if unit not in bases:
                for base in bases:
                    try:
                        scale *= unit._to(base) ** power
                    except UnitsError:
                        pass
                    except ValueError:
                        # on python2, sqrt(negative number) does not
                        # automatically lead to a complex number, but this is
                        # needed for the corner case of mag=-0.4*dex
                        scale *= cmath.exp(power * cmath.log(unit._to(base)))
                        unit = base
                        break
                    else:
                        unit = base
                        break

            if unit in new_parts:
                new_parts[unit] = add_powers(new_parts[unit], power)
            else:
                new_parts[unit] = power
            return scale

        new_parts = {}
        scale = self.scale

        for b, p in zip(self.bases, self.powers):
            if decompose and b not in bases:
                b = b.decompose(bases=bases)

            if isinstance(b, CompositeUnit):
                try:
                    scale *= b._scale ** p
                except ValueError:
                    # on python2, sqrt(negative number) does not
                    # automatically lead to a complex number, but this is
                    # needed for the corner case of mag=-0.4*dex
                    scale *= cmath.exp(p * cmath.log(b._scale))
                for b_sub, p_sub in zip(b._bases, b._powers):
                    scale = add_unit(b_sub, p_sub * p, scale)
            else:
                scale = add_unit(b, p, scale)

        new_parts = [x for x in six.iteritems(new_parts) if x[1] != 0]
        new_parts.sort(key=lambda x: (-x[1], getattr(x[0], 'name', '')))

        self._bases = [x[0] for x in new_parts]
        self._powers = [validate_power(x[1], support_tuples=True)
                        for x in new_parts]
        self._scale = sanitize_scale(scale)

    def __copy__(self):
        """
        For compatibility with python copy module.
        """
        return CompositeUnit(self._scale, self._bases[:], self._powers[:])

    def decompose(self, bases=set()):
        if len(bases) == 0 and self._decomposed_cache is not None:
            return self._decomposed_cache

        for base in self.bases:
            if (not isinstance(base, IrreducibleUnit) or
                    (len(bases) and base not in bases)):
                break
        else:
            if len(bases) == 0:
                self._decomposed_cache = self
            return self

        x = CompositeUnit(self.scale, self.bases, self.powers, decompose=True,
                          decompose_bases=bases)
        if len(bases) == 0:
            self._decomposed_cache = x
        return x

    def is_unity(self):
        unit = self.decompose()
        return len(unit.bases) == 0 and unit.scale == 1.0


si_prefixes = [
    (['Y'], ['yotta'], 1e24),
    (['Z'], ['zetta'], 1e21),
    (['E'], ['exa'], 1e18),
    (['P'], ['peta'], 1e15),
    (['T'], ['tera'], 1e12),
    (['G'], ['giga'], 1e9),
    (['M'], ['mega'], 1e6),
    (['k'], ['kilo'], 1e3),
    (['h'], ['hecto'], 1e2),
    (['da'], ['deka', 'deca'], 1e1),
    (['d'], ['deci'], 1e-1),
    (['c'], ['centi'], 1e-2),
    (['m'], ['milli'], 1e-3),
    (['u'], ['micro'], 1e-6),
    (['n'], ['nano'], 1e-9),
    (['p'], ['pico'], 1e-12),
    (['f'], ['femto'], 1e-15),
    (['a'], ['atto'], 1e-18),
    (['z'], ['zepto'], 1e-21),
    (['y'], ['yocto'], 1e-24)
]


binary_prefixes = [
    (['Ki'], ['kibi'], 2. ** 10),
    (['Mi'], ['mebi'], 2. ** 20),
    (['Gi'], ['gibi'], 2. ** 30),
    (['Ti'], ['tebi'], 2. ** 40),
    (['Pi'], ['pebi'], 2. ** 50),
    (['Ei'], ['exbi'], 2. ** 60)
]


def _add_prefixes(u, excludes=[], namespace=None, prefixes=False):
    """
    Set up all of the standard metric prefixes for a unit.  This
    function should not be used directly, but instead use the
    `prefixes` kwarg on `def_unit`.

    Parameters
    ----------
    excludes : list of str, optional
        Any prefixes to exclude from creation to avoid namespace
        collisions.

    namespace : dict, optional
        When provided, inject the unit (and all of its aliases) into
        the given namespace dictionary.

    prefixes : list, optional
        When provided, it is a list of prefix definitions of the form:

            (short_names, long_tables, factor)
    """
    if prefixes is True:
        prefixes = si_prefixes
    elif prefixes is False:
        prefixes = []

    for short, full, factor in prefixes:
        names = []
        format = {}
        for prefix in short:
            if prefix in excludes:
                continue

            for alias in u.short_names:
                names.append(prefix + alias)

                # This is a hack to use Greek mu as a prefix
                # for some formatters.
                if prefix == 'u':
                    format['latex'] = r'\mu ' + u.get_format_name('latex')
                    format['unicode'] = 'μ' + u.get_format_name('unicode')

                for key, val in six.iteritems(u._format):
                    format.setdefault(key, prefix + val)

        for prefix in full:
            if prefix in excludes:
                continue

            for alias in u.long_names:
                names.append(prefix + alias)

        if len(names):
            PrefixUnit(names, CompositeUnit(factor, [u], [1],
                                            _error_check=False),
                       namespace=namespace, format=format)


def def_unit(s, represents=None, register=None, doc=None,
             format=None, prefixes=False, exclude_prefixes=[],
             namespace=None):
    """
    Factory function for defining new units.

    Parameters
    ----------
    s : str or list of str
        The name of the unit.  If a list, the first element is the
        canonical (short) name, and the rest of the elements are
        aliases.

    represents : UnitBase instance, optional
        The unit that this named unit represents.  If not provided,
        a new `IrreducibleUnit` is created.

    doc : str, optional
        A docstring describing the unit.

    format : dict, optional
        A mapping to format-specific representations of this unit.
        For example, for the ``Ohm`` unit, it might be nice to
        have it displayed as ``\\Omega`` by the ``latex``
        formatter.  In that case, `format` argument should be set
        to::

            {'latex': r'\\Omega'}

    prefixes : bool or list, optional
        When `True`, generate all of the SI prefixed versions of the
        unit as well.  For example, for a given unit ``m``, will
        generate ``mm``, ``cm``, ``km``, etc.  When a list, it is a list of
        prefix definitions of the form:

            (short_names, long_tables, factor)

        Default is `False`.  This function always returns the base
        unit object, even if multiple scaled versions of the unit were
        created.

    exclude_prefixes : list of str, optional
        If any of the SI prefixes need to be excluded, they may be
        listed here.  For example, ``Pa`` can be interpreted either as
        "petaannum" or "Pascal".  Therefore, when defining the
        prefixes for ``a``, ``exclude_prefixes`` should be set to
        ``["P"]``.

    namespace : dict, optional
        When provided, inject the unit (and all of its aliases and
        prefixes), into the given namespace dictionary.

    Returns
    -------
    unit : `UnitBase` object
        The newly-defined unit, or a matching unit that was already
        defined.
    """
    if register is not None:
        warnings.warn(
            "The register kwarg was removed in astropy 0.3. "
            "Use the namespace kwarg to inject the unit into "
            "a namespace and add_enabled_units() to enable it "
            "in the global unit registry.",
            DeprecationWarning)

    if represents is not None:
        result = Unit(s, represents, namespace=namespace, doc=doc,
                      format=format)
    else:
        result = IrreducibleUnit(
            s, namespace=namespace, doc=doc, format=format)

    if prefixes:
        _add_prefixes(result, excludes=exclude_prefixes, namespace=namespace,
                      prefixes=prefixes)
    return result


def _condition_arg(value):
    """
    Validate value is acceptable for conversion purposes.

    Will convert into an array if not a scalar, and can be converted
    into an array

    Parameters
    ----------
    value : int or float value, or sequence of such values

    Returns
    -------
    Scalar value or numpy array

    Raises
    ------
    ValueError
        If value is not as expected
    """
    if isinstance(value, (float, six.integer_types, complex)):
        return value

    if isinstance(value, np.ndarray) and value.dtype.kind in ['i', 'f', 'c']:
        return value

    try:
        avalue = np.array(value)
        assert avalue.dtype.kind in ['i', 'f', 'c']
        return avalue
    except (ValueError, AssertionError):
        raise ValueError("Value not scalar compatible or convertible to "
                         "an int, float, or complex array")


dimensionless_unscaled = CompositeUnit(1, [], [], _error_check=False)
# Abbreviation of the above, see #1980
one = dimensionless_unscaled
