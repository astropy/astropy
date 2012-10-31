# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core units classes and functions
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from fractions import Fraction
import re
import sys
import textwrap
import warnings
import numbers

import numpy as np
from numpy import ma

from . import format as unit_format

# TODO: Support functional units, e.g. log(x), ln(x)

__all__ = [
    'UnitsException', 'UnitsWarning', 'UnitBase', 'NamedUnit',
    'IrreducibleUnit', 'Unit', 'def_unit', 'CompositeUnit',
    'PrefixUnit', 'UnrecognizedUnit']


class UnitsException(Exception):
    """
    The base class for unit-specific exceptions.
    """
    pass


class UnitsWarning(Warning):
    """
    The base class for unit-specific exceptions.
    """
    pass


class UnitBase(object):
    """
    Abstract base class for units.

    Most of the arithmetic operations on units are defined in this
    base class.

    Should not be instantiated by users directly.
    """
    _registry = []
    _namespace = {}

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

    def __str__(self):
        """Return string representation for unit"""
        return unit_format.Generic().to_string(self)

    def __repr__(self):
        return 'Unit("' + str(self) + '")'

    def to_string(self, format='generic'):
        """
        Output the unit in the given format as a string.

        Parameters
        ----------
        format : `astropy.format.Base` instance or str
            The name of a format or a formatter object.  If not
            provided, defaults to the generic format.
        """
        f = unit_format.get_format(format)
        return f.to_string(self)

    @staticmethod
    def _get_namespace():
        """
        Get the namespace that units will be registered to.
        """
        return UnitBase._namespace

    @staticmethod
    def _set_namespace(d):
        """
        Set the namespace that units will be registered to.
        """
        UnitBase._namespace = d

    def __pow__(self, p):
        if isinstance(p, tuple) and len(p) == 2:
            p = Fraction(p[0], p[1])
        else:
            # allow two possible floating point fractions, all others illegal
            if not int(2 * p) == 2 * p:
                raise ValueError(
                    "floating values for unit powers must be integers or "
                    "integers + 0.5")
        return CompositeUnit(1, [self], [p]).simplify()

    def __div__(self, m):
        # Strictly speaking, we should be using old-style division here.
        # However, I think it's less surprising for this to behave the
        # same way whether __future__ division is being used or not
        if isinstance(m, UnitBase):
            return CompositeUnit(1, [self, m], [1, -1]).simplify()
        elif isinstance(m, numbers.Number):
            from .quantity import Quantity
            return Quantity(1./m, self)
        else:
            return CompositeUnit(1.0 / m, [self], [1]).simplify()

    def __rdiv__(self, m):
        if isinstance(m, numbers.Number):
            from .quantity import Quantity
            return Quantity(m, CompositeUnit(1.0, [self], [-1]).simplify())
        else:
            return CompositeUnit(m, [self], [-1]).simplify()

    def __truediv__(self, m):
        return self.__div__(m)

    def __rtruediv__(self, m):
        return self.__rdiv__(m)

    def __mul__(self, m):
        if isinstance(m, UnitBase):
            return CompositeUnit(1, [self, m], [1, 1]).simplify()
        elif isinstance(m, numbers.Number):
            from .quantity import Quantity
            return Quantity(m, self)
        else:
            return CompositeUnit(m, [self], [1]).simplify()

    def __rmul__(self, m):
        if isinstance(m, numbers.Number):
            from .quantity import Quantity
            return Quantity(m, self)
        else:
            return CompositeUnit(m, [self], [1]).simplify()

    if sys.version_info[0] >= 3:
        def __hash__(self):
            # Since this class defines __eq__, it will become unhashable
            # on Python 3.x, so we need to define our own hash.
            return id(self)

    def __eq__(self, other):
        other = Unit(other, parse_strict='silent')
        try:
            return np.allclose(self.to(other, 1), 1.0)
        except UnitsException:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        other = Unit(other)
        return self.to(other, 1) < 1.

    def __gt__(self, other):
        other = Unit(other)
        return self.to(other, 1) > 1.

    def __le__(self, other):
        other = Unit(other)
        return self.to(other, 1) <= 1.

    def __ge__(self, other):
        other = Unit(other)
        return self.to(other, 1) >= 1.

    def __neg__(self):
        return self * -1.

    def simplify(self):
        """
        Compresses a possibly composite unit down to a single
        instance.
        """
        return self

    def is_dimensionless(self):
        """
        Returns `True` if this unit translates into a scalar quantity
        without a unit.

        Examples
        --------
        >>> ((2 * u.m) / (3 * u.m)).is_dimensionless()
        True
        >>> (2 * u.m).is_dimensionless()
        False
        """
        return False

    def is_equivalent(self, other, equivalencies=[]):
        """
        Returns `True` if this unit is equivalent to `other`.

        Parameters
        ----------
        other : unit object or string
           The unit to convert to.

        equivalencies : list of equivalence pairs, optional
           A list of equivalence pairs to try if the units are not
           directly convertible.  See :ref:`unit_equivalencies`.

        Returns
        -------
        bool
        """
        other = Unit(other, parse_strict='silent')

        if isinstance(other, UnrecognizedUnit):
            return False

        try:
            (self / other).dimensionless_constant()
        except UnitsException:
            unit = self.decompose()
            other = other.decompose()
            for equiv in equivalencies:
                a = equiv[0]
                b = equiv[1]
                if (unit.is_equivalent(a) and
                    other.is_equivalent(b)):
                    return True
                elif (unit.is_equivalent(b) and
                      other.is_equivalent(a)):
                    return True
            return False
        else:
            return True

    def _apply_equivalences(self, unit, other, equivalencies):
        """
        Internal function (used from `get_converter`) to apply
        equivalence pairs.
        """
        def make_converter(scale1, func, scale2):
            def convert(v):
                return func(_condition_arg(v) * scale1) * scale2
            return convert

        orig_unit = unit
        orig_other = other

        unit = self.decompose()
        other = other.decompose()

        for equiv in equivalencies:
            if len(equiv) == 2:
                funit, tunit = equiv
                a, b = lambda x: x
            if len(equiv) == 3:
                funit, tunit, a = equiv
                b = a
            elif len(equiv) == 4:
                funit, tunit, a, b = equiv
            else:
                raise ValueError("Invalid equivalence entry")
            if (unit.is_equivalent(funit) and
                other.is_equivalent(tunit)):
                scale1 = (unit / funit).dimensionless_constant()
                scale2 = (tunit / other).dimensionless_constant()
                return make_converter(scale1, a, scale2)
            elif (other.is_equivalent(funit) and
                  unit.is_equivalent(tunit)):
                scale1 = (unit / tunit).dimensionless_constant()
                scale2 = (funit / other).dimensionless_constant()
                return make_converter(scale1, b, scale2)

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

        raise UnitsException(
            "{0} and {1} are not convertible".format(
                unit_str, other_str))

    def get_converter(self, other, equivalencies=[]):
        """
        Return the conversion function to convert values from `self`
        to the specified unit.

        Parameters
        ----------
        other : unit object or string
           The unit to convert to.

        equivalencies : list of equivalence pairs, optional
           A list of equivalence pairs to try if the units are not
           directly convertible.  See :ref:`unit_equivalencies`.

        Returns
        -------
        func : callable
            A callable that normally expects a single argument that is
            a scalar value or an array of values (or anything that may
            be converted to an array).

        Raises
        ------
        UnitsException
            If units are inconsistent
        """
        other = Unit(other)

        try:
            scale = (self / other).dimensionless_constant()
        except UnitsException:
            return self._apply_equivalences(
                self, other, equivalencies)
        return lambda val: scale * _condition_arg(val)

    def to(self, other, value=1.0, equivalencies=[]):
        """
        Return the converted values in the specified unit.

        Parameters
        ----------
        other : unit object or string
            The unit to convert to.

        value : scalar int or float, or sequence that can be converted to array, optional
            Value(s) in the current unit to be converted to the
            specified unit.  If not provided, defaults to 1.0

        equivalencies : list of equivalence pairs, optional
           A list of equivalence pairs to try if the units are not
           directly convertible.  See :ref:`unit_equivalencies`.

        Returns
        -------
        values : scalar or array
            Converted value(s). Input value sequences are returned as
            numpy arrays.

        Raises
        ------
        UnitException
            If units are inconsistent
        """
        other = Unit(other)
        return self.get_converter(
            other, equivalencies=equivalencies)(value)

    def in_units(self, other, value=1.0, equivalencies=[]):
        """
        Alias for `to` for backward compatibility with pynbody.
        """
        return self.to(
            other, value=value, equivalencies=equivalencies)

    def decompose(self):
        """
        Return a unit object composed of only irreducible units.

        Parameters
        ----------
        None

        Returns
        -------
        unit : CompositeUnit object
            New object containing only irreducible unit objects.
        """
        return self

    @property
    def physical_type(self):
        """
        Return the physical type on the unit.

        Examples
        --------
        >>> u.m.physical_type
        'length'
        """
        from . import physical
        return physical.get_physical_type(self)

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

    def find_equivalent_units(self, equivalencies=[]):
        """
        Return a list of all the units that are the same type as the
        specified unit.

        Parameters
        ----------
        u : Unit instance or string
            The `Unit` to find similar units to.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to also list.  See
            :ref:`unit_equivalencies`.

        Returns
        -------
        units : list of `UnitBase`
            A list of unit objects that match `u`.  A subclass of
            `list` (`EquivalentUnitsList`) is returned that
            pretty-prints the list of units when output.
        """
        units = [self]
        for equiv in equivalencies:
            funit, tunit = equiv[:2]
            if self.is_equivalent(funit):
                units.append(tunit)
            elif self.is_equivalent(tunit):
                units.append(funit)

        equivalencies = set()
        for tunit in UnitBase._registry:
            if not isinstance(tunit, PrefixUnit):
                for u in units:
                    try:
                        tunit.get_converter(u)
                    except UnitsException:
                        pass
                    else:
                        equivalencies.add(tunit)

        return self.EquivalentUnitsList(equivalencies)


class NamedUnit(UnitBase):
    """
    The base class of units that have a name.

    Parameters
    ----------
    st : str or list of str
        The name of the unit.  If a list, the first element is the
        canonical (short) name, and the rest of the elements are
        aliases.

    register : boolean, optional
        When `True`, also register the unit in the standard unit
        namespace.  Default is `False`.

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
    def __init__(self, st, register=False, doc=None, format=None):
        UnitBase.__init__(self)

        if isinstance(st, (bytes, unicode)):
            self._names = [st]
        else:
            if len(st) == 0:
                raise ValueError(
                    "st list must have at least one entry")
            self._names = st[:]

        if format is None:
            format = {}
        self._format = format

        if doc is None:
            doc = self._generate_doc()

        doc = textwrap.dedent(doc)
        doc = textwrap.fill(doc)

        self.__doc__ = doc

        self._register_unit(register)

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

    def _register_unit(self, register):
        """
        Registers the unit in the registry, and optionally in another
        namespace.  It is registered under all of the names and
        aliases given to the constructor.

        The namespace used is set with `UnitBase._set_namespace`.

        Parameters
        ----------
        register : bool
            When `True`, register the unit in the external namespace
            as well as the central registry.
        """
        if not self._names:
            raise UnitsException("unit has no string representation")

        for st in self._names:
            if not re.match("^[A-Za-z_]+$", st):
                # will cause problems for simple string parser in
                # unit() factory
                raise ValueError(
                    "Invalid unit name {0!r}".format(st))

            if register:
                if st in self._namespace:
                    raise ValueError(
                        "Object with name {0!r} already exists "
                        "in namespace".format(st))
                self._namespace[st] = self

        self._registry.append(self)


class IrreducibleUnit(NamedUnit):
    """
    Irreducible units are the units that all other units are defined
    in terms of.

    Examples are meters, seconds, kilograms, amperes, etc.  There is
    only once instance of such a unit per type.
    """
    def decompose(self):
        return CompositeUnit(1, [self], [1])
    decompose.__doc__ = UnitBase.decompose.__doc__


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
    def __init__(self, st):
        IrreducibleUnit.__init__(self, st)

    def __repr__(self):
        return "UnrecognizedUnit({0})".format(str(self))

    def __str__(self):
        return self.name

    def to_string(self, format='generic'):
        return self.name

    def _register_unit(self, register):
        pass

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

    def is_equivalent(self, other, equivalencies=[]):
        return self == other

    def get_converter(self, other, equivalencies=[]):
        raise ValueError(
            "The unit {0!r} is unrecognized.  It can not be converted "
            "to other units.".format(self.name))

    def get_format_name(self, format):
        return self.name


class _UnitMetaClass(type):
    """
    This metaclass exists because the Unit constructor should
    sometimes return instances that already exist.  This "overrides"
    the constructor before the new instance is actually created, so we
    can return an existing one.
    """
    def __call__(self, s, represents=None, format=None, register=False,
                 doc=None, parse_strict='raise'):


        from .quantity import Quantity
        if isinstance(represents, Quantity):
            if represents.value == 1:
                represents = represents.unit
            elif isinstance(represents.unit, CompositeUnit):
                represents = CompositeUnit(represents.value, bases=represents.unit.bases, powers=represents.unit.powers)
            else:
                represents = CompositeUnit(represents.value, bases=[represents.unit], powers=[1])

        if isinstance(s, Quantity):
            if s.value == 1:
                s = s.unit
            elif isinstance(s.unit, CompositeUnit):
                s = CompositeUnit(s.value*s.unit.scale, bases=s.unit.bases, powers=s.unit.powers)
            else:
                s = CompositeUnit(s.value, bases=[s.unit], powers=[1])

        if isinstance(represents, UnitBase):
            # This has the effect of calling the real __new__ and
            # __init__ on the Unit class.
            return super(_UnitMetaClass, self).__call__(
                s, represents, format=format, register=register, doc=doc)
            raise TypeError("Can not convert {0!r} to a unit".format(s))

        elif isinstance(s, UnitBase):
            return s

        elif isinstance(s, (bytes, unicode)):
            if len(s.strip()) == 0:
                # Return the NULL unit
                return CompositeUnit(1.0, [], [])

            if format is None:
                format = 'generic'

            f = unit_format.get_format(format)
            try:
                return f.parse(s)
            except ValueError as e:
                if parse_strict == 'raise':
                    raise
                elif parse_strict == 'warn':
                    warnings.warn(
                        "'{0}' did not parse using format '{1}'. {2}".format(
                            s, format, str(e)),
                            UnitsWarning)
                elif parse_strict != 'silent':
                    raise ValueError(
                        "'parse_strict' must be 'warn', 'raise' or 'silent'")
                return UnrecognizedUnit(s)

        elif isinstance(s, (int, float, np.floating, np.integer)):
            return CompositeUnit(s, [], [])

        elif s is None:
            raise ValueError("None is not a valid Unit")

        else:
            raise TypeError("{0} can not be converted to a Unit".format(s))


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

      The optional `parse_strict` keyword controls what happens when an
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

    register : boolean, optional
        When `True`, also register the unit in the standard unit
        namespace.  Default is `False`.

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
    __metaclass__ = _UnitMetaClass

    def __init__(self, st, represents=None, register=False, doc=None,
                 format=None):

        represents = Unit(represents)
        self._represents = represents

        NamedUnit.__init__(self, st, register=register, doc=doc,
                           format=format)

    def decompose(self):
        return self._represents.decompose()
    decompose.__doc__ = UnitBase.decompose.__doc__


class PrefixUnit(Unit):
    """
    A unit that is simply a SI-prefixed version of another unit.

    For example, `mm` is a `PrefixUnit` of ``.001 * m``.

    The constructor is the same as for `Unit`.
    """
    pass


class CompositeUnit(UnitBase):
    """
    Create a composite unit using expressions of previously defined
    units.

    Direct use of this class is not recommended. Instead use the
    factory function `Unit(...)` and arithmetic operators to compose
    units.

    Parameters
    ----------
    scale : number
        A scaling factor for the unit.

    bases : sequence of `UnitBase`
        A sequence of units this unit is composed of.

    powers : sequence of numbers
        A sequence of powers (in parallel with `bases`) for each
        of the base units.
    """
    def __init__(self, scale, bases, powers):
        if scale == 1.:
            scale = 1
        self._scale = scale
        for base in bases:
            if not isinstance(base, UnitBase):
                raise TypeError("bases must be sequence of UnitBase instances")
        self._bases = bases
        self._powers = powers

    def __repr__(self):
        if len(self._bases):
            return super(CompositeUnit, self).__repr__()
        else:
            return 'Unit(dimensionless)'

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

    def _expand_and_gather(self, decompose=False):
        bases = {}
        scale = self.scale

        for i, (b, p) in enumerate(zip(self.bases, self.powers)):
            if decompose:
                b = b.decompose()

            if isinstance(b, CompositeUnit):
                scale *= b.scale ** p
                for b_sub, p_sub in zip(b.bases, b.powers):
                    bases[b_sub] = p_sub * p + bases.get(b_sub, 0)

            else:
                bases[b] = p + bases.get(b, 0)

        bases = [(b, p) for (b, p) in bases.items() if p != 0]
        bases.sort(key=lambda x: x[1], reverse=True)

        self._bases = [x[0] for x in bases]
        self._powers = [x[1] for x in bases]
        self._scale = scale

    def __copy__(self):
        """
        For compatibility with python copy module.
        """
        return CompositeUnit(self._scale, self._bases[:], self._powers[:])

    def simplify(self):
        self._expand_and_gather()
        return self
    simplify.__doc__ = UnitBase.simplify.__doc__

    def decompose(self):
        x = CompositeUnit(self.scale, self.bases, self.powers)
        x._expand_and_gather(True)
        return x
    decompose.__doc__ = UnitBase.decompose.__doc__

    def is_dimensionless(self):
        x = self.decompose()
        return (len(x.powers) == 0)
    is_dimensionless.__doc__ = UnitBase.is_dimensionless.__doc__

    def dimensionless_constant(self):
        """
        If this unit is dimensionless, return its scalar quantity.

        Direct use of this method is not recommended. It is generally
        better to use the `to` or `get_converter` methods
        instead.
        """
        x = self.decompose()
        c = x.scale
        if len(x.bases):
            raise UnitsException(
                "'{0}' is not dimensionless".format(self.to_string()))
        return c


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


def _add_prefixes(u, excludes=[], register=False):
    """
    Set up all of the standard metric prefixes for a unit.  This
    function should not be used directly, but instead use the
    `prefixes` kwarg on `def_unit`.

    Parameters
    ----------
    excludes : list of str, optional
        Any prefixes to exclude from creation to avoid namespace
        collisions.

    register : bool, optional
        When `True`, also register the unit in the standard unit
        namespace.  Default is `False`.
    """
    for short, long, factor in si_prefixes:
        exclude = False
        for prefix in short:
            if prefix in excludes:
                exclude = True
        if exclude:
            continue

        names = []
        format = {}
        for prefix in short:
            for alias in [u.name] + [x for x in u.aliases if len(x) <= 2]:
                names.append(prefix + alias)

                # This is a hack to use Greek mu as a prefix
                # for some formatters.
                if prefix == 'u':
                    format['latex'] = r'\mu ' + u.get_format_name('latex')
                    format['unicode'] = 'μ' + u.get_format_name('unicode')

                for key, val in u._format.items():
                    format.setdefault(key, prefix + val)

        for prefix in long:
            for alias in u.aliases:
                if len(alias) > 2:
                    names.append(prefix + alias)

        PrefixUnit(names, CompositeUnit(factor, [u], [1]), register=register, format=format)

def def_unit(s, represents=None, register=None, doc=None,
             format=None, prefixes=False, exclude_prefixes=[]):
    """
    Factory function for defining new units.

    Parameters
    ----------
    names : str or list of str
        The name of the unit.  If a list, the first element is the
        canonical (short) name, and the rest of the elements are
        aliases.

    represents : UnitBase instance, optional
        The unit that this named unit represents.  If not provided,
        a new `IrreducibleUnit` is created.

    register : boolean, optional
        When `True`, also register the unit in the standard unit
        namespace.  Default is `False`.

    doc : str, optional
        A docstring describing the unit.

    format : dict, optional
        A mapping to format-specific representations of this unit.
        For example, for the ``Ohm`` unit, it might be nice to
        have it displayed as ``\\Omega`` by the ``latex``
        formatter.  In that case, `format` argument should be set
        to::

            {'latex': r'\\Omega'}

    prefixes : bool, optional
        When `True`, generate all of the SI prefixed versions of the
        unit as well.  For example, for a given unit `m`, will generate
        `mm`, `cm`, `km`, etc.  Default is `False`.  This function
        always returns the base unit object, even if multiple scaled
        versions of the unit were created.

    exclude_prefixes : list of str, optional
        If any of the SI prefixes need to be excluded, they may be
        listed here.  For example, `Pa` can be interpreted either as
        "petaannum" or "Pascal".  Therefore, when defining the
        prefixes for `a`, `exclude_prefixes` should be set to
        ``["P"]``.

    Returns
    -------
    unit : `UnitBase` object
        The newly-defined unit, or a matching unit that was already
        defined.
    """

    if register is None:
        register = False
    if represents is not None:
        result = Unit(s, represents, register=register, doc=doc,
                      format=format)
    else:
        result = IrreducibleUnit(s, register=register, doc=doc, format=format)

    if prefixes:
        _add_prefixes(result, excludes=exclude_prefixes, register=register)
    return result


def _condition_arg(value):
    """
    Validate value is acceptable for conversion purposes.

    Will convert into an array if not a scalar, and can be converted
    into an array

    Parameters
    ----------
    value: int or float value, or sequence of such values
        that can be converted into an array if not scalar

    Returns
    -------
    Scalar value or numpy array

    Raises
    ------
    ValueError
        If value is not as expected
    """
    if isinstance(value, float) or isinstance(value, int):
        return value
    else:
        try:
            avalue = np.array(value)
            dt = str(avalue.dtype)
            if not (dt.startswith('int') or dt.startswith('float')):
                raise ValueError("Must be convertable to int or float array")
            if ma.isMaskedArray(value):
                return value
            return avalue
        except ValueError:
            raise ValueError(
                "Value not scalar compatible or convertable into a float or "
                "integer array")
