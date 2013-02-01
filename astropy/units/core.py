# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Core units classes and functions
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import inspect
import re
import sys
import textwrap
import warnings

import numpy as np
from numpy import ma

from ..utils.compat.fractions import Fraction
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
        return CompositeUnit(1, [self], [p])._simplify()

    def __div__(self, m):
        # Strictly speaking, we should be using old-style division here.
        # However, I think it's less surprising for this to behave the
        # same way whether __future__ division is being used or not
        from .quantity import Quantity
        if isinstance(m, UnitBase):
            return CompositeUnit(1, [self, m], [1, -1])._simplify()
        elif isinstance(m, Quantity):
            return Quantity(1, self) / m
        else:
            return Quantity(1. / m, self)

    def __rdiv__(self, m):
        from .quantity import Quantity
        return Quantity(m, CompositeUnit(1.0, [self], [-1])._simplify())

    def __truediv__(self, m):
        return self.__div__(m)

    def __rtruediv__(self, m):
        return self.__rdiv__(m)

    def __mul__(self, m):
        from .quantity import Quantity
        if isinstance(m, UnitBase):
            return CompositeUnit(1, [self, m], [1, 1])._simplify()
        elif isinstance(m, Quantity):
            return Quantity(1, self) * m
        else:
            return Quantity(m, self)

    def __rmul__(self, m):
        from .quantity import Quantity
        return Quantity(m, self)

    if sys.version_info[0] >= 3:
        def __hash__(self):
            # Since this class defines __eq__, it will become unhashable
            # on Python 3.x, so we need to define our own hash.
            return id(self)

    def __eq__(self, other):
        try:
            other = Unit(other, parse_strict='silent')
        except (ValueError, UnitsException):
            return False
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

    def _simplify(self):
        """
        Compresses a possibly composite unit down to a single
        instance.
        """
        return self

    def is_unity(self):
        """
        Returns `True` if this unit translates into a scalar quantity
        without a unit or a scale.

        Examples
        --------
        >>> ((2 * u.m) / (3 * u.m)).is_dimensionless()
        True
        >>> ((2 * u.m) / (3 * u.m)).is_unity()
        False
        >>> ((2 * u.m) / (2 * u.m)).is_unity()
        True
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

    def decompose(self, bases=[]):
        """
        Return a unit object composed of only irreducible units.

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
        unit : CompositeUnit object
            New object containing only irreducible unit objects.
        """
        raise NotImplementedError()

    def _compose(self, equivalencies=[], namespace=[], max_depth=2, depth=0):
        def is_final_result(unit):
            # Returns True if this result contains only the expected
            # units
            for base in unit._bases:
                if base not in namespace:
                    return False
            return True

        def sort_results(results):
            if not len(results):
                return []

            # Sort the results so the simplest ones appear first.
            # Simplest is defined as "the minimum sum of absolute
            # powers" (i.e. the fewest bases), and preference should
            # be given to results where the sum of powers is positive
            # and the scale is exactly equal to 1.0
            results = list(results)
            results.sort(key=lambda x: str(x))
            results.sort(key=lambda x: np.sum(np.abs(x.powers)))
            results.sort(key=lambda x: np.sum(x.powers) < 0.0)
            results.sort(key=lambda x: not np.allclose(x.scale, 1.0))

            last_result = results[0]
            filtered = [last_result]
            for result in results[1:]:
                if str(result) != str(last_result):
                    filtered.append(result)
                last_result = result

            return filtered

        unit = self.decompose()

        # Prevent too many levels of recursion
        if depth >= max_depth:
            return [unit]

        # Special case for dimensionless unit
        if len(unit._bases) == 0:
            return [unit]

        # Make a list including all of the equivalent units
        units = [unit]
        for equiv in equivalencies:
            funit, tunit = equiv[:2]
            if self.is_equivalent(funit):
                units.append(tunit.decompose())
            elif self.is_equivalent(tunit):
                units.append(funit.decompose())

        # Store partial results
        partial_results = []
        # Store final results that reduce to a single unit or pair of
        # units
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
                    for base, power in zip(u._bases, u._powers):
                        if tunit_decomposed.is_equivalent(base):
                            tunit = tunit ** power
                            tunit_decomposed = tunit_decomposed ** power
                            break

                composed = (u / tunit_decomposed).decompose()
                factored = composed * tunit
                len_bases = len(composed._bases)
                if is_final_result(factored) and len_bases <= 1:
                    final_results[len_bases].add(factored)
                else:
                    partial_results.append(
                        (len_bases, composed, tunit))

        # Do we have any minimal results?
        for final_result in final_results:
            if len(final_result):
                return sort_results(final_result)

        partial_results.sort(key=lambda x: x[0])

        # ...we have to recurse and try to further compose
        results = []
        for len_bases, composed, tunit in partial_results:
            try:
                composed_list = composed._compose(
                    equivalencies=equivalencies, namespace=namespace,
                    max_depth=max_depth, depth=depth + 1)
            except UnitsException:
                composed_list = []
            for subcomposed in composed_list:
                results.append(
                    (len(subcomposed._bases), subcomposed, tunit))

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
                        subresults.add(composed * tunit)

            if len(subresults):
                return sort_results(subresults)

        for base in self.bases:
            if base not in namespace:
                raise UnitsException(
                    "Cannot represent unit {0} in terms of the given "
                    "units".format(self))

        return [self]

    def compose(self, equivalencies=[], units=None, max_depth=2):
        """
        Return the simplest possible composite unit(s) that represent
        the given unit.  Since there may be multiple equally simple
        compositions of the unit, a list of units is always returned.

        Parameters
        ----------
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to also list.  See
            :ref:`unit_equivalencies`.

        units : set of units to compose to, optional
            If not provided, any known units may be used to compose
            into.  Otherwise, ``units`` is a dict, module or sequence
            containing the units to compose into.

        max_depth : int, optional
            The maximum recursion depth to use when composing into
            composite units.

        Returns
        -------
        units : list of `CompositeUnit`
            A list of candidate compositions.  These will all be
            equally simple, but it may not be possible to
            automatically determine which of the candidates are
            better.
        """
        def filter_units(units):
            filtered_namespace = set()
            for tunit in units:
                if (isinstance(tunit, UnitBase) and
                    not isinstance(tunit, PrefixUnit)):
                    filtered_namespace.add(tunit)
            return filtered_namespace

        if units is None:
            units = filter_units(UnitBase._registry)
        elif isinstance(units, dict):
            units = set(units.values())
        elif inspect.ismodule(units):
            units = filter_units(units.__dict__.values())
        else:
            units = set(units)

        return self._compose(
            equivalencies=equivalencies, namespace=units,
            max_depth=max_depth, depth=0)

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
            sequence member named `bases` containing the base units of
            the system.

        Returns
        -------
        units : list of `CompositeUnit`
            The list is ranked so that units containing only the base
            units of that system will appear first.
        """
        bases = set(system.bases)

        def score(compose):
            sum = 0
            for base in compose._bases:
                if base in bases:
                    sum += 1
            return sum / float(len(compose._bases))

        x = self.decompose(bases=bases)
        composed = x.compose(units=system)
        composed = sorted(composed, key=score, reverse=True)
        return composed

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

    def find_equivalent_units(self, equivalencies=[], units=None):
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

        units : set of units to search in, optional
            If not provided, all defined units will be searched for
            equivalencies.  Otherwise, may be a dict, module or
            sequence containing the units to search for equivalencies.

        Returns
        -------
        units : list of `UnitBase`
            A list of unit objects that match `u`.  A subclass of
            `list` (`EquivalentUnitsList`) is returned that
            pretty-prints the list of units when output.
        """
        results = self.compose(
            equivalencies=equivalencies, units=units, max_depth=1)
        results = [
            x._bases[0] for x in results if len(x._bases) == 1]
        return self.EquivalentUnitsList(results)


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
        return [1.0]

    def decompose(self, bases=[]):
        bases = set(bases)

        if len(bases) and not self in bases:
            for base in bases:
                if self.is_equivalent(base):
                    return CompositeUnit(self.to(base), [base], [1])

            raise UnitsException(
                "Unit {0} can not be decomposed into the requested "
                "bases".format(self))

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

    def decompose(self, bases=[]):
        return self._represents.decompose(bases=bases)
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
            if self._scale != 1.0:
                return 'Unit(dimensionless with a scale of {0})'.format(
                    self._scale)
            else:
                return 'Unit(dimensionless)'

    def __hash__(self):
        parts = zip((hash(x) for x in self._bases), self._powers)
        parts.sort()
        return hash(tuple([self._scale] + parts))

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

    def _expand_and_gather(self, decompose=False, bases=[]):
        def add_unit(unit, power, scale):
            if len(bases) and unit not in bases:
                for base in bases:
                    if unit.is_equivalent(base):
                        scale *= unit.to(base) ** power
                        unit = base
                        break

            new_parts[unit] = power + new_parts.get(unit, 0)
            return scale

        new_parts = {}
        scale = self.scale

        for b, p in zip(self.bases, self.powers):
            if decompose and (not len(bases) or b not in bases):
                b = b.decompose(bases=bases)

            if isinstance(b, CompositeUnit):
                scale *= b.scale ** p
                for b_sub, p_sub in zip(b.bases, b.powers):
                    scale = add_unit(b_sub, p_sub * p, scale)
            else:
                scale = add_unit(b, p, scale)

        new_parts = [(b, p) for (b, p) in new_parts.items() if p != 0]
        new_parts.sort(key=lambda x: x[1], reverse=True)

        self._bases = [x[0] for x in new_parts]
        self._powers = [x[1] for x in new_parts]
        self._scale = scale

    def __copy__(self):
        """
        For compatibility with python copy module.
        """
        return CompositeUnit(self._scale, self._bases[:], self._powers[:])

    def _simplify(self):
        self._expand_and_gather()
        return self
    _simplify.__doc__ = UnitBase._simplify.__doc__

    def decompose(self, bases=[]):
        x = CompositeUnit(self.scale, self.bases, self.powers)
        x._expand_and_gather(True, bases=bases)
        return x
    decompose.__doc__ = UnitBase.decompose.__doc__

    def is_unity(self):
        x = self.decompose()
        return (len(x.bases) == 0 and x.scale == 1.0)
    is_unity.__doc__ = UnitBase.is_unity.__doc__

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

    Returns
    -------
    Scalar value or numpy array

    Raises
    ------
    ValueError
        If value is not as expected
    """
    if isinstance(value, (float, int, long)):
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
