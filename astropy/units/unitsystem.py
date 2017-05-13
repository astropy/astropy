# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import division, print_function

from collections import OrderedDict

from .core import one, Unit
# from .physical import _physical_unit_mapping

__all__ = ['UnitSystem', 'DimensionlessUnitSystem']

class UnitSystem(object):
    """
    Represents a system of units. At minimum, this consists of a set of
    length, time, mass, and angle units, but may also contain preferred
    representations for composite units. For example, the base unit system
    could be ``{kpc, Myr, Msun, radian}``, but you can also specify a preferred
    speed, such as ``km/s``.

    This class functions like a dictionary with keys set by physical types.
    If a unit for a particular physical type is not specified on creation,
    a composite unit will be created with the base units. See Examples below
    for some demonstrations.

    Parameters
    ----------
    *units
        The units that define the unit system. At minimum, this must
        contain length, time, mass, and angle units.

    Examples
    --------
    If only base units are specified, any physical type specified as a key
    to this object will be composed out of the base units::

        >>> usys = UnitSystem(u.m, u.s, u.kg, u.radian)
        >>> usys['energy']
        Unit("kg m2 / s2")

    However, custom representations for composite units can also be specified
    when initializing::

        >>> usys = UnitSystem(u.m, u.s, u.kg, u.radian, u.erg)
        >>> usys['energy']
        Unit("erg")

    This is useful for Galactic dynamics where lengths and times are usually
    given in terms of ``kpc`` and ``Myr``, but speeds are given in ``km/s``::

        >>> usys = UnitSystem(u.kpc, u.Myr, u.Msun, u.radian, u.km/u.s)
        >>> usys['speed']
        Unit("km / s")

    """

    def __init__(self, units, *args):

        self._required_physical_types = ['length', 'time', 'mass', 'angle']
        self._core_units = []

        if isinstance(units, UnitSystem):
            self._registry = units._registry.copy()
            self._core_units = units._core_units
            return

        if len(args) > 0:
            units = (units,) + tuple(args)

        self._registry = OrderedDict()
        for unit in units:
            typ = unit.physical_type
            if typ in self._registry:
                raise ValueError("Multiple units passed in with type '{0}'".format(typ))
            self._registry[typ] = unit

        for phys_type in self._required_physical_types:
            if phys_type not in self._registry:
                raise ValueError("You must specify a unit with physical type '{0}'".format(phys_type))
            self._core_units.append(self._registry[phys_type])

    def __getitem__(self, key):
        from .physical import _physical_unit_mapping
        if key in self._registry:
            return self._registry[key]

        else:
            unit = None
            for k,v in _physical_unit_mapping.items():
                if v == key:
                    unit = Unit(" ".join(["{}**{}".format(x,y) for x,y in k]))
                    break

            if unit is None:
                raise ValueError("Physical type '{0}' doesn't exist in unit "
                                 "registry.".format(key))

            unit = unit.decompose(self._core_units)
            unit._scale = 1.
            return unit

    def __len__(self):
        return len(self._core_units)

    def __iter__(self):
        for uu in self._core_units:
            yield uu

    def __str__(self):
        unit_str = ",".join([str(uu) for uu in self._core_units])
        return "UnitSystem ({0})".format(unit_str)

    def __repr__(self):
        return "<{0}>".format(self.__str__())

    def __eq__(self, other):
        for k in self._registry:
            if not self[k] == other[k]:
                return False

        for k in other._registry:
            if not self[k] == other[k]:
                return False

        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_dict(self):
        """
        Return a dictionary representation of the unit system with keys
        set by the physical types and values set by the unit objects.
        """
        return self._registry.copy()

    def get_constant(self, name):
        """
        Retrieve a constant with specified name in this unit system.

        Parameters
        ----------
        name : str
            The name of the constant, e.g., G.

        Returns
        -------
        const : `~astropy.units.Quantity`
            The value of the constant represented in this unit system.

        Examples
        --------

            >>> usys = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun)
            >>> usys.get_constant('c')
            <Quantity 306.6013937879527 kpc / Myr>

        """
        from .. import constants as const
        try:
            c = getattr(const, name)
        except AttributeError:
            raise ValueError("Constant name '{}' doesn't exist in "
                             "astropy.constants".format(name))

        return c.decompose(self._core_units)

class DimensionlessUnitSystem(UnitSystem):

    def __init__(self):
        self._core_units = [one]
        self._registry = OrderedDict()
        self._registry['dimensionless'] = one

    def __getitem__(self, key):
        return one

    def __str__(self):
        return "UnitSystem (dimensionless)"

    def to_dict(self):
        raise ValueError("Cannot represent dimensionless unit system as dict!")

    def get_constant(self, name):
        raise ValueError("Cannot get constant in dimensionless units!")
