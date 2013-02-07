# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools
import types
import warnings

from ..units.core import Unit, UnitsException
from ..units.quantity import Quantity
from ..utils import lazyproperty

__all__ = ['Constant']


class ConstantMeta(type):
    def __new__(mcls, name, bases, d):
        def wrap(meth):
            @functools.wraps(meth)
            def wrapper(self, *args, **kwargs):
                name_lower = self.name.lower()
                if not self._checked_units:
                    instances = Constant._registry[name_lower]
                    for inst in instances.values():
                        try:
                            self.unit.to(inst.unit)
                        except UnitsException:
                            Constant._has_incompatible_units.add(name_lower)
                    self._checked_units = True

                if (not self.system and
                        name_lower in Constant._has_incompatible_units):
                    systems = sorted(filter(None, instances))
                    raise TypeError(
                        'Constant {0!r} does not have physically compatible '
                        'units across all systems of units and cannot be '
                        'combined with other values without specifying a '
                        'system (eg. {1}.{2})'.format(self.abbrev, self.abbrev,
                                                      systems[0]))

                return meth(self, *args, **kwargs)

            return wrapper

        # The wrapper applies to so many of the __ methods that it's easier to
        # just exclude the ones it doesn't apply to
        exclude = set(['__init__', '__str__', '__repr__'])
        for attr, value in vars(Quantity).items():
            if (isinstance(value, types.FunctionType) and
                    attr.startswith('__') and attr.endswith('__') and
                    attr not in exclude):
                d[attr] = wrap(value)

        return super(ConstantMeta, mcls).__new__(mcls, name, bases, d)


class Constant(Quantity):
    """A physical or astronomical constant.

    These objects are quantities that are meant to represent physical constants

    Attributes
    ----------
    name : str
        The name of this constant.
    uncertainty : float
        The uncertainty in the value of this constant.
    reference : str
        The source used for the value of this constant.
    units : `astropy.units.UnitBase` instance
        The units of this constant. Can be set either as a string or
        `astropy.units.UnitBase`.
    """

    __metaclass__ = ConstantMeta

    _registry = {}
    _has_incompatible_units = set()

    def __new__(cls, abbrev, name, value, unit, uncertainty, reference,
                system=None):
        name_lower = name.lower()
        instances = Constant._registry.setdefault(name_lower, {})
        if system in instances:
            warnings.warn('Constant {0!r} is already has a definition in the '
                          '{1!r} system'.format(name, system))

        inst = super(Constant, cls).__new__(cls, abbrev, name, value, unit,
                                            uncertainty, reference, system)

        for c in instances.values():
            if system is not None and not hasattr(c.__class__, system):
                setattr(c, system, inst)
            if c.system is not None and not hasattr(inst.__class__, c.system):
                setattr(inst, c.system, c)

        instances[system] = inst

        return inst

    def __init__(self, abbrev, name, value, unit, uncertainty, reference,
                 system=None):
        self._abbrev = abbrev
        self._name = name
        self._value = value
        self._unit = unit
        self._uncertainty = uncertainty
        self._reference = reference
        self._system = system

        self._checked_units = False

    def __repr__(self):
        return ('<Constant name={0!r} value={1} error={2} units={3!r} '
                'reference={4!r}>'.format(self.name, self.value,
                                          self.uncertainty, str(self.unit),
                                          self.reference))

    def __str__(self):
        return ('  Name   = {0}\n'
                '  Value  = {1}\n'
                '  Error  = {2}\n'
                '  Units  = {3}\n'
                '  Reference = {4}'.format(self.name, self.value,
                                           self.uncertainty, self.unit,
                                           self.reference))

    @property
    def abbrev(self):
        return self._abbrev

    @property
    def name(self):
        return self._name

    @lazyproperty
    def unit(self):
        return Unit(self._unit)

    @property
    def uncertainty(self):
        return self._uncertainty

    @property
    def reference(self):
        return self._reference

    @property
    def system(self):
        return self._system

    @property
    def si(self):
        """If the Constant is defined in the SI system return that instance of
        the constant, else convert to a Quantity in the appropriate SI units.
        """

        instances = Constant._registry[self.name.lower()]
        return instances.get('si') or super(Constant, self).si

    @property
    def cgs(self):
        """If the Constant is defined in the CGS system return that instance of
        the constant, else convert to a Quantity in the appropriate CGS units.
        """

        instances = Constant._registry[self.name.lower()]
        return instances.get('cgs') or super(Constant, self).cgs


class EMConstant(Constant):
    """An electromagnetic constant"""

    @property
    def cgs(self):
        raise TypeError("Cannot convert EM constants to cgs because there "
                        "are different systems for E.M constants within the "
                        "c.g.s system (ESU, Gaussian, etc.). Instead, "
                        "directly use the constant with the appropriate "
                        "suffix (e.g. e.esu, e.gauss, etc.).")
