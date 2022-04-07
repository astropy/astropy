# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines colloquially used Imperial units.  They are
available in the `astropy.units.imperial` namespace, but not in the
top-level `astropy.units` namespace, e.g.::

    >>> import astropy.units as u
    >>> mph = u.imperial.mile / u.hour
    >>> mph
    Unit("mi / h")

To include them in `~astropy.units.UnitBase.compose` and the results of
`~astropy.units.UnitBase.find_equivalent_units`, do::

    >>> import astropy.units as u
    >>> u.imperial.enable()  # doctest: +SKIP
"""


from . import si
from .core import UnitBase, def_unit

_ns = globals()

###########################################################################
# LENGTH

def_unit(['inch'], 2.54 * si.cm, namespace=_ns,
         doc="International inch")
def_unit(['ft', 'foot'], 12 * inch, namespace=_ns,
         doc="International foot")
def_unit(['yd', 'yard'], 3 * ft, namespace=_ns,
         doc="International yard")
def_unit(['mi', 'mile'], 5280 * ft, namespace=_ns,
         doc="International mile")
def_unit(['mil', 'thou'], 0.001 * inch, namespace=_ns,
         doc="Thousandth of an inch")
def_unit(['nmi', 'nauticalmile', 'NM'], 1852 * si.m, namespace=_ns,
         doc="Nautical mile")
def_unit(['fur', 'furlong'], 660 * ft, namespace=_ns,
         doc="Furlong")


###########################################################################
# AREAS

def_unit(['ac', 'acre'], 43560 * ft ** 2, namespace=_ns,
         doc="International acre")


###########################################################################
# VOLUMES

def_unit(['gallon'], si.liter / 0.264172052, namespace=_ns,
         doc="U.S. liquid gallon")
def_unit(['quart'], gallon / 4, namespace=_ns,
         doc="U.S. liquid quart")
def_unit(['pint'], quart / 2, namespace=_ns,
         doc="U.S. liquid pint")
def_unit(['cup'], pint / 2, namespace=_ns,
         doc="U.S. customary cup")
def_unit(['foz', 'fluid_oz', 'fluid_ounce'], cup / 8, namespace=_ns,
         doc="U.S. fluid ounce")
def_unit(['tbsp', 'tablespoon'], foz / 2, namespace=_ns,
         doc="U.S. customary tablespoon")
def_unit(['tsp', 'teaspoon'], tbsp / 3, namespace=_ns,
         doc="U.S. customary teaspoon")


###########################################################################
# MASS

def_unit(['oz', 'ounce'], 28.349523125 * si.g, namespace=_ns,
         doc="International avoirdupois ounce: mass")
def_unit(['lb', 'lbm', 'pound'], 16 * oz, namespace=_ns,
         doc="International avoirdupois pound: mass")
def_unit(['st', 'stone'], 14 * lb, namespace=_ns,
         doc="International avoirdupois stone: mass")
def_unit(['ton'], 2000 * lb, namespace=_ns,
         doc="International avoirdupois ton: mass")
def_unit(['slug'], 32.174049 * lb, namespace=_ns,
         doc="slug: mass")


###########################################################################
# SPEED

def_unit(['kn', 'kt', 'knot', 'NMPH'], nmi / si.h, namespace=_ns,
         doc="nautical unit of speed: 1 nmi per hour")


###########################################################################
# FORCE

def_unit('lbf', slug * ft * si.s**-2, namespace=_ns,
         doc="Pound: force")
def_unit(['kip', 'kilopound'], 1000 * lbf, namespace=_ns,
         doc="Kilopound: force")


##########################################################################
# ENERGY

def_unit(['BTU', 'btu'], 1.05505585 * si.kJ, namespace=_ns,
         doc="British thermal unit")
def_unit(['cal', 'calorie'], 4.184 * si.J, namespace=_ns,
         doc="Thermochemical calorie: pre-SI metric unit of energy")
def_unit(['kcal', 'Cal', 'Calorie', 'kilocal', 'kilocalorie'],
         1000 * cal, namespace=_ns,
         doc="Calorie: colloquial definition of Calorie")


##########################################################################
# PRESSURE

def_unit('psi', lbf * inch ** -2, namespace=_ns,
         doc="Pound per square inch: pressure")


###########################################################################
# POWER

# Imperial units
def_unit(['hp', 'horsepower'], si.W / 0.00134102209, namespace=_ns,
         doc="Electrical horsepower")


###########################################################################
# TEMPERATURE

def_unit(['deg_F', 'Fahrenheit'], namespace=_ns, doc='Degrees Fahrenheit',
         format={'latex': r'{}^{\circ}F', 'unicode': '°F'})
def_unit(['deg_R', 'Rankine'], namespace=_ns, doc='Rankine scale: absolute scale of thermodynamic temperature')


###########################################################################
# CLEANUP

del UnitBase
del def_unit


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary

if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())


def enable():
    """
    Enable Imperial units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable Imperial
    units only temporarily.
    """
    # Local import to avoid cyclical import
    # Local import to avoid polluting namespace
    import inspect

    from .core import add_enabled_units
    return add_enabled_units(inspect.getmodule(enable))
