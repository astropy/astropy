# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Colloquially used Imperial units.

These units are available in the `astropy.units.imperial` namespace, but not in the
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

__all__: list[str] = []  #  Units are added at the end

from . import si
from .core import def_unit
from .docgen import generate_dunder_all, generate_unit_summary

###########################################################################
# LENGTH

inch = def_unit(["inch"], 2.54 * si.cm, doc="International inch")
ft = foot = def_unit(["ft", "foot"], 12 * inch, doc="International foot")
yd = yard = def_unit(["yd", "yard"], 3 * ft, doc="International yard")
mi = mile = def_unit(["mi", "mile"], 5280 * ft, doc="International mile")
mil = thou = def_unit(["mil", "thou"], 0.001 * inch, doc="Thousandth of an inch")
nmi = nauticalmile = NM = def_unit(
    ["nmi", "nauticalmile", "NM"], 1852 * si.m, doc="Nautical mile"
)
fur = furlong = def_unit(["fur", "furlong"], 660 * ft, doc="Furlong")


###########################################################################
# AREAS

ac = acre = def_unit(["ac", "acre"], 43560 * ft**2, doc="International acre")


###########################################################################
# VOLUMES

gallon = def_unit(["gallon"], si.liter / 0.264172052, doc="U.S. liquid gallon")
quart = def_unit(["quart"], gallon / 4, doc="U.S. liquid quart")
pint = def_unit(["pint"], quart / 2, doc="U.S. liquid pint")
cup = def_unit(["cup"], pint / 2, doc="U.S. customary cup")
foz = fluid_oz = fluid_ounce = def_unit(
    ["foz", "fluid_oz", "fluid_ounce"], cup / 8, doc="U.S. fluid ounce"
)
tbsp = tablespoon = def_unit(
    ["tbsp", "tablespoon"], foz / 2, doc="U.S. customary tablespoon"
)
tsp = teaspoon = def_unit(["tsp", "teaspoon"], tbsp / 3, doc="U.S. customary teaspoon")


###########################################################################
# MASS

oz = ounce = def_unit(
    ["oz", "ounce"], 28.349523125 * si.g, doc="International avoirdupois ounce: mass"
)
lb = lbm = pound = def_unit(
    ["lb", "lbm", "pound"], 16 * oz, doc="International avoirdupois pound: mass"
)
st = stone = def_unit(
    ["st", "stone"], 14 * lb, doc="International avoirdupois stone: mass"
)
ton = def_unit(["ton"], 2000 * lb, doc="International avoirdupois ton: mass")
slug = def_unit(["slug"], 32.174049 * lb, doc="slug: mass")


###########################################################################
# SPEED

kn = kt = knot = NMPH = def_unit(
    ["kn", "kt", "knot", "NMPH"],
    nmi / si.h,
    doc="nautical unit of speed: 1 nmi per hour",
)


###########################################################################
# FORCE

lbf = def_unit("lbf", slug * ft * si.s**-2, doc="Pound: force")
kip = kilopound = def_unit(["kip", "kilopound"], 1000 * lbf, doc="Kilopound: force")


##########################################################################
# ENERGY

BTU = btu = def_unit(["BTU", "btu"], 1.05505585 * si.kJ, doc="British thermal unit")
cal = calorie = def_unit(
    ["cal", "calorie"],
    4.184 * si.J,
    doc="Thermochemical calorie: pre-SI metric unit of energy",
)
kcal = Cal = Calorie = kilocal = kilocalorie = def_unit(
    ["kcal", "Cal", "Calorie", "kilocal", "kilocalorie"],
    1000 * cal,
    doc="Calorie: colloquial definition of Calorie",
)


##########################################################################
# PRESSURE

psi = def_unit("psi", lbf * inch**-2, doc="Pound per square inch: pressure")


###########################################################################
# POWER

# Imperial units
hp = horsepower = def_unit(
    ["hp", "horsepower"], si.W / 0.00134102209, doc="Electrical horsepower"
)


###########################################################################
# TEMPERATURE

deg_F = Fahrenheit = def_unit(
    ["deg_F", "Fahrenheit"],
    doc="Degrees Fahrenheit",
    format={"latex": r"{}^{\circ}F", "unicode": "°F"},
)
deg_R = Rankine = def_unit(
    ["deg_R", "Rankine"],
    5 / 9 * si.K,
    doc="Rankine scale: absolute scale of thermodynamic temperature",
    format={"latex": r"{}^{\circ}R", "unicode": "°R"},
)

###########################################################################
# ALL & DOCSTRING

__all__ += generate_dunder_all(globals())  # noqa: PLE0605

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    __doc__ += generate_unit_summary(globals())


def enable():
    """
    Enable Imperial units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable Imperial
    units only temporarily.
    """
    from .core import add_enabled_units

    return add_enabled_units(globals())


# Cleanup for the benefit of e.g. IDEs or other tools that don't execute
# generate_dunder_all(). This is simpler than writing out __all__ explicitly.
del si, def_unit, generate_dunder_all, generate_unit_summary
