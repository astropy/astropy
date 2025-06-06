# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines the SI units.  They are also available in
(and should be used through) the `astropy.units` namespace.
"""
# avoid ruff complaints about undefined names defined by def_unit
# ruff: noqa: F821

import numpy as np

from .core import CompositeUnit, UnitBase, def_unit
from .docgen import generate_unit_summary

__all__: list[str] = []  #  Units are added at the end

_ns = globals()


###########################################################################
# DIMENSIONLESS

def_unit(
    ["percent", "pct"],
    CompositeUnit(0.01, [], []),
    namespace=_ns,
    prefixes=False,
    doc="percent: one hundredth of unity, factor 0.01",
    format={"generic": "%", "console": "%", "cds": "%", "latex": r"\%", "unicode": "%"},
)

###########################################################################
# LENGTH

# We exclude c as a prefix so we can define cm as a derived unit,
# not a PrefixUnit. That way it can be used in cgs as a base unit.
def_unit(
    ["m", "meter"],
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["c"],
    doc="meter: base unit of length in SI",
)
def_unit(
    ["cm", "centimeter"],
    0.01 * m,
    namespace=_ns,
    prefixes=False,
    doc="cm: 0.01 m in SI, base unit of length in cgs",
)
def_unit(
    ["micron"],
    um,
    namespace=_ns,
    doc="micron: alias for micrometer (um)",
    format={"latex": r"\mu m", "unicode": "\N{MICRO SIGN}m"},
)
def_unit(
    ["Angstrom", "AA", "angstrom", "Å"],
    0.1 * nm,
    namespace=_ns,
    doc="ångström: 10 ** -10 m",
    prefixes=[(["m", "milli"], ["milli", "m"], 1.0e-3)],
    format={
        "latex": r"\mathring{A}",
        "ogip": "angstrom",
        "unicode": "Å",
        "vounit": "Angstrom",
    },
)


###########################################################################
# AREA

def_unit(
    ["ha", "hectare"],
    1e4 * m**2,
    namespace=_ns,
    prefixes=False,
    doc="hectare: unit of area, used to express area of land",
)


###########################################################################
# VOLUMES

def_unit(
    (["l", "L", "ℓ"], ["liter"]),
    1000 * cm**3.0,
    namespace=_ns,
    prefixes=True,
    format={"latex": r"\mathcal{l}", "unicode": "ℓ"},
    doc="liter: metric unit of volume",
)


###########################################################################
# ANGULAR MEASUREMENTS

def_unit(
    ["rad", "radian"],
    namespace=_ns,
    prefixes=True,
    doc=(
        "radian: angular measurement of the ratio between the length "
        "on an arc and its radius"
    ),
)
def_unit(
    ["deg", "degree"],
    np.pi / 180.0 * rad,
    namespace=_ns,
    prefixes=True,
    doc="degree: angular measurement 1/360 of full rotation",
)
def_unit(
    ["hourangle"],
    15.0 * deg,
    namespace=_ns,
    prefixes=False,
    doc="hour angle: angular measurement with 24 in a full circle",
    format={"latex": r"{}^{h}", "unicode": "ʰ"},
)
def_unit(
    ["arcmin", "arcminute"],
    1.0 / 60.0 * deg,
    namespace=_ns,
    prefixes=True,
    doc="arc minute: angular measurement",
)
def_unit(
    ["arcsec", "arcsecond"],
    1.0 / 3600.0 * deg,
    namespace=_ns,
    prefixes=True,
    doc="arc second: angular measurement",
)
# These special formats should only be used for the non-prefix versions
deg._format = {"latex": r"{}^{\circ}", "unicode": "°"}
arcmin._format = {"latex": r"{}^{\prime}", "unicode": "′"}
arcsec._format = {"latex": r"{}^{\prime\prime}", "unicode": "″"}
def_unit(
    ["mas"],
    0.001 * arcsec,
    namespace=_ns,
    doc="milli arc second: angular measurement",
)
def_unit(
    ["uas", "\N{MICRO SIGN}as", "\N{GREEK SMALL LETTER MU}as"],
    0.000001 * arcsec,
    namespace=_ns,
    doc="micro arc second: angular measurement",
    format={"latex": r"\mu as", "unicode": "μas"},
)
def_unit(
    ["sr", "steradian"],
    rad**2,
    namespace=_ns,
    prefixes=True,
    doc="steradian: unit of solid angle in SI",
)


###########################################################################
# TIME

def_unit(
    ["s", "second"],
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["a"],
    doc="second: base unit of time in SI.",
)
def_unit(
    ["min", "minute"],
    60 * s,
    prefixes=True,
    namespace=_ns,
)
def_unit(
    ["h", "hour", "hr"],
    3600 * s,
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["p"],
)
def_unit(
    ["d", "day"],
    24 * h,
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["c", "y"],
)
def_unit(
    ["sday"],
    86164.09053 * s,
    namespace=_ns,
    doc="Sidereal day (sday) is the time of one rotation of the Earth.",
)
def_unit(
    ["wk", "week"],
    7 * day,
    namespace=_ns,
)
def_unit(
    ["fortnight"],
    2 * wk,
    namespace=_ns,
)
def_unit(
    ["a", "annum"],
    365.25 * d,
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["P", "h"],  # Avoid possible confusion with Pascal and hectare
)
def_unit(
    ["yr", "year"],
    365.25 * d,
    namespace=_ns,
    prefixes=True,
)


###########################################################################
# FREQUENCY

def_unit(
    ["Hz", "Hertz", "hertz"],
    1 / s,
    namespace=_ns,
    prefixes=True,
    doc="Frequency",
)


###########################################################################
# MASS

def_unit(
    ["kg", "kilogram"],
    namespace=_ns,
    doc="kilogram: base unit of mass in SI.",
)
def_unit(
    ["g", "gram"],
    1.0e-3 * kg,
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["k", "kilo"],
)
def_unit(
    ["t", "tonne"],
    1000 * kg,
    namespace=_ns,
    doc="Metric tonne",
)


###########################################################################
# AMOUNT OF SUBSTANCE

def_unit(
    ["mol", "mole"],
    namespace=_ns,
    prefixes=True,
    doc="mole: amount of a chemical substance in SI.",
)
def_unit(
    ["kat", "katal"],
    mol * s**-1,
    namespace=_ns,
    prefixes=True,
    doc="katal: catalytic activity.",
)


###########################################################################
# TEMPERATURE

def_unit(
    ["K", "Kelvin"],
    namespace=_ns,
    prefixes=True,
    doc="Kelvin: temperature with a null point at absolute zero.",
)
def_unit(
    ["deg_C", "Celsius"],
    namespace=_ns,
    doc="Degrees Celsius",
    format={"latex": r"{}^{\circ}C", "unicode": "°C", "fits": "Celsius"},
)


###########################################################################
# FORCE

def_unit(
    ["N", "Newton", "newton"],
    kg * m * s**-2,
    namespace=_ns,
    prefixes=True,
    doc="Newton: force",
)


##########################################################################
# ENERGY

def_unit(
    ["J", "Joule", "joule"],
    N * m,
    namespace=_ns,
    prefixes=True,
    doc="Joule: energy",
)


##########################################################################
# PRESSURE

def_unit(
    ["Pa", "Pascal", "pascal"],
    J * m**-3,
    namespace=_ns,
    prefixes=True,
    doc="Pascal: pressure",
)


###########################################################################
# POWER

def_unit(
    ["W", "Watt", "watt"],
    J / s,
    namespace=_ns,
    prefixes=True,
    doc="Watt: power",
)


###########################################################################
# ELECTRICAL

def_unit(
    ["A", "ampere", "amp"],
    namespace=_ns,
    prefixes=True,
    doc="ampere: base unit of electric current in SI",
)
def_unit(
    ["C", "coulomb"],
    A * s,
    namespace=_ns,
    prefixes=True,
    doc="coulomb: electric charge",
)
def_unit(
    ["V", "Volt", "volt"],
    J * C**-1,
    namespace=_ns,
    prefixes=True,
    doc="Volt: electric potential or electromotive force",
)
def_unit(
    (["Ohm", "ohm", "Ω"], ["Ohm"]),
    V * A**-1,
    namespace=_ns,
    prefixes=True,
    doc="Ohm: electrical resistance",
    format={"latex": r"\Omega", "ogip": "ohm", "unicode": "Ω"},
)
def_unit(
    ["S", "Siemens", "siemens"],
    A * V**-1,
    namespace=_ns,
    prefixes=True,
    doc="Siemens: electrical conductance",
)
def_unit(
    ["F", "Farad", "farad"],
    C * V**-1,
    namespace=_ns,
    prefixes=True,
    doc="Farad: electrical capacitance",
)


###########################################################################
# MAGNETIC

def_unit(
    ["Wb", "Weber", "weber"],
    V * s,
    namespace=_ns,
    prefixes=True,
    doc="Weber: magnetic flux",
)
def_unit(
    ["T", "Tesla", "tesla"],
    Wb * m**-2,
    namespace=_ns,
    prefixes=True,
    doc="Tesla: magnetic flux density",
)
def_unit(
    ["H", "Henry", "henry"],
    Wb * A**-1,
    namespace=_ns,
    prefixes=True,
    doc="Henry: inductance",
)


###########################################################################
# ILLUMINATION

def_unit(
    ["cd", "candela"],
    namespace=_ns,
    prefixes=True,
    doc="candela: base unit of luminous intensity in SI",
)
def_unit(
    ["lm", "lumen"],
    cd * sr,
    namespace=_ns,
    prefixes=True,
    doc="lumen: luminous flux",
)
def_unit(
    ["lx", "lux"],
    lm * m**-2,
    namespace=_ns,
    prefixes=True,
    doc="lux: luminous emittance",
)

###########################################################################
# RADIOACTIVITY

def_unit(
    ["Bq", "becquerel"],
    1 / s,
    namespace=_ns,
    prefixes=True,
    doc="becquerel: unit of radioactivity",
)
def_unit(
    ["Ci", "curie"],
    Bq * 3.7e10,
    namespace=_ns,
    prefixes=True,
    doc="curie: unit of radioactivity",
)
def_unit(
    ["Gy", "gray"],
    J * kg**-1,
    namespace=_ns,
    prefixes=True,
    doc="gray: absorbed dose of ionizing radiation or kinetic energy released per unit mass (kerma)",
)
def_unit(
    ["Sv", "sievert"],
    J * kg**-1,
    namespace=_ns,
    prefixes=True,
    doc="sievert: equivalent dose of ionizing radiation",
)


###########################################################################
# BASES

bases = {m, s, kg, A, cd, rad, K, mol}


###########################################################################
# ALL & DOCSTRING

__all__ += [n for n, v in _ns.items() if isinstance(v, UnitBase)]

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    __doc__ += generate_unit_summary(globals())
