# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from astropy.constants import cgs as _cgs
from astropy.constants import si as _si
from .core import UnitBase, def_unit

import numpy as _numpy

UnitBase._set_namespace(globals())

# TODO: Finish the docstrings

###########################################################################
# LENGTH

def_unit(['m', 'meter'], register=True, prefixes=True,
    doc="""
    The meter (m) is the base unit of length in SI.
    """)

def_unit(['AU'], _si.au * m, register=True, prefixes=True,
     doc="""
     The astronomical unit (AU) is a unit of length equal to
     149,597,870,700 meters, or approximately the mean Earth--Sun
     distance.
     """)

def_unit(['pc', 'parsec'], _si.pc * m, register=True, prefixes=True,
     doc="""
     The parsec (pc) is approximately 3.26 light-years.
     """)

def_unit(['micron'], um, register=True,
     doc="""
     The micron is an alias for the micrometer (um).
     """,
     format={'latex': r'\mu', 'unicode': 'μ'})

def_unit(['AA', 'Angstrom', 'angstrom'], 0.1 * nm, register=True, prefixes=True,
     doc="""
     The ångström is equal to :math:`10^-10` m.
     """,
     format={'latex': r'\AA', 'unicode': 'Å'})

def_unit(['solRad'], _si.R_sun * m, register=True)
def_unit(['lyr'], 9.460730e15 * m, register=True)

# Imperial measurements
def_unit(['inch'], 2.54 * cm, register=True,
     doc="""
     The International inch is a unit of length defined as 25.4
     millimeters.
     """,
     format={'latex': r'\second', 'unicode': '″'})
def_unit(['ft', 'foot'], 12 * inch, register=True,
     doc="""
     The International foot is a unit of length defined as 12 inches,
     or 0.3048 meters.
     """,
     format={'latex': r'\prime', 'unicode': '′'})
def_unit(['yd', 'yard'], 3*ft, register=True,
     doc="""
     The International yard is a unit of length defined as 3 feet or
     0.9144 meters.
     """)
def_unit(['mi', 'mile'], 5280 * ft, register=True,
     doc="""
     The International mile is a unit of length defined as 5280 feet
     or 1,602.344 meters.
     """)


###########################################################################
# AREAS
def_unit(['barn'], 10 ** -28 * m ** 2, register=True, prefixes=True)

def_unit(['acre'], 43560 * ft ** 2, register=True)


###########################################################################
# VOLUMES

def_unit(['l', 'liter'], 1000 * cm ** 3, register=True, prefixes=True)

# Imperial units
def_unit(['gallon'], liter / 0.264172052, register=True)
def_unit(['quart'], gallon / 4, register=True)
def_unit(['pint'], quart / 2, register=True)
def_unit(['cup'], pint / 2, register=True)
def_unit(['foz', 'fluid_oz', 'fluid_ounce'], cup / 8, register=True)
def_unit(['tbsp', 'tablespoon'], foz / 2, register=True)
def_unit(['tsp', 'teaspoon'], tbsp / 3, register=True)


###########################################################################
# ANGULAR MEASUREMENTS

def_unit(['rad', 'radian'], register=True, prefixes=True,
    doc="""
    The radian (rad) is an angular measurement of the ratio between
    the length on an arc and its radius.
    """)
def_unit(['deg', 'degree'], _numpy.pi / 180.0 * rad, register=True,
         format={'latex': r'\degree', 'unicode': '°'})
def_unit(['arcmin', 'arcminute'], 1.0 / 60.0 * deg, register=True,
         format={'latex': r'\prime', 'unicode': '′'})
def_unit(['arcsec', 'arcsecond'], 1.0 / 3600.0 * deg, register=True,
         format={'latex': r'\second', 'unicode': '″'})
def_unit(['mas'], 1.0 / 3600000 * deg, register=True,
         format={'latex': r'\third', 'unicode': '‴'})

def_unit(['sr', 'steradian'], register=True, prefixes=True,
     doc="""
     The steradian (sr) is the SI unit of solid angle.
     """)


###########################################################################
# TIME

def_unit(['s', 'second'], register=True, prefixes=True,
     doc="""
     The second (s) is the base unit of time in SI.
     """)

def_unit(['min', 'minute'], 60 * s, register=True,
     format={'latex': r'^{m}'})
def_unit(['h', 'hour'], 3600 * s, register=True,
     format={'latex': r'^{h}'})
def_unit(['d', 'day'], 24 * h, register=True,
     format={'latex': r'^{d}'})
def_unit(['sday'], 86164.09053 * s, register=True,
     doc="""
     The sidereal day (sday), in a nutshell, is the length of one
     rotation of the Earth.
     """)
def_unit(['wk', 'week'], 7 * day, register=True)
def_unit(['fortnight'], 2 * wk, register=True)

def_unit(['a', 'annum'], 3.1556926e7 * s, register=True,
     prefixes=True, exclude_prefixes=['P'])
def_unit(['yr', 'year'], 3.1556926e7 * s, register=True,
     prefixes=True)



###########################################################################
# FREQUENCY

def_unit(['Hz', 'Hertz', 'hertz'], 1 / s, register=True, prefixes=True)


###########################################################################
# MASS

def_unit(['kg', 'kilogram'], register=True,
         doc="""
         The kilogram (kg) is the base unit of mass in SI.
         """)
def_unit(['g', 'gram'], 1.0e-3 * kg, register=True, prefixes=True,
         exclude_prefixes=['k'])

def_unit(['solMass'], _si.M_sun * kg, register=True, prefixes=True,
     format={'latex': r'M_{\odot}', 'unicode': 'M⊙'})
def_unit(['M_p', 'proton_mass'], _si.m_p * kg, register=True,
         format={'latex': r'M_{p}', 'unicode': 'Mₚ'})
def_unit(['M_e', 'electron_mass'], _si.m_e * kg, register=True,
         format={'latex': r'M_{e}', 'unicode': 'Mₑ'})
# Unified atomic mass unit
def_unit(['u', 'Da', 'Dalton'], 1.6605387e-27 * kg, register=True)
def_unit(['t', 'tonne'], 1000 * kg, register=True)

# Imperial measurements
# well, force actually, but who uses it that way?
def_unit(['oz', 'ounce'], 28.349523125 * g, register=True)
def_unit(['lb', 'pound'], 16 * oz, register=True)
def_unit(['ton'], 2000 * lb, register=True)


###########################################################################
# AMOUNT OF SUBSTANCE

def_unit(['mol', 'mole'], register=True, prefixes=True,
     doc="""
     The mole (mol) is an SI unit to express amounts of a chemical
     substance.
     """)


###########################################################################
# TEMPERATURE

def_unit(
    ['K', 'Kelvin'], register=True, prefixes=True,
    doc="""
    The Kelvin (K) is a unit of measurement for temperature.  It uses
    as its null point absolute zero.
    """)


###########################################################################
# FORCE

def_unit(['N', 'Newton', 'newton'], kg * m * s ** -2, register=True, prefixes=True)


##########################################################################
# ENERGY

def_unit(['J', 'Joule', 'joule'], N * m, register=True, prefixes=True)
def_unit(['eV', 'electronvolt'], _si.e * J, register=True, prefixes=True)
def_unit(['erg'], 1.0e-7 * J, register=True)
def_unit(['Ry', 'rydberg'], 13.605692 * eV, register=True,
         format={'latex': r'R_{\infty}', 'unicode': 'R∞'})

# Imperial units
def_unit(['BTU', 'btu'], 1.05505585e10 * erg, register=True)
def_unit(['cal', 'calorie'], 41840000 * erg, register=True)
def_unit(['kcal', 'Cal', 'Calorie', 'kilocal', 'kilocalorie'],
         1000 * cal, register=True)


###########################################################################
# POWER

def_unit(['W', 'Watt', 'watt'], J / s, register=True, prefixes=True)

# Imperial units
def_unit(['hp', 'horsepower'], W / 0.00134102209, register=True)


###########################################################################
# ELECTRICAL

def_unit(['A', 'ampere', 'amp'], register=True, prefixes=True,
     doc="""
     The ampere (A) is the SI unit of electric current.
     """)
def_unit(['C', 'coulomb'], A * s, register=True, prefixes=True)
def_unit(['V', 'Volt', 'volt'], J * C ** -1, register=True, prefixes=True)
def_unit(['Ohm', 'ohm'], V * A ** -1, register=True, prefixes=True,
     format={'latex': r'\Omega', 'unicode': 'Ω'})
def_unit(['S', 'Siemens', 'siemens'], A * V ** -1, register=True, prefixes=True)
def_unit(['F', 'Farad', 'farad'], C * V ** -1, register=True, prefixes=True)


###########################################################################
# MAGNETIC

def_unit(['Wb', 'Weber', 'weber'], V * s, register=True, prefixes=True)
def_unit(['T', 'Tesla', 'tesla'], Wb * m ** -2, register=True, prefixes=True)
def_unit(['H', 'Henry', 'henry'], Wb * A ** -1, register=True, prefixes=True)
def_unit(['G', 'Gauss', 'gauss'], 1e-4 * T, register=True, prefixes=True)


###########################################################################
# ILLUMINATION

def_unit(['cd', 'candela'], register=True, prefixes=True,
     doc="""
     The candela (cd) is the SI base unit of luminous intensity.
     """)
def_unit(['lm', 'lumen'], cd * sr, register=True, prefixes=True)
def_unit(['lx', 'lux'], lm * m ** -2, register=True, prefixes=True)
def_unit(['solLum'], _si.L_sun * W, register=True, prefixes=True)
def_unit(['Pa', 'Pascal', 'pascal'], J * m ** -3, register=True, prefixes=True)
def_unit(['dyn', 'dyne'], erg * cm ** -3, register=True)


###########################################################################
# EVENTS

def_unit(['ct', 'count'], register=True)
def_unit(['ph', 'photon'], register=True)
def_unit(['pix', 'pixel'], register=True)


###########################################################################
# SPECTRAL DENSITY

def_unit(['Jy', 'Jansky', 'jansky'], 10 ** -23 * erg / Hz / cm ** 2 / s,
         register=True, prefixes=True)
def_unit(['R', 'Rayleigh', 'rayleigh'],
     (1e10 / (4 * _numpy.pi)) * ph * m ** -2 * s ** -1 * sr ** -1,
     register=True, prefixes=True)

def_unit(['mag'], register=True, prefixes=True,
     doc="""
     Stellar magnitude.
     """)


###########################################################################
# MISCELLANEOUS

# Some of these are very FITS-specific and perhaps considered a mistake.
# Maybe they should be moved into the FITS format class?
def_unit(['D', 'Debye', 'debye'], (1.0 / 3.0) * 1e-29 * C * m, register=True)
# TODO: This is defined by the FITS standard as "relative to the sun".
# Is that mass, volume, what?
def_unit(['Sun'], register=True)
def_unit(['chan'], register=True)
def_unit(['bin'], register=True)
def_unit(['vox', 'voxel'], register=True)
def_unit(['bit'], register=True, prefixes=True)
def_unit(['byte'], register=True, prefixes=True)
def_unit(['adu'], register=True)
def_unit(['beam'], register=True)


###########################################################################
# EQUIVALENCIES

def spectral():
    """
    Returns a list of equivalence pairs that handle spectral
    wavelength, frequency, and energy equivalences.

    Allows conversions between wavelength units, frequency units and
    energy units as they relate to light.
    """

    return [
        (m, Hz, lambda x: _si.c / x),
        (m, J, lambda x: (_si.c * _si.h) / x),
        (Hz, J, lambda x: _si.h * x)
    ]
sp = spectral


def spectral_density(sunit, sfactor):
    """
    Returns a list of equivalence pairs that handle spectral density
    with regard to wavelength and frequency.
    """
    c_Aps = _si.c * 10 ** 10

    flambda = erg / angstrom / cm ** 2 / s
    fnu = erg / Hz / cm ** 2 / s

    def converter(x):
        return x * (sunit.to(AA, sfactor, sp()) ** 2 / c_Aps)

    def iconverter(x):
        return x / (sunit.to(AA, sfactor, sp()) ** 2 / c_Aps)

    return [
        (AA, fnu, converter, iconverter),
        (flambda, fnu, converter, iconverter),
        (AA, Hz, converter, iconverter),
        (flambda, Hz, converter, iconverter)
        ]

sd = spectral_density


###########################################################################
# CLEANUP

del UnitBase
del def_unit
