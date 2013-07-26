# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines the SI units.  They are also available in the
`astropy.units` namespace.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..constants import si as _si
from .core import UnitBase, def_unit, _UnitRegistry

import numpy as _numpy

_UnitRegistry().namespace = globals()

###########################################################################
# LENGTH

def_unit(['m', 'meter'], register=True, prefixes=True,
         doc="meter: base unit of length in SI")

def_unit(['micron'], um, register=True,
         doc="micron: alias for micrometer (um)",
         format={'latex': r'\mu', 'unicode': 'μ'})

def_unit(['Angstrom', 'AA', 'angstrom'], 0.1 * nm, register=True,
         doc="ångström: 10 ** -10 m",
         format={'latex': r'\AA', 'unicode': 'Å', 'vounit': 'angstrom'})


###########################################################################
# VOLUMES

def_unit(['l', 'L', 'liter'], 1000 * cm ** 3.0, register=True, prefixes=True,
         doc="liter: metric unit of volume")


###########################################################################
# ANGULAR MEASUREMENTS

def_unit(['rad', 'radian'], register=True, prefixes=True,
         doc="radian: angular measurement of the ratio between the length "
         "on an arc and its radius")
def_unit(['deg', 'degree'], _numpy.pi / 180.0 * rad, register=True,
         doc="degree: angular measurement 1/360 of full rotation",
         format={'latex': r'{}^{\circ}', 'unicode': '°'})
def_unit(['hourangle'], 15.0 * deg, register=True, prefixes=False,
         doc="hour angle: angular measurement with 24 in a full circle",
         format={'latex': r'{}^{h}', 'unicode': 'ʰ'})
def_unit(['arcmin', 'arcminute'], 1.0 / 60.0 * deg, register=True,
         prefixes=True,
         doc="arc minute: angular measurement",
         format={'latex': r'{}^{\prime}', 'unicode': '′'})
def_unit(['arcsec', 'arcsecond'], 1.0 / 3600.0 * deg, register=True,
         prefixes=True,
         doc="arc second: angular measurement")
# These special formats should only be used for the non-prefix versions
arcsec._format = {'latex': r'{}^{\prime\prime}', 'unicode': '″'}
def_unit(['mas'], 0.001 * arcsec, register=True,
         doc="arc second: angular measurement")
def_unit(['uas'], 0.000001 * arcsec, register=True,
         doc="arc second: angular measurement",
         format={'latex': r'\mu as', 'unicode': 'μas'})

def_unit(['sr', 'steradian'], rad ** 2, register=True, prefixes=True,
         doc="steradian: unit of solid angle in SI")


###########################################################################
# TIME

def_unit(['s', 'second'], register=True, prefixes=True,
         exclude_prefixes=['a'],
         doc="second: base unit of time in SI.")

def_unit(['min', 'minute'], 60 * s, register=True)
def_unit(['h', 'hour', 'hr'], 3600 * s, register=True)
def_unit(['d', 'day'], 24 * h, register=True)
def_unit(['sday'], 86164.09053 * s, register=True,
         doc="Sidereal day (sday) is the time of one rotation of the Earth.")
def_unit(['wk', 'week'], 7 * day, register=True)
def_unit(['fortnight'], 2 * wk, register=True)

def_unit(['a', 'annum'], 365.25 * d, register=True, prefixes=True,
         exclude_prefixes=['P'])
def_unit(['yr', 'year'], 365.25 * d, register=True, prefixes=True)


###########################################################################
# FREQUENCY

def_unit(['Hz', 'Hertz', 'hertz'], 1 / s, register=True, prefixes=True,
         doc="Frequency")


###########################################################################
# MASS

def_unit(['kg', 'kilogram'], register=True,
         doc="kilogram: base unit of mass in SI.")
def_unit(['g', 'gram'], 1.0e-3 * kg, register=True, prefixes=True,
         exclude_prefixes=['k', 'kilo'])

def_unit(['t', 'tonne'], 1000 * kg, register=True,
         doc="Metric tonne")


###########################################################################
# AMOUNT OF SUBSTANCE

def_unit(['mol', 'mole'], register=True, prefixes=True,
         doc="mole: amount of a chemical substance in SI.")


###########################################################################
# TEMPERATURE

def_unit(
    ['K', 'Kelvin'], register=True, prefixes=True,
    doc="Kelvin: temperature with a null point at absolute zero.")


###########################################################################
# FORCE

def_unit(['N', 'Newton', 'newton'], kg * m * s ** -2, register=True,
         prefixes=True, doc="Newton: force")


##########################################################################
# ENERGY

def_unit(['J', 'Joule', 'joule'], N * m, register=True, prefixes=True,
         doc="Joule: energy")
def_unit(['eV', 'electronvolt'], _si.e.value * J, register=True, prefixes=True,
         doc="Electron Volt")
def_unit(['Pa', 'Pascal', 'pascal'], J * m ** -3, register=True, prefixes=True,
         doc="Pascal: pressure")


###########################################################################
# POWER

def_unit(['W', 'Watt', 'watt'], J / s, register=True, prefixes=True,
         doc="Watt: power")


###########################################################################
# ELECTRICAL

def_unit(['A', 'ampere', 'amp'], register=True, prefixes=True,
         doc="ampere: base unit of electric current in SI")
def_unit(['C', 'coulomb'], A * s, register=True, prefixes=True,
         doc="coulomb: electric charge")
def_unit(['V', 'Volt', 'volt'], J * C ** -1, register=True, prefixes=True,
         doc="Volt: electric potential or electromotive force")
def_unit(['Ohm', 'ohm'], V * A ** -1, register=True, prefixes=True,
         doc="Ohm: electrical resistance",
         format={'latex': r'\Omega', 'unicode': 'Ω'})
def_unit(['S', 'Siemens', 'siemens'], A * V ** -1, register=True,
         prefixes=True, doc="Siemens: electrical conductance")
def_unit(['F', 'Farad', 'farad'], C * V ** -1, register=True, prefixes=True,
         doc="Farad: electrical capacitance")


###########################################################################
# MAGNETIC

def_unit(['Wb', 'Weber', 'weber'], V * s, register=True, prefixes=True,
         doc="Weber: magnetic flux")
def_unit(['T', 'Tesla', 'tesla'], Wb * m ** -2, register=True, prefixes=True,
         doc="Tesla: magnetic flux density")
def_unit(['H', 'Henry', 'henry'], Wb * A ** -1, register=True, prefixes=True,
         doc="Henry: inductance")


###########################################################################
# ILLUMINATION

def_unit(['cd', 'candela'], register=True, prefixes=True,
         doc="candela: base unit of luminous intensity in SI")
def_unit(['lm', 'lumen'], cd * sr, register=True, prefixes=True,
         doc="lumen: luminous flux")
def_unit(['lx', 'lux'], lm * m ** -2, register=True, prefixes=True,
         doc="lux: luminous emittence")

###########################################################################
# RADIOACTIVITY

def_unit(['Bq', 'becquerel'], Hz, register=True, prefixes=False,
         doc="becquerel: unit of radioactivity")
def_unit(['Ci', 'curie'], Bq / 3.7e10, register=True, prefixes=False,
         doc="curie: unit of radioactivity")


###########################################################################
# BASES

bases = set([m, s, kg, A, cd, rad, K, mol])


###########################################################################
# CLEANUP

del UnitBase
del def_unit


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary
__doc__ += _generate_unit_summary(globals())
