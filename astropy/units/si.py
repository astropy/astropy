# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines the SI units.  They are also available in the
`astropy.units` namespace.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..constants import si as _si
from .core import UnitBase, Unit, def_unit

import numpy as _numpy

_ns = globals()


###########################################################################
# DIMENSIONLESS

def_unit(['percent', 'pct'], Unit(0.01), namespace=_ns, prefixes=False,
         doc="percent: one hundredth of unity, factor 0.01",
         format={'generic': '%', 'console': '%', 'cds': '%',
                 'latex': r'\%', 'unicode': '%'})

###########################################################################
# LENGTH

def_unit(['m', 'meter'], namespace=_ns, prefixes=True,
         doc="meter: base unit of length in SI")

def_unit(['micron'], um, namespace=_ns,
         doc="micron: alias for micrometer (um)",
         format={'latex': r'\mu m', 'unicode': 'μm'})

def_unit(['Angstrom', 'AA', 'angstrom'], 0.1 * nm, namespace=_ns,
         doc="ångström: 10 ** -10 m",
         format={'latex': r'\mathring{A}', 'unicode': 'Å',
                 'vounit': 'Angstrom'})


###########################################################################
# VOLUMES

def_unit((['l', 'L'], ['liter']), 1000 * cm ** 3.0, namespace=_ns, prefixes=True,
         format={'latex': r'\mathcal{l}', 'unicode': 'ℓ'},
         doc="liter: metric unit of volume")


###########################################################################
# ANGULAR MEASUREMENTS

def_unit(['rad', 'radian'], namespace=_ns, prefixes=True,
         doc="radian: angular measurement of the ratio between the length "
         "on an arc and its radius")
def_unit(['deg', 'degree'], _numpy.pi / 180.0 * rad, namespace=_ns,
         prefixes=True,
         doc="degree: angular measurement 1/360 of full rotation",
         format={'latex': r'{}^{\circ}', 'unicode': '°'})
def_unit(['hourangle'], 15.0 * deg, namespace=_ns, prefixes=False,
         doc="hour angle: angular measurement with 24 in a full circle",
         format={'latex': r'{}^{h}', 'unicode': 'ʰ'})
def_unit(['arcmin', 'arcminute'], 1.0 / 60.0 * deg, namespace=_ns,
         prefixes=True,
         doc="arc minute: angular measurement",
         format={'latex': r'{}^{\prime}', 'unicode': '′'})
def_unit(['arcsec', 'arcsecond'], 1.0 / 3600.0 * deg, namespace=_ns,
         prefixes=True,
         doc="arc second: angular measurement")
# These special formats should only be used for the non-prefix versions
arcsec._format = {'latex': r'{}^{\prime\prime}', 'unicode': '″'}
def_unit(['mas'], 0.001 * arcsec, namespace=_ns,
         doc="arc second: angular measurement")
def_unit(['uas'], 0.000001 * arcsec, namespace=_ns,
         doc="arc second: angular measurement",
         format={'latex': r'\mu as', 'unicode': 'μas'})

def_unit(['sr', 'steradian'], rad ** 2, namespace=_ns, prefixes=True,
         doc="steradian: unit of solid angle in SI")


###########################################################################
# TIME

def_unit(['s', 'second'], namespace=_ns, prefixes=True,
         exclude_prefixes=['a'],
         doc="second: base unit of time in SI.")

def_unit(['min', 'minute'], 60 * s, prefixes=True, namespace=_ns)
def_unit(['h', 'hour', 'hr'], 3600 * s, namespace=_ns, prefixes=True,
         exclude_prefixes=['p'])
def_unit(['d', 'day'], 24 * h, namespace=_ns, prefixes=True,
         exclude_prefixes=['c', 'y'])
def_unit(['sday'], 86164.09053 * s, namespace=_ns,
         doc="Sidereal day (sday) is the time of one rotation of the Earth.")
def_unit(['wk', 'week'], 7 * day, namespace=_ns)
def_unit(['fortnight'], 2 * wk, namespace=_ns)

def_unit(['a', 'annum'], 365.25 * d, namespace=_ns, prefixes=True,
         exclude_prefixes=['P'])
def_unit(['yr', 'year'], 365.25 * d, namespace=_ns, prefixes=True)


###########################################################################
# FREQUENCY

def_unit(['Hz', 'Hertz', 'hertz'], 1 / s, namespace=_ns, prefixes=True,
         doc="Frequency")


###########################################################################
# MASS

def_unit(['kg', 'kilogram'], namespace=_ns,
         doc="kilogram: base unit of mass in SI.")
def_unit(['g', 'gram'], 1.0e-3 * kg, namespace=_ns, prefixes=True,
         exclude_prefixes=['k', 'kilo'])

def_unit(['t', 'tonne'], 1000 * kg, namespace=_ns,
         doc="Metric tonne")


###########################################################################
# AMOUNT OF SUBSTANCE

def_unit(['mol', 'mole'], namespace=_ns, prefixes=True,
         doc="mole: amount of a chemical substance in SI.")


###########################################################################
# TEMPERATURE

def_unit(
    ['K', 'Kelvin'], namespace=_ns, prefixes=True,
    doc="Kelvin: temperature with a null point at absolute zero.")
def_unit(
    ['deg_C', 'Celsius'], namespace=_ns, doc='Degrees Celsius',
    format={'latex': r'{}^{\circ}C', 'unicode': '°C'})


###########################################################################
# FORCE

def_unit(['N', 'Newton', 'newton'], kg * m * s ** -2, namespace=_ns,
         prefixes=True, doc="Newton: force")


##########################################################################
# ENERGY

def_unit(['J', 'Joule', 'joule'], N * m, namespace=_ns, prefixes=True,
         doc="Joule: energy")
def_unit(['eV', 'electronvolt'], _si.e.value * J, namespace=_ns, prefixes=True,
         doc="Electron Volt")
def_unit(['Pa', 'Pascal', 'pascal'], J * m ** -3, namespace=_ns, prefixes=True,
         doc="Pascal: pressure")
def_unit(['bar'], 1e5 * Pa, namespace=_ns,
         doc="bar: pressure")


###########################################################################
# POWER

def_unit(['W', 'Watt', 'watt'], J / s, namespace=_ns, prefixes=True,
         doc="Watt: power")


###########################################################################
# ELECTRICAL

def_unit(['A', 'ampere', 'amp'], namespace=_ns, prefixes=True,
         doc="ampere: base unit of electric current in SI")
def_unit(['C', 'coulomb'], A * s, namespace=_ns, prefixes=True,
         doc="coulomb: electric charge")
def_unit(['V', 'Volt', 'volt'], J * C ** -1, namespace=_ns, prefixes=True,
         doc="Volt: electric potential or electromotive force")
def_unit((['Ohm', 'ohm'], ['Ohm']), V * A ** -1, namespace=_ns, prefixes=True,
         doc="Ohm: electrical resistance",
         format={'latex': r'\Omega', 'unicode': 'Ω'})
def_unit(['S', 'Siemens', 'siemens'], A * V ** -1, namespace=_ns,
         prefixes=True, doc="Siemens: electrical conductance")
def_unit(['F', 'Farad', 'farad'], C * V ** -1, namespace=_ns, prefixes=True,
         doc="Farad: electrical capacitance")


###########################################################################
# MAGNETIC

def_unit(['Wb', 'Weber', 'weber'], V * s, namespace=_ns, prefixes=True,
         doc="Weber: magnetic flux")
def_unit(['T', 'Tesla', 'tesla'], Wb * m ** -2, namespace=_ns, prefixes=True,
         doc="Tesla: magnetic flux density")
def_unit(['H', 'Henry', 'henry'], Wb * A ** -1, namespace=_ns, prefixes=True,
         doc="Henry: inductance")


###########################################################################
# ILLUMINATION

def_unit(['cd', 'candela'], namespace=_ns, prefixes=True,
         doc="candela: base unit of luminous intensity in SI")
def_unit(['lm', 'lumen'], cd * sr, namespace=_ns, prefixes=True,
         doc="lumen: luminous flux")
def_unit(['lx', 'lux'], lm * m ** -2, namespace=_ns, prefixes=True,
         doc="lux: luminous emittence")

###########################################################################
# RADIOACTIVITY

def_unit(['Bq', 'becquerel'], Hz, namespace=_ns, prefixes=False,
         doc="becquerel: unit of radioactivity")
def_unit(['Ci', 'curie'], Bq / 3.7e10, namespace=_ns, prefixes=False,
         doc="curie: unit of radioactivity")


###########################################################################
# BASES

bases = set([m, s, kg, A, cd, rad, K, mol])


###########################################################################
# CLEANUP

del UnitBase
del Unit
del def_unit


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary
__doc__ += _generate_unit_summary(globals())
