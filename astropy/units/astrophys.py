# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines the astrophysics-specific units.  They are also
available in the `astropy.units` namespace.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import si
from ..constants import si as _si
from .core import (UnitBase, def_unit, si_prefixes, binary_prefixes,
                   set_enabled_units)

# To ensure si units of the constants can be interpreted.
set_enabled_units([si])

import numpy as _numpy

_ns = globals()

###########################################################################
# LENGTH

def_unit((['AU', 'au'], ['astronomical_unit']), _si.au, namespace=_ns, prefixes=True,
         doc="astronomical unit: approximately the mean Earth--Sun "
         "distance.")

def_unit(['pc', 'parsec'], _si.pc, namespace=_ns, prefixes=True,
         doc="parsec: approximately 3.26 light-years.")

def_unit(['solRad', 'R_sun', 'Rsun'], _si.R_sun, namespace=_ns,
         doc="Solar radius", prefixes=False,
         format={'latex': r'R_{\odot}', 'unicode': 'R⊙'})
def_unit(['jupiterRad', 'R_jup', 'Rjup', 'R_jupiter', 'Rjupiter'],
         _si.R_jup, namespace=_ns, prefixes=False, doc="Jupiter radius",
         # LaTeX jupiter symbol requires wasysym
         format={'latex': r'R_{\rm J}', 'unicode': 'R♃'})
def_unit(['earthRad', 'R_earth', 'Rearth'], _si.R_earth, namespace=_ns,
         prefixes=False, doc="Earth radius",
         # LaTeX earth symbol requires wasysym
         format={'latex': r'R_{\oplus}', 'unicode': 'R⊕'})

def_unit(['lyr', 'lightyear'], (_si.c * si.yr).to(si.m),
         namespace=_ns, prefixes=True, doc="Light year")


###########################################################################
# AREAS

def_unit(['barn', 'barn'], 10 ** -28 * si.m ** 2, namespace=_ns, prefixes=True,
         doc="barn: unit of area used in HEP")


###########################################################################
# ANGULAR MEASUREMENTS

def_unit(['cycle', 'cy'], 2.0 * _numpy.pi * si.rad,
         namespace=_ns, prefixes=False,
         doc="cycle: angular measurement, a full turn or rotation")

###########################################################################
# MASS

def_unit(['solMass', 'M_sun', 'Msun'], _si.M_sun, namespace=_ns,
         prefixes=False, doc="Solar mass",
         format={'latex': r'M_{\odot}', 'unicode': 'M⊙'})
def_unit(['jupiterMass', 'M_jup', 'Mjup', 'M_jupiter', 'Mjupiter'],
         _si.M_jup, namespace=_ns, prefixes=False, doc="Jupiter mass",
         # LaTeX jupiter symbol requires wasysym
         format={'latex': r'M_{\rm J}', 'unicode': 'M♃'})
def_unit(['earthMass', 'M_earth', 'Mearth'], _si.M_earth, namespace=_ns,
         prefixes=False, doc="Earth mass",
         # LaTeX earth symbol requires wasysym
         format={'latex': r'M_{\oplus}', 'unicode': 'M⊕'})
def_unit(['M_p'], _si.m_p, namespace=_ns, doc="Proton mass",
         format={'latex': r'M_{p}', 'unicode': 'Mₚ'})
def_unit(['M_e'], _si.m_e, namespace=_ns, doc="Electron mass",
         format={'latex': r'M_{e}', 'unicode': 'Mₑ'})
# Unified atomic mass unit
def_unit(['u', 'Da', 'Dalton'], _si.u, namespace=_ns,
         prefixes=True, exclude_prefixes=['a', 'da'],
         doc="Unified atomic mass unit")

##########################################################################
# ENERGY

# Here, explicitly convert the planck constant to 'eV s' since the constant
# can override that to give a more precise value that takes into account
# covariances between e and h.  Eventually, this may also be replaced with
# just `_si.Ryd.to(eV)`.
def_unit(['Ry', 'rydberg'],
         (_si.Ryd * _si.c * _si.h.to(si.eV * si.s)).to(si.eV),
         namespace=_ns, prefixes=True,
         doc="Rydberg: Energy of a photon whose wavenumber is the Rydberg "
         "constant",
         format={'latex': r'R_{\infty}', 'unicode': 'R∞'})


###########################################################################
# ILLUMINATION

def_unit(['solLum', 'L_sun', 'Lsun'], _si.L_sun, namespace=_ns,
         prefixes=False, doc="Solar luminance",
         format={'latex': r'L_{\odot}', 'unicode': 'L⊙'})


###########################################################################
# SPECTRAL DENSITY

def_unit((['ph', 'photon'], ['photon']),
         format={'ogip': 'photon', 'vounit': 'photon'},
         namespace=_ns, prefixes=True)
def_unit(['Jy', 'Jansky', 'jansky'], 1e-26 * si.W / si.m ** 2 / si.Hz,
         namespace=_ns, prefixes=True,
         doc="Jansky: spectral flux density")
def_unit(['R', 'Rayleigh', 'rayleigh'],
         (1e10 / (4 * _numpy.pi)) *
         ph * si.m ** -2 * si.s ** -1 * si.sr ** -1,
         namespace=_ns, prefixes=True,
         doc="Rayleigh: photon flux")


###########################################################################
# MISCELLANEOUS

# Some of these are very FITS-specific and perhaps considered a mistake.
# Maybe they should be moved into the FITS format class?
# TODO: This is defined by the FITS standard as "relative to the sun".
# Is that mass, volume, what?
def_unit(['Sun'], namespace=_ns)


###########################################################################
# EVENTS

def_unit((['ct', 'count'], ['count']),
         format={'fits': 'count', 'ogip': 'count', 'vounit': 'count'},
         namespace=_ns, prefixes=True, exclude_prefixes=['p'])
def_unit((['pix', 'pixel'], ['pixel']),
         format={'ogip': 'pixel', 'vounit': 'pixel'},
         namespace=_ns, prefixes=True)


###########################################################################
# MISCELLANEOUS

def_unit(['chan'], namespace=_ns, prefixes=True)
def_unit(['bin'], namespace=_ns, prefixes=True)
def_unit((['vox', 'voxel'], ['voxel']),
         format={'fits': 'voxel', 'ogip': 'voxel', 'vounit': 'voxel'},
         namespace=_ns, prefixes=True)
def_unit((['bit', 'b'], ['bit']), namespace=_ns,
         prefixes=si_prefixes + binary_prefixes)
def_unit((['byte', 'B'], ['byte']), 8 * bit, namespace=_ns,
         format={'vounit': 'byte'},
         prefixes=si_prefixes + binary_prefixes,
         exclude_prefixes=['d'])
def_unit(['adu'], namespace=_ns, prefixes=True)
def_unit(['beam'], namespace=_ns, prefixes=True)
def_unit(['electron'], doc="Number of electrons", namespace=_ns,
         format={'latex': r'e^{-}', 'unicode': 'e⁻'})


###########################################################################
# CLEANUP

del UnitBase
del def_unit
del si


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary
if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())
