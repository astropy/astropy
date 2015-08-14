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
from .core import UnitBase, def_unit, si_prefixes, binary_prefixes

import numpy as _numpy

_ns = globals()

###########################################################################
# LENGTH

def_unit((['AU', 'au'], ['astronomical_unit']), _si.au.value * si.m, namespace=_ns, prefixes=True,
         doc="astronomical unit: approximately the mean Earth--Sun "
         "distance.")

def_unit(['pc', 'parsec'], _si.pc.value * si.m, namespace=_ns, prefixes=True,
         doc="parsec: approximately 3.26 light-years.")

def_unit(['solRad', 'R_sun', 'Rsun'], _si.R_sun.value * si.m, namespace=_ns,
         doc="Solar radius", prefixes=True,
         format={'latex': r'R_{\odot}', 'unicode': 'R⊙'})
def_unit(['lyr', 'lightyear'], _si.c.value * si.yr.to(si.s) * si.m,
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

def_unit(['solMass', 'M_sun', 'Msun'], _si.M_sun.value * si.kg, namespace=_ns,
         prefixes=True, doc="Solar mass",
         format={'latex': r'M_{\odot}', 'unicode': 'M⊙'})
def_unit(['M_p'], _si.m_p.value * si.kg, namespace=_ns,
         doc="Proton mass",
         format={'latex': r'M_{p}', 'unicode': 'Mₚ'})
def_unit(['M_e'], _si.m_e.value * si.kg, namespace=_ns,
         doc="Electron mass",
         format={'latex': r'M_{e}', 'unicode': 'Mₑ'})
# Unified atomic mass unit
def_unit(['u', 'Da', 'Dalton'], 1.6605387e-27 * si.kg, namespace=_ns,
         prefixes=True, exclude_prefixes=['a', 'da'],
         doc="Unified atomic mass unit")


##########################################################################
# ENERGY

def_unit(['Ry', 'rydberg'], 13.605692 * si.eV, namespace=_ns,
         prefixes=True,
         doc="Rydberg: Energy of a photon whose wavenumber is the Rydberg "
         "constant",
         format={'latex': r'R_{\infty}', 'unicode': 'R∞'})


###########################################################################
# ILLUMINATION

def_unit(['solLum', 'L_sun', 'Lsun'], _si.L_sun.value * si.W, namespace=_ns,
         prefixes=True, doc="Solar luminance",
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
