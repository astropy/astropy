# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines the astrophysics-specific units.  They are also
available in the `astropy.units` namespace.

The `mag` unit is provided for compatibility with the FITS unit string
standard.  However, it is not very useful as-is since it is "orphaned"
and can not be converted to any other unit.  A future astropy
magnitudes library is planned to address this shortcoming.

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import si
from ..constants import si as _si
from .core import UnitBase, def_unit

import numpy as _numpy

UnitBase._set_namespace(globals())

###########################################################################
# LENGTH

def_unit(['AU', 'au'], _si.au.value * si.m, register=True, prefixes=True,
         doc="astronomical unit: approximately the mean Earth--Sun "
         "distance.")

def_unit(['pc', 'parsec'], _si.pc.value * si.m, register=True, prefixes=True,
         doc="parsec: approximately 3.26 light-years.")

def_unit(['solRad', 'R_sun'], _si.R_sun.value * si.m, register=True,
         doc="Solar radius",
         format={'latex': r'R_{\odot}', 'unicode': 'R⊙'})
def_unit(['lyr', 'lightyear'], 9.460730e15 * si.m, register=True,
         doc="Light year")


###########################################################################
# AREAS

def_unit(['barn'], 10 ** -28 * si.m ** 2, register=True, prefixes=True,
         doc="barn: unit of area used in HEP")


###########################################################################
# MASS

def_unit(['solMass', 'M_sun'], _si.M_sun.value * si.kg, register=True,
         prefixes=True, doc="Solar mass",
         format={'latex': r'M_{\odot}', 'unicode': 'M⊙'})
def_unit(['M_p'], _si.m_p.value * si.kg, register=True,
         doc="Proton mass",
         format={'latex': r'M_{p}', 'unicode': 'Mₚ'})
def_unit(['M_e'], _si.m_e.value * si.kg, register=True,
         doc="Electron mass",
         format={'latex': r'M_{e}', 'unicode': 'Mₑ'})
# Unified atomic mass unit
def_unit(['u', 'Da', 'Dalton'], 1.6605387e-27 * si.kg, register=True,
         doc="Unified atomic mass unit")


##########################################################################
# ENERGY

def_unit(['Ry', 'rydberg'], 13.605692 * si.eV, register=True,
         doc="Rydberg: Energy of a photon whose wavenumber is the Rydberg "
         "constant",
         format={'latex': r'R_{\infty}', 'unicode': 'R∞'})


###########################################################################
# ILLUMINATION

def_unit(['solLum', 'L_sun'], _si.L_sun.value * si.W, register=True,
         prefixes=True, doc="Solar luminance",
         format={'latex': r'L_{\odot}', 'unicode': 'L⊙'})


###########################################################################
# SPECTRAL DENSITY

def_unit(['ph', 'photon'], register=True)
def_unit(['Jy', 'Jansky', 'jansky'], 1e-26 * si.W / si.m ** 2 / si.Hz,
         register=True, prefixes=True,
         doc="Jansky: spectral flux density")
def_unit(['R', 'Rayleigh', 'rayleigh'],
         (1e10 / (4 * _numpy.pi)) *
         ph * si.m ** -2 * si.s ** -1 * si.sr ** -1,
         register=True, prefixes=True,
         doc="Rayleigh: photon flux")

def_unit(['mag'], register=True, prefixes=True,
         doc="Stellar magnitude.")


###########################################################################
# MISCELLANEOUS

# Some of these are very FITS-specific and perhaps considered a mistake.
# Maybe they should be moved into the FITS format class?
# TODO: This is defined by the FITS standard as "relative to the sun".
# Is that mass, volume, what?
def_unit(['Sun'], register=True)


###########################################################################
# EVENTS

def_unit(['ct', 'count'], register=True)
def_unit(['pix', 'pixel'], register=True)


###########################################################################
# MISCELLANEOUS

def_unit(['chan'], register=True)
def_unit(['bin'], register=True)
def_unit(['vox', 'voxel'], register=True)
def_unit(['bit'], register=True, prefixes=True)
def_unit(['byte'], register=True, prefixes=True)
def_unit(['adu'], register=True)
def_unit(['beam'], register=True)


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
__doc__ += _generate_unit_summary(globals())
