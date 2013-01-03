# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines the CGS units.  They are also available in the
top-level `astropy.units` namespace.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from ..utils.compat.fractions import Fraction

from . import si
from .core import UnitBase, def_unit


UnitBase._set_namespace(globals())

##########################################################################
# ACCELERATION

def_unit(['Gal', 'gal'], si.cm / si.s ** 2, register=True, prefixes=True,
         doc="Gal: CGS unit of acceleration")


##########################################################################
# ENERGY

# Use CGS definition of erg
def_unit(['erg'], si.g * si.cm ** 2 / si.s ** 2, register=True, prefixes=True,
         doc="erg: CGS unit of energy")


##########################################################################
# FORCE

def_unit(['dyn', 'dyne'], si.g * si.cm / si.s ** 2, register=True,
         prefixes=True,
         doc="dyne: CGS unit of force")


##########################################################################
# PRESSURE

def_unit(['Ba', 'Barye', 'barye'], si.g / (si.cm * si.s ** 2), register=True,
         prefixes=True,
         doc="Barye: CGS unit of pressure")


##########################################################################
# DYNAMIC VISCOSITY

def_unit(['P', 'poise'], si.g / (si.cm * si.s), register=True,
         prefixes=True,
         doc="poise: CGS unit of dynamic viscosity")


##########################################################################
# KINEMATIC VISCOSITY

def_unit(['St', 'stokes'], si.cm ** 2 / si.s, register=True,
         prefixes=True,
         doc="stokes: CGS unit of kinematic viscosity")


##########################################################################
# WAVENUMBER

def_unit(['k', 'Kayser', 'kayser'], si.cm ** -1, register=True,
         prefixes=True,
         doc="kayser: CGS unit of wavenumber")


###########################################################################
# ELECTRICAL

def_unit(['D', 'Debye', 'debye'], Fraction(1, 3) * 1e-29 * si.C * si.m,
         register=True,
         doc="Debye: CGS unit of electric dipole moment")

def_unit(['Fr', 'Franklin', 'statcoulomb', 'statC', 'esu'],
         si.g ** Fraction(1, 2) * si.cm ** Fraction(3, 2) * si.s ** -1,
         register=True)

###########################################################################
# MAGNETIC

def_unit(['G', 'Gauss', 'gauss'], 1e-4 * si.T, register=True, prefixes=True,
         doc="Gauss: CGS unit for magnetic field")


# Mapping from SI units to CGS units
_cgs_bases = {si.meter : si.cm,
              si.joule : erg,
              si.newton : dyne,
              si.pascal : barye,
              si.tesla : gauss,
              si.coulomb : statcoulomb,
              si.kilogram : si.g}

###########################################################################
# CLEANUP

del UnitBase
del def_unit
del si
del Fraction


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary
__doc__ += _generate_unit_summary(globals())
