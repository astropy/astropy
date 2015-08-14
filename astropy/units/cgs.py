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


_ns = globals()

def_unit(['cm', 'centimeter'], si.cm, namespace=_ns, prefixes=False)
g = si.g
s = si.s
C = si.C
rad = si.rad
sr = si.sr
cd = si.cd
K = si.K
deg_C = si.deg_C
mol = si.mol


##########################################################################
# ACCELERATION

def_unit(['Gal', 'gal'], cm / s ** 2, namespace=_ns, prefixes=True,
         doc="Gal: CGS unit of acceleration")


##########################################################################
# ENERGY

# Use CGS definition of erg
def_unit(['erg'], g * cm ** 2 / s ** 2, namespace=_ns, prefixes=True,
         doc="erg: CGS unit of energy")


##########################################################################
# FORCE

def_unit(['dyn', 'dyne'], g * cm / s ** 2, namespace=_ns,
         prefixes=True,
         doc="dyne: CGS unit of force")


##########################################################################
# PRESSURE

def_unit(['Ba', 'Barye', 'barye'], g / (cm * s ** 2), namespace=_ns,
         prefixes=True,
         doc="Barye: CGS unit of pressure")


##########################################################################
# DYNAMIC VISCOSITY

def_unit(['P', 'poise'], g / (cm * s), namespace=_ns,
         prefixes=True,
         doc="poise: CGS unit of dynamic viscosity")


##########################################################################
# KINEMATIC VISCOSITY

def_unit(['St', 'stokes'], cm ** 2 / s, namespace=_ns,
         prefixes=True,
         doc="stokes: CGS unit of kinematic viscosity")


##########################################################################
# WAVENUMBER

def_unit(['k', 'Kayser', 'kayser'], cm ** -1, namespace=_ns,
         prefixes=True,
         doc="kayser: CGS unit of wavenumber")


###########################################################################
# ELECTRICAL

def_unit(['D', 'Debye', 'debye'], Fraction(1, 3) * 1e-29 * C * si.m,
         namespace=_ns, prefixes=True,
         doc="Debye: CGS unit of electric dipole moment")

def_unit(['Fr', 'Franklin', 'statcoulomb', 'statC', 'esu'],
         g ** Fraction(1, 2) * cm ** Fraction(3, 2) * s ** -1,
         namespace=_ns,
         doc='Franklin: CGS (ESU) unit of charge')

def_unit(['statA', 'statampere'], Fr * s ** -1, namespace=_ns,
         doc='statampere: CGS (ESU) unit of current')

def_unit(['Bi', 'Biot', 'abA', 'abampere', 'emu'],
         g ** Fraction(1, 2) * cm ** Fraction(1, 2) * s ** -1, namespace=_ns,
         doc='Biot: CGS (EMU) unit of current')

def_unit(['abC', 'abcoulomb'], Bi * s, namespace=_ns,
         doc='abcoulomb: CGS (EMU) of charge')

###########################################################################
# MAGNETIC

def_unit(['G', 'Gauss', 'gauss'], 1e-4 * si.T, namespace=_ns, prefixes=True,
         doc="Gauss: CGS unit for magnetic field")


###########################################################################
# BASES

bases = set([cm, g, s, C, rad, cd, K, mol])


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
if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())
