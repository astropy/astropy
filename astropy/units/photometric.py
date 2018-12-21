# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines magnitude zero points and related photometric quantities.

The corresponding magnitudes are given in the description of each unit
(the actual definitions are in `~astropy.units.function.logarithmic`).
"""


import numpy as _numpy
from .core import UnitBase, def_unit, Unit

from astropy.constants import si as _si
from . import cgs, si, astrophys


_ns = globals()

def_unit(['Bol', 'L_bol'], _si.L_bol0, namespace=_ns, prefixes=False,
         doc="Luminosity corresponding to absolute bolometric magnitude zero "
         "(magnitude ``M_bol``).")
def_unit(['bol', 'f_bol'], _si.L_bol0 / (4 * _numpy.pi * (10.*astrophys.pc)**2),
         namespace=_ns, prefixes=False, doc="Irradiance corresponding to "
         "appparent bolometric magnitude zero (magnitude ``m_bol``).")
def_unit(['AB', 'ABflux'], 10.**(48.6/-2.5) * cgs.erg * cgs.cm**-2 / si.s / si.Hz,
         namespace=_ns, prefixes=False,
         doc="AB magnitude zero flux density (magnitude ``ABmag``).")
def_unit(['ST', 'STflux'], 10.**(21.1/-2.5) * cgs.erg * cgs.cm**-2 / si.s / si.AA,
         namespace=_ns, prefixes=False,
         doc="ST magnitude zero flux density (magnitude ``STmag``).")

def_unit(['mgy', 'maggy'],
         namespace=_ns, prefixes=[(['n'], ['nano'], 1e-9)],
         doc="Maggies - a linear flux unit that is the flux for a mag=0 object."
             "To tie this onto a specific calibrated unit system, the "
             "zero_point_flux equivalency should be used.")


def zero_point_flux(flux0):
    """
    An equivalency for converting linear flux units ("maggys") defined relative
    to a standard source into a standardized system.

    Parameters
    ----------
    flux0 : u.Quantity
        The flux of a magnitude-0 object in the "maggy" system.
    """
    flux_unit0 = Unit(flux0)
    return [(maggy, flux_unit0)]


###########################################################################
# CLEANUP

del UnitBase
del def_unit
del cgs, si, astrophys

###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary
if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())
