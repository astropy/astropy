# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines magnitude zero points and related photometric quantities.

The corresponding magnitudes are given in the description of each unit.
(the actual definitions are in `astropy.units.function.units`).

Both the units and magnitudes are available in (and should be used
through) the `astropy.units` namespace.

"""

import numpy as np

from astropy.constants.si import L_bol0

from . import astrophys, cgs, si
from .core import UnitBase, def_unit
from .docgen import generate_unit_summary

__all__ = []  #  Units are added at the end

_ns = globals()

def_unit(
    ["Bol", "L_bol"],
    L_bol0,
    namespace=_ns,
    prefixes=False,
    doc=(
        "Luminosity corresponding to absolute bolometric magnitude zero "
        "(magnitude ``M_bol``)."
    ),
)
def_unit(
    ["bol", "f_bol"],
    L_bol0 / (4 * np.pi * (10.0 * astrophys.pc) ** 2),
    namespace=_ns,
    prefixes=False,
    doc=(
        "Irradiance corresponding to apparent bolometric magnitude zero "
        "(magnitude ``m_bol``)."
    ),
)
def_unit(
    ["AB", "ABflux"],
    10.0 ** (48.6 / -2.5) * cgs.erg * cgs.cm**-2 / si.s / si.Hz,
    namespace=_ns,
    prefixes=False,
    doc="AB magnitude zero flux density (magnitude ``ABmag``).",
)
def_unit(
    ["ST", "STflux"],
    10.0 ** (21.1 / -2.5) * cgs.erg * cgs.cm**-2 / si.s / si.AA,
    namespace=_ns,
    prefixes=False,
    doc="ST magnitude zero flux density (magnitude ``STmag``).",
)
def_unit(
    ["mgy", "maggy"],
    namespace=_ns,
    prefixes=[(["n"], ["nano"], 1e-9)],
    doc=(
        "Maggies - a linear flux unit that is the flux for a mag=0 object."
        "To tie this onto a specific calibrated unit system, the "
        "zero_point_flux equivalency should be used."
    ),
)


###########################################################################
# ALL & DOCSTRING

__all__ += [n for n, v in _ns.items() if isinstance(v, UnitBase)]


if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    __doc__ += generate_unit_summary(globals())
