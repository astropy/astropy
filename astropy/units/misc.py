# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines miscellaneous units. They are also
available in the `astropy.units` namespace.
"""

import numpy as np

from astropy.constants import si as _si

from . import si
from .core import UnitBase, binary_prefixes, def_unit, set_enabled_units, si_prefixes

# To ensure si units of the constants can be interpreted.
set_enabled_units([si])


__all__: list[str] = []  #  Units are added at the end

_ns = globals()

###########################################################################
# AREAS

def_unit(
    ["barn", "barn"],
    10**-28 * si.m**2,
    namespace=_ns,
    prefixes=True,
    doc="barn: unit of area used in HEP",
)


###########################################################################
# ANGULAR MEASUREMENTS

def_unit(
    ["cycle", "cy"],
    2.0 * np.pi * si.rad,
    namespace=_ns,
    prefixes=False,
    doc="cycle: angular measurement, a full turn or rotation",
)
def_unit(
    ["spat", "sp"],
    4.0 * np.pi * si.sr,
    namespace=_ns,
    prefixes=False,
    doc="spat: the solid angle of the sphere, 4pi sr",
)

##########################################################################
# PRESSURE

def_unit(
    ["bar"],
    1e5 * si.Pa,
    namespace=_ns,
    prefixes=[(["m"], ["milli"], 1.0e-3)],
    doc="bar: pressure",
)
# The torr is almost the same as mmHg but not quite.
# See https://en.wikipedia.org/wiki/Torr
# Define the unit here despite it not being an astrophysical unit.
# It may be moved if more similar units are created later.
def_unit(
    ["Torr", "torr"],
    _si.atm.value / 760.0 * si.Pa,
    namespace=_ns,
    prefixes=[(["m"], ["milli"], 1.0e-3)],
    doc=(
        "Unit of pressure based on an absolute scale, now defined as "
        "exactly 1/760 of a standard atmosphere"
    ),
)

###########################################################################
# MASS

def_unit(
    ["M_p"],
    _si.m_p,
    namespace=_ns,
    doc="Proton mass",
    format={"latex": r"M_{p}", "unicode": "Mₚ"},
)
def_unit(
    ["M_e"],
    _si.m_e,
    namespace=_ns,
    doc="Electron mass",
    format={"latex": r"M_{e}", "unicode": "Mₑ"},
)
# Unified atomic mass unit
def_unit(
    ["u", "Da", "Dalton"],
    _si.u,
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["a", "da"],
    doc="Unified atomic mass unit",
)


###########################################################################
# COMPUTER

def_unit(
    (["bit", "b"], ["bit"]),
    namespace=_ns,
    prefixes=si_prefixes + binary_prefixes,
)
def_unit(
    (["byte", "B"], ["byte"]),
    8 * bit,
    namespace=_ns,
    format={"vounit": "byte"},
    prefixes=si_prefixes + binary_prefixes,
    exclude_prefixes=["d"],
)
def_unit(
    (["pix", "pixel"], ["pixel"]),
    format={"ogip": "pixel", "vounit": "pixel"},
    namespace=_ns,
    prefixes=True,
)
def_unit(
    (["vox", "voxel"], ["voxel"]),
    format={"fits": "voxel", "ogip": "voxel", "vounit": "voxel"},
    namespace=_ns,
    prefixes=True,
)

###########################################################################
# ALL & DOCSTRING

__all__ += [n for n, v in _ns.items() if isinstance(v, UnitBase)]


if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    from .utils import generate_unit_summary as _generate_unit_summary

    __doc__ += _generate_unit_summary(globals())
