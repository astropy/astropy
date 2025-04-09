# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines miscellaneous units. They are also available in
(and should be used through) the `astropy.units` namespace.
"""
# avoid ruff complaints about undefined names defined by def_unit
# ruff: noqa: F821

import numpy as np

from astropy.constants import si as _si

from . import si
from .core import UnitBase, binary_prefixes, def_unit, si_prefixes
from .utils import generate_unit_summary

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
def_unit(
    ["u", "Da", "Dalton"],
    _si.u,
    namespace=_ns,
    prefixes=True,
    exclude_prefixes=["a", "da"],
    doc="Unified atomic mass unit",
)

##########################################################################
# ENERGY

def_unit(
    ["eV", "electronvolt"],
    _si.e.value * si.J,
    namespace=_ns,
    prefixes=True,
    doc="Electron Volt",
)

# Here, explicitly convert the planck constant to 'eV s' since the constant
# can override that to give a more precise value that takes into account
# covariances between e and h.  Eventually, this may also be replaced with
# just `_si.Ryd.to(eV)`.
def_unit(
    ["Ry", "rydberg"],
    (_si.Ryd * _si.c * _si.h.to(eV * si.s)).to(eV),
    namespace=_ns,
    prefixes=True,
    doc="Rydberg: Energy of a photon whose wavenumber is the Rydberg constant",
    format={"latex": r"R_{\infty}", "unicode": "R∞"},
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
    __doc__ += generate_unit_summary(globals())
