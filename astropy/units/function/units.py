# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines units that can also be used as functions of other units.
If called, their arguments are used to initialize the corresponding function
unit (e.g., ``u.mag(u.ct/u.s)``).  Note that the prefixed versions cannot be
called, as it would be unclear what, e.g., ``u.mmag(u.ct/u.s)`` would mean.
"""
from astropy.units.core import _add_prefixes

from .mixin import IrreducibleFunctionUnit, RegularFunctionUnit

_ns = globals()

###########################################################################
# Logarithmic units

# These calls are what core.def_unit would do, but we need to use the callable
# unit versions.  The actual function unit classes get added in logarithmic.

dex = IrreducibleFunctionUnit(
    ["dex"], namespace=_ns, doc="Dex: Base 10 logarithmic unit"
)

dB = RegularFunctionUnit(
    ["dB", "decibel"],
    0.1 * dex,
    namespace=_ns,
    doc="Decibel: ten per base 10 logarithmic unit",
)

mag = RegularFunctionUnit(
    ["mag"],
    -0.4 * dex,
    namespace=_ns,
    doc="Astronomical magnitude: -2.5 per base 10 logarithmic unit",
)

_add_prefixes(mag, namespace=_ns, prefixes=True)

###########################################################################
# CLEANUP

del RegularFunctionUnit
del IrreducibleFunctionUnit

###########################################################################
# DOCSTRING

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    from astropy.units.utils import generate_unit_summary as _generate_unit_summary

    __doc__ += _generate_unit_summary(globals())
