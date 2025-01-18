# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Cosmological units."""

__all__ = ["littleh", "redshift"]

import astropy.units as u

_ns = globals()


# This is not formally a unit, but is used in that way in many contexts, and
# an appropriate equivalency is only possible if it's treated as a unit.
redshift = u.def_unit(
    ["redshift"],
    prefixes=False,
    namespace=_ns,
    doc="Cosmological redshift.",
    format={"latex": r""},
)
u.def_physical_type(redshift, "redshift")

# This is not formally a unit, but is used in that way in many contexts, and
# an appropriate equivalency is only possible if it's treated as a unit (see
# https://arxiv.org/pdf/1308.4150.pdf for more)
# Also note that h or h100 or h_100 would be a better name, but they either
# conflict or have numbers in them, which is disallowed
littleh = u.def_unit(
    ["littleh"],
    namespace=_ns,
    prefixes=False,
    doc='Reduced/"dimensionless" Hubble constant',
    format={"latex": r"h_{100}"},
)
