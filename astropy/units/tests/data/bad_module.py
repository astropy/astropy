# Licensed under a 3-clause BSD style license - see LICENSE.rst

# ruff: noqa: INP001

# The machinery that dynamically creates docstrings for the modules that contain unit
# definitions assumes that all the units are `NamedUnit` instances. If the units are
# created using `def_unit()` then they are, but multiplying or dividing existing units
# creates `CompositeUnit` instances. This module deliberately creates a `CompositeUnit`
# to allow testing that the resulting error message is informative.

from astropy import units as u
from astropy.units.utils import generate_unit_summary

km_per_h = u.km / u.h

__doc__ = generate_unit_summary(globals())
