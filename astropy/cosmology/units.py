# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Cosmological units and equivalencies."""

from __future__ import annotations

from astropy.units.utils import generate_unit_summary as _generate_unit_summary

__all__ = [
    # redshift equivalencies
    "dimensionless_redshift",
    "littleh",
    "redshift",
    "redshift_distance",
    "redshift_hubble",
    "redshift_temperature",
    # other equivalencies
    "with_H0",
    "with_redshift",
]


from astropy.units import add_enabled_equivalencies as _add_enabled_equivalencies

from ._src.units import littleh, redshift
from ._src.units_equivalencies import (
    dimensionless_redshift,
    redshift_distance,
    redshift_hubble,
    redshift_temperature,
    with_H0,
    with_redshift,
)

# ===================================================================
# Enable the set of default equivalencies.
# If the cosmology package is imported, this is added to the list astropy-wide.
_add_enabled_equivalencies(dimensionless_redshift())


# =============================================================================
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
if __doc__ is not None:
    from ._src.units import _ns

    __doc__ += "\n" + _generate_unit_summary(_ns)
