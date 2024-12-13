# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions for `astropy.cosmology`."""

__all__ = ["cosmology_equal", "z_at_value"]

from .comparison import cosmology_equal
from .optimize import z_at_value
