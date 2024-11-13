# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions for `astropy.cosmology`."""

from .comparison import cosmology_equal
from .optimize import z_at_value

__all__ = ["z_at_value", "cosmology_equal"]
