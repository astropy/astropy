# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Cosmology-related constants.

This module and its contents adhere to the `Cosmology API`_.
"""

import astropy.constants as _consts

__all__ = [
    "G",  # Gravitational constant
    "c",  # Speed of light
]


G = _consts.G.to("pc km2 s-2 Msun-1")
c = _consts.c.to("km s-1")
