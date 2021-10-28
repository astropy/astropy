# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Read/Write/Interchange methods for `astropy.cosmology`. **NOT public API**.
"""

# Import the interchange to register them into the io registry.
from . import mapping, table  # noqa: F403
