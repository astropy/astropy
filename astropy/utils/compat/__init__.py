# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The content of this module is solely for internal use of ``astropy``
and subject to changes without deprecations. Do not use it in external
packages or code.

"""

# Importing this module will also install monkey-patches defined in it
from .numpycompat import *


def __getattr__(attr):
    if attr == "COPY_IF_NEEDED":
        from .numpycompat import COPY_IF_NEEDED  # Will emit warning

        return COPY_IF_NEEDED

    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}.")
