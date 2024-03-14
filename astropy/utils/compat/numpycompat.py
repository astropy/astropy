# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""

import warnings

import numpy as np

from astropy.utils import minversion

__all__ = [
    "NUMPY_LT_1_24",
    "NUMPY_LT_1_25",
    "NUMPY_LT_1_26",
    "NUMPY_LT_2_0",
    "COPY_IF_NEEDED",
    "sanitize_copy_arg",
]

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_1_24 = not minversion(np, "1.24")
NUMPY_LT_1_25 = not minversion(np, "1.25")
NUMPY_LT_1_26 = not minversion(np, "1.26")
NUMPY_LT_2_0 = not minversion(np, "2.0.dev")


COPY_IF_NEEDED = False if NUMPY_LT_2_0 else None


def sanitize_copy_arg(copy, /):
    if not NUMPY_LT_2_0 and copy is False:
        warnings.warn(
            "With copy=False, copies are still made where they cannot be avoided. "
            "In the future, an error will be raised in that scenario. "
            "This is done following a change in numpy 2.0.\n"
            "To silence this warning, pass copy=None, "
            "or downgrade numpy<2.0 (not recommended).",
            category=FutureWarning,
            stacklevel=3,
        )
        return None

    return copy
