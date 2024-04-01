# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""

import numpy as np

from astropy.utils import minversion

__all__ = [
    "NUMPY_LT_2_1",
    "None",
    "sanitize_copy_arg",
]

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_2_1 = not minversion(np, "2.1.dev")


def sanitize_copy_arg(copy, /):
    if copy is False:
        return None

    return copy
