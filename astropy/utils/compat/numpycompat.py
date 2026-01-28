# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""

import warnings

import numpy as np

from astropy.utils import minversion
from astropy.utils.exceptions import AstropyPendingDeprecationWarning

__all__ = [
    "NUMPY_LT_2_1",
    "NUMPY_LT_2_2",
    "NUMPY_LT_2_3",
    "NUMPY_LT_2_4",
    "NUMPY_LT_2_4_1",
    "NUMPY_LT_2_5",
    "is_chararray",
]

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_2_1 = not minversion(np, "2.1.0.dev")
NUMPY_LT_2_2 = not minversion(np, "2.2.0.dev0")
NUMPY_LT_2_3 = not minversion(np, "2.3.0.dev0")
NUMPY_LT_2_4 = not minversion(np, "2.4.0.dev0")
NUMPY_LT_2_4_1 = not minversion(np, "2.4.1.dev0")
NUMPY_LT_2_5 = not minversion(np, "2.5.0.dev0")


def is_chararray(obj, /) -> bool:
    if not hasattr(np, "char"):
        return False

    with warnings.catch_warnings():
        # numpy.char.chararray is deprecated in not NUMPY_LT_2_5
        warnings.simplefilter("ignore")
        from numpy.char import chararray

    return isinstance(obj, chararray)


def __getattr__(attr):
    # MHvK: Added in 8.0. Regular deprecation in 9.0, remove in 10.0?
    if attr == "COPY_IF_NEEDED":
        warnings.warn(
            "COPY_IF_NEEDED is no longer needed now that astropy only "
            "supports numpy >= 2. It can be safely replaced with 'None'. ",
            AstropyPendingDeprecationWarning,
        )
        return None

    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}.")
