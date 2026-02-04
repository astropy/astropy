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
    "chararray",
    "get_chararray",
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


with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    class chararray(np.char.chararray):
        pass


def get_chararray(obj, itemsize=None, copy=True, unicode=None, order=None):
    if isinstance(obj, (bytes, str)):
        if unicode is None:
            if isinstance(obj, str):
                unicode = True
            else:
                unicode = False

        if itemsize is None:
            itemsize = len(obj)
        shape = len(obj) // itemsize

        return chararray(
            shape, itemsize=itemsize, unicode=unicode, buffer=obj, order=order
        )

    if isinstance(obj, (list, tuple)):
        obj = np.asarray(obj)

    if isinstance(obj, np.ndarray) and issubclass(obj.dtype.type, np.character):
        # If we just have a vanilla chararray, create a chararray
        # view around it.
        if not isinstance(obj, chararray):
            obj = obj.view(chararray)

        if itemsize is None:
            itemsize = obj.itemsize
            # itemsize is in 8-bit chars, so for Unicode, we need
            # to divide by the size of a single Unicode character,
            # which for NumPy is always 4
            if issubclass(obj.dtype.type, np.str_):
                itemsize //= 4

        if unicode is None:
            if issubclass(obj.dtype.type, np.str_):
                unicode = True
            else:
                unicode = False

        if unicode:
            dtype = np.str_
        else:
            dtype = np.bytes_

        if order is not None:
            obj = np.asarray(obj, order=order)
        if (
            copy
            or (itemsize != obj.itemsize)
            or (not unicode and isinstance(obj, np.str_))
            or (unicode and isinstance(obj, np.bytes_))
        ):
            obj = obj.astype((dtype, int(itemsize)))
        return obj

    if isinstance(obj, np.ndarray) and issubclass(obj.dtype.type, object):
        if itemsize is None:
            # Since no itemsize was specified, convert the input array to
            # a list so the ndarray constructor will automatically
            # determine the itemsize for us.
            obj = obj.tolist()
            # Fall through to the default case

    if unicode:
        dtype = np.str_
    else:
        dtype = np.bytes_

    if itemsize is None:
        val = np.array(obj, dtype=dtype, order=order, subok=True)
    else:
        val = np.array(obj, dtype=(dtype, itemsize), order=order, subok=True)
    return val.view(chararray)
