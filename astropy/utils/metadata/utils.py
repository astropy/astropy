# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains helper functions handling metadata."""

import numpy as np

from astropy.utils.decorators import deprecated

from .exceptions import MergeConflictError

__all__ = ["common_dtype", "result_type"]


def dtype(arr):
    return getattr(arr, "dtype", np.dtype("O"))


@deprecated(
    "8.0",
    message=(
        "The {func} {obj_type} is deprecated and may be removed in a "
        "future version. Use {alternative} instead, but note that this "
        "returns a dtype rather than a string."
    ),
    alternative="astropy.utils.metadata.result_type",
)
def common_dtype(arrs):
    """
    Use numpy to find the common dtype for a list of ndarrays.

    Only allow arrays within the following fundamental numpy data types:
    ``np.bool_``, ``np.object_``, ``np.number``, ``np.character``, ``np.void``

    Parameters
    ----------
    arrs : list of ndarray
        Arrays for which to find the common dtype

    Returns
    -------
    dtype_str : str
        String representation of dytpe (dtype ``str`` attribute)
    """
    dt = result_type(arrs)
    return dt.str if dt.names is None else dt.descr


def result_type(arrs):
    """
    Use numpy to find the common type for a list of ndarray.

    The difference with `numpy.result_type` is that all arrays should
    share the same fundamental numpy data type, one of:
    ``np.bool_``, ``np.object_``, ``np.number``, ``np.character``, ``np.void``
    Hence, a mix like integer and string will raise instead of resulting
    in a string type.

    Parameters
    ----------
    arrs : list of ndarray
        Array likes for which to find the common dtype. Anything not
        an array (e.g, |Time|), will be considered an object array.

    Returns
    -------
    dtype : ~numpy.dtype
        The common type.
    """
    dtypes = [dtype(arr) for arr in arrs]
    np_types = (np.bool_, np.object_, np.number, np.character, np.void)
    uniq_types = {
        tuple(issubclass(dt.type, np_type) for np_type in np_types) for dt in dtypes
    }
    if len(uniq_types) > 1:
        # Embed into the exception the actual list of incompatible types.
        incompat_types = [dt.name for dt in dtypes]
        tme = MergeConflictError(f"Arrays have incompatible types {incompat_types}")
        tme._incompat_types = incompat_types
        raise tme

    return np.result_type(*dtypes)
