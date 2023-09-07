# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains helper functions handling metadata."""

import numpy as np

from astropy.utils.misc import dtype_bytes_or_chars

from .exceptions import MergeConflictError

__all__ = ["common_dtype"]


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

    def dtype(arr):
        return getattr(arr, "dtype", np.dtype("O"))

    np_types = (np.bool_, np.object_, np.number, np.character, np.void)
    uniq_types = {
        tuple(issubclass(dtype(arr).type, np_type) for np_type in np_types)
        for arr in arrs
    }
    if len(uniq_types) > 1:
        # Embed into the exception the actual list of incompatible types.
        incompat_types = [dtype(arr).name for arr in arrs]
        tme = MergeConflictError(f"Arrays have incompatible types {incompat_types}")
        tme._incompat_types = incompat_types
        raise tme

    arrs = [np.empty(1, dtype=dtype(arr)) for arr in arrs]

    # For string-type arrays need to explicitly fill in non-zero
    # values or the final arr_common = .. step is unpredictable.
    for i, arr in enumerate(arrs):
        if arr.dtype.kind in ("S", "U"):
            arrs[i] = [
                ("0" if arr.dtype.kind == "U" else b"0")
                * dtype_bytes_or_chars(arr.dtype)
            ]

    arr_common = np.array([arr[0] for arr in arrs])
    return (
        arr_common.dtype.str
        if arr_common.dtype.names is None
        else arr_common.dtype.descr
    )
