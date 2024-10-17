# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The functions in this module provide faster versions of np.nan*
functions using the optional bottleneck package if it is installed. If
bottleneck is not installed, then the np.nan* functions are used.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np

from astropy.stats.funcs import mad_std
from astropy.utils.exceptions import AstropyDeprecationWarning

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

__all__ = [
    "nansum",  # noqa: F822
    "nanmean",  # noqa: F822
    "nanmedian",  # noqa: F822
    "nanstd",  # noqa: F822
    "nanmadstd",  # noqa: F822
]


def __getattr__(item):
    if item not in __all__:
        raise AttributeError(f"Module {__name__!r} has no attribute {item!r}.")

    if item == "nanmadstd":
        return _nanmadstd

    warnings.warn(
        f"Importing {item} from {__name__} is deprecated, please use numpy directly",
        category=AstropyDeprecationWarning,
        stacklevel=2,
    )
    return getattr(np, item)


def _nanmadstd(
    array: ArrayLike,
    axis: int | tuple[int, ...] | None = None,
) -> float | NDArray:
    """mad_std function that ignores NaNs by default."""
    return mad_std(array, axis=axis, ignore_nan=True)
