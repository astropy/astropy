# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import functools
from typing import TYPE_CHECKING

import numpy as np

from astropy.stats.funcs import mad_std
from astropy.units import Quantity
from astropy.utils.compat.optional_deps import HAS_BOTTLENECK

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import ArrayLike, NDArray


if HAS_BOTTLENECK:
    import bottleneck

    def _move_tuple_axes_last(
        array: ArrayLike,
        axis: tuple[int, ...] | None = None,
    ) -> ArrayLike:
        """
        Move the specified axes of a NumPy array to the last positions
        and combine them.

        Bottleneck can only take integer axis, not tuple, so this
        function takes all the axes to be operated on and combines them
        into the last dimension of the array so that we can then use
        axis=-1.

        Parameters
        ----------
        array : `~numpy.ndarray`
            The input array.

        axis : tuple of int
            The axes on which to move and combine.

        Returns
        -------
        array_new : `~numpy.ndarray`
            Array with the axes being operated on moved into the last
            dimension.
        """
        other_axes = tuple(i for i in range(array.ndim) if i not in axis)

        # Move the specified axes to the last positions
        array_new = np.transpose(array, other_axes + axis)

        # Reshape the array by combining the moved axes
        return array_new.reshape(array_new.shape[: len(other_axes)] + (-1,))

    def _apply_bottleneck(
        function: Callable,
        array: ArrayLike,
        axis: int | tuple[int, ...] | None = None,
        **kwargs,
    ) -> float | NDArray | Quantity:
        """Wrap bottleneck function to handle tuple axis.

        Also takes care to ensure the output is of the expected type,
        i.e., a quantity, numpy array, or numpy scalar.
        """
        if isinstance(axis, tuple):
            array = _move_tuple_axes_last(array, axis=axis)
            axis = -1

        result = function(array, axis=axis, **kwargs)
        if isinstance(array, Quantity):
            return array.__array_wrap__(result)
        elif isinstance(result, float):
            # For compatibility with numpy, always return a numpy scalar.
            return np.float64(result)
        else:
            return result

    nanmean = functools.partial(_apply_bottleneck, bottleneck.nanmean)
    nanmedian = functools.partial(_apply_bottleneck, bottleneck.nanmedian)
    nanstd = functools.partial(_apply_bottleneck, bottleneck.nanstd)

else:
    nanmean = np.nanmean
    nanmedian = np.nanmedian
    nanstd = np.nanstd


def nanmadstd(
    array: ArrayLike,
    axis: int | tuple[int, ...] | None = None,
) -> float | NDArray:
    """mad_std function that ignores NaNs by default."""
    return mad_std(array, axis=axis, ignore_nan=True)
