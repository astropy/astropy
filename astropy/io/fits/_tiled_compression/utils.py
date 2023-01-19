import numpy as np


def _iter_array_tiles(data_shape, tile_shape):
    """
    Given an array shape and a tile shape, iterate over the tiles in the array
    returning at each iteration the slices for the array.
    """

    ndim = len(data_shape)
    istart = np.zeros(ndim, dtype=int)

    while True:

        # In the following, we don't need to special case tiles near the edge
        # as Numpy will automatically ignore parts of the slices that are out
        # of bounds.
        tile_slices = tuple(
            [
                slice(istart[idx], istart[idx] + tile_shape[idx])
                for idx in range(len(istart))
            ]
        )

        yield tile_slices

        istart[-1] += tile_shape[-1]

        for idx in range(ndim - 1, 0, -1):
            if istart[idx] >= data_shape[idx]:
                istart[idx] = 0
                istart[idx - 1] += tile_shape[idx - 1]

        if istart[0] >= data_shape[0]:
            break
