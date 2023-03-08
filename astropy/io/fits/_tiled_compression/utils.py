from math import ceil

import numpy as np


def _n_tiles(data_shape, tile_shape):
    return [int(ceil(d / t)) for d, t in zip(data_shape, tile_shape)]


def _iter_array_tiles(
    data_shape, tile_shape, first_tile_index=None, last_tile_index=None
):
    ndim = len(tile_shape)

    n_tiles = _n_tiles(data_shape, tile_shape)

    if first_tile_index is None:
        first_tile_index = (0,) * ndim

    if last_tile_index is None:
        last_tile_index = tuple(n - 1 for n in n_tiles)

    tile_index = list(first_tile_index)

    while True:
        tile_slices = tuple(
            slice(
                (tile_index[idx] - first_tile_index[idx]) * tile_shape[idx],
                (tile_index[idx] - first_tile_index[idx] + 1) * tile_shape[idx],
            )
            for idx in range(ndim)
        )

        row_index = tile_index[0]
        for dim in range(1, ndim):
            row_index = tile_index[dim] + row_index * n_tiles[dim]

        yield row_index, tile_slices

        tile_index[-1] += 1

        for idx in range(ndim - 1, 0, -1):
            if tile_index[idx] > last_tile_index[idx]:
                tile_index[idx] = first_tile_index[idx]
                tile_index[idx - 1] += 1

        if tile_index[0] > last_tile_index[0]:
            break


def _tile_shape(header):
    return np.array(
        [header[f"ZTILE{idx}"] for idx in range(header["ZNAXIS"], 0, -1)], dtype=int
    )


def _data_shape(header):
    return np.array(
        [header[f"ZNAXIS{idx}"] for idx in range(header["ZNAXIS"], 0, -1)], dtype=int
    )
