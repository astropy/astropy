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


def _tile_index_to_row_index(tile_index, n_tiles):
    row_index = tile_index[0]
    for dim in range(1, len(tile_index)):
        row_index = tile_index[dim] + row_index * n_tiles[dim]
    return row_index


def _tile_index_to_tile_slices(tile_index, tile_shape):
    return tuple(
        [
            slice(
                tile_index[idx] * tile_shape[idx],
                (tile_index[idx] + 1) * tile_shape[idx],
            )
            for idx in range(len(tile_index))
        ]
    )


def _n_tiles(data_shape, tile_shape):
    return tuple(int(np.ceil(d / t)) for d, t in zip(data_shape, tile_shape))


def _iter_array_tiles_subset(data_shape, tile_shape, first_tile_index, last_tile_index):
    ndim = len(tile_shape)

    n_tiles = _n_tiles(data_shape, tile_shape)

    first_tile_index = np.asarray(first_tile_index, dtype=int)
    last_tile_index = np.asarray(last_tile_index, dtype=int)

    tile_index = first_tile_index.copy()

    while True:
        tile_slices = _tile_index_to_tile_slices(
            tile_index - first_tile_index, tile_shape
        )

        yield _tile_index_to_row_index(tile_index, n_tiles), tile_slices

        tile_index[-1] += 1

        for idx in range(ndim - 1, 0, -1):
            if tile_index[idx] > last_tile_index[idx]:
                tile_index[idx] = first_tile_index[idx]
                tile_index[idx - 1] += 1

        if tile_index[0] > last_tile_index[0]:
            break


def _tile_shape(header):
    return tuple(header[f"ZTILE{idx}"] for idx in range(header["ZNAXIS"], 0, -1))


def _data_shape(header):
    return tuple(header[f"ZNAXIS{idx}"] for idx in range(header["ZNAXIS"], 0, -1))
