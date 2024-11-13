# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from astropy.io.fits.hdu.base import BITPIX2DTYPE
from astropy.io.fits.hdu.compressed._tiled_compression import (
    decompress_image_data_section,
)
from astropy.io.fits.hdu.compressed.utils import _n_tiles
from astropy.utils.shapes import simplify_basic_index

__all__ = ["CompImageSection"]


class CompImageSection:
    """
    Class enabling subsets of CompImageHDU data to be loaded lazily via slicing.

    Slices of this object load the corresponding section of an image array from
    the underlying FITS file, and applies any BSCALE/BZERO factors.

    Section slices cannot be assigned to, and modifications to a section are
    not saved back to the underlying file.

    See the :ref:`astropy:data-sections` section of the Astropy documentation
    for more details.
    """

    def __init__(self, hdu):
        self.hdu = hdu

        self._data_shape = self.hdu.shape
        self._tile_shape = self.hdu.tile_shape

        self._n_dim = len(self._data_shape)
        self._n_tiles = np.array(
            _n_tiles(self._data_shape, self._tile_shape), dtype=int
        )

    @property
    def shape(self):
        return tuple(self._data_shape)

    @property
    def ndim(self):
        return len(self.hdu.shape)

    @property
    def dtype(self):
        return np.dtype(BITPIX2DTYPE[self.hdu._bitpix])

    def __getitem__(self, index):
        if self.hdu._bintable is None:
            return self.hdu.data[index]

        # Shortcut if the whole data is requested (this is used by the
        # data property, so we optimize it as it is frequently used)
        if index is Ellipsis:
            first_tile_index = np.zeros(self._n_dim, dtype=int)
            last_tile_index = self._n_tiles - 1
            data = decompress_image_data_section(
                self.hdu._bintable.data,
                self.hdu.compression_type,
                self.hdu._bintable.header,
                self.hdu._bintable,
                self.hdu.header,
                first_tile_index,
                last_tile_index,
            )
            if self.hdu._do_not_scale_image_data:
                return data
            scaled_data = self.hdu._scale_data(data)
            self.hdu._update_header_scale_info(scaled_data.dtype)
            return scaled_data

        index = simplify_basic_index(index, shape=self._data_shape)

        # Determine for each dimension the first and last tile to extract

        first_tile_index = np.zeros(self._n_dim, dtype=int)
        last_tile_index = np.zeros(self._n_dim, dtype=int)

        final_array_index = []

        for dim, idx in enumerate(index):
            if isinstance(idx, slice):
                if idx.step > 0:
                    first_tile_index[dim] = idx.start // self._tile_shape[dim]
                    last_tile_index[dim] = (idx.stop - 1) // self._tile_shape[dim]
                else:
                    stop = 0 if idx.stop is None else max(idx.stop - 1, 0)
                    first_tile_index[dim] = stop // self._tile_shape[dim]
                    last_tile_index[dim] = idx.start // self._tile_shape[dim]

                # Because slices such as slice(5, 0, 1) can exist (which
                # would be empty) we need to make sure last_tile_index is
                # always larger than first_tile_index
                last_tile_index = np.maximum(last_tile_index, first_tile_index)

                if idx.step < 0 and idx.stop is None:
                    final_array_index.append(idx)
                else:
                    final_array_index.append(
                        slice(
                            idx.start - self._tile_shape[dim] * first_tile_index[dim],
                            idx.stop - self._tile_shape[dim] * first_tile_index[dim],
                            idx.step,
                        )
                    )
            else:
                first_tile_index[dim] = idx // self._tile_shape[dim]
                last_tile_index[dim] = first_tile_index[dim]
                final_array_index.append(
                    idx - self._tile_shape[dim] * first_tile_index[dim]
                )

        data = decompress_image_data_section(
            self.hdu._bintable.data,
            self.hdu.compression_type,
            self.hdu._bintable.header,
            self.hdu._bintable,
            self.hdu.header,
            first_tile_index,
            last_tile_index,
        )

        data = data[tuple(final_array_index)]

        if self.hdu._do_not_scale_image_data:
            return data

        scaled_data = self.hdu._scale_data(data)
        self.hdu._update_header_scale_info(scaled_data.dtype)
        return scaled_data
