# Licensed under a 3-clause BSD style license - see PYFITS.rst

import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning

__all__ = ["_validate_tile_shape"]


def _validate_tile_shape(*, tile_shape, compression_type, image_header):
    naxis = image_header["NAXIS"]

    if not tile_shape:
        tile_shape = []
    elif len(tile_shape) != naxis:
        warnings.warn(
            "Provided tile size not appropriate for the data.  "
            "Default tile size will be used.",
            AstropyUserWarning,
        )
        tile_shape = []
    else:
        tile_shape = list(tile_shape)

    # Set default tile dimensions for HCOMPRESS_1

    if compression_type == "HCOMPRESS_1":
        if image_header["NAXIS1"] < 4 or image_header["NAXIS2"] < 4:
            raise ValueError("Hcompress minimum image dimension is 4 pixels")
        elif tile_shape:
            if tile_shape[-1] < 4 or tile_shape[-2] < 4:
                # user specified tile size is too small
                raise ValueError("Hcompress minimum tile dimension is 4 pixels")
            major_dims = len([ts for ts in tile_shape if ts > 1])
            if major_dims > 2:
                raise ValueError(
                    "HCOMPRESS can only support 2-dimensional tile sizes."
                    "All but two of the tile_shape dimensions must be set "
                    "to 1."
                )

        if tile_shape and (tile_shape[-1] == 0 and tile_shape[-2] == 0):
            # compress the whole image as a single tile
            tile_shape[-1] = image_header["NAXIS1"]
            tile_shape[-2] = image_header["NAXIS2"]

            for i in range(2, naxis):
                # set all higher tile dimensions = 1
                tile_shape[i] = 1
        elif not tile_shape:
            # The Hcompress algorithm is inherently 2D in nature, so the
            # row by row tiling that is used for other compression
            # algorithms is not appropriate.  If the image has less than 30
            # rows, then the entire image will be compressed as a single
            # tile.  Otherwise the tiles will consist of 16 rows of the
            # image.  This keeps the tiles to a reasonable size, and it
            # also includes enough rows to allow good compression
            # efficiency.  It the last tile of the image happens to contain
            # less than 4 rows, then find another tile size with between 14
            # and 30 rows (preferably even), so that the last tile has at
            # least 4 rows.

            # 1st tile dimension is the row length of the image
            tile_shape = [image_header["NAXIS1"]]

            if image_header["NAXIS2"] <= 30:
                tile_shape.insert(0, image_header["NAXIS1"])
            else:
                # look for another good tile dimension
                naxis2 = image_header["NAXIS2"]
                for dim in [16, 24, 20, 30, 28, 26, 22, 18, 14]:
                    if naxis2 % dim == 0 or naxis2 % dim > 3:
                        tile_shape.insert(0, dim)
                        break
                else:
                    tile_shape.insert(0, 17)

            for i in range(2, naxis):
                # set all higher tile dimensions = 1
                tile_shape.insert(0, 1)

        # check if requested tile size causes the last tile to have
        # less than 4 pixels

        remain = image_header["NAXIS1"] % tile_shape[-1]  # 1st dimen

        original_tile_shape = tile_shape[:]

        if remain > 0 and remain < 4:
            tile_shape[-1] += 1  # try increasing tile size by 1

            remain = image_header["NAXIS1"] % tile_shape[-1]

            if remain > 0 and remain < 4:
                raise ValueError("Last tile along 1st dimension has less than 4 pixels")

        remain = image_header["NAXIS2"] % tile_shape[-2]  # 2nd dimen

        if remain > 0 and remain < 4:
            tile_shape[-2] += 1  # try increasing tile size by 1

            remain = image_header["NAXIS2"] % tile_shape[-2]

            if remain > 0 and remain < 4:
                raise ValueError("Last tile along 2nd dimension has less than 4 pixels")

        if tile_shape != original_tile_shape:
            warnings.warn(
                f"The tile shape should be such that no tiles have "
                f"fewer than 4 pixels. The tile shape has "
                f"automatically been changed from {original_tile_shape} "
                f"to {tile_shape}, but in future this will raise an "
                f"error and the correct tile shape should be specified "
                f"directly.",
                AstropyDeprecationWarning,
            )

    if len(tile_shape) == 0 and image_header["NAXIS"] > 0:
        tile_shape = [1] * (naxis - 1) + [image_header["NAXIS1"]]

    return tuple(tile_shape)
