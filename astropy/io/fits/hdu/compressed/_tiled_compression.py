"""
This module contains low level helper functions for compressing and
decompressing buffer for the Tiled Table Compression algorithms as specified in
the FITS 4 standard.
"""

import sys
from math import prod

import numpy as np

from astropy.io.fits.hdu.base import BITPIX2DTYPE

from ._codecs import PLIO1, Gzip1, Gzip2, HCompress1, NoCompress, Rice1
from ._quantization import DITHER_METHODS, QuantizationFailedException, Quantize
from .utils import _data_shape, _iter_array_tiles, _tile_shape

ALGORITHMS = {
    "GZIP_1": Gzip1,
    "GZIP_2": Gzip2,
    "RICE_1": Rice1,
    "RICE_ONE": Rice1,
    "PLIO_1": PLIO1,
    "HCOMPRESS_1": HCompress1,
    "NOCOMPRESS": NoCompress,
}

DEFAULT_ZBLANK = -2147483648


__all__ = [
    "compress_image_data",
    "decompress_image_data_section",
]


def _decompress_tile(buf, *, algorithm: str, **settings):
    """
    Decompress the buffer of a tile using the given compression algorithm.

    Parameters
    ----------
    buf
        The compressed buffer to be decompressed.
    algorithm
        A supported decompression algorithm.
    settings
        Any parameters for the given compression algorithm
    """
    return ALGORITHMS[algorithm](**settings).decode(buf)


def _compress_tile(buf, *, algorithm: str, **settings):
    """
    Compress the buffer of a tile using the given compression algorithm.

    Parameters
    ----------
    buf
        The decompressed buffer to be compressed.
    algorithm
        A supported compression algorithm.
    settings
        Any parameters for the given compression algorithm
    """
    return ALGORITHMS[algorithm](**settings).encode(buf)


def _header_to_settings(header):
    """
    Extract the settings which are constant given a header
    """
    settings = {}
    compression_type = header["ZCMPTYPE"]
    if compression_type == "GZIP_2":
        settings["itemsize"] = abs(header["ZBITPIX"]) // 8
    elif compression_type in ("RICE_1", "RICE_ONE"):
        settings["blocksize"] = _get_compression_setting(header, "BLOCKSIZE", 32)
        settings["bytepix"] = _get_compression_setting(header, "BYTEPIX", 4)
    elif compression_type == "HCOMPRESS_1":
        settings["bytepix"] = 8
        settings["scale"] = int(_get_compression_setting(header, "SCALE", 0))
        settings["smooth"] = _get_compression_setting(header, "SMOOTH", 0)

    return settings


def _update_tile_settings(settings, compression_type, actual_tile_shape):
    """
    Update the settings with tile-specific settings
    """
    if compression_type in ("PLIO_1", "RICE_1", "RICE_ONE"):
        # We have to calculate the tilesize from the shape of the tile not the
        # header, so that it's correct for edge tiles etc.
        settings["tilesize"] = prod(actual_tile_shape)
    elif compression_type == "HCOMPRESS_1":
        # HCOMPRESS requires 2D tiles, so to find the shape of the 2D tile we
        # need to ignore all length 1 tile dimensions
        # Also cfitsio expects the tile shape in C order
        shape_2d = tuple(nd for nd in actual_tile_shape if nd != 1)
        if len(shape_2d) != 2:
            raise ValueError(f"HCOMPRESS expects two dimensional tiles, got {shape_2d}")
        settings["nx"] = shape_2d[0]
        settings["ny"] = shape_2d[1]

    return settings


def _finalize_array(tile_buffer, *, bitpix, tile_shape, algorithm, lossless):
    """
    Convert a buffer to an array.

    This is a helper function which takes a raw buffer (as output by .decode)
    and translates it into a numpy array with the correct dtype, endianness and
    shape.
    """
    tile_size = prod(tile_shape)

    if algorithm.startswith("GZIP") or algorithm == "NOCOMPRESS":
        # This algorithm is taken from fitsio
        # https://github.com/astropy/astropy/blob/a8cb1668d4835562b89c0d0b3448ac72ca44db63/cextern/cfitsio/lib/imcompress.c#L6345-L6388
        tile_bytesize = len(tile_buffer)
        if tile_bytesize == tile_size * 2:
            dtype = ">i2"
        elif tile_bytesize == tile_size * 4:
            if bitpix < 0 and lossless:
                dtype = ">f4"
            else:
                dtype = ">i4"
        elif tile_bytesize == tile_size * 8:
            if bitpix < 0 and lossless:
                dtype = ">f8"
            else:
                dtype = ">i8"
        else:
            # Just return the raw bytes
            dtype = ">u1"
        tile_data = np.asarray(tile_buffer).view(dtype).reshape(tile_shape)
    else:
        # For RICE_1 compression the tiles that are on the edge can end up
        # being padded, so we truncate excess values
        if algorithm in ("RICE_1", "RICE_ONE", "PLIO_1") and tile_size < len(
            tile_buffer
        ):
            tile_buffer = tile_buffer[:tile_size]

        if tile_buffer.data.format == "b":
            # NOTE: this feels like a Numpy bug - need to investigate
            tile_data = np.asarray(tile_buffer, dtype=np.uint8).reshape(tile_shape)
        else:
            tile_data = np.asarray(tile_buffer).reshape(tile_shape)

    return tile_data


def _check_compressed_header(header):
    # NOTE: this could potentially be moved up into CompImageHDU, e.g. in a
    # _verify method.

    # Check for overflows which might cause issues when calling C code

    for kw in ["ZNAXIS", "ZVAL1", "ZVAL2", "ZBLANK", "BLANK"]:
        if kw in header:
            if header[kw] > 0 and header[kw] > np.iinfo(np.intc).max:
                raise OverflowError(f"{kw} value {header[kw]} is too large")

    for i in range(1, header["ZNAXIS"] + 1):
        for kw_name in ["ZNAXIS", "ZTILE"]:
            kw = f"{kw_name}{i}"
            if kw in header:
                if header[kw] > 0 and header[kw] > np.iinfo(np.int32).max:
                    raise OverflowError(f"{kw} value {header[kw]} is too large")

    for i in range(1, header["NAXIS"] + 1):
        kw = f"NAXIS{i}"
        if kw in header:
            if header[kw] > 0 and header[kw] > np.iinfo(np.int64).max:
                raise OverflowError(f"{kw} value {header[kw]} is too large")

    for kw in ["TNULL1", "PCOUNT", "THEAP"]:
        if kw in header:
            if header[kw] > 0 and header[kw] > np.iinfo(np.int64).max:
                raise OverflowError(f"{kw} value {header[kw]} is too large")

    for kw in ["ZVAL3"]:
        if kw in header:
            # float() below to avoid a NEP 50 warning about a cast to float32 inf.
            if header[kw] > float(np.finfo(np.float32).max):
                raise OverflowError(f"{kw} value {header[kw]} is too large")

    # Validate data types

    for kw in ["ZSCALE", "ZZERO", "TZERO1", "TSCAL1"]:
        if kw in header:
            if not np.isreal(header[kw]):
                raise TypeError(f"{kw} should be floating-point")

    for kw in ["TTYPE1", "TFORM1", "ZCMPTYPE", "ZNAME1", "ZQUANTIZ"]:
        if kw in header:
            if not isinstance(header[kw], str):
                raise TypeError(f"{kw} should be a string")

    for kw in ["ZDITHER0"]:
        if kw in header:
            if not np.isreal(header[kw]) or not float(header[kw]).is_integer():
                raise TypeError(f"{kw} should be an integer")

    for suffix in range(1, header["TFIELDS"] + 1):
        if header.get(f"TTYPE{suffix}", "").endswith("COMPRESSED_DATA"):
            for valid in ["PB", "PI", "PJ", "QB", "QI", "QJ"]:
                if header[f"TFORM{suffix}"].startswith((valid, f"1{valid}")):
                    break
            else:
                raise RuntimeError(f"Invalid TFORM{suffix}: {header[f'TFORM{suffix}']}")

    # Check values

    for kw in ["TFIELDS", "PCOUNT"] + [
        f"NAXIS{idx + 1}" for idx in range(header["NAXIS"])
    ]:
        if kw in header:
            if header[kw] < 0:
                raise ValueError(f"{kw} should not be negative.")

    for kw in ["ZNAXIS", "TFIELDS"]:
        if kw in header:
            if header[kw] < 0 or header[kw] > 999:
                raise ValueError(f"{kw} should be in the range 0 to 999")

    if header["ZBITPIX"] not in [8, 16, 32, 64, -32, -64]:
        raise ValueError(f"Invalid value for BITPIX: {header['ZBITPIX']}")

    if header["ZCMPTYPE"] not in ALGORITHMS:
        raise ValueError(f"Unrecognized compression type: {header['ZCMPTYPE']}")

    # Check that certain keys are present

    header["ZNAXIS"]
    header["ZBITPIX"]


def _get_compression_setting(header, name, default):
    # Settings for the various compression algorithms are stored in pairs of
    # keywords called ZNAME? and ZVAL? - a given compression setting could be
    # in any ZNAME? so we need to check through all the possible ZNAMEs which
    # one matches the required setting.

    for i in range(1, 1000):
        if f"ZNAME{i}" not in header:
            break
        if header[f"ZNAME{i}"].lower() == name.lower():
            return header[f"ZVAL{i}"]

    return default


def _column_dtype(compressed_coldefs, column_name):
    tform = compressed_coldefs[column_name].format.removeprefix("1")
    if tform[1] == "B":
        dtype = np.uint8
    elif tform[1] == "I":
        dtype = ">i2"
    elif tform[1] == "J":
        dtype = ">i4"
    return np.dtype(dtype)


def _get_data_from_heap(hdu, size, offset, dtype, heap_cache=None):
    if heap_cache is None:
        return hdu._get_raw_data(size, dtype, hdu._data_offset + hdu._theap + offset)
    else:
        itemsize = dtype.itemsize
        data = heap_cache[offset : offset + size * itemsize]
        if itemsize > 1:
            return data.view(dtype)
        else:
            return data


def decompress_image_data_section(
    compressed_data,
    compression_type,
    compressed_header,
    bintable,
    image_header,
    first_tile_index,
    last_tile_index,
):
    """
    Decompress the data in a `~astropy.io.fits.CompImageHDU`.

    Parameters
    ----------
    compressed_data : `~astropy.io.fits.FITS_rec`
        The compressed data
    compression_type : str
        The compression algorithm
    compressed_header : `~astropy.io.fits.Header`
        The header of the compressed binary table
    bintable : `~astropy.io.fits.BinTableHDU`
        The binary table HDU, used to access the raw heap data
    first_tile_index : iterable
        The indices of the first tile to decompress along each dimension
    last_tile_index : iterable
        The indices of the last tile to decompress along each dimension

    Returns
    -------
    data : `numpy.ndarray`
        The decompressed data array.
    """
    compressed_coldefs = compressed_data._coldefs

    _check_compressed_header(compressed_header)

    tile_shape = _tile_shape(compressed_header)
    data_shape = _data_shape(compressed_header)

    first_array_index = first_tile_index * tile_shape
    last_array_index = (last_tile_index + 1) * tile_shape

    last_array_index = np.minimum(data_shape, last_array_index)

    buffer_shape = tuple((last_array_index - first_array_index).astype(int))

    image_data = np.empty(
        buffer_shape, dtype=BITPIX2DTYPE[compressed_header["ZBITPIX"]]
    )

    quantized = "ZSCALE" in compressed_data.dtype.names

    if image_data.size == 0:
        return image_data

    settings = _header_to_settings(compressed_header)
    zbitpix = compressed_header["ZBITPIX"]
    dither_method = DITHER_METHODS[compressed_header.get("ZQUANTIZ", "NO_DITHER")]
    dither_seed = compressed_header.get("ZDITHER0", 0)

    # NOTE: in the following and below we convert the column to a Numpy array
    # for performance reasons, as accessing rows from a FITS_rec column is
    # otherwise slow.
    compressed_data_column = np.array(compressed_data["COMPRESSED_DATA"])
    compressed_data_dtype = _column_dtype(compressed_coldefs, "COMPRESSED_DATA")

    if "ZBLANK" in compressed_coldefs.dtype.names:
        zblank_column = np.array(compressed_data["ZBLANK"])
    else:
        zblank_column = None

    if "ZSCALE" in compressed_coldefs.dtype.names:
        zscale_column = np.array(compressed_data["ZSCALE"])
    else:
        zscale_column = None

    if "ZZERO" in compressed_coldefs.dtype.names:
        zzero_column = np.array(compressed_data["ZZERO"])
    else:
        zzero_column = None

    zblank_header = compressed_header.get("ZBLANK", image_header.get("BLANK", None))

    gzip_compressed_data_column = None
    gzip_compressed_data_dtype = None

    # If all the data is requested, read in all the heap.
    if tuple(buffer_shape) == tuple(data_shape):
        heap_cache = bintable._get_raw_data(
            compressed_header["PCOUNT"],
            np.uint8,
            bintable._data_offset + bintable._theap,
        )
    else:
        heap_cache = None

    for row_index, tile_slices in _iter_array_tiles(
        data_shape, tile_shape, first_tile_index, last_tile_index
    ):
        # For tiles near the edge, the tile shape from the header might not be
        # correct so we have to pass the shape manually.
        actual_tile_shape = image_data[tile_slices].shape

        settings = _update_tile_settings(settings, compression_type, actual_tile_shape)

        if compressed_data_column[row_index][0] == 0:
            if gzip_compressed_data_column is None:
                gzip_compressed_data_column = np.array(
                    compressed_data["GZIP_COMPRESSED_DATA"]
                )
                gzip_compressed_data_dtype = _column_dtype(
                    compressed_coldefs, "GZIP_COMPRESSED_DATA"
                )

            # When quantizing floating point data, sometimes the data will not
            # quantize efficiently. In these cases the raw floating point data can
            # be losslessly GZIP compressed and stored in the `GZIP_COMPRESSED_DATA`
            # column.

            cdata = _get_data_from_heap(
                bintable,
                *gzip_compressed_data_column[row_index],
                gzip_compressed_data_dtype,
                heap_cache=heap_cache,
            )

            tile_buffer = _decompress_tile(cdata, algorithm="GZIP_1")

            tile_data = _finalize_array(
                tile_buffer,
                bitpix=zbitpix,
                tile_shape=actual_tile_shape,
                algorithm="GZIP_1",
                lossless=True,
            )

        else:
            cdata = _get_data_from_heap(
                bintable,
                *compressed_data_column[row_index],
                compressed_data_dtype,
                heap_cache=heap_cache,
            )

            if compression_type == "GZIP_2":
                # Decompress with GZIP_1 just to find the total number of
                # elements in the uncompressed data.
                # TODO: find a way to avoid doing this for all tiles
                tile_data = np.asarray(_decompress_tile(cdata, algorithm="GZIP_1"))
                settings["itemsize"] = tile_data.size // int(prod(actual_tile_shape))

            tile_buffer = _decompress_tile(
                cdata, algorithm=compression_type, **settings
            )

            tile_data = _finalize_array(
                tile_buffer,
                bitpix=zbitpix,
                tile_shape=actual_tile_shape,
                algorithm=compression_type,
                lossless=not quantized,
            )

            if zblank_column is None:
                zblank = zblank_header
            else:
                zblank = zblank_column[row_index]

            if zblank is not None:
                blank_mask = tile_data == zblank

            if quantized:
                q = Quantize(
                    row=(row_index + dither_seed) if dither_method != -1 else 0,
                    dither_method=dither_method,
                    quantize_level=None,
                    bitpix=zbitpix,
                )
                tile_data = np.asarray(
                    q.decode_quantized(
                        tile_data, zscale_column[row_index], zzero_column[row_index]
                    )
                ).reshape(actual_tile_shape)

            if zblank is not None:
                if not tile_data.flags.writeable:
                    tile_data = tile_data.copy()
                tile_data[blank_mask] = zblank_header if zbitpix > 0 else np.nan

        image_data[tile_slices] = tile_data

    return image_data


def compress_image_data(
    image_data,
    compression_type,
    compressed_header,
    compressed_coldefs,
):
    """
    Compress the data in a `~astropy.io.fits.CompImageHDU`.

    The input HDU is expected to have a uncompressed numpy array as it's
    ``.data`` attribute.

    Parameters
    ----------
    image_data : `~numpy.ndarray`
        The image data to compress
    compression_type : str
        The compression algorithm
    compressed_header : `~astropy.io.fits.Header`
        The header of the compressed binary table
    compressed_coldefs : `~astropy.io.fits.ColDefs`
        The ColDefs object for the compressed binary table

    Returns
    -------
    nbytes : `int`
        The number of bytes of the heap.
    heap : `bytes`
        The bytes of the FITS table heap.
    """
    if not isinstance(image_data, np.ndarray):
        raise TypeError("Image data must be a numpy.ndarray")

    _check_compressed_header(compressed_header)

    # TODO: This implementation is memory inefficient as it generates all the
    # compressed bytes before forming them into the heap, leading to 2x the
    # potential memory usage. Directly storing the compressed bytes into an
    # expanding heap would fix this.

    tile_shape = _tile_shape(compressed_header)
    data_shape = _data_shape(compressed_header)

    compressed_bytes = []
    gzip_fallback = []
    scales = []
    zeros = []
    zblank = None

    noisebit = _get_compression_setting(compressed_header, "noisebit", 0)

    settings = _header_to_settings(compressed_header)

    for irow, tile_slices in _iter_array_tiles(data_shape, tile_shape):
        tile_data = image_data[tile_slices]

        if tile_data.dtype.kind == "u":
            if tile_data.dtype.itemsize == 4:
                tile_data = (tile_data.astype(np.int64) - 2**31).astype(np.int32)
            elif tile_data.dtype.itemsize == 2:
                tile_data = (tile_data.astype(np.int32) - 2**15).astype(np.int16)

        settings = _update_tile_settings(settings, compression_type, tile_data.shape)

        quantize = "ZSCALE" in compressed_coldefs.dtype.names

        if tile_data.dtype.kind == "f" and quantize:
            dither_method = DITHER_METHODS[
                compressed_header.get("ZQUANTIZ", "NO_DITHER")
            ]
            dither_seed = compressed_header.get("ZDITHER0", 0)
            q = Quantize(
                row=(irow + dither_seed) if dither_method != -1 else 0,
                dither_method=dither_method,
                quantize_level=noisebit,
                bitpix=compressed_header["ZBITPIX"],
            )
            original_shape = tile_data.shape

            # If there are any NaN values in the data, we should reset them to
            # a value that will not affect the quantization (an already existing
            # data value in the array) and we can then reset this after quantization
            # to ZBLANK and set the appropriate header keyword
            nan_mask = np.isnan(tile_data)
            any_nan = np.any(nan_mask)
            if any_nan:
                # Note that we need to copy here to avoid modifying the input array.
                tile_data = tile_data.copy()
                if np.all(nan_mask):
                    tile_data[nan_mask] = 0
                else:
                    tile_data[nan_mask] = np.nanmin(tile_data)

            try:
                tile_data, scale, zero = q.encode_quantized(tile_data)
            except QuantizationFailedException:
                if any_nan:
                    # reset NaN values since we will losslessly compress.
                    tile_data[nan_mask] = np.nan

                scales.append(0)
                zeros.append(0)
                gzip_fallback.append(True)

            else:
                tile_data = np.asarray(tile_data).reshape(original_shape)

                if any_nan:
                    if not tile_data.flags.writeable:
                        tile_data = tile_data.copy()
                    # For now, we just use the default ZBLANK value and assume
                    # this is the same for all tiles. We could generalize this
                    # to allow different ZBLANK values (for example if the data
                    # includes this value by chance) and to allow different values
                    # per tile, which is allowed by the FITS standard.
                    tile_data[nan_mask] = DEFAULT_ZBLANK
                    zblank = DEFAULT_ZBLANK

                scales.append(scale)
                zeros.append(zero)
                gzip_fallback.append(False)

        else:
            scales.append(0)
            zeros.append(0)
            gzip_fallback.append(False)

        if gzip_fallback[-1]:
            cbytes = _compress_tile(tile_data, algorithm="GZIP_1")
        else:
            cbytes = _compress_tile(tile_data, algorithm=compression_type, **settings)
        compressed_bytes.append(cbytes)

    if zblank is not None:
        compressed_header["ZBLANK"] = zblank

    table = np.zeros(
        len(compressed_bytes), dtype=compressed_coldefs.dtype.newbyteorder(">")
    )

    if "ZSCALE" in table.dtype.names:
        table["ZSCALE"] = np.array(scales)
        table["ZZERO"] = np.array(zeros)

    for irow, cbytes in enumerate(compressed_bytes):
        table["COMPRESSED_DATA"][irow, 0] = len(cbytes)

    table["COMPRESSED_DATA"][:1, 1] = 0
    table["COMPRESSED_DATA"][1:, 1] = np.cumsum(table["COMPRESSED_DATA"][:-1, 0])

    for irow in range(len(compressed_bytes)):
        if gzip_fallback[irow]:
            table["GZIP_COMPRESSED_DATA"][irow] = table["COMPRESSED_DATA"][irow]
            table["COMPRESSED_DATA"][irow] = 0

    # For PLIO_1, the size of each heap element is a factor of two lower than
    # the real size - not clear if this is deliberate or bug somewhere.
    if compression_type == "PLIO_1":
        table["COMPRESSED_DATA"][:, 0] //= 2

    # For PLIO_1, it looks like the compressed data is always stored big endian
    if compression_type == "PLIO_1":
        for irow in range(len(compressed_bytes)):
            if not gzip_fallback[irow]:
                array = np.frombuffer(compressed_bytes[irow], dtype="i2")
                if array.dtype.byteorder == "<" or (
                    array.dtype.byteorder == "=" and sys.byteorder == "little"
                ):
                    compressed_bytes[irow] = array.astype(">i2", copy=False).tobytes()

    compressed_bytes = b"".join(compressed_bytes)

    table_bytes = table.tobytes()

    heap = table.tobytes() + compressed_bytes

    return heap
