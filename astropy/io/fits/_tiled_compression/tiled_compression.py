"""
This module contains low level helper functions for compressing and
decompressing buffer for the Tiled Table Compression algorithms as specified in
the FITS 4 standard.
"""
import numpy as np

from astropy.io.fits.hdu.base import BITPIX2DTYPE

from .codecs import PLIO1, Gzip1, Gzip2, HCompress1, Rice1
from .quantization import DITHER_METHODS, QuantizationFailedException, Quantize

ALGORITHMS = {
    "GZIP_1": Gzip1,
    "GZIP_2": Gzip2,
    "RICE_1": Rice1,
    "RICE_ONE": Rice1,
    "PLIO_1": PLIO1,
    "HCOMPRESS_1": HCompress1,
}

DEFAULT_ZBLANK = -2147483648


__all__ = ["compress_hdu", "decompress_hdu"]


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


def _tile_shape(header):
    return tuple(header[f"ZTILE{idx}"] for idx in range(header["ZNAXIS"], 0, -1))


def _data_shape(header):
    return tuple(header[f"ZNAXIS{idx}"] for idx in range(header["ZNAXIS"], 0, -1))


def _header_to_settings(header, actual_tile_shape):

    settings = {}

    if header["ZCMPTYPE"] == "GZIP_2":
        settings["itemsize"] = abs(header["ZBITPIX"]) // 8
    elif header["ZCMPTYPE"] == "PLIO_1":
        # We have to calculate the tilesize from the shape of the tile not the
        # header, so that it's correct for edge tiles etc.
        settings["tilesize"] = np.product(actual_tile_shape)
    elif header["ZCMPTYPE"] in ("RICE_1", "RICE_ONE"):
        settings["blocksize"] = _get_compression_setting(header, "BLOCKSIZE", 32)
        settings["bytepix"] = _get_compression_setting(header, "BYTEPIX", 4)
        settings["tilesize"] = np.product(actual_tile_shape)
    elif header["ZCMPTYPE"] == "HCOMPRESS_1":
        settings["bytepix"] = 8
        settings["scale"] = int(_get_compression_setting(header, "SCALE", 0))
        settings["smooth"] = _get_compression_setting(header, "SMOOTH", 0)
        # HCOMPRESS requires 2D tiles, so to find the shape of the 2D tile we
        # need to ignore all length 1 tile dimensions
        # Also cfitsio expects the tile shape in C order
        shape_2d = tuple(nd for nd in actual_tile_shape if nd != 1)
        if len(shape_2d) != 2:
            raise ValueError(f"HCOMPRESS expects two dimensional tiles, got {shape_2d}")
        settings["nx"] = shape_2d[0]
        settings["ny"] = shape_2d[1]

    return settings


def _finalize_array(
    tile_buffer, *, bitpix, tile_shape, algorithm, lossless
):
    """
    Convert a buffer to an array.

    This is a helper function which takes a raw buffer (as output by .decode)
    and translates it into a numpy array with the correct dtype, endianness and
    shape.
    """

    if algorithm.startswith("GZIP"):
        # This algorithm is taken from fitsio
        # https://github.com/astropy/astropy/blob/a8cb1668d4835562b89c0d0b3448ac72ca44db63/cextern/cfitsio/lib/imcompress.c#L6345-L6388
        tilelen = np.product(tile_shape)
        tilebytesize = len(tile_buffer)
        if tilebytesize == tilelen * 2:
            dtype = ">i2"
        elif tilebytesize == tilelen * 4:
            if bitpix < 0 and lossless:
                dtype = ">f4"
            else:
                dtype = ">i4"
        elif tilebytesize == tilelen * 8:
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
        if algorithm in ("RICE_1", "RICE_ONE", "PLIO_1"):
            tile_buffer = tile_buffer[: np.product(tile_shape)]

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
            if header[kw] > np.finfo(np.float32).max:
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

    if "TFORM1" in header:
        for valid in ["1PB", "1PI", "1PJ", "1QB", "1QI", "1QJ"]:
            if header["TFORM1"].startswith(valid):
                break
        else:
            raise RuntimeError(f"Invalid TFORM1: {header['TFORM1']}")

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


def decompress_hdu(hdu):
    """
    Decompress the data in a `~astropy.io.fits.CompImageHDU`.

    Parameters
    ----------
    hdu : `astropy.io.fits.CompImageHDU`
        Input HDU to decompress the data for.

    Returns
    -------

    data : `numpy.ndarray`
        The decompressed data array.
    """

    _check_compressed_header(hdu._header)

    tile_shape = _tile_shape(hdu._header)
    data_shape = _data_shape(hdu._header)

    data = np.zeros(data_shape, dtype=BITPIX2DTYPE[hdu._header["ZBITPIX"]])

    quantized = "ZSCALE" in hdu.compressed_data.dtype.names

    if len(hdu.compressed_data) == 0:
        return None

    override_itemsize = None

    istart = np.zeros(data.ndim, dtype=int)
    for irow, row in enumerate(hdu.compressed_data):

        # In the following, we don't need to special case tiles near the edge
        # as Numpy will automatically ignore parts of the slices that are out
        # of bounds.
        tile_slices = tuple(
            [
                slice(istart[idx], istart[idx] + tile_shape[idx])
                for idx in range(len(istart))
            ]
        )

        # For tiles near the edge, the tile shape from the header might not be
        # correct so we have to pass the shape manually.
        actual_tile_shape = data[tile_slices].shape
        settings = _header_to_settings(hdu._header, actual_tile_shape)

        cdata = row["COMPRESSED_DATA"]

        # When quantizing floating point data, sometimes the data will not
        # quantize efficiently. In these cases the raw floating point data can
        # be losslessly GZIP compressed and stored in the `GZIP_COMPRESSED_DATA`
        # column.
        gzip_fallback = len(cdata) == 0

        if gzip_fallback:

            tile_buffer = _decompress_tile(
                row["GZIP_COMPRESSED_DATA"], algorithm="GZIP_1"
            )

            tile_data = _finalize_array(
                tile_buffer,
                bitpix=hdu._header["ZBITPIX"],
                tile_shape=actual_tile_shape,
                algorithm="GZIP_1",
                lossless=True,
            )

        else:

            if hdu._header["ZCMPTYPE"] == "GZIP_2":
                # Decompress with GZIP_1 just to find the total number of
                # elements in the uncompressed data. We just need to do this once
                # as this will be the same for all tiles.
                if override_itemsize is None:
                    tile_data = np.asarray(
                        _decompress_tile(cdata, algorithm="GZIP_1")
                    )
                    override_itemsize = tile_data.size // int(np.product(actual_tile_shape))
                settings["itemsize"] = override_itemsize

            tile_buffer = _decompress_tile(
                cdata, algorithm=hdu._header["ZCMPTYPE"], **settings
            )

            tile_data = _finalize_array(
                tile_buffer,
                bitpix=hdu._header["ZBITPIX"],
                tile_shape=actual_tile_shape,
                algorithm=hdu._header["ZCMPTYPE"],
                lossless=not quantized,
            )

            if "ZBLANK" in row.array.names:
                zblank = row["ZBLANK"]
            elif "ZBLANK" in hdu._header:
                zblank = hdu._header["ZBLANK"]
            else:
                zblank = None

            if zblank is not None:
                blank_mask = tile_data == zblank

            if quantized:
                dither_method = DITHER_METHODS[hdu._header.get("ZQUANTIZ", "NO_DITHER")]
                dither_seed = hdu._header.get("ZDITHER0", 0)
                q = Quantize(
                    row=(irow + dither_seed) if dither_method != -1 else 0,
                    dither_method=dither_method,
                    quantize_level=None,
                    bitpix=hdu._header["ZBITPIX"],
                )
                tile_data = np.asarray(
                    q.decode_quantized(tile_data, row["ZSCALE"], row["ZZERO"])
                ).reshape(actual_tile_shape)

            if zblank is not None:
                if not tile_data.flags.writeable:
                    tile_data = tile_data.copy()
                tile_data[blank_mask] = np.nan

        data[tile_slices] = tile_data
        istart[-1] += tile_shape[-1]
        for idx in range(data.ndim - 1, 0, -1):
            if istart[idx] >= data_shape[idx]:
                istart[idx] = 0
                istart[idx - 1] += tile_shape[idx - 1]

    return data


def compress_hdu(hdu):
    """
    Compress the data in a `~astropy.io.fits.CompImageHDU`.

    The input HDU is expected to have a uncompressed numpy array as it's
    ``.data`` attribute.

    Parameters
    ----------
    hdu : `astropy.io.fits.CompImageHDU`
        Input HDU to compress the data for.

    Returns
    -------
    nbytes : `int`
        The number of bytes of the heap.
    heap : `bytes`
        The bytes of the FITS table heap.
    """

    if not isinstance(hdu.data, np.ndarray):
        raise TypeError("CompImageHDU.data must be a numpy.ndarray")

    _check_compressed_header(hdu._header)

    # TODO: This implementation is memory inefficient as it generates all the
    # compressed bytes before forming them into the heap, leading to 2x the
    # potential memory usage. Directly storing the compressed bytes into an
    # expanding heap would fix this.

    tile_shape = _tile_shape(hdu._header)
    data_shape = _data_shape(hdu._header)

    compressed_bytes = []
    gzip_fallback = []
    scales = []
    zeros = []
    zblank = None

    irow = 0
    istart = np.zeros(len(data_shape), dtype=int)

    noisebit = _get_compression_setting(hdu._header, "noisebit", 0)

    while True:

        if hdu.data is None:
            break

        # In the following, we don't need to special case tiles near the edge
        # as Numpy will automatically ignore parts of the slices that are out
        # of bounds.
        slices = tuple(
            [
                slice(istart[idx], istart[idx] + tile_shape[idx])
                for idx in range(len(istart))
            ]
        )

        data = hdu.data[slices]

        settings = _header_to_settings(hdu._header, data.shape)

        quantize = "ZSCALE" in hdu.columns.dtype.names

        if data.dtype.kind == "f" and quantize:
            noisebit = _get_compression_setting(hdu._header, "noisebit", 0)
            dither_method = DITHER_METHODS[hdu._header.get("ZQUANTIZ", "NO_DITHER")]
            dither_seed = hdu._header.get("ZDITHER0", 0)
            q = Quantize(
                row=(irow + dither_seed) if dither_method != -1 else 0,
                dither_method=dither_method,
                quantize_level=noisebit,
                bitpix=hdu._header["ZBITPIX"],
            )
            original_shape = data.shape

            # If there are any NaN values in the data, we should reset them to
            # a value that will not affect the quantization (an already existing
            # data value in the array) and we can then reset this after quantization
            # to ZBLANK and set the appropriate header keyword
            nan_mask = np.isnan(data)
            any_nan = np.any(nan_mask)
            if any_nan:
                # Note that we need to copy here to avoid modifying the input array.
                data = data.copy()
                if np.all(nan_mask):
                    data[nan_mask] = 0
                else:
                    data[nan_mask] = np.nanmin(data)

            try:
                data, scale, zero = q.encode_quantized(data)
            except QuantizationFailedException:

                if any_nan:
                    # reset NaN values since we will losslessly compress.
                    data[nan_mask] = np.nan

                scales.append(0)
                zeros.append(0)
                gzip_fallback.append(True)

            else:
                data = np.asarray(data).reshape(original_shape)

                if any_nan:
                    if not data.flags.writeable:
                        data = data.copy()
                    # For now, we just use the default ZBLANK value and assume
                    # this is the same for all tiles. We could generalize this
                    # to allow different ZBLANK values (for example if the data
                    # includes this value by chance) and to allow different values
                    # per tile, which is allowed by the FITS standard.
                    data[nan_mask] = DEFAULT_ZBLANK
                    zblank = DEFAULT_ZBLANK

                scales.append(scale)
                zeros.append(zero)
                gzip_fallback.append(False)

        else:
            scales.append(0)
            zeros.append(0)
            gzip_fallback.append(False)

        # The original compress_hdu assumed the data was in native endian, so we
        # change this here:
        if hdu._header["ZCMPTYPE"].startswith("GZIP") or gzip_fallback[-1]:
            # This is apparently needed so that our heap data agrees with
            # the C implementation!?
            data = data.astype(data.dtype.newbyteorder(">"))
        else:
            if not data.dtype.isnative:
                data = data.astype(data.dtype.newbyteorder("="))

        if gzip_fallback[-1]:
            cbytes = _compress_tile(data, algorithm="GZIP_1")
        else:
            cbytes = _compress_tile(data, algorithm=hdu._header["ZCMPTYPE"], **settings)
        compressed_bytes.append(cbytes)

        istart[-1] += tile_shape[-1]

        for idx in range(data.ndim - 1, 0, -1):
            if istart[idx] >= data_shape[idx]:
                istart[idx] = 0
                istart[idx - 1] += tile_shape[idx - 1]

        if istart[0] >= data_shape[0]:
            break

        irow += 1

    if zblank is not None:
        hdu._header["ZBLANK"] = zblank

    table = np.zeros(len(compressed_bytes), dtype=hdu.columns.dtype.newbyteorder(">"))

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
    if hdu._header["ZCMPTYPE"] == "PLIO_1":
        table["COMPRESSED_DATA"][:, 0] //= 2

    # For PLIO_1, it looks like the compressed data is byteswapped
    if hdu._header["ZCMPTYPE"] == "PLIO_1":
        for irow in range(len(compressed_bytes)):
            if not gzip_fallback[irow]:
                compressed_bytes[irow] = (
                    np.frombuffer(compressed_bytes[irow], dtype="<i2")
                    .astype(">i2")
                    .tobytes()
                )

    compressed_bytes = b"".join(compressed_bytes)

    table_bytes = table.tobytes()

    if len(table_bytes) != hdu._theap:
        raise Exception(
            f"Unexpected compressed table size (expected {hdu._theap}, got {len(table_bytes)})"
        )
    heap = table.tobytes() + compressed_bytes

    return len(compressed_bytes), np.frombuffer(heap, dtype=np.uint8)
