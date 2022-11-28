"""
This module contains low level helper functions for compressing and decompressing buffer for the Tiled Table Compression algorithms as specified in the FITS 4 standard.
"""
from abc import abstractmethod
from gzip import compress as gzip_compress
from gzip import decompress as gzip_decompress

import numpy as np

from astropy.io.fits.hdu.base import BITPIX2DTYPE
from astropy.io.fits.tiled_compression._compression import (
    compress_hcompress_1_c,
    compress_plio_1_c,
    compress_rice_1_c,
    decompress_hcompress_1_c,
    decompress_plio_1_c,
    decompress_rice_1_c,
)

__all__ = [
    "Gzip1",
    "Gzip2",
    "Rice1",
    "PLIO1",
    "HCompress1",
    "compress_tile",
    "decompress_tile",
    "compress_hdu",
    "decompress_hdu",
]


# We define our compression classes in the form of a numcodecs class. We make
# the dependency on numcodecs optional as we can use them internally without it
# (for now).
try:
    from numcodecs.abc import Codec
except ImportError:

    class Codec:
        codec_id = None
        """Codec identifier."""

        @abstractmethod
        def encode(self, buf):  # pragma: no cover
            """Encode data in `buf`.
            Parameters
            ----------
            buf : buffer-like
                Data to be encoded. May be any object supporting the new-style
                buffer protocol.
            Returns
            -------
            enc : buffer-like
                Encoded data. May be any object supporting the new-style buffer
                protocol.
            """

        @abstractmethod
        def decode(self, buf, out=None):  # pragma: no cover
            """Decode data in `buf`.
            Parameters
            ----------
            buf : buffer-like
                Encoded data. May be any object supporting the new-style buffer
                protocol.
            out : buffer-like, optional
                Writeable buffer to store decoded data. N.B. if provided, this buffer must
                be exactly the right size to store the decoded data.
            Returns
            -------
            dec : buffer-like
                Decoded data. May be any object supporting the new-style
                buffer protocol.
            """


class Gzip1(Codec):
    """
    The FITS GZIP 1 compression and decompression algorithm.

    The Gzip algorithm is used in the free GNU software compression utility of
    the same name. It was created by J. L. Gailly and M. Adler, based on the
    DEFLATE algorithm (Deutsch 1996), which is a combination of LZ77 (Ziv &
    Lempel 1977) and Huffman coding.
    """

    codec_id = "FITS_GZIP1"

    def decode(self, buf):
        """
        Decompress buffer using the GZIP_1 algorithm.

        {_GZIP_1_DESCRIPTION}

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            The decompressed buffer.
        """
        cbytes = np.frombuffer(buf, dtype=np.uint8).tobytes()
        dbytes = gzip_decompress(cbytes)
        return np.frombuffer(dbytes, dtype=np.uint8).data

    def encode(self, buf):
        """
        Compress the data in the buffer using the GZIP_1 algorithm.

        Parameters
        ----------
        buf
            The buffer to compress.

        Returns
        -------
        buf
            A buffer with compressed data.
        """
        dbytes = np.asarray(buf).tobytes()
        return gzip_compress(dbytes)


class Gzip2(Codec):
    """
    The FTIS GZIP2 compression and decompression algorithm.

    The gzip2 algorithm is a variation on 'GZIP 1'. In this case the buffer in
    the array of data values are shuffled so that they are arranged in order of
    decreasing significance before being compressed.

    For example, a five-element contiguous array of two-byte (16-bit) integer
    values, with an original big-endian byte order of:

    .. math::
        A1 A2 B1 B2 C1 C2 D1 D2 E1 E2

    will have the following byte order after shuffling:

    .. math::
        A1 B1 C1 D1 E1 A2 B2 C2 D2 E2,

    where A1, B1, C1, D1, and E1 are the most-significant buffer from
    each of the integer values.

    Byte shuffling shall only be performed for integer or floating-point
    numeric data types; logical, bit, and character types must not be shuffled.

    Parameters
    ----------
    itemsize
        The number of buffer per value (e.g. 2 for a 16-bit integer)

    """

    codec_id = "FITS_GZIP2"

    def __init__(self, itemsize: int):
        super().__init__()
        self.itemsize = itemsize

    def decode(self, buf):
        """
        Decompress buffer using the GZIP_2 algorithm.

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            The decompressed buffer.
        """
        cbytes = np.frombuffer(buf, dtype=np.uint8).tobytes()
        # Start off by unshuffling buffer
        unshuffled_buffer = gzip_decompress(cbytes)
        array = np.frombuffer(unshuffled_buffer, dtype=np.uint8)
        return array.reshape((self.itemsize, -1)).T.ravel().data

    def encode(self, buf):
        """
        Compress the data in the buffer using the GZIP_2 algorithm.

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            The decompressed buffer.
        """
        # Start off by shuffling buffer
        array = np.asarray(buf).ravel().view(np.uint8)
        shuffled_buffer = array.reshape((-1, self.itemsize)).T.ravel().tobytes()
        return gzip_compress(shuffled_buffer)


class Rice1(Codec):
    """
    The FTIS RICE1 compression and decompression algorithm.

    The Rice algorithm [1]_ is simple and very fast It requires only enough
    memory to hold a single block of 16 or 32 pixels at a time. It codes the
    pixels in small blocks and so is able to adapt very quickly to changes in
    the input image statistics (e.g., Rice has no problem handling cosmic rays,
    bright stars, saturated pixels, etc.).

    Parameters
    ----------
    blocksize
        The blocksize to use, each tile is coded into blocks a number of pixels
        wide. The default value in FITS headers is 32 pixels per block.

    bytepix
        The number of 8-bit buffer in each original integer pixel value.

    References
    ----------
    .. [1] Rice, R. F., Yeh, P.-S., and Miller, W. H. 1993, in Proc. of the 9th
    AIAA Computing in Aerospace Conf., AIAA-93-4541-CP, American Institute of
    Aeronautics and Astronautics [https://doi.org/10.2514/6.1993-4541]
    """

    codec_id = "FITS_RICE1"

    def __init__(self, blocksize: int, bytepix: int, tilesize: int):
        self.blocksize = blocksize
        self.bytepix = bytepix
        self.tilesize = tilesize

    def decode(self, buf):
        """
        Decompress buffer using the RICE_1 algorithm.

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            The decompressed buffer.
        """
        cbytes = np.frombuffer(buf, dtype=np.uint8).tobytes()
        dbytes = decompress_rice_1_c(
            cbytes, self.blocksize, self.bytepix, self.tilesize
        )
        return np.frombuffer(dbytes, dtype=f"i{self.bytepix}").data

    def encode(self, buf):
        """
        Compress the data in the buffer using the RICE_1 algorithm.

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            A buffer with decompressed data.
        """
        dbytes = np.asarray(buf).astype(f"i{self.bytepix}").tobytes()
        return compress_rice_1_c(dbytes, self.blocksize, self.bytepix)


class PLIO1(Codec):
    """
    The FTIS PLIO1 compression and decompression algorithm.

    The IRAF PLIO (pixel list) algorithm was developed to store integer-valued
    image masks in a compressed form. Such masks often have large regions of
    constant value hence are highly compressible. The compression algorithm
    used is based on run-length encoding, with the ability to dynamically
    follow level changes in the image, allowing a 16-bit encoding to be used
    regardless of the image depth.
    """

    codec_id = "FITS_PLIO1"

    def __init__(self, tilesize: int):
        self.tilesize = tilesize

    def decode(self, buf):
        """
        Decompress buffer using the PLIO_1 algorithm.

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            The decompressed buffer.
        """
        cbytes = np.frombuffer(buf, dtype=np.uint8).tobytes()
        dbytes = decompress_plio_1_c(cbytes, self.tilesize)
        return np.frombuffer(dbytes, dtype="i4").data

    def encode(self, buf):
        """
        Compress the data in the buffer using the PLIO_1 algorithm.
        """
        dbytes = np.asarray(buf).astype("i4").tobytes()
        return compress_plio_1_c(dbytes)


class HCompress1(Codec):
    """
    The FTIS PLIO1 compression and decompression algorithm.

    Hcompress is an the image compression package written by Richard L. White
    for use at the Space Telescope Science Institute. Hcompress was used to
    compress the STScI Digitized Sky Survey and has also been used to compress
    the preview images in the Hubble Data Archive.

    The technique gives very good compression for astronomical images and is
    relatively fast. The calculations are carried out using integer arithmetic
    and are entirely reversible. Consequently, the program can be used for
    either lossy or lossless compression, with no special approach needed for
    the lossless case.

    Parameters
    ----------
    scale
        The integer scale parameter determines the amount of compression. Scale
        = 0 or 1 leads to lossless compression, i.e. the decompressed image has
        exactly the same pixel values as the original image. If the scale
        factor is greater than 1 then the compression is lossy: the
        decompressed image will not be exactly the same as the original

    smooth
        At high compressions factors the decompressed image begins to appear
        blocky because of the way information is discarded. This blockiness
        ness is greatly reduced, producing more pleasing images, if the image
        is smoothed slightly during decompression.
    """

    codec_id = "FITS_HCOMPRESS1"

    def __init__(self, scale: int, smooth: bool, bytepix: int, nx: int, ny: int):
        self.scale = scale
        self.smooth = smooth
        self.bytepix = bytepix
        # NOTE: we should probably make this less confusing, but nx is shape[0] and ny is shape[1]
        self.nx = nx
        self.ny = ny

    def decode(self, buf):
        """
        Decompress buffer using the HCOMPRESS_1 algorithm.

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            A buffer with decompressed data.
        """
        cbytes = np.frombuffer(buf, dtype=np.uint8).tobytes()
        dbytes = decompress_hcompress_1_c(
            cbytes, self.nx, self.ny, self.scale, self.smooth, self.bytepix
        )
        return np.frombuffer(dbytes, dtype=f"i{self.bytepix}").data

    def encode(self, buf):
        """
        Compress the data in the buffer using the HCOMPRESS_1 algorithm.

        Parameters
        ----------
        buf
            The buffer to decompress.

        Returns
        -------
        buf
            A buffer with decompressed data.
        """
        dbytes = np.asarray(buf).astype(f"i{self.bytepix}").tobytes()
        return compress_hcompress_1_c(
            dbytes, self.nx, self.ny, self.scale, self.bytepix
        )


ALGORITHMS = {
    "GZIP_1": Gzip1,
    "GZIP_2": Gzip2,
    "RICE_1": Rice1,
    "PLIO_1": PLIO1,
    "HCOMPRESS_1": HCompress1,
}


def decompress_tile(buf, *, algorithm: str, **kwargs):
    """
    Decompress the buffer of a tile using the given compression algorithm.

    Parameters
    ----------
    buf
        The compressed buffer to be decompressed.
    algorithm
        A supported decompression algorithm.
    kwargs
        Any parameters for the given compression algorithm
    """
    return ALGORITHMS[algorithm](**kwargs).decode(buf)


def compress_tile(buf, *, algorithm: str, **kwargs):
    """
    Compress the buffer of a tile using the given compression algorithm.

    Parameters
    ----------
    buf
        The decompressed buffer to be compressed.
    algorithm
        A supported compression algorithm.
    kwargs
        Any parameters for the given compression algorithm
    """
    return ALGORITHMS[algorithm](**kwargs).encode(buf)


def _tile_shape(header):
    return tuple(header[f"ZTILE{idx}"] for idx in range(header["ZNAXIS"], 0, -1))


def _data_shape(header):
    return tuple(header[f"ZNAXIS{idx}"] for idx in range(header["ZNAXIS"], 0, -1))


def _header_to_settings(header):

    tile_shape = _tile_shape(header)

    settings = {}

    if header["ZCMPTYPE"] == "GZIP_2":
        settings["itemsize"] = header["ZBITPIX"] // 8
    elif header["ZCMPTYPE"] == "PLIO_1":
        settings["tilesize"] = np.product(tile_shape)
    elif header["ZCMPTYPE"] == "RICE_1":
        settings["blocksize"] = header.get("ZVAL1", 32)
        settings["bytepix"] = header.get("ZVAL2", 4)
        settings["tilesize"] = np.product(tile_shape)
    elif header["ZCMPTYPE"] == "HCOMPRESS_1":
        settings["bytepix"] = 4
        settings["scale"] = int(header["ZVAL1"])
        settings["smooth"] = header["ZVAL2"]
        settings["nx"] = header["ZTILE2"]
        settings["ny"] = header["ZTILE1"]

    return settings


def _buffer_to_array(tile_buffer, header, tile_shape=None):
    """
    Convert a buffer to an array using the header.

    This is a helper function which takes a raw buffer (as output by .decode)
    and using the FITS header translates it into a numpy array with the correct
    dtype, endianess and shape.
    """

    if tile_shape is None:
        tile_shape = _tile_shape(header)

    if header["ZCMPTYPE"].startswith("GZIP"):
        # This algorithm is taken from fitsio
        # https://github.com/astropy/astropy/blob/a8cb1668d4835562b89c0d0b3448ac72ca44db63/cextern/cfitsio/lib/imcompress.c#L6345-L6388
        tilelen = np.product(tile_shape)
        tilebytesize = len(tile_buffer)
        if tilebytesize == tilelen * 2:
            dtype = ">i2"
        elif tilebytesize == tilelen * 4:
            # TOOD: support float32?
            dtype = ">i4"
        elif tilebytesize == tilelen * 8:
            dtype = ">f8"
        else:
            # Just return the raw bytes
            dtype = ">u1"
        tile_data = np.asarray(tile_buffer).view(dtype).reshape(tile_shape)
    else:

        # For RICE_1 compression the tiles that are on the edge can end up
        # being padded, so we truncate excess values
        if header["ZCMPTYPE"] in ("RICE_1", "PLIO_1"):
            tile_buffer = tile_buffer[: np.product(tile_shape)]

        if tile_buffer.format == "b":
            # NOTE: this feels like a Numpy bug - need to investigate
            tile_data = np.asarray(tile_buffer, dtype=np.uint8).reshape(tile_shape)
        else:
            tile_data = np.asarray(tile_buffer).reshape(tile_shape)

    return tile_data


def _check_compressed_header(header):

    # Check for overflows which might cause issues when calling C code

    for kw in ["ZNAXIS", "ZVAL1", "ZVAL2", "ZBLANK", "BLANK"]:
        if kw in header:
            if header[kw] > 0 and np.intc(header[kw]) < 0:
                raise OverflowError()

    for i in range(1, header["ZNAXIS"] + 1):
        for kw_name in ["ZNAXIS", "ZTILE"]:
            kw = f"{kw_name}{i}"
            if kw in header:
                if header[kw] > 0 and np.int32(header[kw]) < 0:
                    raise OverflowError()

    for i in range(1, header["NAXIS"] + 1):
        kw = f"NAXIS{i}"
        if kw in header:
            if header[kw] > 0 and np.int64(header[kw]) < 0:
                raise OverflowError()

    for kw in ["TNULL1", "PCOUNT", "THEAP"]:
        if kw in header:
            if header[kw] > 0 and np.int64(header[kw]) < 0:
                raise OverflowError()

    for kw in ["ZVAL3"]:
        if kw in header:
            if np.isinf(np.float32(header[kw])):
                raise OverflowError()

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

    if header["ZCMPTYPE"] not in [
        "GZIP_1",
        "GZIP_2",
        "PLIO_1",
        "RICE_1",
        "HCOMPRESS_1",
    ]:
        raise ValueError(f"Unrecognized compression type: {header['ZCMPTYPE']}")

    # Check that certain keys are present

    header["ZNAXIS"]
    header["ZBITPIX"]


def decompress_hdu(hdu):
    """
    Drop-in replacement for decompress_hdu from compressionmodule.c
    """

    _check_compressed_header(hdu._header)

    tile_shape = _tile_shape(hdu._header)
    data_shape = _data_shape(hdu._header)

    settings = _header_to_settings(hdu._header)

    data = np.zeros(data_shape, dtype=BITPIX2DTYPE[hdu._header["ZBITPIX"]])

    istart = np.zeros(data.ndim, dtype=int)
    for cdata in hdu.compressed_data["COMPRESSED_DATA"]:
        tile_buffer = decompress_tile(
            cdata, algorithm=hdu._header["ZCMPTYPE"], **settings
        )

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

        tile_data = _buffer_to_array(
            tile_buffer, hdu._header, tile_shape=actual_tile_shape
        )

        data[tile_slices] = tile_data
        istart[-1] += tile_shape[-1]
        for idx in range(data.ndim - 1, 0, -1):
            if istart[idx] >= data_shape[idx]:
                istart[idx] = 0
                istart[idx - 1] += tile_shape[idx - 1]

    return data


def compress_hdu(hdu):
    """
    Drop-in replacement for compress_hdu from compressionmodule.c
    """

    if not isinstance(hdu.data, np.ndarray):
        raise TypeError("CompImageHDU.data must be a numpy.ndarray")

    _check_compressed_header(hdu._header)

    # For now this is very inefficient, just a proof of concept!

    settings = _header_to_settings(hdu._header)

    tile_shape = _tile_shape(hdu._header)
    data_shape = _data_shape(hdu._header)

    compressed_bytes = []

    istart = np.zeros(len(data_shape), dtype=int)

    while True:

        # In the following, we don't need to special case tiles near the edge
        # as Numpy will automatically ignore parts of the slices that are out
        # of bounds.
        slices = tuple(
            [
                slice(istart[idx], istart[idx] + tile_shape[idx])
                for idx in range(len(istart))
            ]
        )

        # TODO: deal with data not being integer number of tiles
        data = hdu.data[slices]
        # The original compress_hdu assumed the data was in native endian, so we
        # change this here:
        if hdu._header["ZCMPTYPE"].startswith("GZIP"):
            # This is apparently needed so that our heap data agrees with
            # the C implementation!?
            data = data.astype(data.dtype.newbyteorder(">"))
        else:
            if not data.dtype.isnative:
                data = data.astype(data.dtype.newbyteorder("="))
        cbytes = compress_tile(data, algorithm=hdu._header["ZCMPTYPE"], **settings)
        compressed_bytes.append(cbytes)

        istart[-1] += tile_shape[-1]
        for idx in range(data.ndim - 1, 0, -1):
            if istart[idx] >= data_shape[idx]:
                istart[idx] = 0
                istart[idx - 1] += tile_shape[idx - 1]

        if istart[0] >= data_shape[0]:
            break

    header_len = len(compressed_bytes) * 2

    heap_header = np.zeros(header_len, ">i4")
    for i in range(len(compressed_bytes)):
        heap_header[i * 2] = len(compressed_bytes[i])
        heap_header[1 + i * 2] = heap_header[: i * 2 : 2].sum()

    # For PLIO_1, the size of each heap element is a factor of two lower than
    # the real size - not clear if this is deliberate or bug somewhere.
    if hdu._header["ZCMPTYPE"] == "PLIO_1":
        for i in range(len(compressed_bytes)):
            heap_header[i * 2] /= 2

    compressed_bytes = b"".join(compressed_bytes)

    # For PLIO_1, it looks like the compressed data is byteswapped
    if hdu._header["ZCMPTYPE"] == "PLIO_1":
        compressed_bytes = (
            np.frombuffer(compressed_bytes, dtype="<i2").astype(">i2").tobytes()
        )

    heap = heap_header.tobytes() + compressed_bytes

    return len(compressed_bytes), np.frombuffer(heap, dtype=np.uint8)
