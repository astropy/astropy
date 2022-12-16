"""
This file contains the FITS compression algorithms in numcodecs style Codecs.

This module is currently private as these classes have not been tested
externally to Astropy and their behaviour may have to change.
"""
from gzip import compress as gzip_compress
from gzip import decompress as gzip_decompress

import numpy as np

from astropy.io.fits.tiled_compression._compression import (
    compress_hcompress_1_c,
    compress_plio_1_c,
    compress_rice_1_c,
    decompress_hcompress_1_c,
    decompress_plio_1_c,
    decompress_rice_1_c,
)

try:
    from numcodecs.abc import Codec
except ImportError:

    class Codec:
        codec_id = None


__all__ = [
    "Gzip1",
    "Gzip2",
    "Rice1",
    "PLIO1",
    "HCompress1",
]


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
        array = np.asarray(buf).ravel()
        itemsize = array.dtype.itemsize
        array = array.view(np.uint8)
        shuffled_buffer = array.reshape((-1, itemsize)).T.ravel().tobytes()
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
        return compress_plio_1_c(dbytes, self.tilesize)


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

    References
    ----------
    .. [1] White, R. L. 1992, in Proceedings of the NASA Space and Earth Science
    Data Compression Workshop, ed. J. C. Tilton, Snowbird, UT;
    https://archive.org/details/nasa_techdoc_19930016742
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
