"""
This module contains the FITS compression algorithms in numcodecs style Codecs.
"""

from gzip import compress as gzip_compress
from gzip import decompress as gzip_decompress

import numpy as np

try:
    import imagecodecs

    HAS_IMAGECODECS = True
except ImportError:
    HAS_IMAGECODECS = False

from astropy.io.fits.hdu.compressed._compression import (
    compress_hcompress_1_c,
    compress_plio_1_c,
    compress_rice_1_c,
    decompress_hcompress_1_c,
    decompress_plio_1_c,
    decompress_rice_1_c,
)

# If numcodecs is installed, we use Codec as a base class for the codecs below
# so that they can optionally be used as codecs in any package relying on
# numcodecs - however this is optional and if numcodecs is not installed we use
# our own base class. This does not affect any compressed data functionality
# in astropy.io.fits.
try:
    from numcodecs.abc import Codec
except ImportError:

    class Codec:
        codec_id = None


__all__ = [
    "JPEGXL",
    "JPEGLS",
    "JPEG2K",
    "PLIO1",
    "Gzip1",
    "Gzip2",
    "HCompress1",
    "NoCompress",
    "Rice1",
]


def _as_big_endian_array(data):
    return data.astype(np.asarray(data).dtype.newbyteorder(">"), copy=False)


def _as_native_endian_array(data):
    if data.dtype.isnative:
        return data
    else:
        return data.astype(np.asarray(data).dtype.newbyteorder("="), copy=False)


class NoCompress(Codec):
    """
    A dummy compression/decompression algorithm that stores the data as-is.

    While the data is not compressed/decompressed, it is converted to big
    endian during encoding as this is what is expected in FITS files.
    """

    codec_id = "FITS_NOCOMPRESS"

    def decode(self, buf):
        """
        Decompress buffer using the NOCOMPRESS algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to decompress.

        Returns
        -------
        buf : np.ndarray
            The decompressed buffer.
        """
        return np.frombuffer(buf, dtype=np.uint8)

    def encode(self, buf):
        """
        Compress the data in the buffer using the NOCOMPRESS algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to compress.

        Returns
        -------
        bytes
            The compressed bytes.
        """
        return _as_big_endian_array(buf).tobytes()


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
        buf : bytes or array_like
            The buffer to decompress.

        Returns
        -------
        buf : np.ndarray
            The decompressed buffer.
        """
        # In principle we should be able to not have .tobytes() here and avoid
        # the copy but this does not work correctly in Python 3.11.
        cbytes = np.frombuffer(buf, dtype=np.uint8).tobytes()
        dbytes = gzip_decompress(cbytes)
        return np.frombuffer(dbytes, dtype=np.uint8)

    def encode(self, buf):
        """
        Compress the data in the buffer using the GZIP_1 algorithm.

        Parameters
        ----------
        buf _like
            The buffer to compress.

        Returns
        -------
        bytes
            The compressed bytes.
        """
        # Data bytes should be stored as big endian in files
        # In principle we should be able to not have .tobytes() here and avoid
        # the copy but this does not work correctly in Python 3.11.
        dbytes = _as_big_endian_array(buf).tobytes()
        return gzip_compress(dbytes)


class Gzip2(Codec):
    """
    The FITS GZIP2 compression and decompression algorithm.

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

    def __init__(self, *, itemsize: int):
        super().__init__()
        self.itemsize = itemsize

    def decode(self, buf):
        """
        Decompress buffer using the GZIP_2 algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to decompress.

        Returns
        -------
        buf : np.ndarray
            The decompressed buffer.
        """
        cbytes = np.frombuffer(buf, dtype=np.uint8).tobytes()
        # Start off by unshuffling buffer
        unshuffled_buffer = gzip_decompress(cbytes)
        array = np.frombuffer(unshuffled_buffer, dtype=np.uint8)
        return array.reshape((self.itemsize, -1)).T.ravel()

    def encode(self, buf):
        """
        Compress the data in the buffer using the GZIP_2 algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to compress.

        Returns
        -------
        bytes
            The compressed bytes.
        """
        # Data bytes should be stored as big endian in files
        array = _as_big_endian_array(buf).ravel()
        # Shuffle the buffer
        itemsize = array.dtype.itemsize
        array = array.view(np.uint8)
        shuffled_buffer = array.reshape((-1, itemsize)).T.ravel().tobytes()
        return gzip_compress(shuffled_buffer)


class Rice1(Codec):
    """
    The FITS RICE1 compression and decompression algorithm.

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

    def __init__(self, *, blocksize: int, bytepix: int, tilesize: int):
        self.blocksize = blocksize
        self.bytepix = bytepix
        self.tilesize = tilesize

    def decode(self, buf):
        """
        Decompress buffer using the RICE_1 algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to decompress.

        Returns
        -------
        buf : np.ndarray
            The decompressed buffer.
        """
        cbytes = np.frombuffer(_as_native_endian_array(buf), dtype=np.uint8).tobytes()
        dbytes = decompress_rice_1_c(
            cbytes, self.blocksize, self.bytepix, self.tilesize
        )
        return np.frombuffer(dbytes, dtype=f"i{self.bytepix}")

    def encode(self, buf):
        """
        Compress the data in the buffer using the RICE_1 algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to compress.

        Returns
        -------
        bytes
            The compressed bytes.
        """
        # We convert the data to native endian because it is passed to the
        # C compression code which will interpret it as being native endian.
        dbytes = (
            _as_native_endian_array(buf)
            .astype(f"i{self.bytepix}", copy=False)
            .tobytes()
        )
        return compress_rice_1_c(dbytes, self.blocksize, self.bytepix)


class PLIO1(Codec):
    """
    The FITS PLIO1 compression and decompression algorithm.

    The IRAF PLIO (pixel list) algorithm was developed to store integer-valued
    image masks in a compressed form. Such masks often have large regions of
    constant value hence are highly compressible. The compression algorithm
    used is based on run-length encoding, with the ability to dynamically
    follow level changes in the image, allowing a 16-bit encoding to be used
    regardless of the image depth.
    """

    codec_id = "FITS_PLIO1"

    def __init__(self, *, tilesize: int):
        self.tilesize = tilesize

    def decode(self, buf):
        """
        Decompress buffer using the PLIO_1 algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to decompress.

        Returns
        -------
        buf : np.ndarray
            The decompressed buffer.
        """
        cbytes = np.frombuffer(_as_native_endian_array(buf), dtype=np.uint8).tobytes()
        dbytes = decompress_plio_1_c(cbytes, self.tilesize)
        return np.frombuffer(dbytes, dtype="i4")

    def encode(self, buf):
        """
        Compress the data in the buffer using the PLIO_1 algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to compress.

        Returns
        -------
        bytes
            The compressed bytes.
        """
        # We convert the data to native endian because it is passed to the
        # C compression code which will interpret it as being native endian.
        dbytes = _as_native_endian_array(buf).astype("i4", copy=False).tobytes()
        return compress_plio_1_c(dbytes, self.tilesize)


class HCompress1(Codec):
    """
    The FITS HCompress compression and decompression algorithm.

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

    def __init__(self, *, scale: int, smooth: bool, bytepix: int, nx: int, ny: int):
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
        buf : bytes or array_like
            The buffer to decompress.

        Returns
        -------
        buf : np.ndarray
            The decompressed buffer.
        """
        cbytes = np.frombuffer(_as_native_endian_array(buf), dtype=np.uint8).tobytes()
        dbytes = decompress_hcompress_1_c(
            cbytes, self.nx, self.ny, self.scale, self.smooth, self.bytepix
        )
        # fits_hdecompress* always returns 4 byte integers irrespective of bytepix
        return np.frombuffer(dbytes, dtype="i4")

    def encode(self, buf):
        """
        Compress the data in the buffer using the HCOMPRESS_1 algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to compress.

        Returns
        -------
        bytes
            The compressed bytes.
        """
        # We convert the data to native endian because it is passed to the
        # C compression code which will interpret it as being native endian.
        dbytes = (
            _as_native_endian_array(buf)
            .astype(f"i{self.bytepix}", copy=False)
            .tobytes()
        )
        return compress_hcompress_1_c(
            dbytes, self.nx, self.ny, self.scale, self.bytepix
        )


class AbstractJPEGXLJPEGLS(Codec):
    """
    Abstract class for JPEGXL and JPEGLS compression and decompression algorithm.

    All JPEG algorithms will likely always produce better compression performance than RICE
    or HCOMPRESS. https://openreview.net/forum?id=kQCHCkNk7s&nesting=2&sort=date-desc

    Default tile shape is 512 x 512, or as close to it as possible if
    the data is smaller.

    Data types that are supported to be passed into JPEG encode / decode:
    uint8
    int16
    int32
    
    Future support for:
    float32
    int64

    For anything other than uint8:
    Note that astropy passes signed integer data to all codecs
    Since JPEG only supports compressing unsigned integer data,
    We convert the data to unsigned integer during compression
    And convert back to signed integer during decompression
    """

    def __init__(self, bytepix: int, encode_func: callable, decode_func: callable):
        self.encode_func = encode_func
        self.decode_func = decode_func
        self.bytepix = bytepix

    def decode(self, buf):
        """
        Decompress buffer using supplied JPEG algorithm.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to decompress.

        Returns
        -------
        buf : np.ndarray
            The decompressed buffer.
        """
        # Ensure input is bytes
        buf_bytes = buf if isinstance(buf, (bytes, bytearray)) else buf.tobytes()
        if self.bytepix == 1:
            decoded = self.decode_func(buf_bytes)
        elif self.bytepix == 2:
            decoded = self.decode_func(buf_bytes).astype(np.int16)
        elif self.bytepix == 4:
            length = int.from_bytes(buf_bytes[:2], "big")
            upper_enc = buf_bytes[2 : 2 + length]
            lower_enc = buf_bytes[2 + length :]
            upper = self.decode_func(upper_enc).astype(np.uint16)
            lower = self.decode_func(lower_enc).astype(np.uint16)
            dec_uint32 = (upper.astype(np.uint32) << 16) | lower.astype(np.uint32)
            decoded = dec_uint32.astype(np.int32)
        else:
            raise RuntimeError("Unsupported data type for JPEG XL compression.")
        return decoded

    def encode(self, buf):
        """
        Compress the data in the buffer.

        Parameters
        ----------
        buf : bytes or array_like
            The buffer to compress.

        Returns
        -------
        bytes
            The compressed bytes.
        """
        if self.bytepix == 1:
            encoded = self.encode_func(buf)
        elif self.bytepix == 2:
            encoded = self.encode_func(buf.astype(np.uint16))
        elif self.bytepix == 4:
            # Split into upper and lower 16 bits
            upper_bits = (buf >> 16).astype(np.uint16)
            lower_bits = (buf & 0xFFFF).astype(np.uint16)
            # Encode both parts separately
            encoded_upper = self.encode_func(upper_bits)
            encoded_lower = self.encode_func(lower_bits)
            # Store length of first part and concatenate
            length_bytes = len(encoded_upper).to_bytes(2, byteorder="big")
            encoded = length_bytes + encoded_upper + encoded_lower
        else:
            raise RuntimeError("Unsupported data type for JPEG XL compression.")
        return encoded


class JPEGXL(AbstractJPEGXLJPEGLS):
    """
    The JPEG-XL compression and decompression algorithm.

    JPEG XL is a modern image compression format designed for high-quality
    and efficient compression of photographic and other complex images. It
    offers both lossy and lossless compression modes with excellent performance
    characteristics.
    
    This implementation uses the imagecodecs library to provide JPEG XL
    compression and decompression capabilities.

    JPEG XL can process 2D or 3D arrays.

    JPEG XL already uses tiling under the hood for memory efficiency,
    but if you need to take advantage of FITS tiling for any reason,
    you can provide a tile_shape.

    Parameters
    ----------
    bytepix : int
        Number of bytes per pixel (1, 2, 4, or 8)
    effort : int
        The effort level to provide to JPEGXL for compression/decompression (0-9).
    """

    codec_id = "FITS_JPEGXL"

    def __init__(self, *, bytepix: int, effort: int):
        encode_func = lambda x: imagecodecs.jpegxl_encode(x, effort=effort)
        decode_func = imagecodecs.jpegxl_decode
        super().__init__(bytepix=bytepix, encode_func=encode_func, decode_func=decode_func)


class JPEGLS(AbstractJPEGXLJPEGLS):
    """
    The JPEG-LS compression and decompression algorithm.

    JPEG LS is a lossless compression algorithm designed to be
    extremely efficient, while still more performant than RICE1 or HCOMPRESS.
    
    This implementation uses the imagecodecs library to provide JPEGLS
    compression and decompression capabilities.

    JPEG LS can process only 2D arrays.

    Parameters
    ----------
    bytepix : int
        Number of bytes per pixel (1, 2, 4, or 8)
    naxes : int
        Number of axes in the original uncompressed data.
    """

    codec_id = "FITS_JPEGLS"

    def __init__(self, *, bytepix: int, naxes: int):
        self.naxes = naxes

        def encode_func(x):
            print("x.shape:", x.shape)
            squeezed = np.squeeze(x)
            encoded = imagecodecs.jpegls_encode(squeezed)
            return encoded

        def decode_func(x):
            decoded = imagecodecs.jpegls_decode(x)
            if self.naxes > 2:
                # Expand dims to match original shape: (1, 1, ..., H, W)
                for _ in range(self.naxes - 2):
                    decoded = np.expand_dims(decoded, axis=0)
            return decoded

        super().__init__(bytepix=bytepix, encode_func=encode_func, decode_func=decode_func)


class JPEG2K(Codec):
    """
    The JPEG 2000 compression and decompression algorithm.

    JPEG 2000 is an image compression standard and coding system. It was created
    by the Joint Photographic Experts Group committee in 2000 with the intention
    of superseding their original JPEG standard.

    This implementation uses the imagecodecs library to provide JPEG 2000
    compression and decompression capabilities.

    JPEG 2000 can process 2D or 3D arrays.

    Data types that are supported to be passed into JPEG encode / decode:
    uint8
    int16
    int32

    Parameters
    ----------
    bytepix : int
        Number of bytes per pixel (1, 2, 4, or 8)
    naxes : int
        Number of axes in the original uncompressed data.
    """

    codec_id = "FITS_JPEG2K"

    def __init__(self, *, bytepix: int, naxes: int):
        self.bytepix = bytepix
        self.naxes = naxes
    
    def decode(self, buf):
        decoded = imagecodecs.jpeg2k_decode(buf)
        if len(decoded.shape) == 3:
            # Convert from H W C to C H W
            # To undo JPEG2K axis swap during compression
            decoded = np.moveaxis(decoded, -1, 0)
        return decoded
    
    def encode(self, buf):
        encoded = imagecodecs.jpeg2k_encode(buf)
        return encoded