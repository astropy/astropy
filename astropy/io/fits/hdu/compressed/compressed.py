# Licensed under a 3-clause BSD style license - see LICENSE.rst

import ctypes
import math
import time
import warnings
from contextlib import suppress

import numpy as np

from astropy.io.fits import conf
from astropy.io.fits.fitsrec import FITS_rec
from astropy.io.fits.hdu.base import BITPIX2DTYPE, DELAYED, DTYPE2BITPIX, ExtensionHDU
from astropy.io.fits.hdu.compressed._tiled_compression import compress_image_data
from astropy.io.fits.hdu.image import ImageHDU
from astropy.io.fits.hdu.table import BinTableHDU
from astropy.io.fits.util import (
    _get_array_mmap,
    _is_int,
    _is_pseudo_integer,
    _pseudo_zero,
)
from astropy.utils import lazyproperty
from astropy.utils.decorators import deprecated_renamed_argument
from astropy.utils.exceptions import AstropyUserWarning

from .header import (
    CompImageHeader,
    _bintable_header_to_image_header,
    _image_header_to_bintable_header_and_coldefs,
)
from .section import CompImageSection
from .settings import (
    CMTYPE_ALIASES,
    DEFAULT_COMPRESSION_TYPE,
    DEFAULT_DITHER_SEED,
    DEFAULT_HCOMP_SCALE,
    DEFAULT_HCOMP_SMOOTH,
    DEFAULT_QUANTIZE_LEVEL,
    DEFAULT_QUANTIZE_METHOD,
    DITHER_SEED_CHECKSUM,
    DITHER_SEED_CLOCK,
)

# This global variable is used e.g., when calling fits.open with
# disable_image_compression which temporarily changes the global variable to
# False. This should ideally be refactored to avoid relying on global module
# variables.
COMPRESSION_ENABLED = True


# TODO: Fix this class so that it doesn't actually inherit from BinTableHDU,
# but instead has an internal BinTableHDU reference
class CompImageHDU(BinTableHDU):
    """
    Compressed Image HDU class.
    """

    _manages_own_heap = True
    """
    The calls to CFITSIO lay out the heap data in memory, and we write it out
    the same way CFITSIO organizes it.  In principle this would break if a user
    manually changes the underlying compressed data by hand, but there is no
    reason they would want to do that (and if they do that's their
    responsibility).
    """

    _load_variable_length_data = False
    """
    We don't want to always load all the tiles so by setting this option
    we can then access the tiles as needed.
    """

    _default_name = "COMPRESSED_IMAGE"

    @deprecated_renamed_argument(
        "tile_size",
        None,
        since="5.3",
        message="The tile_size argument has been deprecated. Use tile_shape "
        "instead, but note that this should be given in the reverse "
        "order to tile_size (tile_shape should be in Numpy C order).",
    )
    def __init__(
        self,
        data=None,
        header=None,
        name=None,
        compression_type=DEFAULT_COMPRESSION_TYPE,
        tile_shape=None,
        hcomp_scale=DEFAULT_HCOMP_SCALE,
        hcomp_smooth=DEFAULT_HCOMP_SMOOTH,
        quantize_level=DEFAULT_QUANTIZE_LEVEL,
        quantize_method=DEFAULT_QUANTIZE_METHOD,
        dither_seed=DEFAULT_DITHER_SEED,
        do_not_scale_image_data=False,
        uint=False,
        scale_back=False,
        tile_size=None,
    ):
        """
        Parameters
        ----------
        data : array, optional
            Uncompressed image data

        header : `~astropy.io.fits.Header`, optional
            Header to be associated with the image; when reading the HDU from a
            file (data=DELAYED), the header read from the file

        name : str, optional
            The ``EXTNAME`` value; if this value is `None`, then the name from
            the input image header will be used; if there is no name in the
            input image header then the default name ``COMPRESSED_IMAGE`` is
            used.

        compression_type : str, optional
            Compression algorithm: one of
            ``'RICE_1'``, ``'RICE_ONE'``, ``'PLIO_1'``, ``'GZIP_1'``,
            ``'GZIP_2'``, ``'HCOMPRESS_1'``, ``'NOCOMPRESS'``

        tile_shape : tuple, optional
            Compression tile shape, which should be specified using the default
            Numpy convention for array shapes (C order). The default is to
            treat each row of image as a tile.

        hcomp_scale : float, optional
            HCOMPRESS scale parameter

        hcomp_smooth : float, optional
            HCOMPRESS smooth parameter

        quantize_level : float, optional
            Floating point quantization level; see note below

        quantize_method : int, optional
            Floating point quantization dithering method; can be either
            ``NO_DITHER`` (-1; default), ``SUBTRACTIVE_DITHER_1`` (1), or
            ``SUBTRACTIVE_DITHER_2`` (2); see note below

        dither_seed : int, optional
            Random seed to use for dithering; can be either an integer in the
            range 1 to 1000 (inclusive), ``DITHER_SEED_CLOCK`` (0; default), or
            ``DITHER_SEED_CHECKSUM`` (-1); see note below

        Notes
        -----
        The astropy.io.fits package supports 2 methods of image compression:

            1) The entire FITS file may be externally compressed with the gzip
               or pkzip utility programs, producing a ``*.gz`` or ``*.zip``
               file, respectively.  When reading compressed files of this type,
               Astropy first uncompresses the entire file into a temporary file
               before performing the requested read operations.  The
               astropy.io.fits package does not support writing to these types
               of compressed files.  This type of compression is supported in
               the ``_File`` class, not in the `CompImageHDU` class.  The file
               compression type is recognized by the ``.gz`` or ``.zip`` file
               name extension.

            2) The `CompImageHDU` class supports the FITS tiled image
               compression convention in which the image is subdivided into a
               grid of rectangular tiles, and each tile of pixels is
               individually compressed.  The details of this FITS compression
               convention are described at the `FITS Support Office web site
               <https://fits.gsfc.nasa.gov/registry/tilecompression.html>`_.
               Basically, the compressed image tiles are stored in rows of a
               variable length array column in a FITS binary table.  The
               astropy.io.fits recognizes that this binary table extension
               contains an image and treats it as if it were an image
               extension.  Under this tile-compression format, FITS header
               keywords remain uncompressed.  At this time, Astropy does not
               support the ability to extract and uncompress sections of the
               image without having to uncompress the entire image.

        The astropy.io.fits package supports 3 general-purpose compression
        algorithms plus one other special-purpose compression technique that is
        designed for data masks with positive integer pixel values.  The 3
        general purpose algorithms are GZIP, Rice, and HCOMPRESS, and the
        special-purpose technique is the IRAF pixel list compression technique
        (PLIO).  The ``compression_type`` parameter defines the compression
        algorithm to be used.

        The FITS image can be subdivided into any desired rectangular grid of
        compression tiles.  With the GZIP, Rice, and PLIO algorithms, the
        default is to take each row of the image as a tile.  The HCOMPRESS
        algorithm is inherently 2-dimensional in nature, so the default in this
        case is to take 16 rows of the image per tile.  In most cases, it makes
        little difference what tiling pattern is used, so the default tiles are
        usually adequate.  In the case of very small images, it could be more
        efficient to compress the whole image as a single tile.  Note that the
        image dimensions are not required to be an integer multiple of the tile
        dimensions; if not, then the tiles at the edges of the image will be
        smaller than the other tiles.  The ``tile_shape`` parameter may be
        provided as a list of tile sizes, one for each dimension in the image.
        For example a ``tile_shape`` value of ``(100,100)`` would divide a 300 X
        300 image into 9 100 X 100 tiles.

        The 4 supported image compression algorithms are all 'lossless' when
        applied to integer FITS images; the pixel values are preserved exactly
        with no loss of information during the compression and uncompression
        process.  In addition, the HCOMPRESS algorithm supports a 'lossy'
        compression mode that will produce larger amount of image compression.
        This is achieved by specifying a non-zero value for the ``hcomp_scale``
        parameter.  Since the amount of compression that is achieved depends
        directly on the RMS noise in the image, it is usually more convenient
        to specify the ``hcomp_scale`` factor relative to the RMS noise.
        Setting ``hcomp_scale = 2.5`` means use a scale factor that is 2.5
        times the calculated RMS noise in the image tile.  In some cases it may
        be desirable to specify the exact scaling to be used, instead of
        specifying it relative to the calculated noise value.  This may be done
        by specifying the negative of the desired scale value (typically in the
        range -2 to -100).

        Very high compression factors (of 100 or more) can be achieved by using
        large ``hcomp_scale`` values, however, this can produce undesirable
        'blocky' artifacts in the compressed image.  A variation of the
        HCOMPRESS algorithm (called HSCOMPRESS) can be used in this case to
        apply a small amount of smoothing of the image when it is uncompressed
        to help cover up these artifacts.  This smoothing is purely cosmetic
        and does not cause any significant change to the image pixel values.
        Setting the ``hcomp_smooth`` parameter to 1 will engage the smoothing
        algorithm.

        Floating point FITS images (which have ``BITPIX`` = -32 or -64) usually
        contain too much 'noise' in the least significant bits of the mantissa
        of the pixel values to be effectively compressed with any lossless
        algorithm.  Consequently, floating point images are first quantized
        into scaled integer pixel values (and thus throwing away much of the
        noise) before being compressed with the specified algorithm (either
        GZIP, RICE, or HCOMPRESS).  This technique produces much higher
        compression factors than simply using the GZIP utility to externally
        compress the whole FITS file, but it also means that the original
        floating point value pixel values are not exactly preserved.  When done
        properly, this integer scaling technique will only discard the
        insignificant noise while still preserving all the real information in
        the image.  The amount of precision that is retained in the pixel
        values is controlled by the ``quantize_level`` parameter.  Larger
        values will result in compressed images whose pixels more closely match
        the floating point pixel values, but at the same time the amount of
        compression that is achieved will be reduced.  Users should experiment
        with different values for this parameter to determine the optimal value
        that preserves all the useful information in the image, without
        needlessly preserving all the 'noise' which will hurt the compression
        efficiency.

        The default value for the ``quantize_level`` scale factor is 16, which
        means that scaled integer pixel values will be quantized such that the
        difference between adjacent integer values will be 1/16th of the noise
        level in the image background.  An optimized algorithm is used to
        accurately estimate the noise in the image.  As an example, if the RMS
        noise in the background pixels of an image = 32.0, then the spacing
        between adjacent scaled integer pixel values will equal 2.0 by default.
        Note that the RMS noise is independently calculated for each tile of
        the image, so the resulting integer scaling factor may fluctuate
        slightly for each tile.  In some cases, it may be desirable to specify
        the exact quantization level to be used, instead of specifying it
        relative to the calculated noise value.  This may be done by specifying
        the negative of desired quantization level for the value of
        ``quantize_level``.  In the previous example, one could specify
        ``quantize_level = -2.0`` so that the quantized integer levels differ
        by 2.0.  Larger negative values for ``quantize_level`` means that the
        levels are more coarsely-spaced, and will produce higher compression
        factors.

        The quantization algorithm can also apply one of two random dithering
        methods in order to reduce bias in the measured intensity of background
        regions.  The default method, specified with the constant
        ``SUBTRACTIVE_DITHER_1`` adds dithering to the zero-point of the
        quantization array itself rather than adding noise to the actual image.
        The random noise is added on a pixel-by-pixel basis, so in order
        restore each pixel from its integer value to its floating point value
        it is necessary to replay the same sequence of random numbers for each
        pixel (see below).  The other method, ``SUBTRACTIVE_DITHER_2``, is
        exactly like the first except that before dithering any pixel with a
        floating point value of ``0.0`` is replaced with the special integer
        value ``-2147483647``.  When the image is uncompressed, pixels with
        this value are restored back to ``0.0`` exactly.  Finally, a value of
        ``NO_DITHER`` disables dithering entirely.

        As mentioned above, when using the subtractive dithering algorithm it
        is necessary to be able to generate a (pseudo-)random sequence of noise
        for each pixel, and replay that same sequence upon decompressing.  To
        facilitate this, a random seed between 1 and 10000 (inclusive) is used
        to seed a random number generator, and that seed is stored in the
        ``ZDITHER0`` keyword in the header of the compressed HDU.  In order to
        use that seed to generate the same sequence of random numbers the same
        random number generator must be used at compression and decompression
        time; for that reason the tiled image convention provides an
        implementation of a very simple pseudo-random number generator.  The
        seed itself can be provided in one of three ways, controllable by the
        ``dither_seed`` argument:  It may be specified manually, or it may be
        generated arbitrarily based on the system's clock
        (``DITHER_SEED_CLOCK``) or based on a checksum of the pixels in the
        image's first tile (``DITHER_SEED_CHECKSUM``).  The clock-based method
        is the default, and is sufficient to ensure that the value is
        reasonably "arbitrary" and that the same seed is unlikely to be
        generated sequentially.  The checksum method, on the other hand,
        ensures that the same seed is used every time for a specific image.
        This is particularly useful for software testing as it ensures that the
        same image will always use the same seed.
        """
        compression_type = CMTYPE_ALIASES.get(compression_type, compression_type)

        if tile_shape is None and tile_size is not None:
            tile_shape = tuple(tile_size[::-1])
        elif tile_shape is not None and tile_size is not None:
            raise ValueError(
                "Cannot specify both tile_size and tile_shape. "
                "Note that tile_size is deprecated and tile_shape "
                "alone should be used."
            )

        if data is DELAYED:
            # Reading the HDU from a file
            super().__init__(data=data, header=header)
        else:
            # Create at least a skeleton HDU that matches the input
            # header and data (if any were input)
            super().__init__(data=None, header=header)

            # Store the input image data
            self.data = data

            # Update the table header (_header) to the compressed
            # image format and to match the input data (if any);
            # Create the image header (_image_header) from the input
            # image header (if any) and ensure it matches the input
            # data; Create the initially empty table data array to
            # hold the compressed data.
            self._update_header_data(
                header,
                name,
                compression_type=compression_type,
                tile_shape=tile_shape,
                hcomp_scale=hcomp_scale,
                hcomp_smooth=hcomp_smooth,
                quantize_level=quantize_level,
                quantize_method=quantize_method,
                dither_seed=dither_seed,
            )

        # TODO: A lot of this should be passed on to an internal image HDU o
        # something like that, see ticket #88
        self._do_not_scale_image_data = do_not_scale_image_data
        self._uint = uint
        self._scale_back = scale_back

        self._axes = [
            self._header.get("ZNAXIS" + str(axis + 1), 0)
            for axis in range(self._header.get("ZNAXIS", 0))
        ]

        # store any scale factors from the table header
        if do_not_scale_image_data:
            self._bzero = 0
            self._bscale = 1
        else:
            self._bzero = self._header.get("BZERO", 0)
            self._bscale = self._header.get("BSCALE", 1)
        self._bitpix = self._header["ZBITPIX"]

        self._orig_bzero = self._bzero
        self._orig_bscale = self._bscale
        self._orig_bitpix = self._bitpix

        if (
            self._bitpix > 0
            and "BLANK" not in self._header
            and "ZBLANK" not in self._header
        ):
            # check for column named "ZBLANK"
            for i in range(1, self._header["TFIELDS"] + 1):
                if self._header[f"TTYPE{i}"] == "ZBLANK":
                    # required BLANK keyword is missing
                    # use most negative value as default
                    self._header["BLANK"] = -(1 << (self._bitpix - 1))
                    warnings.warn(
                        f"Setting default value {self._header['BLANK']} for missing BLANK keyword in compressed extension",
                        AstropyUserWarning,
                    )
                    break

        if self._bitpix > 0:
            self._blank = self._header.get("BLANK", self._header.get("ZBLANK"))
        else:
            self._blank = None

    def _remove_unnecessary_default_extnames(self, header):
        """Remove default EXTNAME values if they are unnecessary.

        Some data files (eg from CFHT) can have the default EXTNAME and
        an explicit value.  This method removes the default if a more
        specific header exists. It also removes any duplicate default
        values.
        """
        if "EXTNAME" in header:
            indices = header._keyword_indices["EXTNAME"]
            # Only continue if there is more than one found
            n_extname = len(indices)
            if n_extname > 1:
                extnames_to_remove = [
                    index for index in indices if header[index] == self._default_name
                ]
                if len(extnames_to_remove) == n_extname:
                    # Keep the first (they are all the same)
                    extnames_to_remove.pop(0)
                # Remove them all in reverse order to keep the index unchanged.
                for index in sorted(extnames_to_remove, reverse=True):
                    del header[index]

    @property
    def name(self):
        # Convert the value to a string to be flexible in some pathological
        # cases (see ticket #96)
        # Similar to base class but uses .header rather than ._header
        return str(self.header.get("EXTNAME", self._default_name))

    @name.setter
    def name(self, value):
        # This is a copy of the base class but using .header instead
        # of ._header to ensure that the name stays in sync.
        if not isinstance(value, str):
            raise TypeError("'name' attribute must be a string")
        if not conf.extension_name_case_sensitive:
            value = value.upper()
        if "EXTNAME" in self.header:
            self.header["EXTNAME"] = value
        else:
            self.header["EXTNAME"] = (value, "extension name")

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        if card.keyword != "XTENSION":
            return False

        xtension = card.value
        if isinstance(xtension, str):
            xtension = xtension.rstrip()

        if xtension not in ("BINTABLE", "A3DTABLE"):
            return False

        if "ZIMAGE" not in header or not header["ZIMAGE"]:
            return False

        return COMPRESSION_ENABLED

    def _update_header_data(
        self,
        image_header,
        name=None,
        compression_type=None,
        tile_shape=None,
        hcomp_scale=None,
        hcomp_smooth=None,
        quantize_level=None,
        quantize_method=None,
        dither_seed=None,
    ):
        """
        Update the table header (`_header`) to the compressed
        image format and to match the input data (if any).  Create
        the image header (`_image_header`) from the input image
        header (if any) and ensure it matches the input
        data. Create the initially-empty table data array to hold
        the compressed data.

        This method is mainly called internally, but a user may wish to
        call this method after assigning new data to the `CompImageHDU`
        object that is of a different type.

        Parameters
        ----------
        image_header : `~astropy.io.fits.Header`
            header to be associated with the image

        name : str, optional
            the ``EXTNAME`` value; if this value is `None`, then the name from
            the input image header will be used; if there is no name in the
            input image header then the default name 'COMPRESSED_IMAGE' is used

        compression_type : str, optional
            compression algorithm 'RICE_1', 'PLIO_1', 'GZIP_1', 'GZIP_2',
            'HCOMPRESS_1', 'NOCOMPRESS'; if this value is `None`, use value
            already in the header; if no value already in the header, use
            'RICE_1'

        tile_shape : tuple of int, optional
            compression tile shape (in C order); if this value is `None`, use
            value already in the header; if no value already in the header,
            treat each row of image as a tile

        hcomp_scale : float, optional
            HCOMPRESS scale parameter; if this value is `None`, use the value
            already in the header; if no value already in the header, use 1

        hcomp_smooth : float, optional
            HCOMPRESS smooth parameter; if this value is `None`, use the value
            already in the header; if no value already in the header, use 0

        quantize_level : float, optional
            floating point quantization level; if this value is `None`, use the
            value already in the header; if no value already in header, use 16

        quantize_method : int, optional
            floating point quantization dithering method; can be either
            NO_DITHER (-1), SUBTRACTIVE_DITHER_1 (1; default), or
            SUBTRACTIVE_DITHER_2 (2)

        dither_seed : int, optional
            random seed to use for dithering; can be either an integer in the
            range 1 to 1000 (inclusive), DITHER_SEED_CLOCK (0; default), or
            DITHER_SEED_CHECKSUM (-1)
        """
        # Clean up EXTNAME duplicates
        self._remove_unnecessary_default_extnames(self._header)

        image_hdu = ImageHDU(data=self.data, header=self._header)
        self._image_header = CompImageHeader(self._header, image_hdu.header)
        self._axes = image_hdu._axes
        del image_hdu

        # Determine based on the size of the input data whether to use the Q
        # column format to store compressed data or the P format.
        # The Q format is used only if the uncompressed data is larger than
        # 4 GB.  This is not a perfect heuristic, as one can contrive an input
        # array which, when compressed, the entire binary table representing
        # the compressed data is larger than 4GB.  That said, this is the same
        # heuristic used by CFITSIO, so this should give consistent results.
        # And the cases where this heuristic is insufficient are extreme and
        # almost entirely contrived corner cases, so it will do for now
        if self._has_data:
            huge_hdu = self.data.nbytes > 2**32
        else:
            huge_hdu = False

        # NOTE: for now the function below modifies the compressed binary table
        # self._header in-place, but this could be refactored in future to
        # return the compressed header.

        self._header, self.columns = _image_header_to_bintable_header_and_coldefs(
            image_header,
            self._image_header,
            self._header,
            name=name,
            huge_hdu=huge_hdu,
            compression_type=compression_type,
            tile_shape=tile_shape,
            hcomp_scale=hcomp_scale,
            hcomp_smooth=hcomp_smooth,
            quantize_level=quantize_level,
            quantize_method=quantize_method,
            dither_seed=dither_seed,
            axes=self._axes,
            generate_dither_seed=self._generate_dither_seed,
        )

        if name:
            self.name = name

    def _scale_data(self, data):
        if self._orig_bzero != 0 or self._orig_bscale != 1 or self._blank is not None:
            if self._blank is not None:
                blanks = data == np.array(self._blank, dtype=data.dtype)
            else:
                blanks = None

            new_dtype = self._dtype_for_bitpix()
            data = np.array(data, dtype=new_dtype)

            if self._orig_bscale != 1:
                np.multiply(data, self._orig_bscale, data)
            if self._orig_bzero != 0:
                # We have to explicitly cast self._bzero to prevent numpy from
                # raising an error when doing self.data += self._bzero, and we
                # do this instead of self.data = self.data + self._bzero to
                # avoid doubling memory usage.
                np.add(data, self._orig_bzero, out=data, casting="unsafe")

            if blanks is not None:
                # use float32 version of nan to reduce data size for uint conversion
                # result will still be float64 for larger uint sizes (e.g., uint32)
                # for float types np.where retains the type of the data array
                data = np.where(blanks, np.float32(np.nan), data)

        return data

    @lazyproperty
    def data(self):
        """
        The decompressed data array.

        Note that accessing this will cause all the tiles to be loaded,
        decompressed, and combined into a single data array. If you do
        not need to access the whole array, consider instead using the
        :attr:`~astropy.io.fits.CompImageHDU.section` property.
        """
        if len(self.compressed_data) == 0:
            return None

        # Since .section has general code to load any arbitrary part of the
        # data, we can just use this - and the @lazyproperty on the current
        # property will ensure that we do this only once.
        data = self.section[...]

        # Right out of _ImageBaseHDU.data
        self._update_header_scale_info(data.dtype)

        return data

    @data.setter
    def data(self, data):
        if (data is not None) and (
            not isinstance(data, np.ndarray) or data.dtype.fields is not None
        ):
            raise TypeError(
                f"CompImageHDU data has incorrect type:{type(data)}; "
                f"dtype.fields = {data.dtype.fields}"
            )

    @lazyproperty
    def compressed_data(self):
        # First we will get the table data (the compressed
        # data) from the file, if there is any.
        compressed_data = super().data
        if isinstance(compressed_data, np.rec.recarray):
            # Make sure not to use 'del self.data' so we don't accidentally
            # go through the self.data.fdel and close the mmap underlying
            # the compressed_data array
            del self.__dict__["data"]
            return compressed_data
        else:
            # This will actually set self.compressed_data with the
            # pre-allocated space for the compression data; this is something I
            # might do away with in the future
            self._update_compressed_data()

        return self.compressed_data

    @compressed_data.deleter
    def compressed_data(self):
        # Deleting the compressed_data attribute has to be handled
        # with a little care to prevent a reference leak
        # First delete the ._coldefs attributes under it to break a possible
        # reference cycle
        if "compressed_data" in self.__dict__:
            del self.__dict__["compressed_data"]._coldefs

            # Now go ahead and delete from self.__dict__; normally
            # lazyproperty.__delete__ does this for us, but we can prempt it to
            # do some additional cleanup
            del self.__dict__["compressed_data"]

    @property
    def shape(self):
        """
        Shape of the image array--should be equivalent to ``self.data.shape``.
        """
        # Determine from the values read from the header
        return tuple(reversed(self._axes))

    @lazyproperty
    def header(self):
        # The header attribute is the header for the image data.  It
        # is not actually stored in the object dictionary.  Instead,
        # the _image_header is stored.  If the _image_header attribute
        # has already been defined we just return it.  If not, we must
        # create it from the table header (the _header attribute).
        if hasattr(self, "_image_header"):
            return self._image_header

        # Clean up any possible doubled EXTNAME keywords that use
        # the default. Do this on the original header to ensure
        # duplicates are removed cleanly.
        self._remove_unnecessary_default_extnames(self._header)

        # Convert compressed header to image header and save
        # it off to self._image_header so it can be referenced later
        # unambiguously
        self._image_header = _bintable_header_to_image_header(self._header)

        return self._image_header

    def _summary(self):
        """
        Summarize the HDU: name, dimensions, and formats.
        """
        class_name = self.__class__.__name__

        # if data is touched, use data info.
        if self._data_loaded:
            if self.data is None:
                _shape, _format = (), ""
            else:
                # the shape will be in the order of NAXIS's which is the
                # reverse of the numarray shape
                _shape = list(self.data.shape)
                _format = self.data.dtype.name
                _shape.reverse()
                _shape = tuple(_shape)
                _format = _format[_format.rfind(".") + 1 :]

        # if data is not touched yet, use header info.
        else:
            _shape = ()

            for idx in range(self.header["NAXIS"]):
                _shape += (self.header["NAXIS" + str(idx + 1)],)

            _format = BITPIX2DTYPE[self.header["BITPIX"]]

        return (self.name, self.ver, class_name, len(self.header), _shape, _format)

    def _update_compressed_data(self):
        """
        Compress the image data so that it may be written to a file.
        """
        # Check to see that the image_header matches the image data
        image_bitpix = DTYPE2BITPIX[self.data.dtype.name]

        if image_bitpix != self._orig_bitpix or self.data.shape != self.shape:
            self._update_header_data(self.header)

        # TODO: This is copied right out of _ImageBaseHDU._writedata_internal;
        # it would be cool if we could use an internal ImageHDU and use that to
        # write to a buffer for compression or something. See ticket #88
        # deal with unsigned integer 16, 32 and 64 data
        old_data = self.data
        if _is_pseudo_integer(self.data.dtype):
            # Convert the unsigned array to signed
            self.data = np.array(
                self.data - _pseudo_zero(self.data.dtype),
                dtype=f"=i{self.data.dtype.itemsize}",
            )

        try:
            nrows = self._header["NAXIS2"]
            tbsize = self._header["NAXIS1"] * nrows

            self._header["PCOUNT"] = 0
            if "THEAP" in self._header:
                del self._header["THEAP"]
            self._theap = tbsize

            # First delete the original compressed data, if it exists
            del self.compressed_data

            # Compress the data.
            # compress_image_data returns the size of the heap for the written
            # compressed image table
            heapsize, self.compressed_data = compress_image_data(
                self.data, self.compression_type, self._header, self.columns
            )
        finally:
            self.data = old_data

        table_len = len(self.compressed_data) - heapsize
        if table_len != self._theap:
            raise Exception(
                f"Unexpected compressed table size (expected {self._theap}, got {table_len})"
            )

        # CFITSIO will write the compressed data in big-endian order
        dtype = self.columns.dtype.newbyteorder(">")
        buf = self.compressed_data
        compressed_data = buf[: self._theap].view(dtype=dtype, type=np.rec.recarray)
        self.compressed_data = compressed_data.view(FITS_rec)
        self.compressed_data._coldefs = self.columns
        self.compressed_data._heapoffset = self._theap
        self.compressed_data._heapsize = heapsize

    def scale(self, type=None, option="old", bscale=1, bzero=0):
        """
        Scale image data by using ``BSCALE`` and ``BZERO``.

        Calling this method will scale ``self.data`` and update the keywords of
        ``BSCALE`` and ``BZERO`` in ``self._header`` and ``self._image_header``.
        This method should only be used right before writing to the output
        file, as the data will be scaled and is therefore not very usable after
        the call.

        Parameters
        ----------
        type : str, optional
            destination data type, use a string representing a numpy dtype
            name, (e.g. ``'uint8'``, ``'int16'``, ``'float32'`` etc.).  If is
            `None`, use the current data type.

        option : str, optional
            how to scale the data: if ``"old"``, use the original ``BSCALE``
            and ``BZERO`` values when the data was read/created. If
            ``"minmax"``, use the minimum and maximum of the data to scale.
            The option will be overwritten by any user-specified bscale/bzero
            values.

        bscale, bzero : int, optional
            user specified ``BSCALE`` and ``BZERO`` values.
        """
        if self.data is None:
            return

        # Determine the destination (numpy) data type
        if type is None:
            type = BITPIX2DTYPE[self._bitpix]
        _type = getattr(np, type)

        # Determine how to scale the data
        # bscale and bzero takes priority
        if bscale != 1 or bzero != 0:
            _scale = bscale
            _zero = bzero
        else:
            if option == "old":
                _scale = self._orig_bscale
                _zero = self._orig_bzero
            elif option == "minmax":
                if isinstance(_type, np.floating):
                    _scale = 1
                    _zero = 0
                else:
                    _min = np.minimum.reduce(self.data.flat)
                    _max = np.maximum.reduce(self.data.flat)

                    if _type == np.uint8:  # uint8 case
                        _zero = _min
                        _scale = (_max - _min) / (2.0**8 - 1)
                    else:
                        _zero = (_max + _min) / 2.0

                        # throw away -2^N
                        _scale = (_max - _min) / (2.0 ** (8 * _type.bytes) - 2)

        # Do the scaling
        if _zero != 0:
            # We have to explicitly cast self._bzero to prevent numpy from
            # raising an error when doing self.data -= _zero, and we
            # do this instead of self.data = self.data - _zero to
            # avoid doubling memory usage.
            np.subtract(self.data, _zero, out=self.data, casting="unsafe")
            self.header["BZERO"] = _zero
        else:
            # Delete from both headers
            for header in (self.header, self._header):
                with suppress(KeyError):
                    del header["BZERO"]

        if _scale != 1:
            self.data /= _scale
            self.header["BSCALE"] = _scale
        else:
            for header in (self.header, self._header):
                with suppress(KeyError):
                    del header["BSCALE"]

        if self.data.dtype.type != _type:
            self.data = np.array(np.around(self.data), dtype=_type)  # 0.7.7.1

        # Update the BITPIX Card to match the data
        self._bitpix = DTYPE2BITPIX[self.data.dtype.name]
        self._bzero = self.header.get("BZERO", 0)
        self._bscale = self.header.get("BSCALE", 1)
        # Update BITPIX for the image header specifically
        # TODO: Make this more clear by using self._image_header, but only once
        # this has been fixed so that the _image_header attribute is guaranteed
        # to be valid
        self.header["BITPIX"] = self._bitpix

        # Update the table header to match the scaled data
        self._update_header_data(self.header)

        # Since the image has been manually scaled, the current
        # bitpix/bzero/bscale now serve as the 'original' scaling of the image,
        # as though the original image has been completely replaced
        self._orig_bitpix = self._bitpix
        self._orig_bzero = self._bzero
        self._orig_bscale = self._bscale

    def _prewriteto(self, checksum=False, inplace=False):
        if self._scale_back:
            self.scale(BITPIX2DTYPE[self._orig_bitpix])

        if self._has_data:
            self._update_compressed_data()

            # Use methods in the superclass to update the header with
            # scale/checksum keywords based on the data type of the image data
            self._update_pseudo_int_scale_keywords()

            # Shove the image header and data into a new ImageHDU and use that
            # to compute the image checksum
            image_hdu = ImageHDU(data=self.data, header=self.header)
            image_hdu._update_checksum(checksum)
            if "CHECKSUM" in image_hdu.header:
                # This will also pass through to the ZHECKSUM keyword and
                # ZDATASUM keyword
                self._image_header.set(
                    "CHECKSUM",
                    image_hdu.header["CHECKSUM"],
                    image_hdu.header.comments["CHECKSUM"],
                )
            if "DATASUM" in image_hdu.header:
                self._image_header.set(
                    "DATASUM",
                    image_hdu.header["DATASUM"],
                    image_hdu.header.comments["DATASUM"],
                )
            # Store a temporary backup of self.data in a different attribute;
            # see below
            self._imagedata = self.data

            # Now we need to perform an ugly hack to set the compressed data as
            # the .data attribute on the HDU so that the call to _writedata
            # handles it properly
            self.__dict__["data"] = self.compressed_data

        return super()._prewriteto(checksum=checksum, inplace=inplace)

    def _writeheader(self, fileobj):
        """
        Bypasses `BinTableHDU._writeheader()` which updates the header with
        metadata about the data that is meaningless here; another reason
        why this class maybe shouldn't inherit directly from BinTableHDU...
        """
        return ExtensionHDU._writeheader(self, fileobj)

    def _writedata(self, fileobj):
        """
        Wrap the basic ``_writedata`` method to restore the ``.data``
        attribute to the uncompressed image data in the case of an exception.
        """
        try:
            return super()._writedata(fileobj)
        finally:
            # Restore the .data attribute to its rightful value (if any)
            if hasattr(self, "_imagedata"):
                self.__dict__["data"] = self._imagedata
                del self._imagedata
            else:
                del self.data

    def _close(self, closed=True):
        super()._close(closed=closed)

        # Also make sure to close access to the compressed data mmaps
        if (
            closed
            and self._data_loaded
            and _get_array_mmap(self.compressed_data) is not None
        ):
            del self.compressed_data

    # TODO: This was copied right out of _ImageBaseHDU; get rid of it once we
    # find a way to rewrite this class as either a subclass or wrapper for an
    # ImageHDU
    def _dtype_for_bitpix(self):
        """
        Determine the dtype that the data should be converted to depending on
        the BITPIX value in the header, and possibly on the BSCALE value as
        well.  Returns None if there should not be any change.
        """
        bitpix = self._orig_bitpix
        # Handle possible conversion to uints if enabled
        if self._uint and self._orig_bscale == 1:
            for bits, dtype in (
                (16, np.dtype("uint16")),
                (32, np.dtype("uint32")),
                (64, np.dtype("uint64")),
            ):
                if bitpix == bits and self._orig_bzero == 1 << (bits - 1):
                    return dtype

        if bitpix > 16:  # scale integers to Float64
            return np.dtype("float64")
        elif bitpix > 0:  # scale integers to Float32
            return np.dtype("float32")

    def _update_header_scale_info(self, dtype=None):
        if not self._do_not_scale_image_data and not (
            self._orig_bzero == 0 and self._orig_bscale == 1
        ):
            for keyword in ["BSCALE", "BZERO"]:
                # Make sure to delete from both the image header and the table
                # header; later this will be streamlined
                for header in (self.header, self._header):
                    with suppress(KeyError):
                        del header[keyword]
                        # Since _update_header_scale_info can, currently, be
                        # called *after* _prewriteto(), replace these with
                        # blank cards so the header size doesn't change
                        header.append()

            if dtype is None:
                dtype = self._dtype_for_bitpix()
            if dtype is not None:
                self.header["BITPIX"] = DTYPE2BITPIX[dtype.name]

            self._bzero = 0
            self._bscale = 1
            self._bitpix = self.header["BITPIX"]

    def _generate_dither_seed(self, seed):
        if not _is_int(seed):
            raise TypeError("Seed must be an integer")

        if not -1 <= seed <= 10000:
            raise ValueError(
                "Seed for random dithering must be either between 1 and "
                "10000 inclusive, 0 for autogeneration from the system "
                "clock, or -1 for autogeneration from a checksum of the first "
                f"image tile (got {seed})"
            )

        if seed == DITHER_SEED_CHECKSUM:
            # Determine the tile dimensions from the ZTILEn keywords
            tile_dims = self.tile_shape

            # Get the first tile by using the tile dimensions as the end
            # indices of slices (starting from 0)
            first_tile = self.data[tuple(slice(d) for d in tile_dims)]

            # The checksum algorithm used is literally just the sum of the bytes
            # of the tile data (not its actual floating point values).  Integer
            # overflow is irrelevant.
            csum = first_tile.view(dtype="uint8").sum()

            # Since CFITSIO uses an unsigned long (which may be different on
            # different platforms) go ahead and truncate the sum to its
            # unsigned long value and take the result modulo 10000
            return (ctypes.c_ulong(csum).value % 10000) + 1
        elif seed == DITHER_SEED_CLOCK:
            # This isn't exactly the same algorithm as CFITSIO, but that's okay
            # since the result is meant to be arbitrary. The primary difference
            # is that CFITSIO incorporates the HDU number into the result in
            # the hopes of heading off the possibility of the same seed being
            # generated for two HDUs at the same time.  Here instead we just
            # add in the HDU object's id
            return (
                (sum(int(x) for x in math.modf(time.time())) + id(self)) % 10000
            ) + 1
        else:
            return seed

    @property
    def section(self):
        """
        Efficiently access a section of the image array

        This property can be used to access a section of the data without
        loading and decompressing the entire array into memory.

        The :class:`~astropy.io.fits.CompImageSection` object returned by this
        attribute is not meant to be used directly by itself. Rather, slices of
        the section return the appropriate slice of the data, and loads *only*
        that section into memory. Any valid basic Numpy index can be used to
        slice :class:`~astropy.io.fits.CompImageSection`.

        Note that accessing data using :attr:`CompImageHDU.section` will always
        load tiles one at a time from disk, and therefore when accessing a large
        fraction of the data (or slicing it in a way that would cause most tiles
        to be loaded) you may obtain better performance by using
        :attr:`CompImageHDU.data`.
        """
        return CompImageSection(self)

    @property
    def tile_shape(self):
        """
        The tile shape used for the tiled compression.

        This shape is given in Numpy/C order
        """
        return tuple(
            [
                self._header[f"ZTILE{idx + 1}"]
                for idx in range(self._header["ZNAXIS"] - 1, -1, -1)
            ]
        )

    @property
    def compression_type(self):
        """
        The name of the compression algorithm.
        """
        return self._header.get("ZCMPTYPE", DEFAULT_COMPRESSION_TYPE)
