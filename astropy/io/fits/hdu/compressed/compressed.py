# Licensed under a 3-clause BSD style license - see LICENSE.rst

import ctypes
import math
import time
import warnings

import numpy as np

from astropy.io.fits.fitsrec import FITS_rec
from astropy.io.fits.hdu.base import BITPIX2DTYPE, DELAYED
from astropy.io.fits.hdu.compressed._quantization import DITHER_METHODS
from astropy.io.fits.hdu.compressed._tiled_compression import (
    _get_compression_setting,
    compress_image_data,
)
from astropy.io.fits.hdu.compressed.utils import _tile_shape, _validate_tile_shape
from astropy.io.fits.hdu.image import ImageHDU
from astropy.io.fits.header import Header
from astropy.io.fits.util import _is_int
from astropy.io.fits.verify import _ErrList
from astropy.utils.decorators import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning

from .header import (
    _bintable_header_to_image_header,
    _image_header_to_empty_bintable,
)
from .section import CompImageSection
from .settings import (
    CMTYPE_ALIASES,
    COMPRESSION_TYPES,
    DEFAULT_COMPRESSION_TYPE,
    DEFAULT_DITHER_SEED,
    DEFAULT_HCOMP_SCALE,
    DEFAULT_HCOMP_SMOOTH,
    DEFAULT_QUANTIZE_LEVEL,
    DEFAULT_QUANTIZE_METHOD,
    DITHER_SEED_CHECKSUM,
    DITHER_SEED_CLOCK,
)

__all__ = ["CompImageHDU"]


class CompImageHDU(ImageHDU):
    """
    Compressed Image HDU class.
    """

    _default_name = "COMPRESSED_IMAGE"

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
        uint=True,
        scale_back=None,
        bintable=None,
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

        self._bintable = None

        if data is DELAYED or bintable is not None:
            # NOTE: for now we don't ever read in CompImageHDU directly from
            # files, instead we read in BinTableHDU and pass it in here. In
            # future if we do want to read CompImageHDU in directly, we can
            # use the following code.
            # if data is DELAYED:
            #     # Reading the HDU from a file
            #     self._bintable = _CompBinTableHDU(data=data, header=header)
            # else:

            # If bintable is passed in, it should be a BinTableHDU
            self._bintable = bintable
            self._bintable._load_variable_length_data = False
            self._bintable._manages_own_heap = True
            self._bintable._new = False
            self._bitpix = self._bintable.header["ZBITPIX"]

            header = self._bintable_to_image_header()
            header._modified = False
            for card in header._cards:
                card._modified = False

            super().__init__(
                data=DELAYED,
                header=header,
                name=name,
                do_not_scale_image_data=do_not_scale_image_data,
                uint=uint,
                scale_back=scale_back,
            )

            self.compression_type = self._bintable.header.get(
                "ZCMPTYPE", DEFAULT_COMPRESSION_TYPE
            )
            self.tile_shape = tuple(_tile_shape(self._bintable.header))
            self.hcomp_scale = int(
                _get_compression_setting(bintable.header, "SCALE", DEFAULT_HCOMP_SCALE)
            )
            self.hcomp_smooth = _get_compression_setting(
                bintable.header, "SMOOTH", DEFAULT_HCOMP_SMOOTH
            )
            self.quantize_level = _get_compression_setting(
                bintable.header, "noisebit", DEFAULT_QUANTIZE_LEVEL
            )
            self.quantize_method = DITHER_METHODS[
                bintable.header.get("ZQUANTIZ", "NO_DITHER")
            ]
            self.dither_seed = bintable.header.get("ZDITHER0", DEFAULT_DITHER_SEED)

        else:
            # Create at least a skeleton HDU that matches the input
            # header and data (if any were input)
            super().__init__(
                data=data,
                header=header or Header(),
                name=name,
                do_not_scale_image_data=do_not_scale_image_data,
                uint=uint,
                scale_back=scale_back,
            )

            if header is not None and "SIMPLE" in header:
                self.header["SIMPLE"] = header["SIMPLE"]

            self.compression_type = compression_type
            self.tile_shape = _validate_tile_shape(
                tile_shape=tile_shape,
                compression_type=self.compression_type,
                image_header=self.header,
            )
            self.hcomp_scale = hcomp_scale
            self.hcomp_smooth = hcomp_smooth
            self.quantize_level = quantize_level
            self.quantize_method = quantize_method
            self.dither_seed = dither_seed

            # TODO: just for parameter validation, e.g. tile shape - we shouldn't
            # ideally need this and should instead validate the values as they are
            # set above.
            self._get_bintable_without_data()

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

        return True

    @property
    def compression_type(self):
        return self._compression_type

    @compression_type.setter
    def compression_type(self, value):
        value = CMTYPE_ALIASES.get(value, value)
        if value in COMPRESSION_TYPES:
            self._compression_type = value
        else:
            warnings.warn(
                "Unknown compression type provided (supported are {}). "
                "Default ({}) compression will be used.".format(
                    ", ".join(map(repr, COMPRESSION_TYPES)),
                    DEFAULT_COMPRESSION_TYPE,
                ),
                AstropyUserWarning,
            )
            self._compression_type = DEFAULT_COMPRESSION_TYPE

    def _get_bintable_without_data(self):
        """
        Convert the current ImageHDU (excluding the actual data) to a BinTableHDU
        with the correct header.
        """
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
        # TODO: test above

        # NOTE: for now the function below modifies the compressed binary table
        # bintable._header in-place, but this could be refactored in future to
        # return the compressed header.

        bintable = _image_header_to_empty_bintable(
            self.header,
            name=self.name,
            huge_hdu=huge_hdu,
            compression_type=self.compression_type,
            tile_shape=self.tile_shape,
            hcomp_scale=self.hcomp_scale,
            hcomp_smooth=self.hcomp_smooth,
            quantize_level=self.quantize_level,
            quantize_method=self.quantize_method,
            dither_seed=self.dither_seed,
            axes=self._axes,
            generate_dither_seed=self._generate_dither_seed,
        )

        return bintable

    @property
    def _data_loaded(self):
        """
        Whether the data is fully decompressed into self.data - note that is
        a little different to _data_loaded on other HDUs, but it is conceptually
        the same idea in a way.
        """
        return "data" in self.__dict__ and super().data is not None

    @property
    def _data_shape(self):
        if self._data_loaded:
            return self.data.shape
        else:
            return tuple(reversed(self._axes))

    @lazyproperty
    def data(self):
        """
        The decompressed data array.

        Note that accessing this will cause all the tiles to be loaded,
        decompressed, and combined into a single data array. If you do
        not need to access the whole array, consider instead using the
        :attr:`~astropy.io.fits.CompImageHDU.section` property.
        """
        # If there is no internal binary table, the HDU was not created from a
        # file and therefore the data is just the one on the parent ImageHDU
        # class

        if self._data_loaded:
            return super().data
        elif self._bintable is None or len(self._bintable.data) == 0:
            return None

        # Since .section has general code to load any arbitrary part of the
        # data, we can just use this
        data = self.section[...]

        return data

    @data.setter
    def data(self, data):
        ImageHDU.data.fset(self, data)
        if (
            data is not None
            and hasattr(self, "tile_shape")
            and len(self.tile_shape) != data.ndim
        ):
            self.tile_shape = _validate_tile_shape(
                tile_shape=[],
                compression_type=self.compression_type,
                image_header=self.header,
            )

    @property
    def compressed_data(self):
        return None if self._bintable is None else self._bintable.data

    def _bintable_to_image_header(self):
        if self._bintable is None:
            raise ValueError("bintable is not set")

        # Clean up any possible doubled EXTNAME keywords that use
        # the default. Do this on the original header to ensure
        # duplicates are removed cleanly.
        self._remove_unnecessary_default_extnames(self._bintable.header)

        # Convert compressed header to image header and save
        # it off to self._image_header so it can be referenced later
        # unambiguously
        return _bintable_header_to_image_header(self._bintable.header)

    def _add_data_to_bintable(self, bintable):
        """
        Compress the image data so that it may be written to a file.
        """
        if self.data is None:
            return

        heap = compress_image_data(
            self.data, self.compression_type, bintable.header, bintable.columns
        )

        dtype = bintable.columns.dtype.newbyteorder(">")
        buf = np.frombuffer(heap, dtype=np.uint8)
        data = (
            buf[: bintable._theap]
            .view(dtype=dtype, type=np.rec.recarray)
            .view(FITS_rec)
        )
        data._load_variable_length_data = False
        data._coldefs = bintable.columns
        data._heapoffset = bintable._theap
        data._heapsize = len(buf) - bintable._theap

        bintable.data = data

    def _prewriteto(self, inplace=False):
        if (
            self._bintable is not None
            and not self._has_data
            and not self.header._modified
        ):
            self._tmp_bintable = self._bintable
            self._tmp_bintable._output_checksum = self._output_checksum
            return self._tmp_bintable._prewriteto(inplace=inplace)

        if self._scale_back:
            self._scale_internal(
                BITPIX2DTYPE[self._orig_bitpix], blank=self._orig_blank
            )

        self._tmp_bintable = self._get_bintable_without_data()

        self._add_data_to_bintable(self._tmp_bintable)

        # If a bintable already exists internally we should update that instead
        # of using a whole new BinTableHDU so that mode='update' works.

        if self._bintable is not None:
            self._bintable.header = self._tmp_bintable.header
            self._bintable.data = self._tmp_bintable.data
            self._tmp_bintable = self._bintable

        self._tmp_bintable._output_checksum = self._output_checksum
        return self._tmp_bintable._prewriteto(inplace=inplace)

    def _writeto(self, fileobj, inplace=False, copy=False):
        if self._tmp_bintable is not None:
            # Each time we assign the bintable data to the BinTableHDU, some of
            # the blank keywords get removed, so at this point, just before
            # writing, we should make sure that the number of blank cards in
            # the final binary table to be written matches the number of blanks
            # in the image header.

            image_blanks = self.header._countblanks()
            table_blanks = self._tmp_bintable.header._countblanks()

            for _ in range(image_blanks - table_blanks):
                self._tmp_bintable.header.append()

            return self._tmp_bintable._writeto(fileobj, inplace=inplace, copy=copy)

    def _postwriteto(self):
        self._tmp_bintable = None

    def _close(self, closed=True):
        if self._bintable is not None:
            return self._bintable._close(closed=closed)

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

    def _verify(self, *args, **kwargs):
        # The following is the default _verify for ImageHDU
        errs = super()._verify(*args, **kwargs)

        # However in some cases the decompressed header is actually like a
        # PrimaryHDU header rather than an ImageHDU header, in which case
        # there are certain errors we can ignore
        if "SIMPLE" in self.header:
            errs_filtered = []
            for err in errs:
                if len(err) >= 2 and err[1] in (
                    "'XTENSION' card does not exist.",
                    "'PCOUNT' card does not exist.",
                    "'GCOUNT' card does not exist.",
                ):
                    continue
                errs_filtered.append(err)
            return _ErrList(errs_filtered)
        else:
            return errs

    def fileinfo(self):
        if self._bintable is not None:
            return self._bintable.fileinfo()

    @property
    def _data_offset(self):
        if self._bintable is not None:
            return self._bintable._data_offset

    @_data_offset.setter
    def _data_offset(self, value):
        # We should never set _data_offset to a non-None value. We need to
        # implement this setter as one of the parent classes sets _data_offset
        # to None in __init__.
        if value is not None:
            raise RuntimeError("Cannot set CompImageHDU._data_offset")

    @property
    def _header_offset(self):
        if self._bintable is not None:
            return self._bintable._header_offset

    @_header_offset.setter
    def _header_offset(self, value):
        # We should never set _data_offset to a non-None value. We need to
        # implement this setter as one of the parent classes sets _data_offset
        # to None in __init__.
        if value is not None:
            raise RuntimeError("Cannot set CompImageHDU._header_offset")

    @property
    def _data_size(self):
        if self._bintable is not None:
            return self._bintable._data_size

    @_data_size.setter
    def _data_size(self, value):
        # We should never set _data_offset to a non-None value. We need to
        # implement this setter as one of the parent classes sets _data_offset
        # to None in __init__.
        if value is not None:
            raise RuntimeError("Cannot set CompImageHDU._data_size")
