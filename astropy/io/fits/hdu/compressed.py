# Licensed under a 3-clause BSD style license - see PYFITS.rst

import ctypes
import gc
import math
import re
import time
import warnings

import numpy as np

from .base import DELAYED, ExtensionHDU, BITPIX2DTYPE, DTYPE2BITPIX
from .image import ImageHDU
from .table import BinTableHDU
from ..card import Card
from ..column import Column, ColDefs, TDEF_RE
from ..column import KEYWORD_NAMES as TABLE_KEYWORD_NAMES
from ..fitsrec import FITS_rec
from ..header import Header
from ..util import (_is_pseudo_unsigned, _unsigned_zero, _is_int,
                    _get_array_mmap)

from ....extern.six import string_types, iteritems
from ....extern.six.moves import range
from ....utils import lazyproperty
from ....utils.compat import suppress
from ....utils.exceptions import (AstropyPendingDeprecationWarning,
                                  AstropyUserWarning)

try:
    from .. import compression
    COMPRESSION_SUPPORTED = COMPRESSION_ENABLED = True
except ImportError:
    COMPRESSION_SUPPORTED = COMPRESSION_ENABLED = False


# Quantization dithering method constants; these are right out of fitsio.h
NO_DITHER = -1
SUBTRACTIVE_DITHER_1 = 1
SUBTRACTIVE_DITHER_2 = 2
QUANTIZE_METHOD_NAMES = {
    NO_DITHER: 'NO_DITHER',
    SUBTRACTIVE_DITHER_1: 'SUBTRACTIVE_DITHER_1',
    SUBTRACTIVE_DITHER_2: 'SUBTRACTIVE_DITHER_2'
}
DITHER_SEED_CLOCK = 0
DITHER_SEED_CHECKSUM = -1


# Default compression parameter values
DEFAULT_COMPRESSION_TYPE = 'RICE_1'
DEFAULT_QUANTIZE_LEVEL = 16.
DEFAULT_QUANTIZE_METHOD = NO_DITHER
DEFAULT_DITHER_SEED = DITHER_SEED_CLOCK
DEFAULT_HCOMP_SCALE = 0
DEFAULT_HCOMP_SMOOTH = 0
DEFAULT_BLOCK_SIZE = 32
DEFAULT_BYTE_PIX = 4

CMTYPE_ALIASES = {}

# CFITSIO version-specific features
if COMPRESSION_SUPPORTED:
    try:
        CFITSIO_SUPPORTS_GZIPDATA = compression.CFITSIO_VERSION >= 3.28
        CFITSIO_SUPPORTS_Q_FORMAT = compression.CFITSIO_VERSION >= 3.35
        if compression.CFITSIO_VERSION >= 3.35:
            CMTYPE_ALIASES['RICE_ONE'] = 'RICE_1'
    except AttributeError:
        # This generally shouldn't happen unless running setup.py in an
        # environment where an old build of pyfits exists
        CFITSIO_SUPPORTS_GZIPDATA = True
        CFITSIO_SUPPORTS_Q_FORMAT = True


COMPRESSION_KEYWORDS = set(['ZIMAGE', 'ZCMPTYPE', 'ZBITPIX', 'ZNAXIS',
                            'ZMASKCMP', 'ZSIMPLE', 'ZTENSION', 'ZEXTEND'])


class CompImageHeader(Header):
    """
    Header object for compressed image HDUs designed to keep the compression
    header and the underlying image header properly synchronized.

    This essentially wraps the image header, so that all values are read from
    and written to the image header.  However, updates to the image header will
    also update the table header where appropriate.
    """

    # TODO: The difficulty of implementing this screams a need to rewrite this
    # module

    _keyword_remaps = {
        'SIMPLE': 'ZSIMPLE', 'XTENSION': 'ZTENSION', 'BITPIX': 'ZBITPIX',
        'NAXIS': 'ZNAXIS', 'EXTEND': 'ZEXTEND', 'BLOCKED': 'ZBLOCKED',
        'PCOUNT': 'ZPCOUNT', 'GCOUNT': 'ZGCOUNT', 'CHECKSUM': 'ZHECKSUM',
        'DATASUM': 'ZDATASUM'
    }

    _zdef_re = re.compile(r'(?P<label>^[Zz][a-zA-Z]*)(?P<num>[1-9][0-9 ]*$)?')
    _compression_keywords = set(_keyword_remaps.values()).union(
        ['ZIMAGE', 'ZCMPTYPE', 'ZMASKCMP', 'ZQUANTIZ', 'ZDITHER0'])
    _indexed_compression_keywords = set(['ZNAXIS', 'ZTILE', 'ZNAME', 'ZVAL'])
    # TODO: Once it place it should be possible to manage some of this through
    # the schema system, but it's not quite ready for that yet.  Also it still
    # makes more sense to change CompImageHDU to subclass ImageHDU :/

    def __init__(self, table_header, image_header=None):
        if image_header is None:
            image_header = Header()
        self._cards = image_header._cards
        self._keyword_indices = image_header._keyword_indices
        self._rvkc_indices = image_header._rvkc_indices
        self._modified = image_header._modified
        self._table_header = table_header

    # We need to override and Header methods that can modify the header, and
    # ensure that they sync with the underlying _table_header

    def __setitem__(self, key, value):
        # This isn't pretty, but if the `key` is either an int or a tuple we
        # need to figure out what keyword name that maps to before doing
        # anything else; these checks will be repeated later in the
        # super().__setitem__ call but I don't see another way around it
        # without some major refactoring
        if self._set_slice(key, value, self):
            return

        if isinstance(key, int):
            keyword, index = self._keyword_from_index(key)
        elif isinstance(key, tuple):
            keyword, index = key
        else:
            # We don't want to specify and index otherwise, because that will
            # break the behavior for new keywords and for commentary keywords
            keyword, index = key, None

        if self._is_reserved_keyword(keyword):
            return

        super(CompImageHeader, self).__setitem__(key, value)

        if index is not None:
            remapped_keyword = self._remap_keyword(keyword)
            self._table_header[remapped_keyword, index] = value
        # Else this will pass through to ._update

    def __delitem__(self, key):
        if isinstance(key, slice) or self._haswildcard(key):
            # If given a slice pass that on to the superclass and bail out
            # early; we only want to make updates to _table_header when given
            # a key specifying a single keyword
            return super(CompImageHeader, self).__delitem__(key)

        if isinstance(key, int):
            keyword, index = self._keyword_from_index(key)
        elif isinstance(key, tuple):
            keyword, index = key
        else:
            keyword, index = key, None

        if key not in self:
            raise KeyError("Keyword {!r} not found.".format(key))

        super(CompImageHeader, self).__delitem__(key)

        remapped_keyword = self._remap_keyword(keyword)

        if remapped_keyword in self._table_header:
            if index is not None:
                del self._table_header[(remapped_keyword, index)]
            else:
                del self._table_header[remapped_keyword]

    def append(self, card=None, useblanks=True, bottom=False, end=False):
        # This logic unfortunately needs to be duplicated from the base class
        # in order to determine the keyword
        if isinstance(card, string_types):
            card = Card(card)
        elif isinstance(card, tuple):
            card = Card(*card)
        elif card is None:
            card = Card()
        elif not isinstance(card, Card):
            raise ValueError(
                'The value appended to a Header must be either a keyword or '
                '(keyword, value, [comment]) tuple; got: {!r}'.format(card))

        if self._is_reserved_keyword(card.keyword):
            return

        super(CompImageHeader, self).append(card=card, useblanks=useblanks,
                                            bottom=bottom, end=end)

        remapped_keyword = self._remap_keyword(card.keyword)
        card = Card(remapped_keyword, card.value, card.comment)

        # Here we disable the use of blank cards, because the call above to
        # Header.append may have already deleted a blank card in the table
        # header, thanks to inheritance: Header.append calls 'del self[-1]'
        # to delete a blank card, which calls CompImageHeader.__deltitem__,
        # which deletes the blank card both in the image and the table headers!
        self._table_header.append(card=card, useblanks=False,
                                  bottom=bottom, end=end)

    def insert(self, key, card, useblanks=True, after=False):
        if isinstance(key, int):
            # Determine condition to pass through to append
            if after:
                if key == -1:
                    key = len(self._cards)
                else:
                    key += 1

            if key >= len(self._cards):
                self.append(card, end=True)
                return

        if isinstance(card, string_types):
            card = Card(card)
        elif isinstance(card, tuple):
            card = Card(*card)
        elif not isinstance(card, Card):
            raise ValueError(
                'The value inserted into a Header must be either a keyword or '
                '(keyword, value, [comment]) tuple; got: {!r}'.format(card))

        if self._is_reserved_keyword(card.keyword):
            return

        # Now the tricky part is to determine where to insert in the table
        # header.  If given a numerical index we need to map that to the
        # corresponding index in the table header.  Although rare, there may be
        # cases where there is no mapping in which case we just try the same
        # index
        # NOTE: It is crucial that remapped_index in particular is figured out
        # before the image header is modified
        remapped_index = self._remap_index(key)
        remapped_keyword = self._remap_keyword(card.keyword)

        super(CompImageHeader, self).insert(key, card, useblanks=useblanks,
                                            after=after)

        card = Card(remapped_keyword, card.value, card.comment)

        # Here we disable the use of blank cards, because the call above to
        # Header.insert may have already deleted a blank card in the table
        # header, thanks to inheritance: Header.insert calls 'del self[-1]'
        # to delete a blank card, which calls CompImageHeader.__delitem__,
        # which deletes the blank card both in the image and the table headers!
        self._table_header.insert(remapped_index, card, useblanks=False,
                                  after=after)

    def _update(self, card):
        keyword = card[0]

        if self._is_reserved_keyword(keyword):
            return

        super(CompImageHeader, self)._update(card)

        if keyword in Card._commentary_keywords:
            # Otherwise this will result in a duplicate insertion
            return

        remapped_keyword = self._remap_keyword(keyword)
        self._table_header._update((remapped_keyword,) + card[1:])

    # Last piece needed (I think) for synchronizing with the real header
    # This one is tricky since _relativeinsert calls insert
    def _relativeinsert(self, card, before=None, after=None, replace=False):
        keyword = card[0]

        if self._is_reserved_keyword(keyword):
            return

        # Now we have to figure out how to remap 'before' and 'after'
        if before is None:
            if isinstance(after, int):
                remapped_after = self._remap_index(after)
            else:
                remapped_after = self._remap_keyword(after)
            remapped_before = None
        else:
            if isinstance(before, int):
                remapped_before = self._remap_index(before)
            else:
                remapped_before = self._remap_keyword(before)
            remapped_after = None

        super(CompImageHeader, self)._relativeinsert(card, before=before,
                                                     after=after,
                                                     replace=replace)

        remapped_keyword = self._remap_keyword(keyword)

        card = Card(remapped_keyword, card[1], card[2])
        self._table_header._relativeinsert(card, before=remapped_before,
                                           after=remapped_after,
                                           replace=replace)

    @classmethod
    def _is_reserved_keyword(cls, keyword, warn=True):
        msg = ('Keyword {!r} is reserved for use by the FITS Tiled Image '
               'Convention and will not be stored in the header for the '
               'image being compressed.'.format(keyword))

        if keyword == 'TFIELDS':
            if warn:
                warnings.warn(msg)
            return True

        m = TDEF_RE.match(keyword)

        if m and m.group('label').upper() in TABLE_KEYWORD_NAMES:
            if warn:
                warnings.warn(msg)
            return True

        m = cls._zdef_re.match(keyword)

        if m:
            label = m.group('label').upper()
            num = m.group('num')
            if num is not None and label in cls._indexed_compression_keywords:
                if warn:
                    warnings.warn(msg)
                return True
            elif label in cls._compression_keywords:
                if warn:
                    warnings.warn(msg)
                return True

        return False

    @classmethod
    def _remap_keyword(cls, keyword):
        # Given a keyword that one might set on an image, remap that keyword to
        # the name used for it in the COMPRESSED HDU header
        # This is mostly just a lookup in _keyword_remaps, but needs handling
        # for NAXISn keywords

        is_naxisn = False
        if keyword[:5] == 'NAXIS':
            with suppress(ValueError):
                index = int(keyword[5:])
                is_naxisn = index > 0

        if is_naxisn:
            return 'ZNAXIS{}'.format(index)

        # If the keyword does not need to be remapped then just return the
        # original keyword
        return cls._keyword_remaps.get(keyword, keyword)

    def _remap_index(self, idx):
        # Given an integer index into this header, map that to the index in the
        # table header for the same card.  If the card doesn't exist in the
        # table header (generally should *not* be the case) this will just
        # return the same index
        # This *does* also accept a keyword or (keyword, repeat) tuple and
        # obtains the associated numerical index with self._cardindex
        if not isinstance(idx, int):
            idx = self._cardindex(idx)

        keyword, repeat = self._keyword_from_index(idx)
        remapped_insert_keyword = self._remap_keyword(keyword)

        with suppress(IndexError, KeyError):
            idx = self._table_header._cardindex((remapped_insert_keyword,
                                                 repeat))

        return idx


# TODO: Fix this class so that it doesn't actually inherit from BinTableHDU,
# but instead has an internal BinTableHDU reference
class CompImageHDU(BinTableHDU):
    """
    Compressed Image HDU class.
    """

    # Maps deprecated keyword arguments to __init__ to their new names
    DEPRECATED_KWARGS = {
        'compressionType': 'compression_type', 'tileSize': 'tile_size',
        'hcompScale': 'hcomp_scale', 'hcompSmooth': 'hcomp_smooth',
        'quantizeLevel': 'quantize_level'
    }

    _manages_own_heap = True
    """
    The calls to CFITSIO lay out the heap data in memory, and we write it out
    the same way CFITSIO organizes it.  In principle this would break if a user
    manually changes the underlying compressed data by hand, but there is no
    reason they would want to do that (and if they do that's their
    responsibility).
    """

    def __init__(self, data=None, header=None, name=None,
                 compression_type=DEFAULT_COMPRESSION_TYPE,
                 tile_size=None,
                 hcomp_scale=DEFAULT_HCOMP_SCALE,
                 hcomp_smooth=DEFAULT_HCOMP_SMOOTH,
                 quantize_level=DEFAULT_QUANTIZE_LEVEL,
                 quantize_method=DEFAULT_QUANTIZE_METHOD,
                 dither_seed=DEFAULT_DITHER_SEED,
                 do_not_scale_image_data=False,
                 uint=False, scale_back=False, **kwargs):
        """
        Parameters
        ----------
        data : array, optional
            Uncompressed image data

        header : Header instance, optional
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
            ``'GZIP_2'``, ``'HCOMPRESS_1'``

        tile_size : int, optional
            Compression tile sizes.  Default treats each row of image as a
            tile.

        hcomp_scale : float, optional
            HCOMPRESS scale parameter

        hcomp_smooth : float, optional
            HCOMPRESS smooth parameter

        quantize_level : float, optional
            Floating point quantization level; see note below

        quantize_method : int, optional
            Floating point quantization dithering method; can be either
            ``NO_DITHER`` (-1), ``SUBTRACTIVE_DITHER_1`` (1; default), or
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
        smaller than the other tiles.  The ``tile_size`` parameter may be
        provided as a list of tile sizes, one for each dimension in the image.
        For example a ``tile_size`` value of ``[100,100]`` would divide a 300 X
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

        if not COMPRESSION_SUPPORTED:
            # TODO: Raise a more specific Exception type
            raise Exception('The astropy.io.fits.compression module is not '
                            'available.  Creation of compressed image HDUs is '
                            'disabled.')

        compression_type = CMTYPE_ALIASES.get(compression_type, compression_type)

        # Handle deprecated keyword arguments
        compression_opts = {}
        for oldarg, newarg in self.DEPRECATED_KWARGS.items():
            if oldarg in kwargs:
                warnings.warn('Keyword argument {} to {} is pending '
                              'deprecation; use {} instead'.format(
                        oldarg, self.__class__.__name__, newarg),
                              AstropyPendingDeprecationWarning)
                compression_opts[newarg] = kwargs[oldarg]
                del kwargs[oldarg]
            else:
                compression_opts[newarg] = locals()[newarg]
        # Include newer compression options that don't required backwards
        # compatibility with deprecated spellings
        compression_opts['quantize_method'] = quantize_method
        compression_opts['dither_seed'] = dither_seed

        if data is DELAYED:
            # Reading the HDU from a file
            super(CompImageHDU, self).__init__(data=data, header=header)
        else:
            # Create at least a skeleton HDU that matches the input
            # header and data (if any were input)
            super(CompImageHDU, self).__init__(data=None, header=header)

            # Store the input image data
            self.data = data

            # Update the table header (_header) to the compressed
            # image format and to match the input data (if any);
            # Create the image header (_image_header) from the input
            # image header (if any) and ensure it matches the input
            # data; Create the initially empty table data array to
            # hold the compressed data.
            self._update_header_data(header, name, **compression_opts)

        # TODO: A lot of this should be passed on to an internal image HDU o
        # something like that, see ticket #88
        self._do_not_scale_image_data = do_not_scale_image_data
        self._uint = uint
        self._scale_back = scale_back

        self._axes = [self._header.get('ZNAXIS' + str(axis + 1), 0)
                      for axis in range(self._header.get('ZNAXIS', 0))]

        # store any scale factors from the table header
        if do_not_scale_image_data:
            self._bzero = 0
            self._bscale = 1
        else:
            self._bzero = self._header.get('BZERO', 0)
            self._bscale = self._header.get('BSCALE', 1)
        self._bitpix = self._header['ZBITPIX']

        self._orig_bzero = self._bzero
        self._orig_bscale = self._bscale
        self._orig_bitpix = self._bitpix

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        if card.keyword != 'XTENSION':
            return False

        xtension = card.value
        if isinstance(xtension, string_types):
            xtension = xtension.rstrip()

        if xtension not in ('BINTABLE', 'A3DTABLE'):
            return False

        if 'ZIMAGE' not in header or not header['ZIMAGE']:
            return False

        if COMPRESSION_SUPPORTED and COMPRESSION_ENABLED:
            return True
        elif not COMPRESSION_SUPPORTED:
            warnings.warn('Failure matching header to a compressed image '
                          'HDU: The compression module is not available.\n'
                          'The HDU will be treated as a Binary Table HDU.',
                          AstropyUserWarning)
            return False
        else:
            # Compression is supported but disabled; just pass silently (#92)
            return False

    def _update_header_data(self, image_header,
                            name=None,
                            compression_type=None,
                            tile_size=None,
                            hcomp_scale=None,
                            hcomp_smooth=None,
                            quantize_level=None,
                            quantize_method=None,
                            dither_seed=None):
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
        image_header : Header instance
            header to be associated with the image

        name : str, optional
            the ``EXTNAME`` value; if this value is `None`, then the name from
            the input image header will be used; if there is no name in the
            input image header then the default name 'COMPRESSED_IMAGE' is used

        compression_type : str, optional
            compression algorithm 'RICE_1', 'PLIO_1', 'GZIP_1', 'GZIP_2',
            'HCOMPRESS_1'; if this value is `None`, use value already in the
            header; if no value already in the header, use 'RICE_1'

        tile_size : sequence of int, optional
            compression tile sizes as a list; if this value is `None`, use
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
            huge_hdu = self.data.nbytes > 2 ** 32

            if huge_hdu and not CFITSIO_SUPPORTS_Q_FORMAT:
                raise IOError(
                    "Astropy cannot compress images greater than 4 GB in size "
                    "({} is {} bytes) without CFITSIO >= 3.35".format(
                        (self.name, self.ver), self.data.nbytes))
        else:
            huge_hdu = False

        # Update the extension name in the table header
        if not name and 'EXTNAME' not in self._header:
            name = 'COMPRESSED_IMAGE'

        if name:
            self._header.set('EXTNAME', name,
                             'name of this binary table extension',
                             after='TFIELDS')
            self.name = name
        else:
            self.name = self._header['EXTNAME']

        # Set the compression type in the table header.
        if compression_type:
            if compression_type not in ['RICE_1', 'GZIP_1', 'GZIP_2', 'PLIO_1',
                                        'HCOMPRESS_1']:
                warnings.warn('Unknown compression type provided.  Default '
                              '({}) compression used.'.format(
                        DEFAULT_COMPRESSION_TYPE), AstropyUserWarning)
                compression_type = DEFAULT_COMPRESSION_TYPE

            self._header.set('ZCMPTYPE', compression_type,
                             'compression algorithm', after='TFIELDS')
        else:
            compression_type = self._header.get('ZCMPTYPE',
                                                DEFAULT_COMPRESSION_TYPE)
            compression_type = CMTYPE_ALIASES.get(compression_type,
                                                  compression_type)

        # If the input image header had BSCALE/BZERO cards, then insert
        # them in the table header.

        if image_header:
            bzero = image_header.get('BZERO', 0.0)
            bscale = image_header.get('BSCALE', 1.0)
            after_keyword = 'EXTNAME'

            if bscale != 1.0:
                self._header.set('BSCALE', bscale, after=after_keyword)
                after_keyword = 'BSCALE'

            if bzero != 0.0:
                self._header.set('BZERO', bzero, after=after_keyword)

            bitpix_comment = image_header.comments['BITPIX']
            naxis_comment = image_header.comments['NAXIS']
        else:
            bitpix_comment = 'data type of original image'
            naxis_comment = 'dimension of original image'

        # Set the label for the first column in the table

        self._header.set('TTYPE1', 'COMPRESSED_DATA', 'label for field 1',
                         after='TFIELDS')

        # Set the data format for the first column.  It is dependent
        # on the requested compression type.

        if compression_type == 'PLIO_1':
            tform1 = '1QI' if huge_hdu else '1PI'
        else:
            tform1 = '1QB' if huge_hdu else '1PB'

        self._header.set('TFORM1', tform1,
                         'data format of field: variable length array',
                         after='TTYPE1')

        # Create the first column for the table.  This column holds the
        # compressed data.
        col1 = Column(name=self._header['TTYPE1'], format=tform1)

        # Create the additional columns required for floating point
        # data and calculate the width of the output table.

        zbitpix = self._image_header['BITPIX']

        if zbitpix < 0 and quantize_level != 0.0:
            # floating point image has 'COMPRESSED_DATA',
            # 'UNCOMPRESSED_DATA', 'ZSCALE', and 'ZZERO' columns (unless using
            # lossless compression, per CFITSIO)
            ncols = 4

            # CFITSIO 3.28 and up automatically use the GZIP_COMPRESSED_DATA
            # store floating point data that couldn't be quantized, instead
            # of the UNCOMPRESSED_DATA column.  There's no way to control
            # this behavior so the only way to determine which behavior will
            # be employed is via the CFITSIO version

            if CFITSIO_SUPPORTS_GZIPDATA:
                ttype2 = 'GZIP_COMPRESSED_DATA'
                # The required format for the GZIP_COMPRESSED_DATA is actually
                # missing from the standard docs, but CFITSIO suggests it
                # should be 1PB, which is logical.
                tform2 = '1QB' if huge_hdu else '1PB'
            else:
                # Q format is not supported for UNCOMPRESSED_DATA columns.
                ttype2 = 'UNCOMPRESSED_DATA'
                if zbitpix == 8:
                    tform2 = '1QB' if huge_hdu else '1PB'
                elif zbitpix == 16:
                    tform2 = '1QI' if huge_hdu else '1PI'
                elif zbitpix == 32:
                    tform2 = '1QJ' if huge_hdu else '1PJ'
                elif zbitpix == -32:
                    tform2 = '1QE' if huge_hdu else '1PE'
                else:
                    tform2 = '1QD' if huge_hdu else '1PD'

            # Set up the second column for the table that will hold any
            # uncompressable data.
            self._header.set('TTYPE2', ttype2, 'label for field 2',
                             after='TFORM1')

            self._header.set('TFORM2', tform2,
                             'data format of field: variable length array',
                             after='TTYPE2')

            col2 = Column(name=ttype2, format=tform2)

            # Set up the third column for the table that will hold
            # the scale values for quantized data.
            self._header.set('TTYPE3', 'ZSCALE', 'label for field 3',
                             after='TFORM2')
            self._header.set('TFORM3', '1D',
                             'data format of field: 8-byte DOUBLE',
                             after='TTYPE3')
            col3 = Column(name=self._header['TTYPE3'],
                          format=self._header['TFORM3'])

            # Set up the fourth column for the table that will hold
            # the zero values for the quantized data.
            self._header.set('TTYPE4', 'ZZERO', 'label for field 4',
                             after='TFORM3')
            self._header.set('TFORM4', '1D',
                             'data format of field: 8-byte DOUBLE',
                             after='TTYPE4')
            after = 'TFORM4'
            col4 = Column(name=self._header['TTYPE4'],
                          format=self._header['TFORM4'])

            # Create the ColDefs object for the table
            cols = ColDefs([col1, col2, col3, col4])
        else:
            # default table has just one 'COMPRESSED_DATA' column
            ncols = 1
            after = 'TFORM1'

            # remove any header cards for the additional columns that
            # may be left over from the previous data
            to_remove = ['TTYPE2', 'TFORM2', 'TTYPE3', 'TFORM3', 'TTYPE4',
                         'TFORM4']

            for k in to_remove:
                try:
                    del self._header[k]
                except KeyError:
                    pass

            # Create the ColDefs object for the table
            cols = ColDefs([col1])

        # Update the table header with the width of the table, the
        # number of fields in the table, the indicator for a compressed
        # image HDU, the data type of the image data and the number of
        # dimensions in the image data array.
        self._header.set('NAXIS1', cols.dtype.itemsize,
                         'width of table in bytes')
        self._header.set('TFIELDS', ncols, 'number of fields in each row',
                         after='GCOUNT')
        self._header.set('ZIMAGE', True, 'extension contains compressed image',
                         after=after)
        self._header.set('ZBITPIX', zbitpix,
                         bitpix_comment, after='ZIMAGE')
        self._header.set('ZNAXIS', self._image_header['NAXIS'], naxis_comment,
                         after='ZBITPIX')

        # Strip the table header of all the ZNAZISn and ZTILEn keywords
        # that may be left over from the previous data

        idx = 1
        while True:
            try:
                del self._header['ZNAXIS' + str(idx)]
                del self._header['ZTILE' + str(idx)]
                idx += 1
            except KeyError:
                break

        # Verify that any input tile size parameter is the appropriate
        # size to match the HDU's data.

        naxis = self._image_header['NAXIS']

        if not tile_size:
            tile_size = []
        elif len(tile_size) != naxis:
            warnings.warn('Provided tile size not appropriate for the data.  '
                          'Default tile size will be used.', AstropyUserWarning)
            tile_size = []

        # Set default tile dimensions for HCOMPRESS_1

        if compression_type == 'HCOMPRESS_1':
            if (self._image_header['NAXIS1'] < 4 or
                    self._image_header['NAXIS2'] < 4):
                raise ValueError('Hcompress minimum image dimension is '
                                 '4 pixels')
            elif tile_size:
                if tile_size[0] < 4 or tile_size[1] < 4:
                    # user specified tile size is too small
                    raise ValueError('Hcompress minimum tile dimension is '
                                     '4 pixels')
                major_dims = len([ts for ts in tile_size if ts > 1])
                if major_dims > 2:
                    raise ValueError(
                        'HCOMPRESS can only support 2-dimensional tile sizes.'
                        'All but two of the tile_size dimensions must be set '
                        'to 1.')

            if tile_size and (tile_size[0] == 0 and tile_size[1] == 0):
                # compress the whole image as a single tile
                tile_size[0] = self._image_header['NAXIS1']
                tile_size[1] = self._image_header['NAXIS2']

                for i in range(2, naxis):
                    # set all higher tile dimensions = 1
                    tile_size[i] = 1
            elif not tile_size:
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
                tile_size.append(self._image_header['NAXIS1'])

                if self._image_header['NAXIS2'] <= 30:
                    tile_size.append(self._image_header['NAXIS1'])
                else:
                    # look for another good tile dimension
                    naxis2 = self._image_header['NAXIS2']
                    for dim in [16, 24, 20, 30, 28, 26, 22, 18, 14]:
                        if naxis2 % dim == 0 or naxis2 % dim > 3:
                            tile_size.append(dim)
                            break
                    else:
                        tile_size.append(17)

                for i in range(2, naxis):
                    # set all higher tile dimensions = 1
                    tile_size.append(1)

            # check if requested tile size causes the last tile to have
            # less than 4 pixels

            remain = self._image_header['NAXIS1'] % tile_size[0]  # 1st dimen

            if remain > 0 and remain < 4:
                tile_size[0] += 1  # try increasing tile size by 1

                remain = self._image_header['NAXIS1'] % tile_size[0]

                if remain > 0 and remain < 4:
                    raise ValueError('Last tile along 1st dimension has '
                                     'less than 4 pixels')

            remain = self._image_header['NAXIS2'] % tile_size[1]  # 2nd dimen

            if remain > 0 and remain < 4:
                tile_size[1] += 1  # try increasing tile size by 1

                remain = self._image_header['NAXIS2'] % tile_size[1]

                if remain > 0 and remain < 4:
                    raise ValueError('Last tile along 2nd dimension has '
                                     'less than 4 pixels')

        # Set up locations for writing the next cards in the header.
        last_znaxis = 'ZNAXIS'

        if self._image_header['NAXIS'] > 0:
            after1 = 'ZNAXIS1'
        else:
            after1 = 'ZNAXIS'

        # Calculate the number of rows in the output table and
        # write the ZNAXISn and ZTILEn cards to the table header.
        nrows = 0

        for idx, axis in enumerate(self._axes):
            naxis = 'NAXIS' + str(idx + 1)
            znaxis = 'ZNAXIS' + str(idx + 1)
            ztile = 'ZTILE' + str(idx + 1)

            if tile_size and len(tile_size) >= idx + 1:
                ts = tile_size[idx]
            else:
                if ztile not in self._header:
                    # Default tile size
                    if not idx:
                        ts = self._image_header['NAXIS1']
                    else:
                        ts = 1
                else:
                    ts = self._header[ztile]
                tile_size.append(ts)

            if not nrows:
                nrows = (axis - 1) // ts + 1
            else:
                nrows *= ((axis - 1) // ts + 1)

            if image_header and naxis in image_header:
                self._header.set(znaxis, axis, image_header.comments[naxis],
                                 after=last_znaxis)
            else:
                self._header.set(znaxis, axis,
                                 'length of original image axis',
                                 after=last_znaxis)

            self._header.set(ztile, ts, 'size of tiles to be compressed',
                             after=after1)
            last_znaxis = znaxis
            after1 = ztile

        # Set the NAXIS2 header card in the table hdu to the number of
        # rows in the table.
        self._header.set('NAXIS2', nrows, 'number of rows in table')

        self.columns = cols

        # Set the compression parameters in the table header.

        # First, setup the values to be used for the compression parameters
        # in case none were passed in.  This will be either the value
        # already in the table header for that parameter or the default
        # value.
        idx = 1

        while True:
            zname = 'ZNAME' + str(idx)
            if zname not in self._header:
                break
            zval = 'ZVAL' + str(idx)
            if self._header[zname] == 'NOISEBIT':
                if quantize_level is None:
                    quantize_level = self._header[zval]
            if self._header[zname] == 'SCALE   ':
                if hcomp_scale is None:
                    hcomp_scale = self._header[zval]
            if self._header[zname] == 'SMOOTH  ':
                if hcomp_smooth is None:
                    hcomp_smooth = self._header[zval]
            idx += 1

        if quantize_level is None:
            quantize_level = DEFAULT_QUANTIZE_LEVEL

        if hcomp_scale is None:
            hcomp_scale = DEFAULT_HCOMP_SCALE

        if hcomp_smooth is None:
            hcomp_smooth = DEFAULT_HCOMP_SCALE

        # Next, strip the table header of all the ZNAMEn and ZVALn keywords
        # that may be left over from the previous data

        idx = 1

        while True:
            zname = 'ZNAME' + str(idx)
            if zname not in self._header:
                break
            zval = 'ZVAL' + str(idx)
            del self._header[zname]
            del self._header[zval]
            idx += 1

        # Finally, put the appropriate keywords back based on the
        # compression type.

        after_keyword = 'ZCMPTYPE'
        idx = 1

        if compression_type == 'RICE_1':
            self._header.set('ZNAME1', 'BLOCKSIZE', 'compression block size',
                             after=after_keyword)
            self._header.set('ZVAL1', DEFAULT_BLOCK_SIZE, 'pixels per block',
                             after='ZNAME1')

            self._header.set('ZNAME2', 'BYTEPIX',
                             'bytes per pixel (1, 2, 4, or 8)', after='ZVAL1')

            if self._header['ZBITPIX'] == 8:
                bytepix = 1
            elif self._header['ZBITPIX'] == 16:
                bytepix = 2
            else:
                bytepix = DEFAULT_BYTE_PIX

            self._header.set('ZVAL2', bytepix,
                             'bytes per pixel (1, 2, 4, or 8)',
                             after='ZNAME2')
            after_keyword = 'ZVAL2'
            idx = 3
        elif compression_type == 'HCOMPRESS_1':
            self._header.set('ZNAME1', 'SCALE', 'HCOMPRESS scale factor',
                             after=after_keyword)
            self._header.set('ZVAL1', hcomp_scale, 'HCOMPRESS scale factor',
                             after='ZNAME1')
            self._header.set('ZNAME2', 'SMOOTH', 'HCOMPRESS smooth option',
                             after='ZVAL1')
            self._header.set('ZVAL2', hcomp_smooth, 'HCOMPRESS smooth option',
                             after='ZNAME2')
            after_keyword = 'ZVAL2'
            idx = 3

        if self._image_header['BITPIX'] < 0:   # floating point image
            self._header.set('ZNAME' + str(idx), 'NOISEBIT',
                             'floating point quantization level',
                             after=after_keyword)
            self._header.set('ZVAL' + str(idx), quantize_level,
                             'floating point quantization level',
                             after='ZNAME' + str(idx))

            # Add the dither method and seed
            if quantize_method:
                if quantize_method not in [NO_DITHER, SUBTRACTIVE_DITHER_1,
                                           SUBTRACTIVE_DITHER_2]:
                    name = QUANTIZE_METHOD_NAMES[DEFAULT_QUANTIZE_METHOD]
                    warnings.warn('Unknown quantization method provided.  '
                                  'Default method ({}) used.'.format(name))
                    quantize_method = DEFAULT_QUANTIZE_METHOD

                if quantize_method == NO_DITHER:
                    zquantiz_comment = 'No dithering during quantization'
                else:
                    zquantiz_comment = 'Pixel Quantization Algorithm'

                self._header.set('ZQUANTIZ',
                                 QUANTIZE_METHOD_NAMES[quantize_method],
                                 zquantiz_comment,
                                 after='ZVAL' + str(idx))
            else:
                # If the ZQUANTIZ keyword is missing the default is to assume
                # no dithering, rather than whatever DEFAULT_QUANTIZE_METHOD
                # is set to
                quantize_method = self._header.get('ZQUANTIZ', NO_DITHER)

                if isinstance(quantize_method, string_types):
                    for k, v in iteritems(QUANTIZE_METHOD_NAMES):
                        if v.upper() == quantize_method:
                            quantize_method = k
                            break
                    else:
                        quantize_method = NO_DITHER

            if quantize_method == NO_DITHER:
                if 'ZDITHER0' in self._header:
                    # If dithering isn't being used then there's no reason to
                    # keep the ZDITHER0 keyword
                    del self._header['ZDITHER0']
            else:
                if dither_seed:
                    dither_seed = self._generate_dither_seed(dither_seed)
                elif 'ZDITHER0' in self._header:
                    dither_seed = self._header['ZDITHER0']
                else:
                    dither_seed = self._generate_dither_seed(
                            DEFAULT_DITHER_SEED)

                self._header.set('ZDITHER0', dither_seed,
                                 'dithering offset when quantizing floats',
                                 after='ZQUANTIZ')

        if image_header:
            # Move SIMPLE card from the image header to the
            # table header as ZSIMPLE card.

            if 'SIMPLE' in image_header:
                self._header.set('ZSIMPLE', image_header['SIMPLE'],
                                 image_header.comments['SIMPLE'],
                                 before='ZBITPIX')

            # Move EXTEND card from the image header to the
            # table header as ZEXTEND card.

            if 'EXTEND' in image_header:
                self._header.set('ZEXTEND', image_header['EXTEND'],
                                 image_header.comments['EXTEND'])

            # Move BLOCKED card from the image header to the
            # table header as ZBLOCKED card.

            if 'BLOCKED' in image_header:
                self._header.set('ZBLOCKED', image_header['BLOCKED'],
                                 image_header.comments['BLOCKED'])

            # Move XTENSION card from the image header to the
            # table header as ZTENSION card.

            # Since we only handle compressed IMAGEs, ZTENSION should
            # always be IMAGE, even if the caller has passed in a header
            # for some other type of extension.
            if 'XTENSION' in image_header:
                self._header.set('ZTENSION', 'IMAGE',
                                 image_header.comments['XTENSION'],
                                 before='ZBITPIX')

            # Move PCOUNT and GCOUNT cards from image header to the table
            # header as ZPCOUNT and ZGCOUNT cards.

            if 'PCOUNT' in image_header:
                self._header.set('ZPCOUNT', image_header['PCOUNT'],
                                 image_header.comments['PCOUNT'],
                                 after=last_znaxis)

            if 'GCOUNT' in image_header:
                self._header.set('ZGCOUNT', image_header['GCOUNT'],
                                 image_header.comments['GCOUNT'],
                                 after='ZPCOUNT')

            # Move CHECKSUM and DATASUM cards from the image header to the
            # table header as XHECKSUM and XDATASUM cards.

            if 'CHECKSUM' in image_header:
                self._header.set('ZHECKSUM', image_header['CHECKSUM'],
                                 image_header.comments['CHECKSUM'])

            if 'DATASUM' in image_header:
                self._header.set('ZDATASUM', image_header['DATASUM'],
                                 image_header.comments['DATASUM'])
        else:
            # Move XTENSION card from the image header to the
            # table header as ZTENSION card.

            # Since we only handle compressed IMAGEs, ZTENSION should
            # always be IMAGE, even if the caller has passed in a header
            # for some other type of extension.
            if 'XTENSION' in self._image_header:
                self._header.set('ZTENSION', 'IMAGE',
                                 self._image_header.comments['XTENSION'],
                                 before='ZBITPIX')

            # Move PCOUNT and GCOUNT cards from image header to the table
            # header as ZPCOUNT and ZGCOUNT cards.

            if 'PCOUNT' in self._image_header:
                self._header.set('ZPCOUNT', self._image_header['PCOUNT'],
                                 self._image_header.comments['PCOUNT'],
                                 after=last_znaxis)

            if 'GCOUNT' in self._image_header:
                self._header.set('ZGCOUNT', self._image_header['GCOUNT'],
                                 self._image_header.comments['GCOUNT'],
                                 after='ZPCOUNT')

        # When we have an image checksum we need to ensure that the same
        # number of blank cards exist in the table header as there were in
        # the image header.  This allows those blank cards to be carried
        # over to the image header when the hdu is uncompressed.

        if 'ZHECKSUM' in self._header:
            required_blanks = image_header._countblanks()
            image_blanks = self._image_header._countblanks()
            table_blanks = self._header._countblanks()

            for _ in range(required_blanks - image_blanks):
                self._image_header.append()
                table_blanks += 1

            for _ in range(required_blanks - table_blanks):
                self._header.append()

    @lazyproperty
    def data(self):
        # The data attribute is the image data (not the table data).
        data = compression.decompress_hdu(self)

        if data is None:
            return data

        # Scale the data if necessary
        if (self._orig_bzero != 0 or self._orig_bscale != 1):
            new_dtype = self._dtype_for_bitpix()
            data = np.array(data, dtype=new_dtype)

            zblank = None

            if 'ZBLANK' in self.compressed_data.columns.names:
                zblank = self.compressed_data['ZBLANK']
            else:
                if 'ZBLANK' in self._header:
                    zblank = np.array(self._header['ZBLANK'], dtype='int32')
                elif 'BLANK' in self._header:
                    zblank = np.array(self._header['BLANK'], dtype='int32')

            if zblank is not None:
                blanks = (data == zblank)

            if self._bscale != 1:
                np.multiply(data, self._bscale, data)
            if self._bzero != 0:
                # We have to explcitly cast self._bzero to prevent numpy from
                # raising an error when doing self.data += self._bzero, and we
                # do this instead of self.data = self.data + self._bzero to
                # avoid doubling memory usage.
                np.add(data, self._bzero, out=data, casting='unsafe')

            if zblank is not None:
                data = np.where(blanks, np.nan, data)

        # Right out of _ImageBaseHDU.data
        self._update_header_scale_info(data.dtype)

        return data

    @data.setter
    def data(self, data):
        if (data is not None) and (not isinstance(data, np.ndarray) or
                data.dtype.fields is not None):
            raise TypeError('CompImageHDU data has incorrect type:{}; '
                            'dtype.fields = {}'.format(
                    type(data), data.dtype.fields))

    @lazyproperty
    def compressed_data(self):
        # First we will get the table data (the compressed
        # data) from the file, if there is any.
        compressed_data = super(BinTableHDU, self).data
        if isinstance(compressed_data, np.rec.recarray):
            # Make sure not to use 'del self.data' so we don't accidentally
            # go through the self.data.fdel and close the mmap underlying
            # the compressed_data array
            del self.__dict__['data']
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
        if 'compressed_data' in self.__dict__:
            del self.__dict__['compressed_data']._coldefs

            # Now go ahead and delete from self.__dict__; normally
            # lazyproperty.__delete__ does this for us, but we can prempt it to
            # do some additional cleanup
            del self.__dict__['compressed_data']

            # If this file was mmap'd, numpy.memmap will hold open a file
            # handle until the underlying mmap object is garbage-collected;
            # since this reference leak can sometimes hang around longer than
            # welcome go ahead and force a garbage collection
            gc.collect()

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
        if hasattr(self, '_image_header'):
            return self._image_header

        # Start with a copy of the table header.
        image_header = self._header.copy()

        # Delete cards that are related to the table.  And move
        # the values of those cards that relate to the image from
        # their corresponding table cards.  These include
        # ZBITPIX -> BITPIX, ZNAXIS -> NAXIS, and ZNAXISn -> NAXISn.
        # (Note: Used set here instead of list in case there are any duplicate
        # keywords, which there may be in some pathological cases:
        # https://github.com/astropy/astropy/issues/2750
        for keyword in set(image_header):
            if CompImageHeader._is_reserved_keyword(keyword, warn=False):
                del image_header[keyword]

        if 'ZSIMPLE' in self._header:
            image_header.set('SIMPLE', self._header['ZSIMPLE'],
                             self._header.comments['ZSIMPLE'], before=0)
        elif 'ZTENSION' in self._header:
            if self._header['ZTENSION'] != 'IMAGE':
                warnings.warn("ZTENSION keyword in compressed "
                              "extension != 'IMAGE'", AstropyUserWarning)
            image_header.set('XTENSION', 'IMAGE',
                             self._header.comments['ZTENSION'], before=0)
        else:
            image_header.set('XTENSION', 'IMAGE', before=0)

        image_header.set('BITPIX', self._header['ZBITPIX'],
                         self._header.comments['ZBITPIX'], before=1)

        image_header.set('NAXIS', self._header['ZNAXIS'],
                         self._header.comments['ZNAXIS'], before=2)

        last_naxis = 'NAXIS'
        for idx in range(image_header['NAXIS']):
            znaxis = 'ZNAXIS' + str(idx + 1)
            naxis = znaxis[1:]
            image_header.set(naxis, self._header[znaxis],
                             self._header.comments[znaxis],
                             after=last_naxis)
            last_naxis = naxis

        # Delete any other spurious NAXISn keywords:
        naxis = image_header['NAXIS']
        for keyword in list(image_header['NAXIS?*']):
            try:
                n = int(keyword[5:])
            except Exception:
                continue

            if n > naxis:
                del image_header[keyword]

        # Although PCOUNT and GCOUNT are considered mandatory for IMAGE HDUs,
        # ZPCOUNT and ZGCOUNT are optional, probably because for IMAGE HDUs
        # their values are always 0 and 1 respectively
        if 'ZPCOUNT' in self._header:
            image_header.set('PCOUNT', self._header['ZPCOUNT'],
                             self._header.comments['ZPCOUNT'],
                             after=last_naxis)
        else:
            image_header.set('PCOUNT', 0, after=last_naxis)

        if 'ZGCOUNT' in self._header:
            image_header.set('GCOUNT', self._header['ZGCOUNT'],
                             self._header.comments['ZGCOUNT'],
                             after='PCOUNT')
        else:
            image_header.set('GCOUNT', 1, after='PCOUNT')

        if 'ZEXTEND' in self._header:
            image_header.set('EXTEND', self._header['ZEXTEND'],
                             self._header.comments['ZEXTEND'])

        if 'ZBLOCKED' in self._header:
            image_header.set('BLOCKED', self._header['ZBLOCKED'],
                             self._header.comments['ZBLOCKED'])

        # Move the ZHECKSUM and ZDATASUM cards to the image header
        # as CHECKSUM and DATASUM
        if 'ZHECKSUM' in self._header:
            image_header.set('CHECKSUM', self._header['ZHECKSUM'],
                             self._header.comments['ZHECKSUM'])

        if 'ZDATASUM' in self._header:
            image_header.set('DATASUM', self._header['ZDATASUM'],
                             self._header.comments['ZDATASUM'])

        # Remove the EXTNAME card if the value in the table header
        # is the default value of COMPRESSED_IMAGE.
        if ('EXTNAME' in self._header and
                self._header['EXTNAME'] == 'COMPRESSED_IMAGE'):
            del image_header['EXTNAME']

        # Look to see if there are any blank cards in the table
        # header.  If there are, there should be the same number
        # of blank cards in the image header.  Add blank cards to
        # the image header to make it so.
        table_blanks = self._header._countblanks()
        image_blanks = image_header._countblanks()

        for _ in range(table_blanks - image_blanks):
            image_header.append()

        # Create the CompImageHeader that syncs with the table header, and save
        # it off to self._image_header so it can be referenced later
        # unambiguously
        self._image_header = CompImageHeader(self._header, image_header)

        return self._image_header

    def _summary(self):
        """
        Summarize the HDU: name, dimensions, and formats.
        """
        class_name = self.__class__.__name__

        # if data is touched, use data info.
        if self._data_loaded:
            if self.data is None:
                _shape, _format = (), ''
            else:

                # the shape will be in the order of NAXIS's which is the
                # reverse of the numarray shape
                _shape = list(self.data.shape)
                _format = self.data.dtype.name
                _shape.reverse()
                _shape = tuple(_shape)
                _format = _format[_format.rfind('.') + 1:]

        # if data is not touched yet, use header info.
        else:
            _shape = ()

            for idx in range(self.header['NAXIS']):
                _shape += (self.header['NAXIS' + str(idx + 1)],)

            _format = BITPIX2DTYPE[self.header['BITPIX']]

        return (self.name, self.ver, class_name, len(self.header), _shape,
                _format)

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
        if _is_pseudo_unsigned(self.data.dtype):
            # Convert the unsigned array to signed
            self.data = np.array(
                self.data - _unsigned_zero(self.data.dtype),
                dtype='=i{}'.format(self.data.dtype.itemsize))
            should_swap = False
        else:
            should_swap = not self.data.dtype.isnative

        if should_swap:

            if self.data.flags.writeable:
                self.data.byteswap(True)
            else:
                # For read-only arrays, there is no way around making
                # a byteswapped copy of the data.
                self.data = self.data.byteswap(False)

        try:
            nrows = self._header['NAXIS2']
            tbsize = self._header['NAXIS1'] * nrows

            self._header['PCOUNT'] = 0
            if 'THEAP' in self._header:
                del self._header['THEAP']
            self._theap = tbsize

            # First delete the original compressed data, if it exists
            del self.compressed_data

            # Compress the data.
            # The current implementation of compress_hdu assumes the empty
            # compressed data table has already been initialized in
            # self.compressed_data, and writes directly to it
            # compress_hdu returns the size of the heap for the written
            # compressed image table
            heapsize, self.compressed_data = compression.compress_hdu(self)
        finally:
            # if data was byteswapped return it to its original order
            if should_swap:
                self.data.byteswap(True)
            self.data = old_data

        # CFITSIO will write the compressed data in big-endian order
        dtype = self.columns.dtype.newbyteorder('>')
        buf = self.compressed_data
        compressed_data = buf[:self._theap].view(dtype=dtype,
                                                 type=np.rec.recarray)
        self.compressed_data = compressed_data.view(FITS_rec)
        self.compressed_data._coldefs = self.columns
        self.compressed_data._heapoffset = self._theap
        self.compressed_data._heapsize = heapsize

    def scale(self, type=None, option='old', bscale=1, bzero=0):
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
        if (bscale != 1 or bzero != 0):
            _scale = bscale
            _zero = bzero
        else:
            if option == 'old':
                _scale = self._orig_bscale
                _zero = self._orig_bzero
            elif option == 'minmax':
                if isinstance(_type, np.floating):
                    _scale = 1
                    _zero = 0
                else:
                    _min = np.minimum.reduce(self.data.flat)
                    _max = np.maximum.reduce(self.data.flat)

                    if _type == np.uint8:  # uint8 case
                        _zero = _min
                        _scale = (_max - _min) / (2. ** 8 - 1)
                    else:
                        _zero = (_max + _min) / 2.

                        # throw away -2^N
                        _scale = (_max - _min) / (2. ** (8 * _type.bytes) - 2)

        # Do the scaling
        if _zero != 0:
            # We have to explicitly cast self._bzero to prevent numpy from
            # raising an error when doing self.data -= _zero, and we
            # do this instead of self.data = self.data - _zero to
            # avoid doubling memory usage.
            np.subtract(self.data, _zero, out=self.data, casting='unsafe')
            self.header['BZERO'] = _zero
        else:
            # Delete from both headers
            for header in (self.header, self._header):
                with suppress(KeyError):
                    del header['BZERO']

        if _scale != 1:
            self.data /= _scale
            self.header['BSCALE'] = _scale
        else:
            for header in (self.header, self._header):
                with suppress(KeyError):
                    del header['BSCALE']

        if self.data.dtype.type != _type:
            self.data = np.array(np.around(self.data), dtype=_type)  # 0.7.7.1

        # Update the BITPIX Card to match the data
        self._bitpix = DTYPE2BITPIX[self.data.dtype.name]
        self._bzero = self.header.get('BZERO', 0)
        self._bscale = self.header.get('BSCALE', 1)
        # Update BITPIX for the image header specifically
        # TODO: Make this more clear by using self._image_header, but only once
        # this has been fixed so that the _image_header attribute is guaranteed
        # to be valid
        self.header['BITPIX'] = self._bitpix

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
            self._update_uint_scale_keywords()

            # Shove the image header and data into a new ImageHDU and use that
            # to compute the image checksum
            image_hdu = ImageHDU(data=self.data, header=self.header)
            image_hdu._update_checksum(checksum)
            if 'CHECKSUM' in image_hdu.header:
                # This will also pass through to the ZHECKSUM keyword and
                # ZDATASUM keyword
                self._image_header.set('CHECKSUM',
                                       image_hdu.header['CHECKSUM'],
                                       image_hdu.header.comments['CHECKSUM'])
            if 'DATASUM' in image_hdu.header:
                self._image_header.set('DATASUM', image_hdu.header['DATASUM'],
                                       image_hdu.header.comments['DATASUM'])
            # Store a temporary backup of self.data in a different attribute;
            # see below
            self._imagedata = self.data

            # Now we need to perform an ugly hack to set the compressed data as
            # the .data attribute on the HDU so that the call to _writedata
            # handles it properly
            self.__dict__['data'] = self.compressed_data

        return super(CompImageHDU, self)._prewriteto(checksum=checksum,
                                                     inplace=inplace)

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
            return super(CompImageHDU, self)._writedata(fileobj)
        finally:
            # Restore the .data attribute to its rightful value (if any)
            if hasattr(self, '_imagedata'):
                self.__dict__['data'] = self._imagedata
                del self._imagedata
            else:
                del self.data

    def _close(self, closed=True):
        super(CompImageHDU, self)._close(closed=closed)

        # Also make sure to close access to the compressed data mmaps
        if (closed and self._data_loaded and
                _get_array_mmap(self.compressed_data) is not None):
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
            for bits, dtype in ((16, np.dtype('uint16')),
                                (32, np.dtype('uint32')),
                                (64, np.dtype('uint64'))):
                if bitpix == bits and self._orig_bzero == 1 << (bits - 1):
                    return dtype

        if bitpix > 16:  # scale integers to Float64
            return np.dtype('float64')
        elif bitpix > 0:  # scale integers to Float32
            return np.dtype('float32')

    def _update_header_scale_info(self, dtype=None):
        if (not self._do_not_scale_image_data and
                not (self._orig_bzero == 0 and self._orig_bscale == 1)):
            for keyword in ['BSCALE', 'BZERO']:
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
                self.header['BITPIX'] = DTYPE2BITPIX[dtype.name]

            self._bzero = 0
            self._bscale = 1
            self._bitpix = self.header['BITPIX']

    def _generate_dither_seed(self, seed):
        if not _is_int(seed):
            raise TypeError("Seed must be an integer")

        if not -1 <= seed <= 10000:
            raise ValueError(
                "Seed for random dithering must be either between 1 and "
                "10000 inclusive, 0 for autogeneration from the system "
                "clock, or -1 for autogeneration from a checksum of the first "
                "image tile (got {})".format(seed))

        if seed == DITHER_SEED_CHECKSUM:
            # Determine the tile dimensions from the ZTILEn keywords
            naxis = self._header['ZNAXIS']
            tile_dims = [self._header['ZTILE{}'.format(idx + 1)]
                         for idx in range(naxis)]
            tile_dims.reverse()

            # Get the first tile by using the tile dimensions as the end
            # indices of slices (starting from 0)
            first_tile = self.data[tuple(slice(d) for d in tile_dims)]

            # The checksum algorithm used is literally just the sum of the bytes
            # of the tile data (not its actual floating point values).  Integer
            # overflow is irrelevant.
            csum = first_tile.view(dtype='uint8').sum()

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
            return ((sum(int(x) for x in math.modf(time.time())) + id(self)) %
                    10000) + 1
        else:
            return seed
