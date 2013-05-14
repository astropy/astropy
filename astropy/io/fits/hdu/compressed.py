# Licensed under a 3-clause BSD style license - see PYFITS.rst

import operator
import sys
import warnings

import numpy as np

from functools import reduce

from .base import DELAYED, ExtensionHDU
from .image import _ImageBaseHDU, ImageHDU
from .table import BinTableHDU
from ..column import Column, ColDefs, _FormatP, _makep
from ..fitsrec import FITS_rec
from ..header import Header, BLOCK_SIZE
from ..util import _is_pseudo_unsigned, _unsigned_zero

from ....utils import lazyproperty, deprecated

try:
    from .. import compression
    COMPRESSION_SUPPORTED = COMPRESSION_ENABLED = True
except ImportError:
    COMPRESSION_SUPPORTED = COMPRESSION_ENABLED = False


# Default compression parameter values

DEFAULT_COMPRESSION_TYPE = 'RICE_1'
DEFAULT_QUANTIZE_LEVEL = 16.
DEFAULT_HCOMP_SCALE = 0
DEFAULT_HCOMP_SMOOTH = 0
DEFAULT_BLOCK_SIZE = 32
DEFAULT_BYTE_PIX = 4


# CFITSIO version-specific features
if COMPRESSION_SUPPORTED:
    try:
        CFITSIO_SUPPORTS_GZIPDATA = compression.CFITSIO_VERSION >= 3.28
    except AttributeError:
        # This generally shouldn't happen unless running setup.py in an
        # environment where an old build of pyfits exists
        CFITSIO_SUPPORTS_GZIPDATA = True


class CompImageHeader(Header):
    """
    Header object for compressed image HDUs designed to keep the compression
    header and the underlying image header properly synchronized.

    This essentially wraps the image header, so that all values are read from
    and written to the image header.  However, updates to the image header will
    also update the table header where appropriate.
    """

    def __init__(self, table_header, image_header=None):
        if image_header is None:
            image_header = Header()
        self._cards = image_header._cards
        self._keyword_indices = image_header._keyword_indices
        self._modified = image_header._modified
        self._table_header = table_header

    def set(self, keyword, value=None, comment=None, before=None, after=None):
        super(CompImageHeader, self).set(keyword, value, comment, before,
                                         after)

        # update the underlying header (_table_header) unless the update
        # was made to a card that describes the data.

        if (keyword not in ('SIMPLE', 'XTENSION', 'BITPIX', 'PCOUNT', 'GCOUNT',
                            'TFIELDS', 'EXTEND', 'ZIMAGE', 'ZBITPIX',
                            'ZCMPTYPE') and
            keyword[:4] not in ('ZVAL') and
            keyword[:5] not in ('NAXIS', 'TTYPE', 'TFORM', 'ZTILE', 'ZNAME')
                and keyword[:6] not in ('ZNAXIS')):
            self._table_header.set(keyword, value, comment, before, after)

    def add_history(self, value, before=None, after=None):
        super(CompImageHeader, self).add_history(value, before, after)
        self._table_header.add_history(value, before, after)

    def add_comment(self, value, before=None, after=None):
        super(CompImageHeader, self).add_comment(value, before, after)
        self._table_header.add_comment(value, before, after)

    def add_blank(self, value='', before=None, after=None):
        super(CompImageHeader, self).add_blank(value, before, after)
        self._table_header.add_blank(value, before, after)


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

    def __init__(self, data=None, header=None, name=None,
                 compression_type=DEFAULT_COMPRESSION_TYPE,
                 tile_size=None,
                 hcomp_scale=DEFAULT_HCOMP_SCALE,
                 hcomp_smooth=DEFAULT_HCOMP_SMOOTH,
                 quantize_level=DEFAULT_QUANTIZE_LEVEL,
                 do_not_scale_image_data=False,
                 uint=False, scale_back=False, **kwargs):
        """
        Parameters
        ----------
        data : array, optional
            data of the image

        header : Header instance, optional
            header to be associated with the image; when reading the HDU from a
            file (data=DELAYED), the header read from the file

        name : str, optional
            the ``EXTNAME`` value; if this value is `None`, then the name from
            the input image header will be used; if there is no name in the
            input image header then the default name ``COMPRESSED_IMAGE`` is
            used.

        compression_type : str, optional
            compression algorithm 'RICE_1', 'PLIO_1', 'GZIP_1', 'HCOMPRESS_1'

        tile_size : int, optional
            compression tile sizes.  Default treats each row of image as a
            tile.

        hcomp_scale : float, optional
            HCOMPRESS scale parameter

        hcomp_smooth : float, optional
            HCOMPRESS smooth parameter

        quantize_level : float, optional
            floating point quantization level; see note below

        Notes
        -----
        The astropy.io.fits package supports 2 methods of image compression:

            1) The entire FITS file may be externally compressed with the gzip
               or pkzip utility programs, producing a ``*.gz`` or ``*.zip``
               file, respectively.  When reading compressed files of this type,
               Astropy first uncompresses the entire file into a temporary file
               before performing the requested read operations.
               The astropy.io.fits package does not support writing to these
               types of compressed files.  This type of compression is
               supported in the `_File` class, not in the `CompImageHDU` class.
               The file compression type is recognized by the ``.gz`` or
               ``.zip`` file name extension.

            2) The `CompImageHDU` class supports the FITS tiled image
               compression convention in which the image is subdivided into a
               grid of rectangular tiles, and each tile of pixels is
               individually compressed.  The details of this FITS compression
               convention are described at the `FITS Support Office web site
               <http://fits.gsfc.nasa.gov/registry/tilecompression.html>`_.
               Basically, the compressed image tiles are stored in rows of a
               variable length arrray column in a FITS binary table.  The
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
        (PLIO).  The `compression_type` parameter defines the compression
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
        large ``hcomp_scale`` values, however, this can produce undesireable
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
        floating point value pixel values are not exactly perserved.  When done
        properly, this integer scaling technique will only discard the
        insignificant noise while still preserving all the real imformation in
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
        slightly for each tile.  In some cases, it may be desireable to specify
        the exact quantization level to be used, instead of specifying it
        relative to the calculated noise value.  This may be done by specifying
        the negative of desired quantization level for the value of
        ``quantize_level``.  In the previous example, one could specify
        ``quantize_level = -2.0`` so that the quantized integer levels differ
        by 2.0.  Larger negative values for ``quantize_level`` means that the
        levels are more coarsely-spaced, and will produce higher compression
        factors.
        """

        if not COMPRESSION_SUPPORTED:
            # TODO: Raise a more specific Exception type
            raise Exception('The astropy.io.fits.compression module is not '
                            'available.  Creation of compressed image HDUs is '
                            'disabled.')

        # Handle deprecated keyword arguments
        compression_opts = {}
        for oldarg, newarg in self.DEPRECATED_KWARGS.items():
            if oldarg in kwargs:
                warnings.warn('Keyword argument %s to %s is pending '
                              'deprecation; use %s instead' %
                              (oldarg, self.__class__.__name__, newarg),
                              PendingDeprecationWarning)
                compression_opts[newarg] = kwargs[oldarg]
                del kwargs[oldarg]
            else:
                compression_opts[newarg] = locals()[newarg]

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
                      for axis in xrange(self._header.get('ZNAXIS', 0))]

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
        if isinstance(xtension, basestring):
            xtension = xtension.rstrip()

        if xtension not in ('BINTABLE', 'A3DTABLE'):
            return False

        if 'ZIMAGE' not in header or header['ZIMAGE'] != True:
            return False

        if COMPRESSION_SUPPORTED and COMPRESSION_ENABLED:
            return True
        elif not COMPRESSION_SUPPORTED:
            warnings.warn('Failure matching header to a compressed image '
                          'HDU: The compression module is not available.\n'
                          'The HDU will be treated as a Binary Table HDU.')
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
                            quantize_level=None):
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
            compression algorithm 'RICE_1', 'PLIO_1', 'GZIP_1', 'HCOMPRESS_1';
            if this value is `None`, use value already in the header; if no
            value already in the header, use 'RICE_1'

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
        """

        image_hdu = ImageHDU(data=self.data, header=self._header)
        self._image_header = CompImageHeader(self._header, image_hdu.header)
        self._axes = image_hdu._axes
        del image_hdu

        # Update the extension name in the table header
        if not name and not 'EXTNAME' in self._header:
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
            if compression_type not in ['RICE_1', 'GZIP_1', 'PLIO_1',
                                        'HCOMPRESS_1']:
                warnings.warn('Unknown compression type provided.  Default '
                              '(%s) compression used.' %
                              DEFAULT_COMPRESSION_TYPE)
                compression_type = DEFAULT_COMPRESSION_TYPE

            self._header.set('ZCMPTYPE', compression_type,
                             'compression algorithm', after='TFIELDS')
        else:
            compression_type = self._header.get('ZCMPTYPE', 'RICE_1')

        # If the input image header had BSCALE/BZERO cards, then insert
        # them in the table header.

        if image_header:
            bzero = image_header.get('BZERO', 0.0)
            bscale = image_header.get('BSCALE', 1.0)
            afterCard = 'EXTNAME'

            if bscale != 1.0:
                self._header.set('BSCALE', bscale, after=afterCard)
                afterCard = 'BSCALE'

            if bzero != 0.0:
                self._header.set('BZERO', bzero, after=afterCard)

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
            tform1 = '1PI'
        else:
            tform1 = '1PB'

        self._header.set('TFORM1', tform1,
                         'data format of field: variable length array',
                         after='TTYPE1')

        # Create the first column for the table.  This column holds the
        # compressed data.
        col1 = Column(name=self._header['TTYPE1'], format=tform1)

        # Create the additional columns required for floating point
        # data and calculate the width of the output table.

        if self._image_header['BITPIX'] < 0 and quantize_level != 0.0:
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
                tform2 = '1PB'
            else:
                ttype2 = 'UNCOMPRESSED_DATA'
                if self._image_header['BITPIX'] == -32:
                    tform2 = '1PE'
                else:
                    tform2 = '1PD'

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
            keyList = ['TTYPE2', 'TFORM2', 'TTYPE3', 'TFORM3', 'TTYPE4',
                       'TFORM4']

            for k in keyList:
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
        self._header.set('NAXIS1', ncols * 8, 'width of table in bytes')
        self._header.set('TFIELDS', ncols, 'number of fields in each row')
        self._header.set('ZIMAGE', True, 'extension contains compressed image',
                         after=after)
        self._header.set('ZBITPIX', self._image_header['BITPIX'],
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
                          'Default tile size will be used.')
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
                major_dims = len(filter(lambda x: x > 1, tile_size))
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
        nrows = 1

        for idx, axis in enumerate(self._axes):
            naxis = 'NAXIS' + str(idx + 1)
            znaxis = 'ZNAXIS' + str(idx + 1)
            ztile = 'ZTILE' + str(idx + 1)

            if tile_size and len(tile_size) >= idx + 1:
                ts = tile_size[idx]
            else:
                if not ztile in self._header:
                    # Default tile size
                    if not idx:
                        ts = self._image_header['NAXIS1']
                    else:
                        ts = 1
                else:
                    ts = self._header[ztile]
                tile_size.append(ts)

            nrows = nrows * ((axis - 1) // ts + 1)

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

        afterCard = 'ZCMPTYPE'
        idx = 1

        if compression_type == 'RICE_1':
            self._header.set('ZNAME1', 'BLOCKSIZE', 'compression block size',
                             after=afterCard)
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
            afterCard = 'ZVAL2'
            idx = 3
        elif compression_type == 'HCOMPRESS_1':
            self._header.set('ZNAME1', 'SCALE', 'HCOMPRESS scale factor',
                             after=afterCard)
            self._header.set('ZVAL1', hcomp_scale, 'HCOMPRESS scale factor',
                             after='ZNAME1')
            self._header.set('ZNAME2', 'SMOOTH', 'HCOMPRESS smooth option',
                             after='ZVAL1')
            self._header.set('ZVAL2', hcomp_smooth, 'HCOMPRESS smooth option',
                             after='ZNAME2')
            afterCard = 'ZVAL2'
            idx = 3

        if self._image_header['BITPIX'] < 0:   # floating point image
            self._header.set('ZNAME' + str(idx), 'NOISEBIT',
                             'floating point quantization level',
                             after=afterCard)
            self._header.set('ZVAL' + str(idx), quantize_level,
                             'floating point quantization level',
                             after='ZNAME' + str(idx))

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

    @deprecated('3.2', alternative='(refactor your code)', pending=True)
    def updateHeaderData(self, image_header,
                         name=None,
                         compressionType=None,
                         tileSize=None,
                         hcompScale=None,
                         hcompSmooth=None,
                         quantizeLevel=None):
        self._update_header_data(image_header, name=name,
                                 compression_type=compressionType,
                                 tile_size=tileSize,
                                 hcomp_scale=hcompScale,
                                 hcomp_smooth=hcompSmooth,
                                 quantize_level=quantizeLevel)

    @lazyproperty
    def data(self):
        # The data attribute is the image data (not the table data).
        data = compression.decompress_hdu(self)

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
                data += self._bzero

            if zblank is not None:
                data = np.where(blanks, np.nan, data)

        # Right out of _ImageBaseHDU.data
        self._update_header_scale_info(data.dtype)

        return data

    @data.setter
    def data(self, data):
        if (data is not None) and (not isinstance(data, np.ndarray) or
           data.dtype.fields is not None):
                raise TypeError('CompImageHDU data has incorrect type:%s; '
                                'dtype.fields = %s' %
                                (type(data), data.dtype.fields))

    @lazyproperty
    def compressed_data(self):
        # First we will get the table data (the compressed
        # data) from the file, if there is any.
        compressed_data = super(BinTableHDU, self).data
        if isinstance(compressed_data, np.rec.recarray):
            del self.data
            return compressed_data
        else:
            # This will actually set self.compressed_data with the
            # pre-allocated space for the compression data; this is something I
            # might do away with in the future
            self._update_compressed_data()

        return self.compressed_data

    @lazyproperty
    @deprecated('3.2', alternative='the `.compressed_data attribute`',
                pending=True)
    def compData(self):
        return self.compressed_data

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
        # has already been defined we just return it.  If not, we nust
        # create it from the table header (the _header attribute).
        if hasattr(self, '_image_header'):
            return self._image_header

        # Start with a copy of the table header.
        self._image_header = CompImageHeader(self._header, self._header.copy())

        if 'XTENSION' in self._image_header:
            self._image_header['XTENSION'] = ('IMAGE', 'extension type')

        # Delete cards that are related to the table.  And move
        # the values of those cards that relate to the image from
        # their corresponding table cards.  These include
        # nnnZBITPIX -> BITPIX, ZNAXIS -> NAXIS, and ZNAXISn -> NAXISn.
        try:
            del self._image_header['ZIMAGE']
        except KeyError:
            pass

        try:
            del self._image_header['ZCMPTYPE']
        except KeyError:
            pass

        try:
            del self._image_header['ZBITPIX']
            _bitpix = self._header['ZBITPIX']
            self._image_header['BITPIX'] = (_bitpix,
                                            self._header.comments['ZBITPIX'])
        except KeyError:
            pass

        try:
            del self._image_header['ZNAXIS']
            self._image_header['NAXIS'] = (self._header['ZNAXIS'],
                                           self._header.comments['ZNAXIS'])

            last_naxis = 'NAXIS'
            for idx in range(self._image_header['NAXIS']):
                znaxis = 'ZNAXIS' + str(idx + 1)
                naxis = znaxis[1:]
                del self._image_header[znaxis]
                self._image_header.set(naxis, self._header[znaxis],
                                       self._header.comments[znaxis],
                                       after=last_naxis)
                last_naxis = naxis

            if last_naxis == 'NAXIS1':
                # There is only one axis in the image data so we
                # need to delete the extra NAXIS2 card.
                del self._image_header['NAXIS2']
        except KeyError:
            pass

        try:
            for idx in range(self._header['ZNAXIS']):
                del self._image_header['ZTILE' + str(idx + 1)]

        except KeyError:
            pass

        try:
            del self._image_header['ZPCOUNT']
            self._image_header.set('PCOUNT', self._header['ZPCOUNT'],
                                   self._header.comments['ZPCOUNT'])
        except KeyError:
            try:
                del self._image_header['PCOUNT']
            except KeyError:
                pass

        try:
            del self._image_header['ZGCOUNT']
            self._image_header.set('GCOUNT', self._header['ZGCOUNT'],
                                   self._header.comments['ZGCOUNT'])
        except KeyError:
            try:
                del self._image_header['GCOUNT']
            except KeyError:
                pass

        # Add the appropriate BSCALE and BZERO keywords if the data is scaled;
        # though these will be removed again as soon as the data is read
        # (unless do_not_scale_image_data=True)
        if 'GCOUNT' in self._image_header:
            after = 'GCOUNT'
        else:
            after = None
        if 'BSCALE' in self._header:
            self._image_header.set('BSCALE', self._header['BSCALE'],
                                   self._header.comments['BSCALE'],
                                   after=after)
            after = 'BSCALE'

        if 'BZERO' in self._header:
            self._image_header.set('BZERO', self._header['BZERO'],
                                   self._header.comments['BZERO'],
                                   after=after)

        try:
            del self._image_header['ZEXTEND']
            self._image_header.set('EXTEND', self._header['ZEXTEND'],
                                   self._header.comments['ZEXTEND'],
                                   after=last_naxis)
        except KeyError:
            pass

        try:
            del self._image_header['ZBLOCKED']
            self._image_header.set('BLOCKED', self._header['ZBLOCKED'],
                                   self._header.comments['ZBLOCKED'])
        except KeyError:
            pass

        try:
            del self._image_header['TFIELDS']

            for idx in range(self._header['TFIELDS']):
                del self._image_header['TFORM' + str(idx + 1)]
                ttype = 'TTYPE' + str(idx + 1)
                if ttype in self._image_header:
                    del self._image_header[ttype]

        except KeyError:
            pass

        idx = 1

        while True:
            try:
                del self._image_header['ZNAME' + str(idx)]
                del self._image_header['ZVAL' + str(idx)]
                idx += 1
            except KeyError:
                break

        # Move the ZHECKSUM and ZDATASUM cards to the image header
        # as CHECKSUM and DATASUM
        try:
            del self._image_header['ZHECKSUM']
            self._image_header.set('CHECKSUM', self._header['ZHECKSUM'],
                                   self._header.comments['ZHECKSUM'])
        except KeyError:
            pass

        try:
            del self._image_header['ZDATASUM']
            self._image_header.set('DATASUM', self._header['ZDATASUM'],
                                   self._header.comments['ZDATASUM'])
        except KeyError:
            pass

        try:
            del self._image_header['ZSIMPLE']
            self._image_header.set('SIMPLE', self._header['ZSIMPLE'],
                                   self._header.comments['ZSIMPLE'],
                                   before=1)
            del self._image_header['XTENSION']
        except KeyError:
            pass

        try:
            del self._image_header['ZTENSION']
            if self._header['ZTENSION'] != 'IMAGE':
                warnings.warn("ZTENSION keyword in compressed "
                              "extension != 'IMAGE'")
            self._image_header.set('XTENSION', 'IMAGE',
                                   self._header.comments['ZTENSION'])
        except KeyError:
            pass

        # Remove the EXTNAME card if the value in the table header
        # is the default value of COMPRESSED_IMAGE.

        if ('EXTNAME' in self._header and
                self._header['EXTNAME'] == 'COMPRESSED_IMAGE'):
            del self._image_header['EXTNAME']

        # Look to see if there are any blank cards in the table
        # header.  If there are, there should be the same number
        # of blank cards in the image header.  Add blank cards to
        # the image header to make it so.
        table_blanks = self._header._countblanks()
        image_blanks = self._image_header._countblanks()

        for _ in range(table_blanks - image_blanks):
            self._image_header.append()

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

            _format = _ImageBaseHDU.NumCode[self.header['BITPIX']]

        return (self.name, class_name, len(self.header), _shape,
                _format)

    def _update_compressed_data(self):
        """
        Compress the image data so that it may be written to a file.
        """

        # Check to see that the image_header matches the image data
        image_bitpix = _ImageBaseHDU.ImgCode[self.data.dtype.name]

        if (self.header.get('NAXIS', 0) != len(self.data.shape) or
                self.header.get('BITPIX', 0) != image_bitpix or
                self._header.get('ZNAXIS', 0) != len(self.data.shape) or
                self._header.get('ZBITPIX', 0) != image_bitpix or
                self.shape != self.data.shape):
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
                dtype='=i%d' % self.data.dtype.itemsize)
            should_swap = False
        else:
            should_swap = not self.data.dtype.isnative

        if should_swap:
            self.data.byteswap(True)

        try:
            nrows = self._header['NAXIS2']
            tbsize = self._header['NAXIS1'] * nrows

            self._header['PCOUNT'] = 0
            if 'THEAP' in self._header:
                del self._header['THEAP']
            self._theap = tbsize

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

        # Chances are not all the space allocated for the compressed data was
        # needed.  If not, go ahead and truncate the array:
        dataspan = tbsize + heapsize
        if len(self.compressed_data) > dataspan:
            if self.compressed_data.flags.owndata:
                self.compressed_data.resize(dataspan)
            else:
                # Need to copy to a new array; this generally shouldn't happen
                # at all though there are some contrived cases (such as in one
                # of the regression tests) where it can happen.
                self.compressed_data = np.resize(self.compressed_data,
                                                 (dataspan,))

        dtype = np.rec.format_parser(','.join(self.columns._recformats),
                                     self.columns.names, None).dtype
        # CFITSIO will write the compressed data in big-endian order
        dtype = dtype.newbyteorder('>')
        buf = self.compressed_data
        compressed_data = buf[:self._theap].view(dtype=dtype,
                                                 type=np.rec.recarray)
        self.compressed_data = compressed_data.view(FITS_rec)
        self.compressed_data._coldefs = self.columns
        self.compressed_data._heapoffset = self._theap
        self.compressed_data._heapsize = heapsize
        self.compressed_data._buffer = buf
        self.compressed_data.formats = self.columns.formats

        # Update the table header cards to match the compressed data.
        self._update_header()

    @deprecated('3.2', alternative='(refactor your code)', pending=True)
    def updateCompressedData(self):
        self._update_compressed_data()

    def _update_header(self):
        """
        Update the table header cards to match the compressed data.
        """

        # Get the _heapsize attribute to match the data.
        self.compressed_data._scale_back()

        # Check that TFIELDS and NAXIS2 match the data.
        self._header['TFIELDS'] = self.compressed_data._nfields
        self._header['NAXIS2'] = self.compressed_data.shape[0]

        # Calculate PCOUNT, for variable length tables.
        _tbsize = self._header['NAXIS1'] * self._header['NAXIS2']
        _heapstart = self._header.get('THEAP', _tbsize)
        self.compressed_data._gap = _heapstart - _tbsize
        _pcount = self.compressed_data._heapsize + self.compressed_data._gap

        if _pcount > 0:
            self._header['PCOUNT'] = _pcount

        # Update TFORM for variable length columns.
        for idx in range(self.compressed_data._nfields):
            format = self.compressed_data._coldefs._recformats[idx]
            if isinstance(format, _FormatP):
                _max = self.compressed_data.field(idx).max
                format = _FormatP(format.dtype, repeat=format.repeat, max=_max)
                self._header['TFORM' + str(idx + 1)] = format.tform
        # Insure that for RICE_1 that the BLOCKSIZE and BYTEPIX cards
        # are present and set to the hard coded values used by the
        # compression algorithm.
        if self._header['ZCMPTYPE'] == 'RICE_1':
            self._header.set('ZNAME1', 'BLOCKSIZE', 'compression block size',
                             after='ZCMPTYPE')
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

    @deprecated('3.2', alternative='(refactor your code)', pending=True)
    def updateHeader(self):
        self._update_header()

    def scale(self, type=None, option='old', bscale=1, bzero=0):
        """
        Scale image data by using ``BSCALE`` and ``BZERO``.

        Calling this method will scale `self.data` and update the keywords of
        ``BSCALE`` and ``BZERO`` in `self._header` and `self._image_header`.
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
            type = _ImageBaseHDU.NumCode[self._bitpix]
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
                    self.data.shape = dims

                    if _type == np.uint8:  # uint8 case
                        _zero = _min
                        _scale = (_max - _min) / (2. ** 8 - 1)
                    else:
                        _zero = (_max + _min) / 2.

                        # throw away -2^N
                        _scale = (_max - _min) / (2. ** (8 * _type.bytes) - 2)

        # Do the scaling
        if _zero != 0:
            self.data += -_zero
            self.header['BZERO'] = _zero
        else:
            # Delete from both headers
            for header in (self.header, self._header):
                try:
                    del header['BZERO']
                except KeyError:
                    pass

        if _scale != 1:
            self.data /= _scale
            self.header['BSCALE'] = _scale
        else:
            for header in (self.header, self._header):
                try:
                    del header['BSCALE']
                except KeyError:
                    pass

        if self.data.dtype.type != _type:
            self.data = np.array(np.around(self.data), dtype=_type)  # 0.7.7.1

        # Update the BITPIX Card to match the data
        self._bitpix = _ImageBaseHDU.ImgCode[self.data.dtype.name]
        self._bzero = self.header.get('BZERO', 0)
        self._bscale = self.header.get('BSCALE', 1)
        # Update BITPIX for the image header specificially
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

    # TODO: Fix this class so that it doesn't actually inherit from
    # BinTableHDU, but instead has an internal BinTableHDU reference
    def _prewriteto(self, checksum=False, inplace=False):
        if self._scale_back:
            self.scale(_ImageBaseHDU.NumCode[self._orig_bitpix])
        if self._data_loaded and self.data is not None:
            self._update_compressed_data()
        # Doesn't call the super's _prewriteto, since it calls
        # self.data._scale_back(), which is meaningless here.
        return ExtensionHDU._prewriteto(self, checksum=checksum,
                                        inplace=inplace)

    def _writeheader(self, fileobj):
        """
        Bypasses `BinTableHDU._writeheader()` which updates the header with
        metadata about the data that is meaningless here; another reason
        why this class maybe shouldn't inherit directly from BinTableHDU...
        """

        return ExtensionHDU._writeheader(self, fileobj)

    def _writedata_internal(self, fileobj):
        """
        Like the normal `BinTableHDU._writedata_internal`(), but we need to
        make sure the byte swap is done on the compressed data and not the
        image data, which requires a little messing with attributes.
        """

        size = 0

        if self.data is not None:
            imagedata = self.data
            # TODO: Ick; have to assign to __dict__ to bypass _setdata; need to
            # find a way to fix this
            self.__dict__['data'] = self.compressed_data
            # self.data = self.compressed_data
            try:
                size += self._binary_table_byte_swap(fileobj)
            finally:
                self.data = imagedata
            size += self.compressed_data.size * self.compressed_data.itemsize

        return size

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
                    try:
                        del header[keyword]
                        # Since _update_header_scale_info can, currently, be
                        # called *after* _prewriteto(), replace these with
                        # blank cards so the header size doesn't change
                        header.append()
                    except KeyError:
                        pass

            if dtype is None:
                dtype = self._dtype_for_bitpix()
            if dtype is not None:
                self.header['BITPIX'] = _ImageBaseHDU.ImgCode[dtype.name]

            self._bzero = 0
            self._bscale = 1
            self._bitpix = self.header['BITPIX']

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._data_loaded and self.data is not None:
            # We have the data to be used.
            return self._calculate_datasum_from_data(self.compressed_data,
                                                     blocking)
        else:
            # This is the case where the data has not been read from the
            # file yet.  We can handle that in a generic manner so we do
            # it in the base class.  The other possibility is that there
            # is no data at all.  This can also be handled in a generic
            # manner.
            return super(CompImageHDU, self)._calculate_datasum(blocking)
