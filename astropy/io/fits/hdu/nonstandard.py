# Licensed under a 3-clause BSD style license - see PYFITS.rst

import gzip
import io

from ..file import _File
from .base import NonstandardExtHDU
from .hdulist import HDUList
from ..header import Header, _pad_length
from ..util import fileobj_name

from ....extern.six import string_types
from ....utils import lazyproperty


class FitsHDU(NonstandardExtHDU):
    """
    A non-standard extension HDU for encapsulating entire FITS files within a
    single HDU of a container FITS file.  These HDUs have an extension (that is
    an XTENSION keyword) of FITS.

    The FITS file contained in the HDU's data can be accessed by the `hdulist`
    attribute which returns the contained FITS file as an `HDUList` object.
    """

    _extension = 'FITS'

    @lazyproperty
    def hdulist(self):
        self._file.seek(self._data_offset)
        fileobj = io.BytesIO()
        # Read the data into a BytesIO--reading directly from the file
        # won't work (at least for gzipped files) due to problems deep
        # within the gzip module that make it difficult to read gzip files
        # embedded in another file
        fileobj.write(self._file.read(self.size))
        fileobj.seek(0)
        if self._header['COMPRESS']:
            fileobj = gzip.GzipFile(fileobj=fileobj)
        return HDUList.fromfile(fileobj, mode='readonly')

    @classmethod
    def fromfile(cls, filename, compress=False):
        """
        Like `FitsHDU.fromhdulist()`, but creates a FitsHDU from a file on
        disk.

        Parameters
        ----------
        filename : str
            The path to the file to read into a FitsHDU
        compress : bool, optional
            Gzip compress the FITS file
        """

        return cls.fromhdulist(HDUList.fromfile(filename), compress=compress)

    @classmethod
    def fromhdulist(cls, hdulist, compress=False):
        """
        Creates a new FitsHDU from a given HDUList object.

        Parameters
        ----------
        hdulist : HDUList
            A valid Headerlet object.
        compress : bool, optional
            Gzip compress the FITS file
        """

        fileobj = bs = io.BytesIO()
        if compress:
            if hasattr(hdulist, '_file'):
                name = fileobj_name(hdulist._file)
            else:
                name = None
            fileobj = gzip.GzipFile(name, mode='wb', fileobj=bs)

        hdulist.writeto(fileobj)

        if compress:
            fileobj.close()

        # A proper HDUList should still be padded out to a multiple of 2880
        # technically speaking
        padding = (_pad_length(bs.tell()) * cls._padding_byte).encode('ascii')
        bs.write(padding)

        bs.seek(0)

        cards = [
            ('XTENSION', cls._extension, 'FITS extension'),
            ('BITPIX', 8, 'array data type'),
            ('NAXIS', 1, 'number of array dimensions'),
            ('NAXIS1', len(bs.getvalue()), 'Axis length'),
            ('PCOUNT', 0, 'number of parameters'),
            ('GCOUNT', 1, 'number of groups'),
        ]

        # Add the XINDn keywords proposed by Perry, though nothing is done with
        # these at the moment
        if len(hdulist) > 1:
            for idx, hdu in enumerate(hdulist[1:]):
                cards.append(('XIND' + str(idx + 1), hdu._header_offset,
                              'byte offset of extension {}'.format(idx + 1)))

        cards.append(('COMPRESS', compress, 'Uses gzip compression'))
        header = Header(cards)
        return cls._readfrom_internal(_File(bs), header=header)

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        if card.keyword != 'XTENSION':
            return False
        xtension = card.value
        if isinstance(xtension, string_types):
            xtension = xtension.rstrip()
        return xtension == cls._extension

    # TODO: Add header verification

    def _summary(self):
        # TODO: Perhaps make this more descriptive...
        return (self.name, self.ver, self.__class__.__name__, len(self._header))
