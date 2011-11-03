import gzip

from pyfits.card import Card, CardList
from pyfits.file import _File
from pyfits.hdu.base import NonstandardExtHDU
from pyfits.hdu.hdulist import HDUList
from pyfits.header import Header
from pyfits.util import lazyproperty, BytesIO, fileobj_name


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
    def data(self):
        self._file.seek(self._datLoc)
        return self._file.readarray(self.size)

    @lazyproperty
    def hdulist(self):
        self._file.seek(self._datLoc)
        fileobj = BytesIO()
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
        compress : bool (optional)
            Gzip compress the FITS file
        """

        return cls.fromhdulist(HDUList(filename), compress=compress)

    @classmethod
    def fromhdulist(cls, hdulist, compress=False):
        """
        Creates a new FitsHDU from a given HDUList object.

        Parameters
        ----------
        hdulist : HDUList
            A valid Headerlet object.
        compress : bool (optional)
            Gzip compress the FITS file
        """

        fileobj = bs = BytesIO()
        if compress:
            if hasattr(hdulist, '_file'):
                name = fileobj_name(hdulist._file)
            else:
                name = None
            fileobj = gzip.GzipFile(name, mode='wb', fileobj=bs)
        hdulist.writeto(fileobj)
        if compress:
            fileobj.close()
        bs.seek(0)

        cards = [
            Card('XTENSION',  cls._extension, 'FITS extension'),
            Card('BITPIX',    8, 'array data type'),
            Card('NAXIS',     1, 'number of array dimensions'),
            Card('NAXIS1',    len(bs.getvalue()), 'Axis length'),
            Card('PCOUNT',    0, 'number of parameters'),
            Card('GCOUNT',    1, 'number of groups'),
        ]

        # Add the XINDn keywords proposed by Perry, though nothing is done with
        # these at the moment
        if len(hdulist) > 1:
            for idx, hdu in enumerate(hdulist[1:]):
                cards.append(Card('XIND' + str(idx + 1), hdu._hdrLoc,
                                  'byte offset of extension %d' % (idx + 1)))

        cards.append(Card('COMPRESS',  compress, 'Uses gzip compression'))
        header = Header(CardList(cards))
        # TODO: This wrapping of the fileobj should probably be handled by
        # cls.fromstring, though cls.fromstring itself has a strange
        # implementation that I probably need to fix.  For example, it
        # shouldn't care about fileobjs.  There should be a _BaseHDU.fromfile
        # for that (there is _BaseHDU.readfrom which plays that role, but its
        # semantics are also a little unclear...)
        return cls.fromstring(str(header), fileobj=_File(bs))

    @classmethod
    def match_header(cls, header):
        """
        This is a class method used in the pyfits refactoring branch to
        recognize whether or not this class should be used for instantiating
        an HDU object based on values in the header.

        It is included here for forward-compatibility.
        """

        card = header.ascard[0]
        if card.key != 'XTENSION':
            return False
        xtension = card.value
        if isinstance(xtension, basestring):
            xtension = xtension.rstrip()
        return xtension == cls._extension

    # TODO: Add header verification

    def _summary(self):
        # TODO: Perhaps make this more descriptive...
        return (self.name, self.__class__.__name__, len(self._header.ascard))

