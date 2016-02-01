# Licensed under a 3-clause BSD style license - see PYFITS.rst

import gzip
import os

from ..file import _File
from ..header import _pad_length
from .base import _BaseHDU, BITPIX2DTYPE
from .hdulist import HDUList
from .image import PrimaryHDU
from ..util import fileobj_name


class StreamingHDU(object):
    """
    A class that provides the capability to stream data to a FITS file
    instead of requiring data to all be written at once.

    The following pseudocode illustrates its use::

        header = astropy.io.fits.Header()

        for all the cards you need in the header:
            header[key] = (value, comment)

        shdu = astropy.io.fits.StreamingHDU('filename.fits', header)

        for each piece of data:
            shdu.write(data)

        shdu.close()
    """

    def __init__(self, name, header):
        """
        Construct a `StreamingHDU` object given a file name and a header.

        Parameters
        ----------
        name : file path, file object, or file like object
            The file to which the header and data will be streamed.  If opened,
            the file object must be opened in a writeable binary mode such as
            'wb' or 'ab+'.

        header : `Header` instance
            The header object associated with the data to be written
            to the file.

        Notes
        -----
        The file will be opened and the header appended to the end of
        the file.  If the file does not already exist, it will be
        created, and if the header represents a Primary header, it
        will be written to the beginning of the file.  If the file
        does not exist and the provided header is not a Primary
        header, a default Primary HDU will be inserted at the
        beginning of the file and the provided header will be added as
        the first extension.  If the file does already exist, but the
        provided header represents a Primary header, the header will
        be modified to an image extension header and appended to the
        end of the file.
        """

        if isinstance(name, gzip.GzipFile):
            raise TypeError('StreamingHDU not supported for GzipFile objects.')

        self._header = header.copy()

        # handle a file object instead of a file name
        filename = fileobj_name(name) or ''

        # Check if the file already exists.  If it does not, check to see
        # if we were provided with a Primary Header.  If not we will need
        # to prepend a default PrimaryHDU to the file before writing the
        # given header.

        newfile = False

        if filename:
            if not os.path.exists(filename) or os.path.getsize(filename) == 0:
                newfile = True
        elif (hasattr(name, 'len') and name.len == 0):
            newfile = True

        if newfile:
            if 'SIMPLE' not in self._header:
                hdulist = HDUList([PrimaryHDU()])
                hdulist.writeto(name, 'exception')
        else:

            # This will not be the first extension in the file so we
            # must change the Primary header provided into an image
            # extension header.

            if 'SIMPLE' in self._header:
                self._header.set('XTENSION', 'IMAGE', 'Image extension',
                                 after='SIMPLE')
                del self._header['SIMPLE']

                if 'PCOUNT' not in self._header:
                    dim = self._header['NAXIS']

                    if dim == 0:
                        dim = ''
                    else:
                        dim = str(dim)

                    self._header.set('PCOUNT', 0, 'number of parameters',
                                     after='NAXIS' + dim)

                if 'GCOUNT' not in self._header:
                    self._header.set('GCOUNT', 1, 'number of groups',
                                     after='PCOUNT')

        self._ffo = _File(name, 'append')

        # TODO : Fix this once the HDU writing API is cleaned up
        tmp_hdu = _BaseHDU()
        # Passing self._header as an argument to _BaseHDU() will cause its
        # values to be modified in undesired ways...need to have a better way
        # of doing this
        tmp_hdu._header = self._header
        self._header_offset = tmp_hdu._writeheader(self._ffo)[0]
        self._data_offset = self._ffo.tell()
        self._size = self.size

        if self._size != 0:
            self.writecomplete = False
        else:
            self.writecomplete = True

    # Support the 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def write(self, data):
        """
        Write the given data to the stream.

        Parameters
        ----------
        data : ndarray
            Data to stream to the file.

        Returns
        -------
        writecomplete : int
            Flag that when `True` indicates that all of the required
            data has been written to the stream.

        Notes
        -----
        Only the amount of data specified in the header provided to the class
        constructor may be written to the stream.  If the provided data would
        cause the stream to overflow, an `~.exceptions.IOError` exception is
        raised and the data is not written. Once sufficient data has been
        written to the stream to satisfy the amount specified in the header,
        the stream is padded to fill a complete FITS block and no more data
        will be accepted. An attempt to write more data after the stream has
        been filled will raise an `~.exceptions.IOError` exception. If the
        dtype of the input data does not match what is expected by the header,
        a `.exceptions.TypeError` exception is raised.
        """

        size = self._ffo.tell() - self._data_offset

        if self.writecomplete or size + data.nbytes > self._size:
            raise IOError('Attempt to write more data to the stream than the '
                          'header specified.')

        if BITPIX2DTYPE[self._header['BITPIX']] != data.dtype.name:
            raise TypeError('Supplied data does not match the type specified '
                            'in the header.')

        if data.dtype.str[0] != '>':
            # byteswap little endian arrays before writing
            output = data.byteswap()
        else:
            output = data

        self._ffo.writearray(output)

        if self._ffo.tell() - self._data_offset == self._size:
            # the stream is full so pad the data to the next FITS block
            self._ffo.write(_pad_length(self._size) * '\0')
            self.writecomplete = True

        self._ffo.flush()

        return self.writecomplete

    @property
    def size(self):
        """
        Return the size (in bytes) of the data portion of the HDU.
        """

        size = 0
        naxis = self._header.get('NAXIS', 0)

        if naxis > 0:
            simple = self._header.get('SIMPLE', 'F')
            random_groups = self._header.get('GROUPS', 'F')

            if simple == 'T' and random_groups == 'T':
                groups = 1
            else:
                groups = 0

            size = 1

            for idx in range(groups, naxis):
                size = size * self._header['NAXIS' + str(idx + 1)]
            bitpix = self._header['BITPIX']
            gcount = self._header.get('GCOUNT', 1)
            pcount = self._header.get('PCOUNT', 0)
            size = abs(bitpix) * gcount * (pcount + size) // 8
        return size

    def close(self):
        """
        Close the physical FITS file.
        """

        self._ffo.close()
