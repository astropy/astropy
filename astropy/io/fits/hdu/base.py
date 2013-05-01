# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division


import datetime
import inspect
import os
import re
import warnings

import numpy as np

from .. import EXTENSION_NAME_CASE_SENSITIVE
from ..card import Card
from ..file import _File
from ..header import Header, HEADER_END_RE, BLOCK_SIZE, _pad_length
from ..util import (_is_int, _is_pseudo_unsigned, _unsigned_zero,
                    itersubclasses, encode_ascii, decode_ascii,
                    _get_array_mmap)
from ..verify import _Verify, _ErrList

from ....utils import lazyproperty, deprecated


class _Delayed(object):
    pass
DELAYED = _Delayed()


class InvalidHDUException(Exception):
    """
    A custom exception class used mainly to signal to _BaseHDU.__new__ that
    an HDU cannot possibly be considered valid, and must be assumed to be
    corrupted.
    """


def _hdu_class_from_header(cls, header):
    """
    Used primarily by _BaseHDU.__new__ to find an appropriate HDU class to use
    based on values in the header.  See the _BaseHDU.__new__ docstring.
    """

    klass = cls  # By default, if no subclasses are defined
    if header:
        for c in reversed(list(itersubclasses(cls))):
            try:
                # HDU classes built into astropy.io.fits are always considered,
                # but extension HDUs must be explicitly registered
                if not (c.__module__.startswith('astropy.io.fits.') or
                        c in cls._hdu_registry):
                    continue
                if c.match_header(header):
                    klass = c
                    break
            except NotImplementedError:
                continue
            except Exception as e:
                warnings.warn(
                    'An exception occurred matching an HDU header to the '
                    'appropriate HDU type: %s' % unicode(e))
                warnings.warn('The HDU will be treated as corrupted.')
                klass = _CorruptedHDU
                break

    return klass


# TODO: Come up with a better __repr__ for HDUs (and for HDULists, for that
# matter)
class _BaseHDU(object):
    """Base class for all HDU (header data unit) classes."""

    _hdu_registry = set()

    # This HDU type is part of the FITS standard
    _standard = True

    # Byte to use for padding out blocks
    _padding_byte = '\x00'

    def __new__(cls, data=None, header=None, *args, **kwargs):
        """
        Iterates through the subclasses of _BaseHDU and uses that class's
        match_header() method to determine which subclass to instantiate.

        It's important to be aware that the class hierarchy is traversed in a
        depth-last order.  Each match_header() should identify an HDU type as
        uniquely as possible.  Abstract types may choose to simply return False
        or raise NotImplementedError to be skipped.

        If any unexpected exceptions are raised while evaluating
        match_header(), the type is taken to be _CorruptedHDU.
        """

        klass = _hdu_class_from_header(cls, header)
        return super(_BaseHDU, cls).__new__(klass)

    def __init__(self, data=None, header=None, *args, **kwargs):
        self._header = header
        self._file = None
        self._buffer = None
        self._header_offset = None
        self._data_offset = None
        self._data_size = None
        self._new = True
        self._output_checksum = False

        if self._header:
            if 'DATASUM' in self._header and 'CHECKSUM' not in self._header:
                self._output_checksum = 'datasum'
            elif 'CHECKSUM' in self._header:
                self._output_checksum = True

        self._name = ''

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, value):
        self._header = value

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, basestring):
            raise TypeError("'name' attribute must be a string")
        self._name = value

    @property
    def is_image(self):
        return (
            self.name == 'PRIMARY' or
            ('XTENSION' in self._header and
             (self._header['XTENSION'] == 'IMAGE' or
              (self._header['XTENSION'] == 'BINTABLE' and
               'ZIMAGE' in self._header and self._header['ZIMAGE'] == True))))

    @property
    def _data_loaded(self):
        return ('data' in self.__dict__ and self.data is not None and
                self.data is not DELAYED)

    @property
    @deprecated('3.2', alternative='the `._header_offset` attribute',
                pending=True)
    def _hdrLoc(self):
        """The byte offset of this HDU's header in the file it came from;
        available for backwards compatibility--use ._header_offset instead.
        """

        return self._header_offset

    @_hdrLoc.setter
    @deprecated('3.2', alternative='the `._header_offset` attribute',
                pending=True)
    def _hdrLoc(self, value):
        self._header_offset = value

    @property
    @deprecated('3.2', alternative='the `._data_offset` attribute',
                pending=True)
    def _datLoc(self):
        """The byte offset of this HDU's data portion in the file it came from;
        available for backwards compatibility--use ._data_offset instead.
        """

        return self._data_offset

    @_datLoc.setter
    @deprecated('3.2', alternative='the `._data_offset` attribute',
                pending=True)
    def _datLoc(self, value):
        self._data_offset = value

    @property
    @deprecated('3.2', alternative='the `._data_size` attribute',
                pending=True)
    def _datSpan(self):
        """The byte size of this HDU's data portion in the file it came from;
        available for backwards compatibility--use ._data_size instead.
        """

        return self._data_size

    @_datSpan.setter
    @deprecated('3.2', alternative='the `._data_size` attribute',
                pending=True)
    def _datSpan(self, value):
        self._data_size = value

    @classmethod
    def register_hdu(cls, hducls):
        cls._hdu_registry.add(hducls)

    @classmethod
    def unregister_hdu(cls, hducls):
        if hducls in cls._hdu_registry:
            cls._hdu_registry.remove(hducls)

    @classmethod
    def match_header(cls, header):
        raise NotImplementedError

    # TODO: This method is a bit kludgy, especially in how it still usually
    # works with a fileobj, and how the first argument may or may not contain
    # HDU data.  This should be rethought.
    @classmethod
    def fromstring(cls, data, fileobj=None, offset=0, checksum=False,
                   ignore_missing_end=False, **kwargs):
        """
        Creates a new HDU object of the appropriate type from a string
        containing the HDU's entire header and, optionally, its data.

        Note: When creating a new HDU from a string without a backing file
        object, the data of that HDU may be read-only.  It depends on whether
        the underlying string was an immutable Python str/bytes object, or some
        kind of read-write memory buffer such as a `memoryview`.

        Parameters
        ----------
        data : str, bytearray, memoryview, ndarray
           A byte string contining the HDU's header and, optionally, its data.
           If `fileobj` is not specified, and the length of `data` extends
           beyond the header, then the trailing data is taken to be the HDU's
           data.  If `fileobj` is specified then the trailing data is ignored.

        fileobj : file (optional)
           The file-like object that this HDU was read from.

        offset : int (optional)
           If `fileobj` is specified, the offset into the file-like object at
           which this HDU begins.

        checksum : bool (optional)
           Check the HDU's checksum and/or datasum.

        ignore_missing_end : bool (optional)
           Ignore a missing end card in the header data.  Note that without
           the end card the end of the header can't be found, so the entire
           data is just assumed to be the header.

        kwargs : (optional)
           May contain additional keyword arguments specific to an HDU type.
           Any unrecognized kwargs are simply ignored.
        """

        if isinstance(data, Header):
            header = data
            if (not len(header) or
                    header.keys()[0] not in ('SIMPLE', 'XTENSION')):
                raise ValueError('Block does not begin with SIMPLE or '
                                 'XTENSION')
        else:
            try:
                # Test that the given object supports the buffer interface by
                # ensuring an ndarray can be created from it
                np.ndarray((), dtype='ubyte', buffer=data)
            except TypeError:
                raise TypeError(
                    'The provided object %r does not contain an underlying '
                    'memory buffer.  fromstring() requires an object that '
                    'supports the buffer interface such as bytes, str '
                    '(in Python 2.x but not in 3.x), buffer, memoryview, '
                    'ndarray, etc.' % data)

            if data[:8] not in [encode_ascii('SIMPLE  '),
                                encode_ascii('XTENSION')]:
                raise ValueError('Block does not begin with SIMPLE or '
                                 'XTENSION')

            # Make sure the end card is present
            for match in HEADER_END_RE.finditer(data):
                endpos = match.start()
                if endpos % Card.length == 0:
                    hdrlen = endpos + len(match.group())
                    hdrlen += _pad_length(hdrlen)
                    break
            else:
                if ignore_missing_end:
                    hdrlen = len(data)
                else:
                    raise ValueError('Header missing END card.')

            header = Header.fromstring(decode_ascii(data[:hdrlen]))
        # Determine the appropriate arguments to pass to the constructor from
        # self._kwargs.  self._kwargs contains any number of optional arguments
        # that may or may not be valid depending on the HDU type
        cls = _hdu_class_from_header(cls, header)
        args, varargs, varkwargs, defaults = inspect.getargspec(cls.__init__)
        new_kwargs = kwargs.copy()
        if not varkwargs:
            # If __init__ accepts arbitrary keyword arguments, then we can go
            # ahead and pass all keyword arguments; otherwise we need to delete
            # any that are invalid
            for key in kwargs:
                if key not in args:
                    del new_kwargs[key]

        hdu = cls(data=DELAYED, header=header, **new_kwargs)

        hdu._file = fileobj

        if not fileobj and len(data) > hdrlen:
            # Provide an underlying buffer to read the data from
            hdu._buffer = data

        hdu._header_offset = offset            # beginning of the header area
        if fileobj:
            hdu._data_offset = fileobj.tell()  # beginning of the data area
        else:
            hdu._data_offset = hdrlen

        # data area size, including padding
        size = hdu.size
        hdu._data_size = size + _pad_length(size)

        # Checksums are not checked on invalid HDU types
        if checksum and checksum != 'remove' and isinstance(hdu, _ValidHDU):
            hdu._verify_checksum_datasum(checksum)

        return hdu

    @classmethod
    def readfrom(cls, fileobj, checksum=False, ignore_missing_end=False,
                 **kwargs):
        """
        Read the HDU from a file.  Normally an HDU should be opened with
        `fitsopen()` which reads the entire HDU list in a FITS file.  But this
        method is still provided for symmetry with `writeto()`.

        Parameters
        ----------
        fileobj : file object or file-like object
            Input FITS file.  The file's seek pointer is assumed to be at the
            beginning of the HDU.

        checksum : bool
            If `True`, verifies that both ``DATASUM`` and
            ``CHECKSUM`` card values (when present in the HDU header)
            match the header and data of all HDU's in the file.

        ignore_missing_end : bool
            Do not issue an exception when opening a file that is
            missing an ``END`` card in the last header.
        """

        # TODO: Figure out a way to make it possible for the _File
        # constructor to be a noop if the argument is already a _File
        if not isinstance(fileobj, _File):
            fileobj = _File(fileobj)

        hdr_offset = fileobj.tell()
        hdr = Header.fromfile(fileobj, endcard=not ignore_missing_end)

        hdu = cls.fromstring(hdr, fileobj=fileobj, offset=hdr_offset,
                             checksum=checksum,
                             ignore_missing_end=ignore_missing_end, **kwargs)

        # If the checksum had to be checked the data may have already been read
        # from the file, in which case we don't want to see relative
        fileobj.seek(hdu._data_offset + hdu._data_size, os.SEEK_SET)
        return hdu

    def writeto(self, name, output_verify='exception', clobber=False,
                checksum=False):
        """
        Write the HDU to a new file.  This is a convenience method to
        provide a user easier output interface if only one HDU needs
        to be written to a file.

        Parameters
        ----------
        name : file path, file object or file-like object
            Output FITS file.  If opened, must be opened for append
            ("ab+")).

        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  See :ref:`verify` for more info.

        clobber : bool
            Overwrite the output file if exists.

        checksum : bool
            When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards
            to the header of the HDU when written to the file.
        """

        from .hdulist import HDUList

        hdulist = HDUList([self])
        hdulist.writeto(name, output_verify, clobber=clobber,
                        checksum=checksum)

    def _get_raw_data(self, shape, code, offset):
        """
        Return raw array from either the HDU's memory buffer or underlying
        file.
        """

        if isinstance(shape, int):
            shape = (shape,)

        if self._buffer:
            return np.ndarray(shape, dtype=code, buffer=self._buffer,
                              offset=offset)
        else:
            return self._file.readarray(offset=offset, dtype=code, shape=shape)

    # TODO: Rework checksum handling so that it's not necessary to add a
    # checksum argument here
    def _prewriteto(self, checksum=False, inplace=False):
        # If the data is unsigned int 16, 32, or 64 add BSCALE/BZERO
        # cards to header
        if self._data_loaded and self.data is not None and \
           self._standard and _is_pseudo_unsigned(self.data.dtype):
            if 'GCOUNT' in self._header:
                self._header.set('BSCALE', 1, after='GCOUNT')
            else:
                self._header.set('BSCALE', 1)
            self._header.set('BZERO', _unsigned_zero(self.data.dtype),
                             after='BSCALE')

        # Handle checksum

        # If the data is loaded it isn't necessarily 'modified', but we have no
        # way of knowing for sure
        modified = self._header._modified or self._data_loaded

        if checksum == 'remove':
            if 'CHECKSUM' in self._header:
                del self._header['CHECKSUM']

            if 'DATASUM' in self._header:
                del self._header['DATASUM']
        elif modified or self._new:
            if checksum == 'datasum':
                self.add_datasum()
            elif checksum == 'nonstandard_datasum':
                self.add_datasum(blocking='nonstandard')
            elif checksum == 'test':
                self.add_datasum(self._datasum_comment)
                self.add_checksum(self._checksum_comment, True)
            elif checksum == 'nonstandard':
                self.add_checksum(blocking='nonstandard')
            elif checksum:
                self.add_checksum(blocking='standard')

    def _postwriteto(self):
        # If data is unsigned integer 16, 32 or 64, remove the
        # BSCALE/BZERO cards
        if (self._data_loaded and self.data is not None and
                self._standard and _is_pseudo_unsigned(self.data.dtype)):
            for keyword in ('BSCALE', 'BZERO'):
                try:
                    del self._header[keyword]
                except KeyError:
                    pass

    def _writeheader(self, fileobj):
        offset = 0
        if not fileobj.simulateonly:
            try:
                offset = fileobj.tell()
            except (AttributeError, IOError):
                pass

            self._header.tofile(fileobj)

            try:
                size = fileobj.tell() - offset
            except (AttributeError, IOError):
                size = len(str(self._header))
        else:
            size = len(str(self._header))

        # Update hdrLoc with the new offset
        self._header_offset = offset

        return offset, size

    def _writedata(self, fileobj):
        # TODO: A lot of the simulateonly stuff should be moved back into the
        # _File class--basically it should turn write and flush into a noop
        offset = 0
        size = 0

        if not fileobj.simulateonly:
            fileobj.flush()
            try:
                offset = fileobj.tell()
            except IOError:
                offset = 0

        if self.data is not None:
            size += self._writedata_internal(fileobj)
            # pad the FITS data block
            if size > 0 and not fileobj.simulateonly:
                padding = _pad_length(size) * self._padding_byte
                # TODO: Not that this is ever likely, but if for some odd
                # reason _padding_byte is > 0x80 this will fail; but really if
                # somebody's custom fits format is doing that, they're doing it
                # wrong and should be reprimanded harshly.
                fileobj.write(padding.encode('ascii'))

        # flush, to make sure the content is written
        if not fileobj.simulateonly:
            fileobj.flush()

        # Update datLoc with the new offset
        self._data_offset = offset
        self._data_size = size = size + _pad_length(size)

        # return both the location and the size of the data area
        return offset, size

    def _writedata_internal(self, fileobj):
        """
        The beginning and end of most _writedata() implementations are the
        same, but the details of writing the data array itself can vary between
        HDU types, so that should be implemented in this method.

        Should return the size in bytes of the data written.
        """

        if not fileobj.simulateonly:
            fileobj.writearray(self.data)
        return self.data.size * self.data.itemsize

    # TODO: This is the start of moving HDU writing out of the _File class;
    # Though right now this is an internal private method (though still used by
    # HDUList, eventually the plan is to have this be moved into writeto()
    # somehow...
    def _writeto(self, fileobj, inplace=False, copy=False):
        # For now fileobj is assumed to be a _File object
        if not inplace or self._new:
            self._writeheader(fileobj)
            self._writedata(fileobj)
            return

        hdrloc = self._header_offset
        hdrsize = self._data_offset - self._header_offset
        datloc = self._data_offset
        datsize = self._data_size

        if self._header._modified:
            # Seek to the original header location in the file
            self._file.seek(hdrloc)
            # This should update hdrloc with he header location in the new file
            hdrloc, hdrsize = self._writeheader(fileobj)

            # If the data is to be written below with self._writedata, that
            # will also properly update the data location; but it should be
            # updated here too
            datloc = hdrloc + hdrsize
        elif copy:
            # Seek to the original header location in the file
            self._file.seek(hdrloc)
            # Before writing, update the hdrloc with the current file position,
            # which is the hdrloc for the new file
            hdrloc = fileobj.tell()
            fileobj.write(self._file.read(hdrsize))
            # The header size is unchanged, but the data location may be
            # different from before depending on if previous HDUs were resized
            datloc = fileobj.tell()
        if self._data_loaded:
            if self.data is not None:
                # Seek through the array's bases for an memmap'd array; we
                # can't rely on the _File object to give us this info since the
                # user may have replaced the previous mmap'd array
                if copy:
                    # Of course, if we're copying the data to a new file we
                    # don't care about flushing the original mmap; instead just
                    # read it into the new file
                    array_mmap = None
                else:
                    array_mmap = _get_array_mmap(self.data)

                if array_mmap is not None:
                    array_mmap.flush()
                else:
                    self._file.seek(self._data_offset)
                    datloc, datsize = self._writedata(fileobj)
        elif copy:
            # Seek to the data location in the original file
            self._file.seek(self._data_offset)
            fileobj.write(self._file.read(datsize))

        self._header_offset = hdrloc
        self._data_offset = datloc
        self._data_size = datsize + _pad_length(datsize)

_AllHDU = _BaseHDU  # For backwards-compatibility, though nobody should have
                    # been using this directly


# For convenience...
# TODO: register_hdu could be made into a class decorator which would be pretty
# cool, but only once 2.6 support is dropped.
register_hdu = _BaseHDU.register_hdu
unregister_hdu = _BaseHDU.unregister_hdu


class _CorruptedHDU(_BaseHDU):
    """
    A Corrupted HDU class.

    This class is used when one or more mandatory `Card`s are
    corrupted (unparsable), such as the ``BITPIX``, ``NAXIS``, or
    ``END`` cards.  A corrupted HDU usually means that the data size
    cannot be calculated or the ``END`` card is not found.  In the case
    of a missing ``END`` card, the `Header` may also contain the binary
    data

    .. note::
       In future, it may be possible to decipher where the last block
       of the `Header` ends, but this task may be difficult when the
       extension is a `TableHDU` containing ASCII data.
    """

    @property
    def size(self):
        """
        Returns the size (in bytes) of the HDU's data part.
        """

        # Note: On compressed files this might report a negative size; but the
        # file is corrupt anyways so I'm not too worried about it.
        if self._buffer is not None:
            return len(self._buffer) - self._data_offset

        return self._file.size - self._data_offset

    def _summary(self):
        return (self.name, 'CorruptedHDU')

    def verify(self):
        pass


class _NonstandardHDU(_BaseHDU, _Verify):
    """
    A Non-standard HDU class.

    This class is used for a Primary HDU when the ``SIMPLE`` Card has
    a value of `False`.  A non-standard HDU comes from a file that
    resembles a FITS file but departs from the standards in some
    significant way.  One example would be files where the numbers are
    in the DEC VAX internal storage format rather than the standard
    FITS most significant byte first.  The header for this HDU should
    be valid.  The data for this HDU is read from the file as a byte
    stream that begins at the first byte after the header ``END`` card
    and continues until the end of the file.
    """

    _standard = False

    @classmethod
    def match_header(cls, header):
        """
        Matches any HDU that has the 'SIMPLE' keyword but is not a standard
        Primary or Groups HDU.
        """

        # The SIMPLE keyword must be in the first card
        card = header.cards[0]

        # The check that 'GROUPS' is missing is a bit redundant, since the
        # match_header for GroupsHDU will always be called before this one.
        if card.keyword == 'SIMPLE':
            if 'GROUPS' not in header and card.value == False:
                return True
            else:
                raise InvalidHDUException
        else:
            return False

    @property
    def size(self):
        """
        Returns the size (in bytes) of the HDU's data part.
        """

        if self._buffer is not None:
            return len(self._buffer) - self._data_offset

        return self._file.size - self._data_offset

    def _writedata(self, fileobj):
        """
        Differs from the base class `_writedata()` in that it doesn't
        automatically add padding, and treats the data as a string of raw bytes
        instead of an array.
        """

        offset = 0
        size = 0

        if not fileobj.simulateonly:
            fileobj.flush()
            try:
                offset = fileobj.tell()
            except IOError:
                offset = 0

        if self.data is not None:
            if not fileobj.simulateonly:
                fileobj.write(self.data)
                # flush, to make sure the content is written
                fileobj.flush()
                size = len(self.data)

        # return both the location and the size of the data area
        return offset, size

    def _summary(self):
        return (self.name, 'NonstandardHDU', len(self._header))

    @lazyproperty
    def data(self):
        """
        Return the file data.
        """

        return self._get_raw_data(self.size, 'ubyte', self._data_offset)

    def _verify(self, option='warn'):
        errs = _ErrList([], unit='Card')

        # verify each card
        for card in self._header.cards:
            errs.append(card._verify(option))

        return errs


class _ValidHDU(_BaseHDU, _Verify):
    """
    Base class for all HDUs which are not corrupted.
    """

    def __init__(self, data=None, header=None, name=None, **kwargs):
        super(_ValidHDU, self).__init__(data=data, header=header)
        if name is not None:
            self.name = name

        if header and not hasattr(self, '_extver'):
            self._extver = header.get('EXTVER', 1)

    @property
    def name(self):
        # Convert the value to a string to be flexible in some pathological
        # cases (see ticket #96)
        if self._header and 'EXTNAME' in self._header:
            self._name = str(self._header['EXTNAME'])
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, basestring):
            raise TypeError("'name' attribute must be a string")
        if not EXTENSION_NAME_CASE_SENSITIVE():
            value = value.upper()
        if self._header:
            if 'EXTNAME' in self._header:
                self._header['EXTNAME'] = value
            else:
                self._header.append(('EXTNAME', value, 'extension name'))

    @classmethod
    def match_header(cls, header):
        """
        Matches any HDU that is not recognized as having either the SIMPLE or
        XTENSION keyword in its header's first card, but is nonetheless not
        corrupted.

        TODO: Maybe it would make more sense to use _NonstandardHDU in this
        case?  Not sure...
        """

        return header.keys()[0] not in ('SIMPLE', 'XTENSION')

    @property
    def size(self):
        """
        Size (in bytes) of the data portion of the HDU.
        """

        size = 0
        naxis = self._header.get('NAXIS', 0)
        if naxis > 0:
            size = 1
            for idx in range(naxis):
                size = size * self._header['NAXIS' + str(idx + 1)]
            bitpix = self._header['BITPIX']
            gcount = self._header.get('GCOUNT', 1)
            pcount = self._header.get('PCOUNT', 0)
            size = abs(bitpix) * gcount * (pcount + size) // 8
        return size

    def filebytes(self):
        """
        Calculates and returns the number of bytes that this HDU will write to
        a file.

        Parameters
        ----------
        None

        Returns
        -------
        Number of bytes
        """

        f = _File()
        # TODO: Fix this once new HDU writing API is settled on
        return self._writeheader(f)[1] + self._writedata(f)[1]

    def fileinfo(self):
        """
        Returns a dictionary detailing information about the locations
        of this HDU within any associated file.  The values are only
        valid after a read or write of the associated file with no
        intervening changes to the `HDUList`.

        Parameters
        ----------
        None

        Returns
        -------
        dictionary or None

           The dictionary details information about the locations of
           this HDU within an associated file.  Returns `None` when
           the HDU is not associated with a file.

           Dictionary contents:

           ========== ================================================
           Key        Value
           ========== ================================================
           file       File object associated with the HDU
           filemode   Mode in which the file was opened (readonly, copyonwrite,
                      update, append, ostream)
           hdrLoc     Starting byte location of header in file
           datLoc     Starting byte location of data block in file
           datSpan    Data size including padding
           ========== ================================================
        """

        if hasattr(self, '_file') and self._file:
            return {'file': self._file, 'filemode': self._file.mode,
                    'hdrLoc': self._header_offset, 'datLoc': self._data_offset,
                    'datSpan': self._data_size}
        else:
            return None

    def copy(self):
        """
        Make a copy of the HDU, both header and data are copied.
        """

        if self.data is not None:
            data = self.data.copy()
        else:
            data = None
        return self.__class__(data=data, header=self._header.copy())

    # TODO: self.name should be a property that updates the EXTNAME keyword
    # automatically; likewise for self.version.  Likewise it should
    # automatically get its value from the Header keyword.  This method should
    # just be deprecated, as should update_ext_version
    def update_ext_name(self, value, comment=None, before=None,
                        after=None, savecomment=False):
        """
        Update the extension name associated with the HDU.

        If the keyword already exists in the Header, it's value and/or comment
        will be updated.  If it does not exist, a new card will be created
        and it will be placed before or after the specified location.
        If no `before` or `after` is specified, it will be appended at
        the end.

        Parameters
        ----------
        value : str
            value to be used for the new extension name

        comment : str, optional
            to be used for updating, default=None.

        before : str or int, optional
            name of the keyword, or index of the `Card` before which
            the new card will be placed in the Header.  The argument
            `before` takes precedence over `after` if both specified.

        after : str or int, optional
            name of the keyword, or index of the `Card` after which
            the new card will be placed in the Header.

        savecomment : bool, optional
            When `True`, preserve the current comment for an existing
            keyword.  The argument `savecomment` takes precedence over
            `comment` if both specified.  If `comment` is not
            specified then the current comment will automatically be
            preserved.
        """

        if 'extname' in self._header and savecomment:
            comment = None

        self._header.set('EXTNAME', value, comment, before, after)
        self.name = value

    def update_ext_version(self, value, comment=None, before=None,
                           after=None, savecomment=False):
        """
        Update the extension version associated with the HDU.

        If the keyword already exists in the Header, it's value and/or comment
        will be updated.  If it does not exist, a new card will be created
        and it will be placed before or after the specified location.
        If no `before` or `after` is specified, it will be appended at
        the end.

        Parameters
        ----------
        value : str
            value to be used for the new extension version

        comment : str, optional
            to be used for updating, default=None.

        before : str or int, optional
            name of the keyword, or index of the `Card` before which
            the new card will be placed in the Header.  The argument
            `before` takes precedence over `after` if both specified.

        after : str or int, optional
            name of the keyword, or index of the `Card` after which
            the new card will be placed in the Header.

        savecomment : bool, optional
            When `True`, preserve the current comment for an existing
            keyword.  The argument `savecomment` takes precedence over
            `comment` if both specified.  If `comment` is not
            specified then the current comment will automatically be
            preserved.
        """

        if 'extver' in self._header and savecomment:
            comment = None

        self._header.set('extver', value, comment, before, after)
        self._extver = value

    def _verify(self, option='warn'):
        errs = _ErrList([], unit='Card')

        is_valid = lambda v: v in [8, 16, 32, 64, -32, -64]

        # Verify location and value of mandatory keywords.
        # Do the first card here, instead of in the respective HDU classes, so
        # the checking is in order, in case of required cards in wrong order.
        if isinstance(self, ExtensionHDU):
            firstkey = 'XTENSION'
            firstval = self._extension
        else:
            firstkey = 'SIMPLE'
            firstval = True

        self.req_cards(firstkey, 0, None, firstval, option, errs)
        self.req_cards('BITPIX', 1, lambda v: (_is_int(v) and is_valid(v)), 8,
                       option, errs)
        self.req_cards('NAXIS', 2,
                       lambda v: (_is_int(v) and v >= 0 and v <= 999), 0,
                       option, errs)

        naxis = self._header.get('NAXIS', 0)
        if naxis < 1000:
            for ax in range(3, naxis + 3):
                self.req_cards('NAXIS' + str(ax - 2), ax,
                               lambda v: (_is_int(v) and v >= 0), 1, option,
                               errs)

            # Remove NAXISj cards where j is not in range 1, naxis inclusive.
            for keyword in self._header:
                if keyword.startswith('NAXIS') and len(keyword) > 5:
                    try:
                        number = int(keyword[5:])
                        if number <= 0 or number > naxis:
                            raise ValueError
                    except ValueError:
                        err_text = ("NAXISj keyword out of range ('%s' when "
                                    "NAXIS == %d)" % (keyword, naxis))

                        def fix(self=self, keyword=keyword):
                            del self._header[keyword]

                        errs.append(
                            self.run_option(option=option, err_text=err_text,
                                            fix=fix, fix_text="Deleted."))

        # Verify that the EXTNAME keyword exists and is a string
        if 'EXTNAME' in self._header:
            if not isinstance(self._header['EXTNAME'], basestring):
                err_text = 'The EXTNAME keyword must have a string value.'
                fix_text = 'Converted the EXTNAME keyword to a string value.'

                def fix(header=self._header):
                    header['EXTNAME'] = str(header['EXTNAME'])

                errs.append(self.run_option(option, err_text=err_text,
                                            fix_text=fix_text, fix=fix))

        # verify each card
        for card in self._header.cards:
            errs.append(card._verify(option))

        return errs

    # TODO: Improve this API a little bit--for one, most of these arguments
    # could be optional
    def req_cards(self, keyword, pos, test, fix_value, option, errlist):
        """
        Check the existence, location, and value of a required `Card`.

        TODO: Write about parameters

        If `pos` = `None`, it can be anywhere.  If the card does not exist,
        the new card will have the `fix_value` as its value when created.
        Also check the card's value by using the `test` argument.
        """

        errs = errlist
        fix = None

        try:
            index = self._header.index(keyword)
        except ValueError:
            index = None

        fixable = fix_value is not None

        insert_pos = len(self._header) + 1

        # If pos is an int, insert at the given position (and convert it to a
        # lambda)
        if _is_int(pos):
            insert_pos = pos
            pos = lambda x: x == insert_pos

        # if the card does not exist
        if index is None:
            err_text = "'%s' card does not exist." % keyword
            fix_text = "Fixed by inserting a new '%s' card." % keyword
            if fixable:
                # use repr to accomodate both string and non-string types
                # Boolean is also OK in this constructor
                card = (keyword, fix_value)

                def fix(self=self, insert_pos=insert_pos, card=card):
                    self._header.insert(insert_pos, card)

            errs.append(self.run_option(option, err_text=err_text,
                        fix_text=fix_text, fix=fix, fixable=fixable))
        else:
            # if the supposed location is specified
            if pos is not None:
                if not pos(index):
                    err_text = ("'%s' card at the wrong place (card %d)." %
                                (keyword, index))

                    fix_text = ("Fixed by moving it to the right place "
                                "(card %d)." % insert_pos)

                    def fix(self=self, index=index, insert_pos=insert_pos):
                        card = self._header.cards[index]
                        del self._header[index]
                        self._header.insert(insert_pos, card)

                    errs.append(self.run_option(option, err_text=err_text,
                                fix_text=fix_text, fix=fix))

            # if value checking is specified
            if test:
                val = self._header[keyword]
                if not test(val):
                    err_text = ("'%s' card has invalid value '%s'." %
                                (keyword, val))
                    fix_text = "Fixed by setting a new value '%s'." % fix_value

                    if fixable:
                        def fix(self=self, keyword=keyword, val=fix_value):
                            self._header[keyword] = fix_value

                    errs.append(self.run_option(option, err_text=err_text,
                                fix_text=fix_text, fix=fix, fixable=fixable))

        return errs

    def add_datasum(self, when=None, blocking='standard'):
        """
        Add the ``DATASUM`` card to this HDU with the value set to the
        checksum calculated for the data.

        Parameters
        ----------
        when : str, optional
            Comment string for the card that by default represents the
            time when the checksum was calculated

        blocking: str, optional
            "standard" or "nonstandard", compute sum 2880 bytes at a time, or
            not

        Returns
        -------
        checksum : int
            The calculated datasum

        Notes
        -----
        For testing purposes, provide a `when` argument to enable the
        comment value in the card to remain consistent.  This will
        enable the generation of a ``CHECKSUM`` card with a consistent
        value.
        """

        cs = self._calculate_datasum(blocking)

        if when is None:
            when = 'data unit checksum updated %s' % self._get_timestamp()

        self._header['DATASUM'] = (str(cs), when)
        return cs

    def add_checksum(self, when=None, override_datasum=False,
                     blocking='standard'):
        """
        Add the ``CHECKSUM`` and ``DATASUM`` cards to this HDU with
        the values set to the checksum calculated for the HDU and the
        data respectively.  The addition of the ``DATASUM`` card may
        be overridden.

        Parameters
        ----------
        when : str, optional
           comment string for the cards; by default the comments
           will represent the time when the checksum was calculated

        override_datasum : bool, optional
           add the ``CHECKSUM`` card only

        blocking: str, optional
            "standard" or "nonstandard", compute sum 2880 bytes at a time, or
            not

        Notes
        -----
        For testing purposes, first call `add_datasum` with a `when`
        argument, then call `add_checksum` with a `when` argument and
        `override_datasum` set to `True`.  This will provide
        consistent comments for both cards and enable the generation
        of a ``CHECKSUM`` card with a consistent value.
        """

        if not override_datasum:
            # Calculate and add the data checksum to the header.
            data_cs = self.add_datasum(when, blocking)
        else:
            # Just calculate the data checksum
            data_cs = self._calculate_datasum(blocking)

        if when is None:
            when = 'HDU checksum updated %s' % self._get_timestamp()

        # Add the CHECKSUM card to the header with a value of all zeros.
        if 'DATASUM' in self._header:
            self._header.set('CHECKSUM', '0' * 16, when, before='DATASUM')
        else:
            self._header.set('CHECKSUM', '0' * 16, when)

        self._header['CHECKSUM'] = self._calculate_checksum(data_cs, blocking)

    def verify_datasum(self, blocking='standard'):
        """
        Verify that the value in the ``DATASUM`` keyword matches the value
        calculated for the ``DATASUM`` of the current HDU data.

        blocking: str, optional
            "standard" or "nonstandard", compute sum 2880 bytes at a time, or
            not

        Returns
        -------
        valid : int
           - 0 - failure
           - 1 - success
           - 2 - no ``DATASUM`` keyword present
        """

        if 'DATASUM' in self._header:
            datasum = self._calculate_datasum(blocking)
            if datasum == int(self._header['DATASUM']):
                return 1
            elif blocking == 'either':
                # i.e. standard failed,  try nonstandard
                return self.verify_datasum(blocking='nonstandard')
            else:
                # Failed with all permitted blocking kinds
                return 0
        else:
            return 2

    def verify_checksum(self, blocking='standard'):
        """
        Verify that the value in the ``CHECKSUM`` keyword matches the
        value calculated for the current HDU CHECKSUM.

        blocking: str, optional
            "standard" or "nonstandard", compute sum 2880 bytes at a time, or
            not

        Returns
        -------
        valid : int
           - 0 - failure
           - 1 - success
           - 2 - no ``CHECKSUM`` keyword present
        """

        if 'CHECKSUM' in self._header:
            if 'DATASUM' in self._header:
                datasum = self._calculate_datasum(blocking)
            else:
                datasum = 0
            checksum = self._calculate_checksum(datasum, blocking)
            if checksum == self._header['CHECKSUM']:
                return 1
            elif blocking == 'either':
                # i.e. standard failed,  try nonstandard
                return self.verify_checksum(blocking='nonstandard')
            else:
                # Failed with all permitted blocking kinds
                return 0
        else:
            return 2

    def _verify_checksum_datasum(self, blocking):
        """
        Verify the checksum/datasum values if the cards exist in the header.
        Simply displays warnings if either the checksum or datasum don't match.
        """

        # NOTE:  private data members _checksum and _datasum are
        # used by the utility script "fitscheck" to detect missing
        # checksums.

        if 'CHECKSUM' in self._header:
            self._checksum = self._header['CHECKSUM']
            self._checksum_comment = self._header.comments['CHECKSUM']
            if not self.verify_checksum(blocking):
                warnings.warn('Checksum verification failed for HDU %s.\n' %
                              ((self.name, self._extver),))
            del self._header['CHECKSUM']
        else:
            self._checksum = None
            self._checksum_comment = None

        if 'DATASUM' in self._header:
            self._datasum = self._header['DATASUM']
            self._datasum_comment = self._header.comments['DATASUM']

            if not self.verify_datasum(blocking):
                warnings.warn('Datasum verification failed for HDU %s.\n' %
                              ((self.name, self._extver),))
            del self._header['DATASUM']
        else:
            self._checksum = None
            self._checksum_comment = None
            self._datasum = None
            self._datasum_comment = None

    def _get_timestamp(self):
        """
        Return the current timestamp in ISO 8601 format, with microseconds
        stripped off.

        Ex.: 2007-05-30T19:05:11
        """

        return datetime.datetime.now().isoformat()[:19]

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if not self._data_loaded:
            # This is the case where the data has not been read from the file
            # yet.  We find the data in the file, read it, and calculate the
            # datasum.
            if self.size > 0:
                raw_data = self._get_raw_data(self._data_size, 'ubyte',
                                              self._data_offset)
                return self._compute_checksum(raw_data, blocking=blocking)
            else:
                return 0
        elif self.data is not None:
            return self._compute_checksum(self.data.view('ubyte'),
                                          blocking=blocking)
        else:
            return 0

    def _calculate_checksum(self, datasum, blocking):
        """
        Calculate the value of the ``CHECKSUM`` card in the HDU.
        """

        oldChecksum = self._header['CHECKSUM']
        self._header['CHECKSUM'] = '0' * 16

        # Convert the header to a string.
        s = str(self._header)

        # Calculate the checksum of the Header and data.
        cs = self._compute_checksum(np.fromstring(s, dtype='ubyte'), datasum,
                                    blocking=blocking)

        # Encode the checksum into a string.
        s = self._char_encode(~cs)

        # Return the header card value.
        self._header['CHECKSUM'] = oldChecksum

        return s

    def _compute_checksum(self, bytes, sum32=0, blocking="standard"):
        """
        Compute the ones-complement checksum of a sequence of bytes.

        Parameters
        ----------
        bytes
            a memory region to checksum

        sum32
            incremental checksum value from another region

        blocking
            "standard", "nonstandard", or "either"
            selects the block size on which to perform checksumming,
            originally the blocksize was chosen incorrectly.  "nonstandard"
            selects the original approach,  "standard" selects the
            interoperable blocking size of 2880 bytes.  In the context of
            _compute_checksum, "either" is synonymous with "standard".

        Returns
        -------
        ones complement checksum
        """

        blocklen = {'standard': 2880,
                    'nonstandard': len(bytes),
                    'either': 2880,  # do standard first
                    True: 2880}[blocking]

        sum32 = np.uint32(sum32)
        for i in range(0, len(bytes), blocklen):
            length = min(blocklen, len(bytes) - i)   # ????
            sum32 = self._compute_hdu_checksum(bytes[i:i + length], sum32)
        return sum32

    def _compute_hdu_checksum(self, bytes, sum32=0):
        """
        Translated from FITS Checksum Proposal by Seaman, Pence, and Rots.
        Use uint32 literals as a hedge against type promotion to int64.

        This code should only be called with blocks of 2880 bytes
        Longer blocks result in non-standard checksums with carry overflow
        Historically,  this code *was* called with larger blocks and for that
        reason still needs to be for backward compatibility.
        """

        u8 = np.uint32(8)
        u16 = np.uint32(16)
        uFFFF = np.uint32(0xFFFF)

        if bytes.nbytes % 2:
            last = bytes[-1]
            bytes = bytes[:-1]
        else:
            last = np.uint32(0)

        bytes = bytes.view('>u2')

        hi = sum32 >> u16
        lo = sum32 & uFFFF
        hi += np.add.reduce(bytes[0::2])
        lo += np.add.reduce(bytes[1::2])

        if (bytes.nbytes // 2) % 2:
            lo += last << u8
        else:
            hi += last << u8

        hicarry = hi >> u16
        locarry = lo >> u16

        while hicarry or locarry:
            hi = (hi & uFFFF) + locarry
            lo = (lo & uFFFF) + hicarry
            hicarry = hi >> u16
            locarry = lo >> u16

        return (hi << u16) + lo

    # _MASK and _EXCLUDE used for encoding the checksum value into a character
    # string.
    _MASK = [0xFF000000,
             0x00FF0000,
             0x0000FF00,
             0x000000FF]

    _EXCLUDE = [0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 0x40,
                0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60]

    def _encode_byte(self, byte):
        """
        Encode a single byte.
        """

        quotient = byte // 4 + ord('0')
        remainder = byte % 4

        ch = np.array(
            [(quotient + remainder), quotient, quotient, quotient],
            dtype='int32')

        check = True
        while check:
            check = False
            for x in self._EXCLUDE:
                for j in [0, 2]:
                    if ch[j] == x or ch[j + 1] == x:
                        ch[j] += 1
                        ch[j + 1] -= 1
                        check = True
        return ch

    def _char_encode(self, value):
        """
        Encodes the checksum `value` using the algorithm described
        in SPR section A.7.2 and returns it as a 16 character string.

        Parameters
        ----------
        value
            a checksum

        Returns
        -------
        ascii encoded checksum
        """

        value = np.uint32(value)

        asc = np.zeros((16,), dtype='byte')
        ascii = np.zeros((16,), dtype='byte')

        for i in range(4):
            byte = (value & self._MASK[i]) >> ((3 - i) * 8)
            ch = self._encode_byte(byte)
            for j in range(4):
                asc[4 * j + i] = ch[j]

        for i in range(16):
            ascii[i] = asc[(i + 15) % 16]

        return decode_ascii(ascii.tostring())


class ExtensionHDU(_ValidHDU):
    """
    An extension HDU class.

    This class is the base class for the `TableHDU`, `ImageHDU`, and
    `BinTableHDU` classes.
    """

    _extension = ''

    @classmethod
    def match_header(cls, header):
        """
        This class should never be instantiated directly.  Either a standard
        extension HDU type should be used for a specific extension, or
        NonstandardExtHDU should be used.
        """

        raise NotImplementedError

    def writeto(self, name, output_verify='exception', clobber=False,
                checksum=False):
        """
        Works similarly to the normal writeto(), but prepends a default
        `PrimaryHDU` are required by extension HDUs (which cannot stand on
        their own).
        """

        from .hdulist import HDUList
        from .image import PrimaryHDU

        hdulist = HDUList([PrimaryHDU(), self])
        hdulist.writeto(name, output_verify, clobber=clobber,
                        checksum=checksum)

    def _verify(self, option='warn'):

        errs = super(ExtensionHDU, self)._verify(option=option)

        # Verify location and value of mandatory keywords.
        naxis = self._header.get('NAXIS', 0)
        self.req_cards('PCOUNT', naxis + 3, lambda v: (_is_int(v) and v >= 0),
                       0, option, errs)
        self.req_cards('GCOUNT', naxis + 4, lambda v: (_is_int(v) and v == 1),
                       1, option, errs)

        return errs
# For backwards compatilibity, though this needs to be deprecated
# TODO: Mark this as deprecated
_ExtensionHDU = ExtensionHDU


class NonstandardExtHDU(ExtensionHDU):
    """
    A Non-standard Extension HDU class.

    This class is used for an Extension HDU when the ``XTENSION``
    `Card` has a non-standard value.  In this case, Astropy can figure
    out how big the data is but not what it is.  The data for this HDU
    is read from the file as a byte stream that begins at the first
    byte after the header ``END`` card and continues until the
    beginning of the next header or the end of the file.
    """

    _standard = False

    @classmethod
    def match_header(cls, header):
        """
        Matches any extension HDU that is not one of the standard extension HDU
        types.
        """

        card = header.cards[0]
        xtension = card.value
        if isinstance(xtension, basestring):
            xtension = xtension.rstrip()
        # A3DTABLE is not really considered a 'standard' extension, as it was
        # sort of the prototype for BINTABLE; however, since our BINTABLE
        # implementation handles A3DTABLE HDUs it is listed here.
        standard_xtensions = ('IMAGE', 'TABLE', 'BINTABLE', 'A3DTABLE')
        # The check that xtension is not one of the standard types should be
        # redundant.
        return (card.keyword == 'XTENSION' and
                xtension not in standard_xtensions)

    def _summary(self):
        return (self.name, 'NonstandardExtHDU', len(self._header))

    @lazyproperty
    def data(self):
        """
        Return the file data.
        """

        return self._get_raw_data(self.size, 'ubyte', self._data_offset)

# TODO: Mark this as deprecated
_NonstandardExtHDU = NonstandardExtHDU
