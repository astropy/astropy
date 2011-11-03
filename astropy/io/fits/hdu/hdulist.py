import gzip
import os
import signal
import sys
import threading
import warnings

import numpy as np
from numpy import memmap as Memmap

import pyfits
from pyfits.card import Card
from pyfits.column import _FormatP
from pyfits.file import PYTHON_MODES, _File
from pyfits.hdu import compressed
from pyfits.hdu.base import _BaseHDU, _ValidHDU, _NonstandardHDU, ExtensionHDU
from pyfits.hdu.compressed import CompImageHDU
from pyfits.hdu.groups import GroupsHDU
from pyfits.hdu.image import PrimaryHDU, ImageHDU
from pyfits.hdu.table import _TableBaseHDU
from pyfits.util import (Extendable, _is_int, _tmp_name, _with_extensions,
                         _pad_length, BLOCK_SIZE, isfile, fileobj_name,
                         fileobj_closed, fileobj_mode)
from pyfits.verify import _Verify, _ErrList


@_with_extensions
def fitsopen(name, mode="copyonwrite", memmap=None, classExtensions={},
             **kwargs):
    """Factory function to open a FITS file and return an `HDUList` object.

    Parameters
    ----------
    name : file path, file object or file-like object
        File to be opened.

    mode : str
        Open mode, 'copyonwrite' (default), 'readonly', 'update',
        'append', or 'ostream'.

        If `name` is a file object that is already opened, `mode` must
        match the mode the file was opened with, copyonwrite (rb),
        readonly (rb), update (rb+), append (ab+), ostream (w)).

    memmap : bool
        Is memory mapping to be used?

    classExtensions : dict (''Deprecated'')
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    kwargs : dict
        optional keyword arguments, possible values are:

        - **uint** : bool

            Interpret signed integer data where ``BZERO`` is the
            central value and ``BSCALE == 1`` as unsigned integer
            data.  For example, `int16` data with ``BZERO = 32768``
            and ``BSCALE = 1`` would be treated as `uint16` data.

            Note, for backward compatibility, the kwarg **uint16** may
            be used instead.  The kwarg was renamed when support was
            added for integers of any size.

        - **ignore_missing_end** : bool

            Do not issue an exception when opening a file that is
            missing an ``END`` card in the last header.

        - **checksum** : bool

            If `True`, verifies that both ``DATASUM`` and
            ``CHECKSUM`` card values (when present in the HDU header)
            match the header and data of all HDU's in the file.

        - **disable_image_compression** : bool

            If `True`, treates compressed image HDU's like normal
            binary table HDU's.

        - **do_not_scale_image_data** : bool

            If `True`, image data is not scaled using BSCALE/BZERO values
            when read.

    Returns
    -------
        hdulist : an `HDUList` object
            `HDUList` containing all of the header data units in the
            file.

    """

    if memmap is None:
        from pyfits.core import USE_MEMMAP
        memmap = USE_MEMMAP

    if 'uint16' in kwargs and 'uint' not in kwargs:
        kwargs['uint'] = kwargs['uint16']
        del kwargs['uint16']

    if not name:
        raise ValueError('Empty filename: %s' % repr(name))

    return HDUList.fromfile(name, mode, memmap, **kwargs)


class HDUList(list, _Verify):
    """
    HDU list class.  This is the top-level FITS object.  When a FITS
    file is opened, a `HDUList` object is returned.
    """

    __metaclass__ = Extendable

    def __init__(self, hdus=[], file=None):
        """
        Construct a `HDUList` object.

        Parameters
        ----------
        hdus : sequence of HDU objects or single HDU, optional
            The HDU object(s) to comprise the `HDUList`.  Should be
            instances of `_BaseHDU`.

        file : file object, optional
            The opened physical file associated with the `HDUList`.
        """

        self.__file = file
        if hdus is None:
            hdus = []

        # can take one HDU, as well as a list of HDU's as input
        if isinstance(hdus, _ValidHDU):
            hdus = [hdus]
        elif not isinstance(hdus, (HDUList, list)):
            raise TypeError("Invalid input for HDUList.")

        for idx, hdu in enumerate(hdus):
            if not isinstance(hdu, _BaseHDU):
                raise TypeError(
                      "Element %d in the HDUList input is not an HDU." % idx)
        super(HDUList, self).__init__(hdus)

    def __iter__(self):
        for idx in range(len(self)):
            yield self[idx]

    @_with_extensions
    def __getitem__(self, key, classExtensions={}):
        """
        Get an HDU from the `HDUList`, indexed by number or name.
        """

        if isinstance(key, slice):
            hdus = super(HDUList, self).__getitem__(key)
            return HDUList(hdus)

        idx = self.index_of(key)
        return super(HDUList, self).__getitem__(idx)

    def __setitem__(self, key, hdu):
        """
        Set an HDU to the `HDUList`, indexed by number or name.
        """

        _key = self.index_of(key)
        if isinstance(hdu, (slice, list)):
            if _is_int(_key):
                raise ValueError('An element in the HDUList must be an HDU.')
            for item in hdu:
                if not isinstance(item, _BaseHDU):
                    raise ValueError('%s is not an HDU.' % item)
        else:
            if not isinstance(hdu, _BaseHDU):
                raise ValueError('%s is not an HDU.' % hdu)

        try:
            super(HDUList, self).__setitem__(_key, hdu)
        except IndexError:
            raise IndexError('Extension %s is out of bound or not found.'
                             % key)
        self._resize = True
        self._truncate = False

    def __delitem__(self, key):
        """
        Delete an HDU from the `HDUList`, indexed by number or name.
        """

        if isinstance(key, slice):
            end_index = len(self)
        else:
            key = self.index_of(key)
            end_index = len(self) - 1

        super(HDUList, self).__delitem__(key)

        if (key == end_index or key == -1 and not self._resize):
            self._truncate = True
        else:
            self._truncate = False
            self._resize = True

    def __getslice__(self, start, end):
        return self[slice(start, end)]

    def __delslice__(self, start, stop):
        """
        Delete a slice of HDUs from the `HDUList`, indexed by number only.
        """

        del self[slice(start, stop)]

    # Support the 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    @classmethod
    def fromfile(cls, fileobj, mode="copyonwrite", memmap=False, **kwargs):
        """
        Creates an HDUList instance from a file-like object.

        The actual implementation of `fitsopen()`.
        """

        # instantiate a FITS file object (ffo)
        ffo = _File(fileobj, mode=mode, memmap=memmap)
        hdulist = cls(file=ffo)

        saved_compression_supported = compressed.COMPRESSION_SUPPORTED

        try:
            if 'disable_image_compression' in kwargs and \
               kwargs['disable_image_compression']:
                compressed.COMPRESSION_ENABLED = False

            if mode == 'ostream':
                # Output stream--not interested in reading/parsing the HDUs--just
                # writing to the output file
                return hdulist

            # read all HDUs
            while True:
                try:
                    hdu = _BaseHDU.readfrom(ffo, **kwargs)
                    hdulist.append(hdu)
                    hdu._new = False
                except EOFError:
                    break
                # check in the case there is extra space after the last HDU or
                # corrupted HDU
                except ValueError, err:
                    warnings.warn(
                        'Required keywords missing when trying to read '
                        'HDU #%d (note: PyFITS uses zero-based indexing.\n'
                        '          %s\n          There may be extra '
                        'bytes after the last HDU or the file is corrupted.' %
                        (len(hdulist), err))
                    break
                except IOError, err:
                    if ffo.writeonly:
                        break
                    else:
                        raise

            # If we're trying to read only and no header units were found,
            # raise and exception
            if mode == 'readonly' and len(hdulist) == 0:
                raise IOError('Empty FITS file')

            # initialize/reset attributes to be used in "update/append" mode
            # CardList needs its own _mod attribute since it has methods to change
            # the content of header without being able to pass it to the header
            # object
            hdulist._resize = False
            hdulist._truncate = False

        finally:
            compressed.COMPRESSION_SUPPORTED = saved_compression_supported

        return hdulist

    def fileinfo(self, index):
        """
        Returns a dictionary detailing information about the locations
        of the indexed HDU within any associated file.  The values are
        only valid after a read or write of the associated file with
        no intervening changes to the `HDUList`.

        Parameters
        ----------
        index : int
            Index of HDU for which info is to be returned.

        Returns
        -------
        dictionary or None

            The dictionary details information about the locations of
            the indexed HDU within an associated file.  Returns `None`
            when the HDU is not associated with a file.

            Dictionary contents:

            ========== =========================================================
            Key        Value
            ========== =========================================================
            file       File object associated with the HDU
            filename   Name of associated file object
            filemode   Mode in which the file was opened (readonly, copyonwrite,
                       update, append, ostream)
            resized    Flag that when `True` indicates that the data has been
                       resized since the last read/write so the returned values
                       may not be valid.
            hdrLoc     Starting byte location of header in file
            datLoc     Starting byte location of data block in file
            datSpan    Data size including padding
            ========== =========================================================

        """

        if self.__file is not None:
            output = self[index].fileinfo()

            if not output:
                # OK, the HDU associated with this index is not yet
                # tied to the file associated with the HDUList.  The only way
                # to get the file object is to check each of the HDU's in the
                # list until we find the one associated with the file.
                f = None

                for hdu in self:
                   info = hdu.fileinfo()

                   if info:
                      f = info['file']
                      fm = info['filemode']
                      break

                output = {'file': f, 'filemode': fm, 'hdrLoc': None,
                          'datLoc': None, 'datSpan': None}

            output['filename'] = self.__file.name
            output['resized'] = self._wasresized()
        else:
            output = None

        return output

    @_with_extensions
    def insert(self, index, hdu, classExtensions={}):
        """
        Insert an HDU into the `HDUList` at the given `index`.

        Parameters
        ----------
        index : int
            Index before which to insert the new HDU.

        hdu : _BaseHDU instance
            The HDU object to insert

        classExtensions : dict
            A dictionary that maps pyfits classes to extensions of those
            classes.  When present in the dictionary, the extension class
            will be constructed in place of the pyfits class.
        """

        if not isinstance(hdu, _BaseHDU):
            raise ValueError('%s is not an HDU.' % hdu)

        num_hdus = len(self)

        if index == 0 or num_hdus == 0:
            if num_hdus != 0:
                # We are inserting a new Primary HDU so we need to
                # make the current Primary HDU into an extension HDU.
                if isinstance(self[0], GroupsHDU):
                   raise ValueError(
                         "The current Primary HDU is a GroupsHDU.  "
                         "It can't be made into an extension HDU, "
                         "so you can't insert another HDU in front of it.")

                hdu1= ImageHDU(self[0].data, self[0].header)

                # Insert it into position 1, then delete HDU at position 0.
                super(HDUList, self).insert(1, hdu1)
                super(HDUList, self).__delitem__(0)

            if not isinstance(hdu, PrimaryHDU):
                # You passed in an Extension HDU but we need a Primary HDU.
                # If you provided an ImageHDU then we can convert it to
                # a primary HDU and use that.
                if isinstance(hdu, ImageHDU):
                    hdu = PrimaryHDU(hdu.data, hdu.header)
                else:
                    # You didn't provide an ImageHDU so we create a
                    # simple Primary HDU and append that first before
                    # we append the new Extension HDU.
                    phdu = PrimaryHDU()

                    super(HDUList, self).insert(0, phdu)
                    index = 1
        else:
            if isinstance(hdu, GroupsHDU):
               raise ValueError('A GroupsHDU must be inserted as a '
                                'Primary HDU.')

            if isinstance(hdu, PrimaryHDU):
                # You passed a Primary HDU but we need an Extension HDU
                # so create an Extension HDU from the input Primary HDU.
                hdu = ImageHDU(hdu.data, hdu.header)

        super(HDUList, self).insert(index, hdu)
        hdu._new = True
        self._resize = True
        self._truncate = False
        # make sure the EXTEND keyword is in primary HDU if there is extension
        if len(self) > 1:
            self.update_extend()

    @_with_extensions
    def append(self, hdu, classExtensions={}):
        """
        Append a new HDU to the `HDUList`.

        Parameters
        ----------
        hdu : instance of _BaseHDU
            HDU to add to the `HDUList`.

        classExtensions : dict
            A dictionary that maps pyfits classes to extensions of those
            classes.  When present in the dictionary, the extension class
            will be constructed in place of the pyfits class.
        """

        if not isinstance(hdu, _BaseHDU):
            raise ValueError('HDUList can only append an HDU.')

        if len(self) > 0:
            if isinstance(hdu, GroupsHDU):
               raise ValueError(
                     "Can't append a GroupsHDU to a non-empty HDUList")

            if isinstance(hdu, PrimaryHDU):
                # You passed a Primary HDU but we need an Extension HDU
                # so create an Extension HDU from the input Primary HDU.
                # TODO: This isn't necessarily sufficient to copy the HDU;
                # _hdrLoc and friends need to be copied too.
                hdu = ImageHDU(hdu.data, hdu.header)
        else:
            if not isinstance(hdu, PrimaryHDU):
                # You passed in an Extension HDU but we need a Primary
                # HDU.
                # If you provided an ImageHDU then we can convert it to
                # a primary HDU and use that.
                if isinstance(hdu, ImageHDU):
                    hdu = PrimaryHDU(hdu.data, hdu.header)
                else:
                    # You didn't provide an ImageHDU so we create a
                    # simple Primary HDU and append that first before
                    # we append the new Extension HDU.
                    phdu = PrimaryHDU()
                    super(HDUList, self).append(phdu)

        super(HDUList, self).append(hdu)
        hdu._new = True
        self._resize = True
        self._truncate = False

        # make sure the EXTEND keyword is in primary HDU if there is extension
        if len(self) > 1:
            self.update_extend()

    def index_of(self, key):
        """
        Get the index of an HDU from the `HDUList`.

        Parameters
        ----------
        key : int, str or tuple of (string, int)
           The key identifying the HDU.  If `key` is a tuple, it is of
           the form (`key`, `ver`) where `ver` is an ``EXTVER`` value
           that must match the HDU being searched for.

        Returns
        -------
        index : int
           The index of the HDU in the `HDUList`.
        """

        if _is_int(key):
            return key
        elif isinstance(key, tuple):
            _key, _ver = key
        else:
            _key = key
            _ver = None

        if not isinstance(_key, str):
            raise KeyError(key)
        _key = (_key.strip()).upper()

        nfound = 0
        found = None
        for idx, hdu in enumerate(self):
            name = hdu.name
            if isinstance(name, str):
                name = name.strip().upper()
            if name == _key and (_ver is None or _ver == hdu._extver):
                found = idx
                nfound += 1

        if (nfound == 0):
            raise KeyError('Extension %s not found.' % repr(key))
        elif (nfound > 1):
            raise KeyError('There are %d extensions of %s.'
                           % (nfound, repr(key)))
        else:
            return found

    def readall(self):
        """
        Read data of all HDUs into memory.
        """

        for hdu in self:
            if hdu.data is not None:
                continue

    @_with_extensions
    def flush(self, output_verify='exception', verbose=False,
              classExtensions={}):
        """
        Force a write of the `HDUList` back to the file (for append and
        update modes only).

        Parameters
        ----------
        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  See :ref:`verify` for more info.

        verbose : bool
            When `True`, print verbose messages

        classExtensions : dict
            A dictionary that maps pyfits classes to extensions of
            those classes.  When present in the dictionary, the
            extension class will be constructed in place of the pyfits
            class.
        """

        # Get the name of the current thread and determine if this is a single treaded application
        curr_thread = threading.currentThread()
        single_thread = (threading.activeCount() == 1) and \
                        (curr_thread.getName() == 'MainThread')

        # Define new signal interput handler
        if single_thread:
            keyboard_interrupt_sent = False
            def new_sigint(*args):
                warnings.warn('KeyboardInterrupt ignored until flush is '
                              'complete!')
                keyboard_interrupt_sent = True

            # Install new handler
            old_handler = signal.signal(signal.SIGINT, new_sigint)

        if self.__file.mode not in ('append', 'update', 'ostream'):
            warnings.warn("Flush for '%s' mode is not supported."
                          % self.__file.mode)
            return

        self.verify(option=output_verify)

        if self.__file.mode in ('append', 'ostream'):
            for hdu in self:
                if verbose:
                    try:
                        extver = str(hdu.header['extver'])
                    except KeyError:
                        extver = ''

                # only append HDU's which are "new"
                if not hasattr(hdu, '_new') or hdu._new:
                    # only output the checksum if flagged to do so
                    if hasattr(hdu, '_output_checksum'):
                        checksum = hdu._output_checksum
                    else:
                        checksum = False

                    # TODO: Fix this once new HDU writing API is settled on
                    hdu._writeto(self.__file, checksum=checksum)
                    if verbose:
                        print 'append HDU', hdu.name, extver
                    hdu._new = False

        elif self.__file.mode == 'update':
            self._wasresized(verbose)

            # TODO: Much of this section should probably be handled in _File

            # if the HDUList is resized, need to write out the entire contents
            # of the hdulist to the file.
            if self._resize or self.__file.compression:
                old_name = self.__file.name
                old_memmap = self.__file.memmap
                name = _tmp_name(old_name)

                if not self.__file.file_like:
                    old_mode = os.stat(old_name).st_mode
                    #
                    # The underlying file is an acutal file object.
                    # The HDUList is resized, so we need to write it to a tmp
                    # file, delete the original file, and rename the tmp
                    # file to the original file.
                    #
                    if self.__file.compression == 'gzip':
                        new_file = gzip.GzipFile(name, mode='ab+')
                    else:
                        new_file = name

                    hdulist = self.fromfile(new_file, mode='append')

                    if verbose:
                        print 'open a temp file', name

                    for hdu in self:
                        # only output the checksum if flagged to do so
                        if hasattr(hdu, '_output_checksum'):
                            checksum = hdu._output_checksum
                        else:
                            checksum = False

                        # TODO: Fix this once new HDU writing API is settled on
                        (hdu._hdrLoc, hdu._datLoc, hdu._datSpan) = \
                               hdu._writeto(hdulist.__file, checksum=checksum)
                    hdulist.__file.close()
                    self.__file.close()
                    os.remove(self.__file.name)

                    if verbose:
                        print 'delete the original file', old_name

                    # reopen the renamed new file with "update" mode
                    os.rename(name, old_name)
                    os.chmod(old_name, old_mode)

                    if isinstance(new_file, gzip.GzipFile):
                        old_file = gzip.GzipFile(old_name, mode='rb+')
                    else:
                        old_file = old_name

                    ffo = _File(old_file, mode='update', memmap=old_memmap)

                    self.__file = ffo
                    if verbose:
                        print 'reopen the newly renamed file', old_name
                else:
                    #
                    # The underlying file is not a file object, it is a file
                    # like object.  We can't write out to a file, we must
                    # update the file like object in place.  To do this,
                    # we write out to a temporary file, then delete the
                    # contents in our file like object, then write the
                    # contents of the temporary file to the now empty file
                    # like object.
                    #
                    self.writeto(name)
                    hdulist = self.fromfile(name)
                    ffo = self.__file

                    ffo.truncate(0)
                    ffo.seek(0)

                    for hdu in hdulist:
                        # only output the checksum if flagged to do so
                        if hasattr(hdu, '_output_checksum'):
                            checksum = hdu._output_checksum
                        else:
                            checksum = False

                        # TODO: Fix this once new HDU writing API is settled on
                        (hdu._hdrLoc, hdu._datLoc, hdu._datSpan) = \
                                hdu._writeto(ffo, checksum=checksum)

                    # Close the temporary file and delete it.
                    hdulist.close()
                    os.remove(hdulist.__file.name)

                # reset the resize attributes after updating
                self._resize = False
                self._truncate = False
                for hdu in self:
                    hdu.header._mod = False
                    hdu.header.ascard._mod = False
                    hdu._new = False
                    hdu._file = ffo

            # if not resized, update in place
            else:
                for hdu in self:
                    if verbose:
                        try:
                            extver = str(hdu.header['extver'])
                        except KeyError:
                            extver = ''

                    if hdu.header._mod or hdu.header.ascard._mod:
                        # only output the checksum if flagged to do so
                        if hasattr(hdu, '_output_checksum'):
                            checksum = hdu._output_checksum
                        else:
                            checksum = False

                        hdu._file.seek(hdu._hdrLoc)
                        # TODO: Fix this once new HDU writing API is settled on
                        hdu._writeheader(self.__file, checksum=checksum)
                        if verbose:
                            print 'update header in place: Name =', \
                                  hdu.name, extver
                    if hdu._data_loaded:
                        if hdu.data is not None:
                            if isinstance(hdu.data, Memmap):
                                hdu.data.flush()
                            else:
                                hdu._file.seek(hdu._datLoc)
                                # TODO: Fix this once new HDU writing API is settled on
                                hdu._writedata(self.__file)

                            if verbose:
                                print 'update data in place: Name =', \
                                      hdu.name, extver

                # reset the modification attributes after updating
                for hdu in self:
                    hdu.header._mod = False
                    hdu.header.ascard._mod = False

        if single_thread:
            if keyboard_interrupt_sent:
                raise KeyboardInterrupt

            if old_handler is not None:
                signal.signal(signal.SIGINT, old_handler)
            else:
                signal.signal(signal.SIGINT, signal.SIG_DFL)

    def update_extend(self):
        """
        Make sure that if the primary header needs the keyword
        ``EXTEND`` that it has it and it is correct.
        """
        hdr = self[0].header
        if 'extend' in hdr:
            if (hdr['extend'] == False):
                hdr['extend'] = True
        else:
            if hdr['naxis'] == 0:
                hdr.update('extend', True, after='naxis')
            else:
                n = hdr['naxis']
                hdr.update('extend', True, after='naxis' + str(n))

    @_with_extensions
    def writeto(self, fileobj, output_verify='exception', clobber=False,
                classExtensions={}, checksum=False):
        """
        Write the `HDUList` to a new file.

        Parameters
        ----------
        fileobj : file path, file object or file-like object
            File to write to.  If a file object, must be opened for
            append (ab+).

        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  See :ref:`verify` for more info.

        clobber : bool
            When `True`, overwrite the output file if exists.

        classExtensions : dict
            A dictionary that maps pyfits classes to extensions of
            those classes.  When present in the dictionary, the
            extension class will be constructed in place of the pyfits
            class.

        checksum : bool
            When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards
            to the headers of all HDU's written to the file.
        """

        if (len(self) == 0):
            warnings.warn("There is nothing to write.")
            return

        if output_verify == 'warn':
            output_verify = 'exception'
        self.verify(option=output_verify)

        # check if the file object is closed
        closed = fileobj_closed(fileobj)
        fmode = fileobj_mode(fileobj) or 'ab+'
        filename = fileobj_name(fileobj)

        # check if the output file already exists
        if (isfile(fileobj) or
            isinstance(fileobj, (basestring, gzip.GzipFile))):
            if (os.path.exists(filename) and os.path.getsize(filename) != 0):
                if clobber:
                    warnings.warn("Overwriting existing file '%s'." % filename)
                    if not closed:
                        fileobj.close()
                    os.remove(filename)
                else:
                    raise IOError("File '%s' already exists." % filename)
        elif (hasattr(fileobj, 'len') and fileobj.len > 0):
            if clobber:
                warnings.warn("Overwriting existing file '%s'." % filename)
                name.truncate(0)
            else:
                raise IOError("File '%s' already exists." % filename)

        # make sure the EXTEND keyword is there if there is extension
        if len(self) > 1:
            self.update_extend()

        mode = 'copyonwrite'
        for key, val in PYTHON_MODES.iteritems():
            if val == fmode:
                mode = key
                break

        hdulist = fitsopen(fileobj, mode=mode)

        for hdu in self:
            # TODO: Fix this once new HDU writing API is settled on
            hdu._writeto(hdulist.__file, checksum)
        hdulist.close(output_verify=output_verify, closed=closed)


    def close(self, output_verify='exception', verbose=False, closed=True):
        """
        Close the associated FITS file and memmap object, if any.

        Parameters
        ----------
        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  See :ref:`verify` for more info.

        verbose : bool
            When `True`, print out verbose messages.

        closed : bool
            When `True`, close the underlying file object.
        """

        if self.__file:
            if self.__file.mode in ['append', 'update']:
                self.flush(output_verify=output_verify, verbose=verbose)

            if closed and hasattr(self.__file, 'close'):
                self.__file.close()

    def info(self, output=None):
        """
        Summarize the info of the HDUs in this `HDUList`.

        Note that this function prints its results to the console---it
        does not return a value.

        Parameters
        ----------
        output : file, optional
            A file-like object to write the output to.  If False, does not
            output to a file and instead returns a list of tuples representing
            the HDU info.  Writes to sys.stdout by default.
        """

        if output is None:
            output = sys.stdout

        if self.__file is None:
            name = '(No file associated with this HDUList)'
        else:
            name = self.__file.name

        results = ['Filename: %s' % name,
                   'No.    Name         Type      Cards   Dimensions   Format']

        format = '%-3d  %-10s  %-11s  %5d   %-10s   %s%s'
        default = ('', '', 0, '()', '', '')
        for idx, hdu in enumerate(self):
            summary = hdu._summary()
            if len(summary) < len(default):
                summary += default[len(summary):]
            summary = (idx,) + summary
            if output:
                results.append(format % summary)
            else:
                results.append(summary)

        if output:
            output.write('\n'.join(results))
            output.write('\n')
            output.flush()
        else:
            return results[2:]

    def filename(self):
        """
        Return the file name associated with the HDUList object if one exists.
        Otherwise returns None.

        Returns
        -------
        filename : a string containing the file name associated with the
                   HDUList object if an association exists.  Otherwise returns
                   None.
        """
        if self.__file is not None:
           if hasattr(self.__file, 'name'):
              return self.__file.name
        return None

    def _verify(self, option='warn'):
        text = ''
        errs = _ErrList([], unit='HDU')

        # the first (0th) element must be a primary HDU
        if len(self) > 0 and (not isinstance(self[0], PrimaryHDU)) and \
                             (not isinstance(self[0], _NonstandardHDU)):
            err_text = "HDUList's 0th element is not a primary HDU."
            fix_text = 'Fixed by inserting one as 0th HDU.'

            def fix(self=self):
                self.insert(0, PrimaryHDU())

            err = self.run_option(option, err_text=err_text,
                                  fix_text=fix_text, fix=fix)
            errs.append(err)

        if len(self) > 1 and ('EXTEND' not in self[0].header or
                              self[0].header['EXTEND'] is not True):
            err_text = ('Primary HDU does not contain an EXTEND keyword '
                        'equal to T even though there are extension HDUs.')
            fix_text = 'Fixed by inserting or updating the EXTEND keyword.'

            def fix(header=self[0].header):
                naxis = header['NAXIS']
                if naxis == 0:
                    after = 'NAXIS'
                else:
                    after = 'NAXIS' + str(naxis)
                header.update('EXTEND', value=True, after=after)

            errs.append(self.run_option(option, err_text=err_text,
                                        fix_text=fix_text, fix=fix))

        # each element calls their own verify
        for idx, hdu in enumerate(self):
            if idx > 0 and (not isinstance(hdu, ExtensionHDU)):
                err_text = "HDUList's element %s is not an extension HDU." \
                           % str(idx)
                err = self.run_option(option, err_text=err_text, fixable=True)
                errs.append(errs)

            else:
                result = hdu._verify(option)
                if result:
                    errs.append(result)
        return errs

    def _wasresized(self, verbose=False):
        """
        Determine if any changes to the HDUList will require a file resize
        when flushing the file.

        Side effect of setting the objects _resize attribute.
        """

        if not self._resize:

            # determine if any of the HDU is resized
            for hdu in self:

                # Header:
                # Add 1 to .ascard to include the END card
                _nch80 = sum([card._ncards() for card in hdu.header.ascard])
                _bytes = (_nch80+1) * Card.length
                _bytes = _bytes + _pad_length(_bytes)
                if _bytes != (hdu._datLoc-hdu._hdrLoc):
                    self._resize = True
                    self._truncate = False
                    if verbose:
                        print 'One or more header is resized.'
                    break

                # Data:
                if not hdu._data_loaded or hdu.data is None:
                    continue
                _bytes = hdu.data.nbytes
                _bytes = _bytes + _pad_length(_bytes)
                if _bytes != hdu._datSpan:
                    self._resize = True
                    self._truncate = False
                    if verbose:
                        print 'One or more data area is resized.'
                    break

            if self._truncate:
               try:
                   self.__file.truncate(hdu._datLoc + hdu._datSpan)
               except IOError:
                   self._resize = True
               self._truncate = False

        return self._resize
