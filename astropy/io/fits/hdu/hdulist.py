# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import print_function

import gzip
import os
import shutil
import sys
import warnings

from . import compressed
from .base import _BaseHDU, _ValidHDU, _NonstandardHDU, ExtensionHDU
from .groups import GroupsHDU
from .image import PrimaryHDU, ImageHDU
from ..file import _File
from ..header import _pad_length
from ..util import (_is_int, _tmp_name, fileobj_closed, ignore_sigint,
                    _get_array_mmap)
from ..verify import _Verify, _ErrList, VerifyError, VerifyWarning
from ....extern.six import string_types
from ....utils import indent
from ....utils.exceptions import AstropyUserWarning


def fitsopen(name, mode='readonly', memmap=None, save_backup=False, cached=True, **kwargs):
    """Factory function to open a FITS file and return an `HDUList` object.

    Parameters
    ----------
    name : file path, file object or file-like object
        File to be opened.

    mode : str
        Open mode, 'readonly' (default), 'update', 'append', 'denywrite', or
        'ostream'.

        If ``name`` is a file object that is already opened, ``mode`` must
        match the mode the file was opened with, readonly (rb), update (rb+),
        append (ab+), ostream (w), denywrite (rb)).

    memmap : bool
        Is memory mapping to be used?

    save_backup : bool
        If the file was opened in update or append mode, this ensures that a
        backup of the original file is saved before any changes are flushed.
        The backup has the same name as the original file with ".bak" appended.
        If "file.bak" already exists then "file.bak.1" is used, and so on.

    kwargs : dict
        optional keyword arguments, possible values are:

        - **uint** : bool

            Interpret signed integer data where ``BZERO`` is the
            central value and ``BSCALE == 1`` as unsigned integer
            data.  For example, ``int16`` data with ``BZERO = 32768``
            and ``BSCALE = 1`` would be treated as ``uint16`` data.

            Note, for backward compatibility, the kwarg **uint16** may
            be used instead.  The kwarg was renamed when support was
            added for integers of any size.

        - **ignore_missing_end** : bool

            Do not issue an exception when opening a file that is
            missing an ``END`` card in the last header.

        - **checksum** : bool, str

            If `True`, verifies that both ``DATASUM`` and
            ``CHECKSUM`` card values (when present in the HDU header)
            match the header and data of all HDU's in the file.  Updates to a
            file that already has a checksum will preserve and update the
            existing checksums unless this argument is given a value of
            'remove', in which case the CHECKSUM and DATASUM values are not
            checked, and are removed when saving changes to the file.

        - **disable_image_compression** : bool

            If `True`, treats compressed image HDU's like normal
            binary table HDU's.

        - **do_not_scale_image_data** : bool

            If `True`, image data is not scaled using BSCALE/BZERO values
            when read.

        - **ignore_blank** : bool
           If `True`, the BLANK keyword is ignored if present.

        - **scale_back** : bool

            If `True`, when saving changes to a file that contained scaled
            image data, restore the data to the original type and reapply the
            original BSCALE/BZERO values.  This could lead to loss of accuracy
            if scaling back to integer values after performing floating point
            operations on the data.

    Returns
    -------
        hdulist : an `HDUList` object
            `HDUList` containing all of the header data units in the
            file.

    """

    if memmap is None:
        from .. import conf
        # distinguish between True (kwarg explicitly set)
        # and None (preference for memmap in config, might be ignored)
        memmap = None if conf.use_memmap else False
    else:
        memmap = bool(memmap)

    if 'uint16' in kwargs and 'uint' not in kwargs:
        kwargs['uint'] = kwargs['uint16']
        del kwargs['uint16']

    if not name:
        raise ValueError('Empty filename: %s' % repr(name))

    return HDUList.fromfile(name, mode, memmap, save_backup, cached, **kwargs)


class HDUList(list, _Verify):
    """
    HDU list class.  This is the top-level FITS object.  When a FITS
    file is opened, a `HDUList` object is returned.
    """

    def __init__(self, hdus=[], file=None):
        """
        Construct a `HDUList` object.

        Parameters
        ----------
        hdus : sequence of HDU objects or single HDU, optional
            The HDU object(s) to comprise the `HDUList`.  Should be
            instances of HDU classes like `ImageHDU` or `BinTableHDU`.

        file : file object, optional
            The opened physical file associated with the `HDUList`.
        """

        self.__file = file
        self._save_backup = False

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

        self.update_extend()

    def __iter__(self):
        for idx in range(len(self)):
            yield self[idx]

    def __getitem__(self, key):
        """
        Get an HDU from the `HDUList`, indexed by number or name.
        """

        if isinstance(key, slice):
            hdus = super(HDUList, self).__getitem__(key)
            return HDUList(hdus)

        idx = self.index_of(key)
        return super(HDUList, self).__getitem__(idx)

    def __contains__(self, item):
        """
        Returns `True` if ``HDUList.index_of(item)`` succeeds.
        """
        try:
            self.index_of(item)
            return True
        except KeyError:
            return False

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
    def fromfile(cls, fileobj, mode=None, memmap=None,
                 save_backup=False, cached=True, **kwargs):
        """
        Creates an `HDUList` instance from a file-like object.

        The actual implementation of ``fitsopen()``, and generally shouldn't
        be used directly.  Use :func:`open` instead (and see its
        documentation for details of the parameters accepted by this method).
        """

        return cls._readfrom(fileobj=fileobj, mode=mode, memmap=memmap,
                             save_backup=save_backup, cached=cached, **kwargs)

    @classmethod
    def fromstring(cls, data, **kwargs):
        """
        Creates an `HDUList` instance from a string or other in-memory data
        buffer containing an entire FITS file.  Similar to
        :meth:`HDUList.fromfile`, but does not accept the mode or memmap
        arguments, as they are only relevant to reading from a file on disk.

        This is useful for interfacing with other libraries such as CFITSIO,
        and may also be useful for streaming applications.

        Parameters
        ----------
        data : str, buffer, memoryview, etc.
            A string or other memory buffer containing an entire FITS file.  It
            should be noted that if that memory is read-only (such as a Python
            string) the returned :class:`HDUList`'s data portions will also be
            read-only.

        kwargs : dict
            Optional keyword arguments.  See
            :func:`astropy.io.fits.open` for details.

        Returns
        -------
        hdul : HDUList
            An :class:`HDUList` object representing the in-memory FITS file.
        """

        return cls._readfrom(data=data, **kwargs)

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
        fileinfo : dict or None

            The dictionary details information about the locations of
            the indexed HDU within an associated file.  Returns `None`
            when the HDU is not associated with a file.

            Dictionary contents:

            ========== ========================================================
            Key        Value
            ========== ========================================================
            file       File object associated with the HDU
            filename   Name of associated file object
            filemode   Mode in which the file was opened (readonly,
                       update, append, denywrite, ostream)
            resized    Flag that when `True` indicates that the data has been
                       resized since the last read/write so the returned values
                       may not be valid.
            hdrLoc     Starting byte location of header in file
            datLoc     Starting byte location of data block in file
            datSpan    Data size including padding
            ========== ========================================================

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

    def insert(self, index, hdu):
        """
        Insert an HDU into the `HDUList` at the given ``index``.

        Parameters
        ----------
        index : int
            Index before which to insert the new HDU.

        hdu : HDU object
            The HDU object to insert
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
                        "so another HDU cannot be inserted before it.")

                hdu1 = ImageHDU(self[0].data, self[0].header)

                # Insert it into position 1, then delete HDU at position 0.
                super(HDUList, self).insert(1, hdu1)
                super(HDUList, self).__delitem__(0)

            if not isinstance(hdu, (PrimaryHDU, _NonstandardHDU)):
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
        self.update_extend()

    def append(self, hdu):
        """
        Append a new HDU to the `HDUList`.

        Parameters
        ----------
        hdu : HDU object
            HDU to add to the `HDUList`.
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
                # _header_offset and friends need to be copied too.
                hdu = ImageHDU(hdu.data, hdu.header)
        else:
            if not isinstance(hdu, (PrimaryHDU, _NonstandardHDU)):
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
        self.update_extend()

    def index_of(self, key):
        """
        Get the index of an HDU from the `HDUList`.

        Parameters
        ----------
        key : int, str or tuple of (string, int)
           The key identifying the HDU.  If ``key`` is a tuple, it is of the
           form ``(key, ver)`` where ``ver`` is an ``EXTVER`` value that must
           match the HDU being searched for.

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

        if not isinstance(_key, string_types):
            raise KeyError(key)
        _key = (_key.strip()).upper()

        nfound = 0
        found = None
        for idx, hdu in enumerate(self):
            name = hdu.name
            if isinstance(name, string_types):
                name = name.strip().upper()
            # 'PRIMARY' should always work as a reference to the first HDU
            if ((name == _key or (_key == 'PRIMARY' and idx == 0)) and
                (_ver is None or _ver == hdu.ver)):
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

    @ignore_sigint
    def flush(self, output_verify='fix', verbose=False):
        """
        Force a write of the `HDUList` back to the file (for append and
        update modes only).

        Parameters
        ----------
        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  May also be any combination of ``"fix"`` or
            ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
            (e.g. ``"fix+warn"``).  See :ref:`verify` for more info.

        verbose : bool
            When `True`, print verbose messages
        """

        if self.__file.mode not in ('append', 'update', 'ostream'):
            warnings.warn("Flush for '%s' mode is not supported."
                          % self.__file.mode, AstropyUserWarning)
            return

        if self._save_backup and self.__file.mode in ('append', 'update'):
            filename = self.__file.name
            if os.path.exists(filename):
                # The the file doesn't actually exist anymore for some reason
                # then there's no point in trying to make a backup
                backup = filename + '.bak'
                idx = 1
                while os.path.exists(backup):
                    backup = filename + '.bak.' + str(idx)
                    idx += 1
                warnings.warn('Saving a backup of %s to %s.' %
                              (filename, backup), AstropyUserWarning)
                try:
                    shutil.copy(filename, backup)
                except IOError as exc:
                    raise IOError('Failed to save backup to destination %s: '
                                  '%s' % (filename, exc))

        self.verify(option=output_verify)

        if self.__file.mode in ('append', 'ostream'):
            for hdu in self:
                if verbose:
                    try:
                        extver = str(hdu._header['extver'])
                    except KeyError:
                        extver = ''

                # only append HDU's which are "new"
                if hdu._new:
                    hdu._prewriteto(checksum=hdu._output_checksum)
                    try:
                        hdu._writeto(self.__file)
                        if verbose:
                            print('append HDU', hdu.name, extver)
                        hdu._new = False
                    finally:
                        hdu._postwriteto()
        elif self.__file.mode == 'update':
            self._flush_update()

    def update_extend(self):
        """
        Make sure that if the primary header needs the keyword ``EXTEND`` that
        it has it and it is correct.
        """

        if not len(self):
            return

        if not isinstance(self[0], PrimaryHDU):
            # A PrimaryHDU will be automatically inserted at some point, but it
            # might not have been added yet
            return

        hdr = self[0].header

        if 'EXTEND' in hdr:
            if len(self) > 1 and hdr['EXTEND'] == False:
                hdr['EXTEND'] = True
        elif len(self) > 1:
            if hdr['NAXIS'] == 0:
                hdr.set('EXTEND', True, after='NAXIS')
            else:
                n = hdr['NAXIS']
                hdr.set('EXTEND', True, after='NAXIS' + str(n))

    def writeto(self, fileobj, output_verify='exception', clobber=False,
                checksum=False):
        """
        Write the `HDUList` to a new file.

        Parameters
        ----------
        fileobj : file path, file object or file-like object
            File to write to.  If a file object, must be opened in a
            writeable mode.

        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  May also be any combination of ``"fix"`` or
            ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
            (e.g. ``"fix+warn"``).  See :ref:`verify` for more info.

        clobber : bool
            When `True`, overwrite the output file if exists.

        checksum : bool
            When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards
            to the headers of all HDU's written to the file.
        """

        if (len(self) == 0):
            warnings.warn("There is nothing to write.", AstropyUserWarning)
            return

        self.verify(option=output_verify)

        # make sure the EXTEND keyword is there if there is extension
        self.update_extend()

        # make note of whether the input file object is already open, in which
        # case we should not close it after writing (that should be the job
        # of the caller)
        closed = isinstance(fileobj, string_types) or fileobj_closed(fileobj)

        # writeto is only for writing a new file from scratch, so the most
        # sensible mode to require is 'ostream'.  This can accept an open
        # file object that's open to write only, or in append/update modes
        # but only if the file doesn't exist.
        fileobj = _File(fileobj, mode='ostream', clobber=clobber)
        hdulist = self.fromfile(fileobj)

        for hdu in self:
            hdu._prewriteto(checksum=checksum)
            try:
                hdu._writeto(hdulist.__file)
            finally:
                hdu._postwriteto()

        hdulist.close(output_verify=output_verify, closed=closed)

    def close(self, output_verify='exception', verbose=False, closed=True):
        """
        Close the associated FITS file and memmap object, if any.

        Parameters
        ----------
        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  May also be any combination of ``"fix"`` or
            ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
            (e.g. ``"fix+warn"``).  See :ref:`verify` for more info.

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
        output : file, bool, optional
            A file-like object to write the output to.  If `False`, does not
            output to a file and instead returns a list of tuples representing
            the HDU info.  Writes to ``sys.stdout`` by default.
        """

        if output is None:
            output = sys.stdout

        if self.__file is None:
            name = '(No file associated with this HDUList)'
        else:
            name = self.__file.name

        results = ['Filename: %s' % name,
                   'No.    Name         Type      Cards   Dimensions   Format']

        format = '%-3d  %-10s  %-11s  %5d   %-10s   %s   %s'
        default = ('', '', 0, (), '', '')
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

    @classmethod
    def _readfrom(cls, fileobj=None, data=None, mode=None,
                  memmap=None, save_backup=False, cached=True, **kwargs):
        """
        Provides the implementations from HDUList.fromfile and
        HDUList.fromstring, both of which wrap this method, as their
        implementations are largely the same.
        """

        if fileobj is not None:
            if not isinstance(fileobj, _File):
                # instantiate a FITS file object (ffo)
                ffo = _File(fileobj, mode=mode, memmap=memmap, cached=cached)
            else:
                ffo = fileobj
            # The pyfits mode is determined by the _File initializer if the
            # supplied mode was None
            mode = ffo.mode
            hdulist = cls(file=ffo)
        else:
            if mode is None:
                # The default mode
                mode = 'readonly'

            hdulist = cls()
            # This method is currently only called from HDUList.fromstring and
            # HDUList.fromfile.  If fileobj is None then this must be the
            # fromstring case; the data type of ``data`` will be checked in the
            # _BaseHDU.fromstring call.

        hdulist._save_backup = save_backup

        saved_compression_enabled = compressed.COMPRESSION_ENABLED

        try:
            if ('disable_image_compression' in kwargs and
                kwargs['disable_image_compression']):
                compressed.COMPRESSION_ENABLED = False

            # read all HDUs
            while True:
                try:
                    if fileobj is not None:
                        if ffo.writeonly:
                            # Output stream--not interested in reading/parsing
                            # the HDUs--just writing to the output file
                            return hdulist

                        try:
                            hdu = _BaseHDU.readfrom(ffo, **kwargs)
                        except EOFError:
                            break
                        except IOError:
                            if ffo.writeonly:
                                break
                            else:
                                raise
                    else:
                        if not data:
                            break
                        hdu = _BaseHDU.fromstring(data)
                        data = data[hdu._data_offset + hdu._data_size:]
                    hdulist.append(hdu)
                    hdu._new = False
                    if 'checksum' in kwargs:
                        hdu._output_checksum = kwargs['checksum']
                # check in the case there is extra space after the last HDU or
                # corrupted HDU
                except (VerifyError, ValueError) as exc:
                    warnings.warn(
                        'Error validating header for HDU #%d (note: Astropy '
                        'uses zero-based indexing).\n%s\n'
                        'There may be extra bytes after the last HDU or the '
                        'file is corrupted.' %
                        (len(hdulist), indent(str(exc))), VerifyWarning)
                    del exc
                    break

            # If we're trying to read only and no header units were found,
            # raise and exception
            if mode in ('readonly', 'denywrite') and len(hdulist) == 0:
                raise IOError('Empty or corrupt FITS file')

            # initialize/reset attributes to be used in "update/append" mode
            hdulist._resize = False
            hdulist._truncate = False

        finally:
            compressed.COMPRESSION_ENABLED = saved_compression_enabled

        return hdulist

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
                header.set('EXTEND', value=True, after=after)

            errs.append(self.run_option(option, err_text=err_text,
                                        fix_text=fix_text, fix=fix))

        # each element calls their own verify
        for idx, hdu in enumerate(self):
            if idx > 0 and (not isinstance(hdu, ExtensionHDU)):
                err_text = ("HDUList's element %s is not an extension HDU." %
                            str(idx))
                err = self.run_option(option, err_text=err_text, fixable=False)
                errs.append(err)

            else:
                result = hdu._verify(option)
                if result:
                    errs.append(result)
        return errs

    def _flush_update(self):
        """Implements flushing changes to a file in update mode."""

        for hdu in self:
            # Need to all _prewriteto() for each HDU first to determine if
            # resizing will be necessary
            hdu._prewriteto(checksum=hdu._output_checksum, inplace=True)

        try:
            self._wasresized()

            # if the HDUList is resized, need to write out the entire contents of
            # the hdulist to the file.
            if self._resize or self.__file.compression:
                self._flush_resize()
            else:
                # if not resized, update in place
                for hdu in self:
                    hdu._writeto(self.__file, inplace=True)

            # reset the modification attributes after updating
            for hdu in self:
                hdu._header._modified = False
        finally:
            for hdu in self:
                hdu._postwriteto()

    def _flush_resize(self):
        """
        Implements flushing changes in update mode when parts of one or more HDU
        need to be resized.
        """

        old_name = self.__file.name
        old_memmap = self.__file.memmap
        name = _tmp_name(old_name)

        if not self.__file.file_like:
            old_mode = os.stat(old_name).st_mode
            # The underlying file is an actual file object.  The HDUList is
            # resized, so we need to write it to a tmp file, delete the
            # original file, and rename the tmp file to the original file.
            if self.__file.compression == 'gzip':
                new_file = gzip.GzipFile(name, mode='ab+')
            else:
                new_file = name

            hdulist = self.fromfile(new_file, mode='append')

            for hdu in self:
                hdu._writeto(hdulist.__file, inplace=True, copy=True)

            if sys.platform.startswith('win'):
                # Collect a list of open mmaps to the data; this well be used
                # later.  See below.
                mmaps = [(idx, _get_array_mmap(hdu.data), hdu.data)
                         for idx, hdu in enumerate(self) if hdu._has_data]

            hdulist.__file.close()
            self.__file.close()

            if sys.platform.startswith('win'):
                # Close all open mmaps to the data.  This is only necessary on
                # Windows, which will not allow a file to be renamed or deleted
                # until all handles to that file have been closed.
                for idx, mmap, arr in mmaps:
                    if mmap is not None:
                        mmap.close()

            os.remove(self.__file.name)

            # reopen the renamed new file with "update" mode
            os.rename(name, old_name)
            os.chmod(old_name, old_mode)

            if isinstance(new_file, gzip.GzipFile):
                old_file = gzip.GzipFile(old_name, mode='rb+')
            else:
                old_file = old_name

            ffo = _File(old_file, mode='update', memmap=old_memmap)

            self.__file = ffo

            for hdu in self:
                # Need to update the _file attribute and close any open mmaps
                # on each HDU
                if hdu._has_data and _get_array_mmap(hdu.data) is not None:
                    del hdu.data
                hdu._file = ffo

            if sys.platform.startswith('win'):
                # On Windows, all the original data mmaps were closed above.
                # However, it's possible that the user still has references to
                # the old data which would no longer work (possibly even cause
                # a segfault if they try to access it).  This replaces the
                # buffers used by the original arrays with the buffers of mmap
                # arrays created from the new file.  This seems to work, but
                # it's a flaming hack and carries no guarantees that it won't
                # lead to odd behavior in practice.  Better to just not keep
                # references to data from files that had to be resized upon
                # flushing (on Windows--again, this is no problem on Linux).
                for idx, mmap, arr in mmaps:
                    if mmap is not None:
                        arr.data = self[idx].data.data
                del mmaps  # Just to be sure

        else:
            # The underlying file is not a file object, it is a file like
            # object.  We can't write out to a file, we must update the file
            # like object in place.  To do this, we write out to a temporary
            # file, then delete the contents in our file like object, then
            # write the contents of the temporary file to the now empty file
            # like object.
            self.writeto(name)
            hdulist = self.fromfile(name)
            ffo = self.__file

            ffo.truncate(0)
            ffo.seek(0)

            for hdu in hdulist:
                hdu._writeto(ffo, inplace=True, copy=True)

            # Close the temporary file and delete it.
            hdulist.close()
            os.remove(hdulist.__file.name)

        # reset the resize attributes after updating
        self._resize = False
        self._truncate = False
        for hdu in self:
            hdu._header._modified = False
            hdu._new = False
            hdu._file = ffo

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
                nbytes = len(str(hdu._header))
                if nbytes != (hdu._data_offset - hdu._header_offset):
                    self._resize = True
                    self._truncate = False
                    if verbose:
                        print('One or more header is resized.')
                    break

                # Data:
                if not hdu._has_data:
                    continue

                nbytes = hdu.size
                nbytes = nbytes + _pad_length(nbytes)
                if nbytes != hdu._data_size:
                    self._resize = True
                    self._truncate = False
                    if verbose:
                        print('One or more data area is resized.')
                    break

            if self._truncate:
                try:
                    self.__file.truncate(hdu._data_offset + hdu._data_size)
                except IOError:
                    self._resize = True
                self._truncate = False

        return self._resize
