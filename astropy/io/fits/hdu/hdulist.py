# Licensed under a 3-clause BSD style license - see PYFITS.rst

import gzip
import itertools
import os
import re
import shutil
import sys
import warnings

import numpy as np

from astropy.io.fits.file import FILE_MODES, _File
from astropy.io.fits.header import _pad_length
from astropy.io.fits.util import (
    _free_space_check,
    _get_array_mmap,
    _is_int,
    _tmp_name,
    fileobj_closed,
    fileobj_mode,
    ignore_sigint,
    isfile,
)
from astropy.io.fits.verify import VerifyError, VerifyWarning, _ErrList, _Verify
from astropy.utils import indent

# NOTE: Python can be built without bz2.
from astropy.utils.compat.optional_deps import HAS_BZ2
from astropy.utils.exceptions import AstropyUserWarning

from . import compressed
from .base import ExtensionHDU, _BaseHDU, _NonstandardHDU, _ValidHDU
from .groups import GroupsHDU
from .image import ImageHDU, PrimaryHDU

if HAS_BZ2:
    import bz2

__all__ = ["HDUList", "fitsopen"]

# FITS file signature as per RFC 4047
FITS_SIGNATURE = b"SIMPLE  =                    T"


def fitsopen(
    name,
    mode="readonly",
    memmap=None,
    save_backup=False,
    cache=True,
    lazy_load_hdus=None,
    ignore_missing_simple=False,
    *,
    use_fsspec=None,
    fsspec_kwargs=None,
    **kwargs,
):
    """Factory function to open a FITS file and return an `HDUList` object.

    Parameters
    ----------
    name : str, file-like or `pathlib.Path`
        File to be opened.

    mode : str, optional
        Open mode, 'readonly', 'update', 'append', 'denywrite', or
        'ostream'. Default is 'readonly'.

        If ``name`` is a file object that is already opened, ``mode`` must
        match the mode the file was opened with, readonly (rb), update (rb+),
        append (ab+), ostream (w), denywrite (rb)).

    memmap : bool, optional
        Is memory mapping to be used? This value is obtained from the
        configuration item ``astropy.io.fits.Conf.use_memmap``.
        Default is `True`.

    save_backup : bool, optional
        If the file was opened in update or append mode, this ensures that
        a backup of the original file is saved before any changes are flushed.
        The backup has the same name as the original file with ".bak" appended.
        If "file.bak" already exists then "file.bak.1" is used, and so on.
        Default is `False`.

    cache : bool, optional
        If the file name is a URL, `~astropy.utils.data.download_file` is used
        to open the file.  This specifies whether or not to save the file
        locally in Astropy's download cache. Default is `True`.

    lazy_load_hdus : bool, optional
        To avoid reading all the HDUs and headers in a FITS file immediately
        upon opening.  This is an optimization especially useful for large
        files, as FITS has no way of determining the number and offsets of all
        the HDUs in a file without scanning through the file and reading all
        the headers. Default is `True`.

        To disable lazy loading and read all HDUs immediately (the old
        behavior) use ``lazy_load_hdus=False``.  This can lead to fewer
        surprises--for example with lazy loading enabled, ``len(hdul)``
        can be slow, as it means the entire FITS file needs to be read in
        order to determine the number of HDUs.  ``lazy_load_hdus=False``
        ensures that all HDUs have already been loaded after the file has
        been opened.

        .. versionadded:: 1.3

    uint : bool, optional
        Interpret signed integer data where ``BZERO`` is the central value and
        ``BSCALE == 1`` as unsigned integer data.  For example, ``int16`` data
        with ``BZERO = 32768`` and ``BSCALE = 1`` would be treated as
        ``uint16`` data. Default is `True` so that the pseudo-unsigned
        integer convention is assumed.

    ignore_missing_end : bool, optional
        Do not raise an exception when opening a file that is missing an
        ``END`` card in the last header. Default is `False`.

    ignore_missing_simple : bool, optional
        Do not raise an exception when the SIMPLE keyword is missing. Note
        that io.fits will raise a warning if a SIMPLE card is present but
        written in a way that does not follow the FITS Standard.
        Default is `False`.

        .. versionadded:: 4.2

    checksum : bool, str, optional
        If `True`, verifies that both ``DATASUM`` and ``CHECKSUM`` card values
        (when present in the HDU header) match the header and data of all HDU's
        in the file.  Updates to a file that already has a checksum will
        preserve and update the existing checksums unless this argument is
        given a value of 'remove', in which case the CHECKSUM and DATASUM
        values are not checked, and are removed when saving changes to the
        file. Default is `False`.

    disable_image_compression : bool, optional
        If `True`, treats compressed image HDU's like normal binary table
        HDU's.  Default is `False`.

    do_not_scale_image_data : bool, optional
        If `True`, image data is not scaled using BSCALE/BZERO values
        when read.  Default is `False`.

    character_as_bytes : bool, optional
        Whether to return bytes for string columns, otherwise unicode strings
        are returned, but this does not respect memory mapping and loads the
        whole column in memory when accessed. Default is `False`.

    ignore_blank : bool, optional
        If `True`, the BLANK keyword is ignored if present.
        Default is `False`.

    scale_back : bool, optional
        If `True`, when saving changes to a file that contained scaled image
        data, restore the data to the original type and reapply the original
        BSCALE/BZERO values. This could lead to loss of accuracy if scaling
        back to integer values after performing floating point operations on
        the data. Default is `False`.

    output_verify : str
        Output verification option.  Must be one of ``"fix"``,
        ``"silentfix"``, ``"ignore"``, ``"warn"``, or
        ``"exception"``.  May also be any combination of ``"fix"`` or
        ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
        (e.g. ``"fix+warn"``).  See :ref:`astropy:verify` for more info.

    use_fsspec : bool, optional
        Use `fsspec.open` to open the file? Defaults to `False` unless
        ``name`` starts with the Amazon S3 storage prefix ``s3://`` or the
        Google Cloud Storage prefix ``gs://``.  Can also be used for paths
        with other prefixes (e.g., ``http://``) but in this case you must
        explicitly pass ``use_fsspec=True``.
        Use of this feature requires the optional ``fsspec`` package.
        A ``ModuleNotFoundError`` will be raised if the dependency is missing.

        .. versionadded:: 5.2

    fsspec_kwargs : dict, optional
        Keyword arguments passed on to `fsspec.open`. This can be used to
        configure cloud storage credentials and caching behavior.
        For example, pass ``fsspec_kwargs={"anon": True}`` to enable
        anonymous access to Amazon S3 open data buckets.
        See ``fsspec``'s documentation for available parameters.

        .. versionadded:: 5.2

    Returns
    -------
    hdulist : `HDUList`
        `HDUList` containing all of the header data units in the file.

    """
    from astropy.io.fits import conf

    if memmap is None:
        # distinguish between True (kwarg explicitly set)
        # and None (preference for memmap in config, might be ignored)
        memmap = None if conf.use_memmap else False
    else:
        memmap = bool(memmap)

    if lazy_load_hdus is None:
        lazy_load_hdus = conf.lazy_load_hdus
    else:
        lazy_load_hdus = bool(lazy_load_hdus)

    if "uint" not in kwargs:
        kwargs["uint"] = conf.enable_uint

    if not name:
        raise ValueError(f"Empty filename: {name!r}")

    return HDUList.fromfile(
        name,
        mode,
        memmap,
        save_backup,
        cache,
        lazy_load_hdus,
        ignore_missing_simple,
        use_fsspec=use_fsspec,
        fsspec_kwargs=fsspec_kwargs,
        **kwargs,
    )


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
        hdus : BaseHDU or sequence thereof, optional
            The HDU object(s) to comprise the `HDUList`.  Should be
            instances of HDU classes like `ImageHDU` or `BinTableHDU`.

        file : file-like, bytes, optional
            The opened physical file associated with the `HDUList`
            or a bytes object containing the contents of the FITS
            file.
        """
        if isinstance(file, bytes):
            self._data = file
            self._file = None
        else:
            self._file = file
            self._data = None

        # For internal use only--the keyword args passed to fitsopen /
        # HDUList.fromfile/string when opening the file
        self._open_kwargs = {}
        self._in_read_next_hdu = False

        # If we have read all the HDUs from the file or not
        # The assumes that all HDUs have been written when we first opened the
        # file; we do not currently support loading additional HDUs from a file
        # while it is being streamed to.  In the future that might be supported
        # but for now this is only used for the purpose of lazy-loading of
        # existing HDUs.
        if file is None:
            self._read_all = True
        elif self._file is not None:
            # Should never attempt to read HDUs in ostream mode
            self._read_all = self._file.mode == "ostream"
        else:
            self._read_all = False

        if hdus is None:
            hdus = []

        # can take one HDU, as well as a list of HDU's as input
        if isinstance(hdus, _ValidHDU):
            hdus = [hdus]
        elif not isinstance(hdus, (HDUList, list)):
            raise TypeError("Invalid input for HDUList.")

        for idx, hdu in enumerate(hdus):
            if not isinstance(hdu, _BaseHDU):
                raise TypeError(f"Element {idx} in the HDUList input is not an HDU.")

        super().__init__(hdus)

        if file is None:
            # Only do this when initializing from an existing list of HDUs
            # When initializing from a file, this will be handled by the
            # append method after the first HDU is read
            self.update_extend()

    def __len__(self):
        if not self._in_read_next_hdu:
            self.readall()

        return super().__len__()

    def __repr__(self):
        # Special case: if the FITS file is located on a remote file system
        # and has not been fully read yet, we return a simplified repr to
        # avoid downloading the entire file.  We can tell that a file is remote
        # from the fact that the ``fsspec`` package was used to open it.
        is_fsspec_file = self._file and "fsspec" in str(
            self._file._file.__class__.__bases__
        )
        if not self._read_all and is_fsspec_file:
            return f"{type(self)} (partially read)"

        # In order to correctly repr an HDUList we need to load all the
        # HDUs as well
        self.readall()

        return super().__repr__()

    def __iter__(self):
        # While effectively this does the same as:
        # for idx in range(len(self)):
        #     yield self[idx]
        # the more complicated structure is here to prevent the use of len(),
        # which would break the lazy loading
        for idx in itertools.count():
            try:
                yield self[idx]
            except IndexError:
                break

    def __getitem__(self, key):
        """
        Get an HDU from the `HDUList`, indexed by number or name.
        """
        # If the key is a slice we need to make sure the necessary HDUs
        # have been loaded before passing the slice on to super.
        if isinstance(key, slice):
            max_idx = key.stop
            # Check for and handle the case when no maximum was
            # specified (e.g. [1:]).
            if max_idx is None:
                # We need all of the HDUs, so load them
                # and reset the maximum to the actual length.
                max_idx = len(self)

            # Just in case the max_idx is negative...
            max_idx = self._positive_index_of(max_idx)

            number_loaded = super().__len__()

            if max_idx >= number_loaded:
                # We need more than we have, try loading up to and including
                # max_idx. Note we do not try to be clever about skipping HDUs
                # even though key.step might conceivably allow it.
                for i in range(number_loaded, max_idx):
                    # Read until max_idx or to the end of the file, whichever
                    # comes first.
                    if not self._read_next_hdu():
                        break

            try:
                hdus = super().__getitem__(key)
            except IndexError as e:
                # Raise a more helpful IndexError if the file was not fully read.
                if self._read_all:
                    raise e
                else:
                    raise IndexError(
                        "HDU not found, possibly because the index "
                        "is out of range, or because the file was "
                        "closed before all HDUs were read"
                    )
            else:
                return HDUList(hdus)

        # Originally this used recursion, but hypothetically an HDU with
        # a very large number of HDUs could blow the stack, so use a loop
        # instead
        try:
            return self._try_while_unread_hdus(
                super().__getitem__, self._positive_index_of(key)
            )
        except IndexError as e:
            # Raise a more helpful IndexError if the file was not fully read.
            if self._read_all:
                raise e
            else:
                raise IndexError(
                    "HDU not found, possibly because the index "
                    "is out of range, or because the file was "
                    "closed before all HDUs were read"
                )

    def __contains__(self, item):
        """
        Returns `True` if ``item`` is an ``HDU`` _in_ ``self`` or a valid
        extension specification (e.g., integer extension number, extension
        name, or a tuple of extension name and an extension version)
        of a ``HDU`` in ``self``.

        """
        try:
            self._try_while_unread_hdus(self.index_of, item)
        except (KeyError, ValueError):
            return False

        return True

    def __setitem__(self, key, hdu):
        """
        Set an HDU to the `HDUList`, indexed by number or name.
        """
        _key = self._positive_index_of(key)
        if isinstance(hdu, (slice, list)):
            if _is_int(_key):
                raise ValueError("An element in the HDUList must be an HDU.")
            for item in hdu:
                if not isinstance(item, _BaseHDU):
                    raise ValueError(f"{item} is not an HDU.")
        else:
            if not isinstance(hdu, _BaseHDU):
                raise ValueError(f"{hdu} is not an HDU.")

        try:
            self._try_while_unread_hdus(super().__setitem__, _key, hdu)
        except IndexError:
            raise IndexError(f"Extension {key} is out of bound or not found.")

        self._resize = True
        self._truncate = False

    def __delitem__(self, key):
        """
        Delete an HDU from the `HDUList`, indexed by number or name.
        """
        if isinstance(key, slice):
            end_index = len(self)
        else:
            key = self._positive_index_of(key)
            end_index = len(self) - 1

        self._try_while_unread_hdus(super().__delitem__, key)

        if key == end_index or key == -1 and not self._resize:
            self._truncate = True
        else:
            self._truncate = False
            self._resize = True

    # Support the 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        output_verify = self._open_kwargs.get("output_verify", "exception")
        self.close(output_verify=output_verify)

    @classmethod
    def fromfile(
        cls,
        fileobj,
        mode=None,
        memmap=None,
        save_backup=False,
        cache=True,
        lazy_load_hdus=True,
        ignore_missing_simple=False,
        **kwargs,
    ):
        """
        Creates an `HDUList` instance from a file-like object.

        The actual implementation of ``fitsopen()``, and generally shouldn't
        be used directly.  Use :func:`open` instead (and see its
        documentation for details of the parameters accepted by this method).
        """
        return cls._readfrom(
            fileobj=fileobj,
            mode=mode,
            memmap=memmap,
            save_backup=save_backup,
            cache=cache,
            ignore_missing_simple=ignore_missing_simple,
            lazy_load_hdus=lazy_load_hdus,
            **kwargs,
        )

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
        data : str, buffer-like, etc.
            A string or other memory buffer containing an entire FITS file.
            Buffer-like objects include :class:`~bytes`, :class:`~bytearray`,
            :class:`~memoryview`, and :class:`~numpy.ndarray`.
            It should be noted that if that memory is read-only (such as a
            Python string) the returned :class:`HDUList`'s data portions will
            also be read-only.
        **kwargs : dict
            Optional keyword arguments.  See
            :func:`astropy.io.fits.open` for details.

        Returns
        -------
        hdul : HDUList
            An :class:`HDUList` object representing the in-memory FITS file.
        """
        try:
            # Test that the given object supports the buffer interface by
            # ensuring an ndarray can be created from it
            np.ndarray((), dtype="ubyte", buffer=data)
        except TypeError:
            raise TypeError(
                f"The provided object {data} does not contain an underlying "
                "memory buffer.  fromstring() requires an object that "
                "supports the buffer interface such as bytes, buffer, "
                "memoryview, ndarray, etc.  This restriction is to ensure "
                "that efficient access to the array/table data is possible."
            )

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
        if self._file is not None:
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
                        f = info["file"]
                        fm = info["filemode"]
                        break

                output = {
                    "file": f,
                    "filemode": fm,
                    "hdrLoc": None,
                    "datLoc": None,
                    "datSpan": None,
                }

            output["filename"] = self._file.name
            output["resized"] = self._wasresized()
        else:
            output = None

        return output

    def __copy__(self):
        """
        Return a shallow copy of an HDUList.

        Returns
        -------
        copy : `HDUList`
            A shallow copy of this `HDUList` object.

        """
        return self[:]

    # Syntactic sugar for `__copy__()` magic method
    copy = __copy__

    def __deepcopy__(self, memo=None):
        return HDUList([hdu.copy() for hdu in self])

    def pop(self, index=-1):
        """Remove an item from the list and return it.

        Parameters
        ----------
        index : int, str, tuple of (string, int), optional
            An integer value of ``index`` indicates the position from which
            ``pop()`` removes and returns an HDU. A string value or a tuple
            of ``(string, int)`` functions as a key for identifying the
            HDU to be removed and returned. If ``key`` is a tuple, it is
            of the form ``(key, ver)`` where ``ver`` is an ``EXTVER``
            value that must match the HDU being searched for.

            If the key is ambiguous (e.g. there are multiple 'SCI' extensions)
            the first match is returned.  For a more precise match use the
            ``(name, ver)`` pair.

            If even the ``(name, ver)`` pair is ambiguous the numeric index
            must be used to index the duplicate HDU.

        Returns
        -------
        hdu : BaseHDU
            The HDU object at position indicated by ``index`` or having name
            and version specified by ``index``.
        """
        # Make sure that HDUs are loaded before attempting to pop
        self.readall()
        list_index = self.index_of(index)
        return super().pop(list_index)

    def insert(self, index, hdu):
        """
        Insert an HDU into the `HDUList` at the given ``index``.

        Parameters
        ----------
        index : int
            Index before which to insert the new HDU.

        hdu : BaseHDU
            The HDU object to insert
        """
        if not isinstance(hdu, _BaseHDU):
            raise ValueError(f"{hdu} is not an HDU.")

        num_hdus = len(self)

        if index == 0 or num_hdus == 0:
            if num_hdus != 0:
                # We are inserting a new Primary HDU so we need to
                # make the current Primary HDU into an extension HDU.
                if isinstance(self[0], GroupsHDU):
                    raise ValueError(
                        "The current Primary HDU is a GroupsHDU.  "
                        "It can't be made into an extension HDU, "
                        "so another HDU cannot be inserted before it."
                    )

                hdu1 = ImageHDU(self[0].data, self[0].header)

                # Insert it into position 1, then delete HDU at position 0.
                super().insert(1, hdu1)
                super().__delitem__(0)

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

                    super().insert(0, phdu)
                    index = 1
        else:
            if isinstance(hdu, GroupsHDU):
                raise ValueError("A GroupsHDU must be inserted as a Primary HDU.")

            if isinstance(hdu, PrimaryHDU):
                # You passed a Primary HDU but we need an Extension HDU
                # so create an Extension HDU from the input Primary HDU.
                hdu = ImageHDU(hdu.data, hdu.header)

        super().insert(index, hdu)
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
        hdu : BaseHDU
            HDU to add to the `HDUList`.
        """
        if not isinstance(hdu, _BaseHDU):
            raise ValueError("HDUList can only append an HDU.")

        if len(self) > 0:
            if isinstance(hdu, GroupsHDU):
                raise ValueError("Can't append a GroupsHDU to a non-empty HDUList")

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
                    super().append(phdu)

        super().append(hdu)
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
        key : int, str, tuple of (string, int) or BaseHDU
            The key identifying the HDU.  If ``key`` is a tuple, it is of the
            form ``(name, ver)`` where ``ver`` is an ``EXTVER`` value that must
            match the HDU being searched for.

            If the key is ambiguous (e.g. there are multiple 'SCI' extensions)
            the first match is returned.  For a more precise match use the
            ``(name, ver)`` pair.

            If even the ``(name, ver)`` pair is ambiguous (it shouldn't be
            but it's not impossible) the numeric index must be used to index
            the duplicate HDU.

            When ``key`` is an HDU object, this function returns the
            index of that HDU object in the ``HDUList``.

        Returns
        -------
        index : int
            The index of the HDU in the `HDUList`.

        Raises
        ------
        ValueError
            If ``key`` is an HDU object and it is not found in the ``HDUList``.
        KeyError
            If an HDU specified by the ``key`` that is an extension number,
            extension name, or a tuple of extension name and version is not
            found in the ``HDUList``.

        """
        if _is_int(key):
            return key
        elif isinstance(key, tuple):
            _key, _ver = key
        elif isinstance(key, _BaseHDU):
            return self.index(key)
        else:
            _key = key
            _ver = None

        if not isinstance(_key, str):
            raise KeyError(
                "{} indices must be integers, extension names as strings, "
                "or (extname, version) tuples; got {}"
                "".format(self.__class__.__name__, _key)
            )

        _key = (_key.strip()).upper()

        found = None
        for idx, hdu in enumerate(self):
            name = hdu.name
            if isinstance(name, str):
                name = name.strip().upper()
            # 'PRIMARY' should always work as a reference to the first HDU
            if (name == _key or (_key == "PRIMARY" and idx == 0)) and (
                _ver is None or _ver == hdu.ver
            ):
                found = idx
                break

        if found is None:
            raise KeyError(f"Extension {key!r} not found.")
        else:
            return found

    def _positive_index_of(self, key):
        """
        Same as index_of, but ensures always returning a positive index
        or zero.

        (Really this should be called non_negative_index_of but it felt
        too long.)

        This means that if the key is a negative integer, we have to
        convert it to the corresponding positive index.  This means
        knowing the length of the HDUList, which in turn means loading
        all HDUs.  Therefore using negative indices on HDULists is inherently
        inefficient.
        """
        index = self.index_of(key)

        if index >= 0:
            return index

        if abs(index) > len(self):
            raise IndexError(f"Extension {index} is out of bound or not found.")

        return len(self) + index

    def readall(self):
        """
        Read data of all HDUs into memory.
        """
        while self._read_next_hdu():
            pass

    @ignore_sigint
    def flush(self, output_verify="fix", verbose=False):
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
            (e.g. ``"fix+warn"``).  See :ref:`astropy:verify` for more info.

        verbose : bool
            When `True`, print verbose messages
        """
        if self._file.mode not in ("append", "update", "ostream"):
            warnings.warn(
                f"Flush for '{self._file.mode}' mode is not supported.",
                AstropyUserWarning,
            )
            return

        save_backup = self._open_kwargs.get("save_backup", False)
        if save_backup and self._file.mode in ("append", "update"):
            filename = self._file.name
            if os.path.exists(filename):
                # The the file doesn't actually exist anymore for some reason
                # then there's no point in trying to make a backup
                backup = filename + ".bak"
                idx = 1
                while os.path.exists(backup):
                    backup = filename + ".bak." + str(idx)
                    idx += 1
                warnings.warn(
                    f"Saving a backup of {filename} to {backup}.", AstropyUserWarning
                )
                try:
                    shutil.copy(filename, backup)
                except OSError as exc:
                    raise OSError(
                        f"Failed to save backup to destination {filename}"
                    ) from exc

        self.verify(option=output_verify)

        if self._file.mode in ("append", "ostream"):
            for hdu in self:
                if verbose:
                    try:
                        extver = str(hdu._header["extver"])
                    except KeyError:
                        extver = ""

                # only append HDU's which are "new"
                if hdu._new:
                    hdu._prewriteto(checksum=hdu._output_checksum)
                    with _free_space_check(self):
                        hdu._writeto(self._file)
                        if verbose:
                            print("append HDU", hdu.name, extver)
                        hdu._new = False
                    hdu._postwriteto()

        elif self._file.mode == "update":
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

        def get_first_ext():
            try:
                return self[1]
            except IndexError:
                return None

        if "EXTEND" in hdr:
            if not hdr["EXTEND"] and get_first_ext() is not None:
                hdr["EXTEND"] = True
        elif get_first_ext() is not None:
            if hdr["NAXIS"] == 0:
                hdr.set("EXTEND", True, after="NAXIS")
            else:
                n = hdr["NAXIS"]
                hdr.set("EXTEND", True, after="NAXIS" + str(n))

    def writeto(
        self, fileobj, output_verify="exception", overwrite=False, checksum=False
    ):
        """
        Write the `HDUList` to a new file.

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.

        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  May also be any combination of ``"fix"`` or
            ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
            (e.g. ``"fix+warn"``).  See :ref:`astropy:verify` for more info.

        overwrite : bool, optional
            If ``True``, overwrite the output file if it exists. Raises an
            ``OSError`` if ``False`` and the output file exists. Default is
            ``False``.

        checksum : bool
            When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards
            to the headers of all HDU's written to the file.
        """
        if len(self) == 0:
            warnings.warn("There is nothing to write.", AstropyUserWarning)
            return

        self.verify(option=output_verify)

        # make sure the EXTEND keyword is there if there is extension
        self.update_extend()

        # make note of whether the input file object is already open, in which
        # case we should not close it after writing (that should be the job
        # of the caller)
        closed = isinstance(fileobj, str) or fileobj_closed(fileobj)

        mode = FILE_MODES[fileobj_mode(fileobj)] if isfile(fileobj) else "ostream"

        # This can accept an open file object that's open to write only, or in
        # append/update modes but only if the file doesn't exist.
        fileobj = _File(fileobj, mode=mode, overwrite=overwrite)
        hdulist = self.fromfile(fileobj)
        try:
            dirname = os.path.dirname(hdulist._file.name)
        except (AttributeError, TypeError):
            dirname = None

        try:
            with _free_space_check(self, dirname=dirname):
                for hdu in self:
                    hdu._prewriteto(checksum=checksum)
                    hdu._writeto(hdulist._file)
                    hdu._postwriteto()
        finally:
            hdulist.close(output_verify=output_verify, closed=closed)

    def close(self, output_verify="exception", verbose=False, closed=True):
        """
        Close the associated FITS file and memmap object, if any.

        Parameters
        ----------
        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  May also be any combination of ``"fix"`` or
            ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
            (e.g. ``"fix+warn"``).  See :ref:`astropy:verify` for more info.

        verbose : bool
            When `True`, print out verbose messages.

        closed : bool
            When `True`, close the underlying file object.
        """
        try:
            if (
                self._file
                and self._file.mode in ("append", "update")
                and not self._file.closed
            ):
                self.flush(output_verify=output_verify, verbose=verbose)
        finally:
            if self._file and closed and hasattr(self._file, "close"):
                self._file.close()

            # Give individual HDUs an opportunity to do on-close cleanup
            for hdu in self:
                hdu._close(closed=closed)

    def info(self, output=None):
        """
        Summarize the info of the HDUs in this `HDUList`.

        Note that this function prints its results to the console---it
        does not return a value.

        Parameters
        ----------
        output : file-like or bool, optional
            A file-like object to write the output to.  If `False`, does not
            output to a file and instead returns a list of tuples representing
            the HDU info.  Writes to ``sys.stdout`` by default.
        """
        if output is None:
            output = sys.stdout

        if self._file is None:
            name = "(No file associated with this HDUList)"
        else:
            name = self._file.name

        results = [
            f"Filename: {name}",
            "No.    Name      Ver    Type      Cards   Dimensions   Format",
        ]

        format = "{:3d}  {:10}  {:3} {:11}  {:5d}   {}   {}   {}"
        default = ("", "", "", 0, (), "", "")
        for idx, hdu in enumerate(self):
            summary = hdu._summary()
            if len(summary) < len(default):
                summary += default[len(summary) :]
            summary = (idx,) + summary
            if output:
                results.append(format.format(*summary))
            else:
                results.append(summary)

        if output:
            output.write("\n".join(results))
            output.write("\n")
            output.flush()
        else:
            return results[2:]

    def filename(self):
        """
        Return the file name associated with the HDUList object if one exists.
        Otherwise returns None.

        Returns
        -------
        filename : str
            A string containing the file name associated with the HDUList
            object if an association exists.  Otherwise returns None.

        """
        if self._file is not None:
            if hasattr(self._file, "name"):
                return self._file.name
        return None

    @classmethod
    def _readfrom(
        cls,
        fileobj=None,
        data=None,
        mode=None,
        memmap=None,
        cache=True,
        lazy_load_hdus=True,
        ignore_missing_simple=False,
        *,
        use_fsspec=None,
        fsspec_kwargs=None,
        **kwargs,
    ):
        """
        Provides the implementations from HDUList.fromfile and
        HDUList.fromstring, both of which wrap this method, as their
        implementations are largely the same.
        """
        if fileobj is not None:
            if not isinstance(fileobj, _File):
                # instantiate a FITS file object (ffo)
                fileobj = _File(
                    fileobj,
                    mode=mode,
                    memmap=memmap,
                    cache=cache,
                    use_fsspec=use_fsspec,
                    fsspec_kwargs=fsspec_kwargs,
                )
            # The Astropy mode is determined by the _File initializer if the
            # supplied mode was None
            mode = fileobj.mode
            hdulist = cls(file=fileobj)
        else:
            if mode is None:
                # The default mode
                mode = "readonly"

            hdulist = cls(file=data)
            # This method is currently only called from HDUList.fromstring and
            # HDUList.fromfile.  If fileobj is None then this must be the
            # fromstring case; the data type of ``data`` will be checked in the
            # _BaseHDU.fromstring call.

        if (
            not ignore_missing_simple
            and hdulist._file
            and hdulist._file.mode != "ostream"
            and hdulist._file.size > 0
        ):
            pos = hdulist._file.tell()
            # FITS signature is supposed to be in the first 30 bytes, but to
            # allow reading various invalid files we will check in the first
            # card (80 bytes).
            simple = hdulist._file.read(80)
            match_sig = simple[:29] == FITS_SIGNATURE[:-1] and simple[29:30] in (
                b"T",
                b"F",
            )

            if not match_sig:
                # Check the SIMPLE card is there but not written correctly
                match_sig_relaxed = re.match(rb"SIMPLE\s*=\s*[T|F]", simple)

                if match_sig_relaxed:
                    warnings.warn(
                        "Found a SIMPLE card but its format doesn't"
                        " respect the FITS Standard",
                        VerifyWarning,
                    )
                else:
                    if hdulist._file.close_on_error:
                        hdulist._file.close()
                    raise OSError(
                        "No SIMPLE card found, this file does not appear to "
                        "be a valid FITS file. If this is really a FITS file, "
                        "try with ignore_missing_simple=True"
                    )

            hdulist._file.seek(pos)

        # Store additional keyword args that were passed to fits.open
        hdulist._open_kwargs = kwargs

        if fileobj is not None and fileobj.writeonly:
            # Output stream--not interested in reading/parsing
            # the HDUs--just writing to the output file
            return hdulist

        # Make sure at least the PRIMARY HDU can be read
        read_one = hdulist._read_next_hdu()

        # If we're trying to read only and no header units were found,
        # raise an exception
        if not read_one and mode in ("readonly", "denywrite"):
            # Close the file if necessary (issue #6168)
            if hdulist._file.close_on_error:
                hdulist._file.close()

            raise OSError("Empty or corrupt FITS file")

        if not lazy_load_hdus or kwargs.get("checksum") is True:
            # Go ahead and load all HDUs
            while hdulist._read_next_hdu():
                pass

        # initialize/reset attributes to be used in "update/append" mode
        hdulist._resize = False
        hdulist._truncate = False

        return hdulist

    def _try_while_unread_hdus(self, func, *args, **kwargs):
        """
        Attempt an operation that accesses an HDU by index/name
        that can fail if not all HDUs have been read yet.  Keep
        reading HDUs until the operation succeeds or there are no
        more HDUs to read.
        """
        while True:
            try:
                return func(*args, **kwargs)
            except Exception:
                if self._read_next_hdu():
                    continue
                else:
                    raise

    def _read_next_hdu(self):
        """
        Lazily load a single HDU from the fileobj or data string the `HDUList`
        was opened from, unless no further HDUs are found.

        Returns True if a new HDU was loaded, or False otherwise.
        """
        if self._read_all:
            return False

        saved_compression_enabled = compressed.COMPRESSION_ENABLED
        fileobj, data, kwargs = self._file, self._data, self._open_kwargs

        if fileobj is not None and fileobj.closed:
            return False

        try:
            self._in_read_next_hdu = True

            if (
                "disable_image_compression" in kwargs
                and kwargs["disable_image_compression"]
            ):
                compressed.COMPRESSION_ENABLED = False

            # read all HDUs
            try:
                if fileobj is not None:
                    try:
                        # Make sure we're back to the end of the last read
                        # HDU
                        if len(self) > 0:
                            last = self[len(self) - 1]
                            if last._data_offset is not None:
                                offset = last._data_offset + last._data_size
                                fileobj.seek(offset, os.SEEK_SET)

                        hdu = _BaseHDU.readfrom(fileobj, **kwargs)
                    except EOFError:
                        self._read_all = True
                        return False
                    except OSError:
                        # Close the file: see
                        # https://github.com/astropy/astropy/issues/6168
                        #
                        if self._file.close_on_error:
                            self._file.close()

                        if fileobj.writeonly:
                            self._read_all = True
                            return False
                        else:
                            raise
                else:
                    if not data:
                        self._read_all = True
                        return False
                    hdu = _BaseHDU.fromstring(data, **kwargs)
                    self._data = data[hdu._data_offset + hdu._data_size :]

                super().append(hdu)
                if len(self) == 1:
                    # Check for an extension HDU and update the EXTEND
                    # keyword of the primary HDU accordingly
                    self.update_extend()

                hdu._new = False
                if "checksum" in kwargs:
                    hdu._output_checksum = kwargs["checksum"]
            # check in the case there is extra space after the last HDU or
            # corrupted HDU
            except (VerifyError, ValueError) as exc:
                warnings.warn(
                    "Error validating header for HDU #{} (note: Astropy "
                    "uses zero-based indexing).\n{}\n"
                    "There may be extra bytes after the last HDU or the "
                    "file is corrupted.".format(len(self), indent(str(exc))),
                    VerifyWarning,
                )
                del exc
                self._read_all = True
                return False
        finally:
            compressed.COMPRESSION_ENABLED = saved_compression_enabled
            self._in_read_next_hdu = False

        return True

    def _verify(self, option="warn"):
        errs = _ErrList([], unit="HDU")

        # the first (0th) element must be a primary HDU
        if (
            len(self) > 0
            and (not isinstance(self[0], PrimaryHDU))
            and (not isinstance(self[0], _NonstandardHDU))
        ):
            err_text = "HDUList's 0th element is not a primary HDU."
            fix_text = "Fixed by inserting one as 0th HDU."

            def fix(self=self):
                self.insert(0, PrimaryHDU())

            err = self.run_option(option, err_text=err_text, fix_text=fix_text, fix=fix)
            errs.append(err)

        if len(self) > 1 and (
            "EXTEND" not in self[0].header or self[0].header["EXTEND"] is not True
        ):
            err_text = (
                "Primary HDU does not contain an EXTEND keyword "
                "equal to T even though there are extension HDUs."
            )
            fix_text = "Fixed by inserting or updating the EXTEND keyword."

            def fix(header=self[0].header):
                naxis = header["NAXIS"]
                if naxis == 0:
                    after = "NAXIS"
                else:
                    after = "NAXIS" + str(naxis)
                header.set("EXTEND", value=True, after=after)

            errs.append(
                self.run_option(option, err_text=err_text, fix_text=fix_text, fix=fix)
            )

        # each element calls their own verify
        for idx, hdu in enumerate(self):
            if idx > 0 and (not isinstance(hdu, ExtensionHDU)):
                err_text = f"HDUList's element {str(idx)} is not an extension HDU."

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
            if self._resize or self._file.compression:
                self._flush_resize()
            else:
                # if not resized, update in place
                for hdu in self:
                    hdu._writeto(self._file, inplace=True)

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
        old_name = self._file.name
        old_memmap = self._file.memmap
        name = _tmp_name(old_name)

        if not self._file.file_like:
            old_mode = os.stat(old_name).st_mode
            # The underlying file is an actual file object.  The HDUList is
            # resized, so we need to write it to a tmp file, delete the
            # original file, and rename the tmp file to the original file.
            if self._file.compression == "gzip":
                new_file = gzip.GzipFile(name, mode="ab+")
            elif self._file.compression == "bzip2":
                if not HAS_BZ2:
                    raise ModuleNotFoundError(
                        "This Python installation does not provide the bz2 module."
                    )
                new_file = bz2.BZ2File(name, mode="w")
            else:
                new_file = name

            with self.fromfile(new_file, mode="append") as hdulist:
                for hdu in self:
                    hdu._writeto(hdulist._file, inplace=True, copy=True)
                if sys.platform.startswith("win"):
                    # Collect a list of open mmaps to the data; this well be
                    # used later.  See below.
                    mmaps = [
                        (idx, _get_array_mmap(hdu.data), hdu.data)
                        for idx, hdu in enumerate(self)
                        if hdu._has_data
                    ]

                hdulist._file.close()
                self._file.close()
            if sys.platform.startswith("win"):
                # Close all open mmaps to the data.  This is only necessary on
                # Windows, which will not allow a file to be renamed or deleted
                # until all handles to that file have been closed.
                for idx, mmap, arr in mmaps:
                    if mmap is not None:
                        mmap.close()

            os.remove(self._file.name)

            # reopen the renamed new file with "update" mode
            os.rename(name, old_name)
            os.chmod(old_name, old_mode)

            if isinstance(new_file, gzip.GzipFile):
                old_file = gzip.GzipFile(old_name, mode="rb+")
            else:
                old_file = old_name

            ffo = _File(old_file, mode="update", memmap=old_memmap)

            self._file = ffo

            for hdu in self:
                # Need to update the _file attribute and close any open mmaps
                # on each HDU
                if hdu._has_data and _get_array_mmap(hdu.data) is not None:
                    del hdu.data
                hdu._file = ffo

            if sys.platform.startswith("win"):
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
                        # https://github.com/numpy/numpy/issues/8628
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=DeprecationWarning)
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
            ffo = self._file

            ffo.truncate(0)
            ffo.seek(0)

            for hdu in hdulist:
                hdu._writeto(ffo, inplace=True, copy=True)

            # Close the temporary file and delete it.
            hdulist.close()
            os.remove(hdulist._file.name)

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
                        print("One or more header is resized.")
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
                        print("One or more data area is resized.")
                    break

            if self._truncate:
                try:
                    self._file.truncate(hdu._data_offset + hdu._data_size)
                except OSError:
                    self._resize = True
                self._truncate = False

        return self._resize
