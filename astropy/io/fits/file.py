# Licensed under a 3-clause BSD style license - see PYFITS.rst


import bz2
import gzip
import errno
import http.client
import mmap
import operator
import pathlib
import io
import os
import sys
import tempfile
import warnings
import zipfile
import re

from functools import reduce

import numpy as np

from .util import (isreadable, iswritable, isfile, fileobj_open, fileobj_name,
                   fileobj_closed, fileobj_mode, _array_from_file,
                   _array_to_file, _write_string)
from ...utils.data import download_file, _is_url
from ...utils.decorators import classproperty, deprecated_renamed_argument
from ...utils.exceptions import AstropyUserWarning


# Maps astropy.io.fits-specific file mode names to the appropriate file
# modes to use for the underlying raw files
IO_FITS_MODES = {
    'readonly': 'rb',
    'copyonwrite': 'rb',
    'update': 'rb+',
    'append': 'ab+',
    'ostream': 'wb',
    'denywrite': 'rb'}

# Maps OS-level file modes to the appropriate astropy.io.fits specific mode
# to use when given file objects but no mode specified; obviously in
# IO_FITS_MODES there are overlaps; for example 'readonly' and 'denywrite'
# both require the file to be opened in 'rb' mode.  But 'readonly' is the
# default behavior for such files if not otherwise specified.
# Note: 'ab' is only supported for 'ostream' which is output-only.
FILE_MODES = {
    'rb': 'readonly', 'rb+': 'update',
    'wb': 'ostream', 'wb+': 'update',
    'ab': 'ostream', 'ab+': 'append'}

# A match indicates the file was opened in text mode, which is not allowed
TEXT_RE = re.compile(r'^[rwa]((t?\+?)|(\+?t?))$')


# readonly actually uses copyonwrite for mmap so that readonly without mmap and
# with mmap still have to same behavior with regard to updating the array.  To
# get a truly readonly mmap use denywrite
# the name 'denywrite' comes from a deprecated flag to mmap() on Linux--it
# should be clarified that 'denywrite' mode is not directly analogous to the
# use of that flag; it was just taken, for lack of anything better, as a name
# that means something like "read only" but isn't readonly.
MEMMAP_MODES = {'readonly': mmap.ACCESS_COPY,
                'copyonwrite': mmap.ACCESS_COPY,
                'update': mmap.ACCESS_WRITE,
                'append': mmap.ACCESS_COPY,
                'denywrite': mmap.ACCESS_READ}

# TODO: Eventually raise a warning, and maybe even later disable the use of
# 'copyonwrite' and 'denywrite' modes unless memmap=True.  For now, however,
# that would generate too many warnings for too many users.  If nothing else,
# wait until the new logging system is in place.

GZIP_MAGIC = b'\x1f\x8b\x08'
PKZIP_MAGIC = b'\x50\x4b\x03\x04'
BZIP2_MAGIC = b'\x42\x5a'

def _normalize_fits_mode(mode):
    if mode is not None and mode not in IO_FITS_MODES:
        if TEXT_RE.match(mode):
            raise ValueError(
                "Text mode '{}' not supported: "
                "files must be opened in binary mode".format(mode))
        new_mode = FILE_MODES.get(mode)
        if new_mode not in IO_FITS_MODES:
            raise ValueError("Mode '{}' not recognized".format(mode))
        mode = new_mode
    return mode

class _File:
    """
    Represents a FITS file on disk (or in some other file-like object).
    """

    @deprecated_renamed_argument('clobber', 'overwrite', '2.0')
    def __init__(self, fileobj=None, mode=None, memmap=None, overwrite=False,
                 cache=True):
        self.strict_memmap = bool(memmap)
        memmap = True if memmap is None else memmap

        if fileobj is None:
            self._file = None
            self.closed = False
            self.binary = True
            self.mode = mode
            self.memmap = memmap
            self.compression = None
            self.readonly = False
            self.writeonly = False
            self.simulateonly = True
            self.close_on_error = False
            return
        else:
            self.simulateonly = False
            # If fileobj is of type pathlib.Path
            if isinstance(fileobj, pathlib.Path):
                fileobj = str(fileobj)
            elif isinstance(fileobj, bytes):
                # Using bytes as filename is tricky, it's deprecated for Windows
                # in Python 3.5 (because it could lead to false-positives) but
                # was fixed and un-deprecated in Python 3.6.
                # However it requires that the bytes object is encoded with the
                # file system encoding.
                # Probably better to error out and ask for a str object instead.
                # TODO: This could be revised when Python 3.5 support is dropped
                # See also: https://github.com/astropy/astropy/issues/6789
                raise TypeError("names should be `str` not `bytes`.")

        # Holds mmap instance for files that use mmap
        self._mmap = None

        if mode is not None and mode not in IO_FITS_MODES:
            raise ValueError("Mode '{}' not recognized".format(mode))
        if isfile(fileobj):
            objmode = _normalize_fits_mode(fileobj_mode(fileobj))
            if mode is not None and mode != objmode:
                raise ValueError(
                    "Requested FITS mode '{}' not compatible with open file "
                    "handle mode '{}'".format(mode, objmode))
            mode = objmode
        if mode is None:
            mode = 'readonly'

        # Handle raw URLs
        if (isinstance(fileobj, str) and
            mode not in ('ostream', 'append', 'update') and _is_url(fileobj)):
            self.name = download_file(fileobj, cache=cache)
        # Handle responses from URL requests that have already been opened
        elif isinstance(fileobj, http.client.HTTPResponse):
            if mode in ('ostream', 'append', 'update'):
                raise ValueError(
                    "Mode {} not supported for HTTPResponse".format(mode))
            fileobj = io.BytesIO(fileobj.read())
        else:
            self.name = fileobj_name(fileobj)

        self.closed = False
        self.binary = True
        self.mode = mode
        self.memmap = memmap

        # Underlying fileobj is a file-like object, but an actual file object
        self.file_like = False

        # Should the object be closed on error: see
        # https://github.com/astropy/astropy/issues/6168
        self.close_on_error = False

        # More defaults to be adjusted below as necessary
        self.compression = None
        self.readonly = False
        self.writeonly = False

        # Initialize the internal self._file object
        if isfile(fileobj):
            self._open_fileobj(fileobj, mode, overwrite)
        elif isinstance(fileobj, str):
            self._open_filename(fileobj, mode, overwrite)
        else:
            self._open_filelike(fileobj, mode, overwrite)

        self.fileobj_mode = fileobj_mode(self._file)

        if isinstance(fileobj, gzip.GzipFile):
            self.compression = 'gzip'
        elif isinstance(fileobj, zipfile.ZipFile):
            # Reading from zip files is supported but not writing (yet)
            self.compression = 'zip'
        elif isinstance(fileobj, bz2.BZ2File):
            self.compression = 'bzip2'

        if (mode in ('readonly', 'copyonwrite', 'denywrite') or
                (self.compression and mode == 'update')):
            self.readonly = True
        elif (mode == 'ostream' or
                (self.compression and mode == 'append')):
            self.writeonly = True

        # For 'ab+' mode, the pointer is at the end after the open in
        # Linux, but is at the beginning in Solaris.
        if (mode == 'ostream' or self.compression or
            not hasattr(self._file, 'seek')):
            # For output stream start with a truncated file.
            # For compressed files we can't really guess at the size
            self.size = 0
        else:
            pos = self._file.tell()
            self._file.seek(0, 2)
            self.size = self._file.tell()
            self._file.seek(pos)

        if self.memmap:
            if not isfile(self._file):
                self.memmap = False
            elif not self.readonly and not self._mmap_available:
                # Test mmap.flush--see
                # https://github.com/astropy/astropy/issues/968
                self.memmap = False

    def __repr__(self):
        return '<{}.{} {}>'.format(self.__module__, self.__class__.__name__,
                                   self._file)

    # Support the 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def readable(self):
        if self.writeonly:
            return False
        return isreadable(self._file)

    def read(self, size=None):
        if not hasattr(self._file, 'read'):
            raise EOFError
        try:
            return self._file.read(size)
        except OSError:
            # On some versions of Python, it appears, GzipFile will raise an
            # OSError if you try to read past its end (as opposed to just
            # returning '')
            if self.compression == 'gzip':
                return ''
            raise

    def readarray(self, size=None, offset=0, dtype=np.uint8, shape=None):
        """
        Similar to file.read(), but returns the contents of the underlying
        file as a numpy array (or mmap'd array if memmap=True) rather than a
        string.

        Usually it's best not to use the `size` argument with this method, but
        it's provided for compatibility.
        """

        if not hasattr(self._file, 'read'):
            raise EOFError

        if not isinstance(dtype, np.dtype):
            dtype = np.dtype(dtype)

        if size and size % dtype.itemsize != 0:
            raise ValueError('size {} not a multiple of {}'.format(size, dtype))

        if isinstance(shape, int):
            shape = (shape,)

        if not (size or shape):
            warnings.warn('No size or shape given to readarray(); assuming a '
                          'shape of (1,)', AstropyUserWarning)
            shape = (1,)

        if size and not shape:
            shape = (size // dtype.itemsize,)

        if size and shape:
            actualsize = np.prod(shape) * dtype.itemsize

            if actualsize > size:
                raise ValueError('size {} is too few bytes for a {} array of '
                                 '{}'.format(size, shape, dtype))
            elif actualsize < size:
                raise ValueError('size {} is too many bytes for a {} array of '
                                 '{}'.format(size, shape, dtype))

        filepos = self._file.tell()

        try:
            if self.memmap:
                if self._mmap is None:
                    # Instantiate Memmap array of the file offset at 0 (so we
                    # can return slices of it to offset anywhere else into the
                    # file)
                    access_mode = MEMMAP_MODES[self.mode]

                    # For reasons unknown the file needs to point to (near)
                    # the beginning or end of the file. No idea how close to
                    # the beginning or end.
                    # If I had to guess there is some bug in the mmap module
                    # of CPython or perhaps in microsoft's underlying code
                    # for generating the mmap.
                    self._file.seek(0, 0)
                    # This would also work:
                    # self._file.seek(0, 2)   # moves to the end
                    try:
                        self._mmap = mmap.mmap(self._file.fileno(), 0,
                                               access=access_mode,
                                               offset=0)
                    except OSError as exc:
                        # NOTE: mode='readonly' results in the memory-mapping
                        # using the ACCESS_COPY mode in mmap so that users can
                        # modify arrays. However, on some systems, the OS raises
                        # a '[Errno 12] Cannot allocate memory' OSError if the
                        # address space is smaller than the file. The solution
                        # is to open the file in mode='denywrite', which at
                        # least allows the file to be opened even if the
                        # resulting arrays will be truly read-only.
                        if exc.errno == errno.ENOMEM and self.mode == 'readonly':
                            warnings.warn("Could not memory map array with "
                                          "mode='readonly', falling back to "
                                          "mode='denywrite', which means that "
                                          "the array will be read-only",
                                          AstropyUserWarning)
                            self._mmap = mmap.mmap(self._file.fileno(), 0,
                                                   access=MEMMAP_MODES['denywrite'],
                                                   offset=0)
                        else:
                            raise

                return np.ndarray(shape=shape, dtype=dtype, offset=offset,
                                  buffer=self._mmap)
            else:
                count = reduce(operator.mul, shape)
                self._file.seek(offset)
                data = _array_from_file(self._file, dtype, count)
                data.shape = shape
                return data
        finally:
            # Make sure we leave the file in the position we found it; on
            # some platforms (e.g. Windows) mmaping a file handle can also
            # reset its file pointer
            self._file.seek(filepos)

    def writable(self):
        if self.readonly:
            return False
        return iswritable(self._file)

    def write(self, string):
        if hasattr(self._file, 'write'):
            _write_string(self._file, string)

    def writearray(self, array):
        """
        Similar to file.write(), but writes a numpy array instead of a string.

        Also like file.write(), a flush() or close() may be needed before
        the file on disk reflects the data written.
        """

        if hasattr(self._file, 'write'):
            _array_to_file(array, self._file)

    def flush(self):
        if hasattr(self._file, 'flush'):
            self._file.flush()

    def seek(self, offset, whence=0):
        if not hasattr(self._file, 'seek'):
            return
        self._file.seek(offset, whence)
        pos = self._file.tell()
        if self.size and pos > self.size:
            warnings.warn('File may have been truncated: actual file length '
                          '({}) is smaller than the expected size ({})'
                          .format(self.size, pos), AstropyUserWarning)

    def tell(self):
        if not hasattr(self._file, 'tell'):
            raise EOFError
        return self._file.tell()

    def truncate(self, size=None):
        if hasattr(self._file, 'truncate'):
            self._file.truncate(size)

    def close(self):
        """
        Close the 'physical' FITS file.
        """

        if hasattr(self._file, 'close'):
            self._file.close()

        self._maybe_close_mmap()
        # Set self._memmap to None anyways since no new .data attributes can be
        # loaded after the file is closed
        self._mmap = None

        self.closed = True
        self.close_on_error = False

    def _maybe_close_mmap(self, refcount_delta=0):
        """
        When mmap is in use these objects hold a reference to the mmap of the
        file (so there is only one, shared by all HDUs that reference this
        file).

        This will close the mmap if there are no arrays referencing it.
        """

        if (self._mmap is not None and
                sys.getrefcount(self._mmap) == 2 + refcount_delta):
            self._mmap.close()
            self._mmap = None

    def _overwrite_existing(self, overwrite, fileobj, closed):
        """Overwrite an existing file if ``overwrite`` is ``True``, otherwise
        raise an OSError.  The exact behavior of this method depends on the
        _File object state and is only meant for use within the ``_open_*``
        internal methods.
        """

        # The file will be overwritten...
        if ((self.file_like and hasattr(fileobj, 'len') and fileobj.len > 0) or
            (os.path.exists(self.name) and os.path.getsize(self.name) != 0)):
            if overwrite:
                if self.file_like and hasattr(fileobj, 'truncate'):
                    fileobj.truncate(0)
                else:
                    if not closed:
                        fileobj.close()
                    os.remove(self.name)
            else:
                raise OSError("File {!r} already exists.".format(self.name))

    def _try_read_compressed(self, obj_or_name, magic, mode, ext=''):
        """Attempt to determine if the given file is compressed"""
        if ext == '.gz' or magic.startswith(GZIP_MAGIC):
            if mode == 'append':
                raise OSError("'append' mode is not supported with gzip files."
                              "Use 'update' mode instead")
            # Handle gzip files
            kwargs = dict(mode=IO_FITS_MODES[mode])
            if isinstance(obj_or_name, str):
                kwargs['filename'] = obj_or_name
            else:
                kwargs['fileobj'] = obj_or_name
            self._file = gzip.GzipFile(**kwargs)
            self.compression = 'gzip'
        elif ext == '.zip' or magic.startswith(PKZIP_MAGIC):
            # Handle zip files
            self._open_zipfile(self.name, mode)
            self.compression = 'zip'
        elif ext == '.bz2' or magic.startswith(BZIP2_MAGIC):
            # Handle bzip2 files
            if mode in ['update', 'append']:
                raise OSError("update and append modes are not supported "
                              "with bzip2 files")
            # bzip2 only supports 'w' and 'r' modes
            bzip2_mode = 'w' if mode == 'ostream' else 'r'
            self._file = bz2.BZ2File(obj_or_name, mode=bzip2_mode)
            self.compression = 'bzip2'
        return self.compression is not None

    def _open_fileobj(self, fileobj, mode, overwrite):
        """Open a FITS file from a file object (including compressed files)."""

        closed = fileobj_closed(fileobj)
        fmode = fileobj_mode(fileobj) or IO_FITS_MODES[mode]

        if mode == 'ostream':
            self._overwrite_existing(overwrite, fileobj, closed)

        if not closed:
            self._file = fileobj
        elif isfile(fileobj):
            self._file = fileobj_open(self.name, IO_FITS_MODES[mode])

        # Attempt to determine if the file represented by the open file object
        # is compressed
        try:
            # We need to account for the possibility that the underlying file
            # handle may have been opened with either 'ab' or 'ab+', which
            # means that the current file position is at the end of the file.
            if mode in ['ostream', 'append']:
                self._file.seek(0)
            magic = self._file.read(4)
            # No matter whether the underlying file was opened with 'ab' or
            # 'ab+', we need to return to the beginning of the file in order
            # to properly process the FITS header (and handle the possibility
            # of a compressed file).
            self._file.seek(0)
        except (OSError,OSError):
            return

        self._try_read_compressed(fileobj, magic, mode)

    def _open_filelike(self, fileobj, mode, overwrite):
        """Open a FITS file from a file-like object, i.e. one that has
        read and/or write methods.
        """

        self.file_like = True
        self._file = fileobj

        if fileobj_closed(fileobj):
            raise OSError("Cannot read from/write to a closed file-like "
                          "object ({!r}).".format(fileobj))

        if isinstance(fileobj, zipfile.ZipFile):
            self._open_zipfile(fileobj, mode)
            # We can bypass any additional checks at this point since now
            # self._file points to the temp file extracted from the zip
            return

        # If there is not seek or tell methods then set the mode to
        # output streaming.
        if (not hasattr(self._file, 'seek') or
            not hasattr(self._file, 'tell')):
            self.mode = mode = 'ostream'

        if mode == 'ostream':
            self._overwrite_existing(overwrite, fileobj, False)

        # Any "writeable" mode requires a write() method on the file object
        if (self.mode in ('update', 'append', 'ostream') and
            not hasattr(self._file, 'write')):
            raise OSError("File-like object does not have a 'write' "
                          "method, required for mode '{}'.".format(self.mode))

        # Any mode except for 'ostream' requires readability
        if self.mode != 'ostream' and not hasattr(self._file, 'read'):
            raise OSError("File-like object does not have a 'read' "
                          "method, required for mode {!r}.".format(self.mode))

    def _open_filename(self, filename, mode, overwrite):
        """Open a FITS file from a filename string."""

        if mode == 'ostream':
            self._overwrite_existing(overwrite, None, True)

        if os.path.exists(self.name):
            with fileobj_open(self.name, 'rb') as f:
                magic = f.read(4)
        else:
            magic = b''

        ext = os.path.splitext(self.name)[1]

        if not self._try_read_compressed(self.name, magic, mode, ext=ext):
            self._file = fileobj_open(self.name, IO_FITS_MODES[mode])
            self.close_on_error = True

        # Make certain we're back at the beginning of the file
        # BZ2File does not support seek when the file is open for writing, but
        # when opening a file for write, bz2.BZ2File always truncates anyway.
        if not (isinstance(self._file, bz2.BZ2File) and mode == 'ostream'):
            self._file.seek(0)

    @classproperty(lazy=True)
    def _mmap_available(cls):
        """Tests that mmap, and specifically mmap.flush works.  This may
        be the case on some uncommon platforms (see
        https://github.com/astropy/astropy/issues/968).

        If mmap.flush is found not to work, ``self.memmap = False`` is
        set and a warning is issued.
        """

        tmpfd, tmpname = tempfile.mkstemp()
        try:
            # Windows does not allow mappings on empty files
            os.write(tmpfd, b' ')
            os.fsync(tmpfd)
            try:
                mm = mmap.mmap(tmpfd, 1, access=mmap.ACCESS_WRITE)
            except OSError as exc:
                warnings.warn('Failed to create mmap: {}; mmap use will be '
                              'disabled'.format(str(exc)), AstropyUserWarning)
                del exc
                return False
            try:
                mm.flush()
            except OSError:
                warnings.warn('mmap.flush is unavailable on this platform; '
                              'using mmap in writeable mode will be disabled',
                              AstropyUserWarning)
                return False
            finally:
                mm.close()
        finally:
            os.close(tmpfd)
            os.remove(tmpname)

        return True

    def _open_zipfile(self, fileobj, mode):
        """Limited support for zipfile.ZipFile objects containing a single
        a file.  Allows reading only for now by extracting the file to a
        tempfile.
        """

        if mode in ('update', 'append'):
            raise OSError(
                  "Writing to zipped fits files is not currently "
                  "supported")

        if not isinstance(fileobj, zipfile.ZipFile):
            zfile = zipfile.ZipFile(fileobj)
            close = True
        else:
            zfile = fileobj
            close = False

        namelist = zfile.namelist()
        if len(namelist) != 1:
            raise OSError(
              "Zip files with multiple members are not supported.")
        self._file = tempfile.NamedTemporaryFile(suffix='.fits')
        self._file.write(zfile.read(namelist[0]))

        if close:
            zfile.close()
        # We just wrote the contents of the first file in the archive to a new
        # temp file, which now serves as our underlying file object. So it's
        # necessary to reset the position back to the beginning
        self._file.seek(0)
