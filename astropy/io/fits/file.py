# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division, with_statement

import mmap
import os
import sys
import tempfile
import warnings
import zipfile
import bz2

from functools import reduce

import numpy as np
from numpy import memmap as Memmap

from .util import (isreadable, iswritable, isfile, fileobj_open, fileobj_name,
                   fileobj_closed, fileobj_mode, _array_from_file,
                   _array_to_file, _write_string, _GZIP_FILE_TYPES)
from ...extern.six import b, string_types
from ...extern.six.moves import urllib
from ...utils.compat import gzip
from ...utils.data import download_file, _is_url
from ...utils.decorators import classproperty
from ...utils.exceptions import AstropyUserWarning


# Maps PyFITS-specific file mode names to the appropriate file modes to use
# for the underlying raw files
# TODO: This should probably renamed IO_FITS_MODES or something, but since it's
# used primarily internally I'm going to leave PYFITS in the name for now for
# in the off chance any third-party software is trying to do anything with this
# object.
PYFITS_MODES = {
    'readonly': 'rb',
    'copyonwrite': 'rb',
    'update': 'rb+',
    'append': 'ab+',
    'ostream': 'wb',
    'denywrite': 'rb'}

# This is the old name of the PYFITS_MODES dict; it is maintained here for
# backwards compatibility and should be removed no sooner than PyFITS 3.4
PYTHON_MODES = PYFITS_MODES

# Maps OS-level file modes to the appropriate PyFITS specific mode to use
# when given file objects but no mode specified; obviously in PYFITS_MODES
# there are overlaps; for example 'readonly' and 'denywrite' both require
# the file to be opened in 'rb' mode.  But 'readonly' is the default
# behavior for such files if not otherwise specified.
# Note: 'ab' is only supported for 'ostream' which is output-only.
FILE_MODES = {
    'rb': 'readonly', 'rb+': 'update',
    'wb': 'ostream', 'wb+': 'update',
    'ab': 'ostream', 'ab+': 'append'}


# readonly actually uses copyonwrite for mmap so that readonly without mmap and
# with mmap still have to same behavior with regard to updating the array.  To
# get a truly readonly mmap use denywrite
# the name 'denywrite' comes from a deprecated flag to mmap() on Linux--it
# should be clarified that 'denywrite' mode is not directly analogous to the
# use of that flag; it was just taken, for lack of anything better, as a name
# that means something like "read only" but isn't readonly.
MEMMAP_MODES = {'readonly': 'c', 'copyonwrite': 'c', 'update': 'r+',
                'append': 'c', 'denywrite': 'r'}

# TODO: Eventually raise a warning, and maybe even later disable the use of
# 'copyonwrite' and 'denywrite' modes unless memmap=True.  For now, however,
# that would generate too many warnings for too many users.  If nothing else,
# wait until the new logging system is in place.

GZIP_MAGIC = b('\x1f\x8b\x08')
PKZIP_MAGIC = b('\x50\x4b\x03\x04')
BZIP2_MAGIC = b('\x42\x5a')

class _File(object):
    """
    Represents a FITS file on disk (or in some other file-like object).
    """

    def __init__(self, fileobj=None, mode=None, memmap=None, clobber=False,
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
            return
        else:
            self.simulateonly = False

        # Holds mmap instance for files that use mmap
        self._mmap = None

        if mode is None:
            if _is_random_access_file_backed(fileobj):
                fmode = fileobj_mode(fileobj)
                # If the mode is unsupported just leave it as None; we'll
                # catch this case below
                mode = FILE_MODES.get(fmode)
            else:
                mode = 'readonly'  # The default

        if mode not in PYFITS_MODES:
            raise ValueError("Mode '%s' not recognized" % mode)

        if (isinstance(fileobj, string_types) and
            mode not in ('ostream', 'append') and
            _is_url(fileobj)): # This is an URL.
                self.name = download_file(fileobj, cache=cache)
        else:
            self.name = fileobj_name(fileobj)

        self.closed = False
        self.binary = True
        self.mode = mode
        self.memmap = memmap

        # Underlying fileobj is a file-like object, but an actual file object
        self.file_like = False

        # More defaults to be adjusted below as necessary
        self.compression = None
        self.readonly = False
        self.writeonly = False

        # Initialize the internal self._file object
        if _is_random_access_file_backed(fileobj):
            self._open_fileobj(fileobj, mode, clobber)
        elif isinstance(fileobj, string_types):
            self._open_filename(fileobj, mode, clobber)
        else:
            self._open_filelike(fileobj, mode, clobber)

        self.fileobj_mode = fileobj_mode(self._file)

        if isinstance(fileobj, _GZIP_FILE_TYPES):
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
        return '<%s.%s %s>' % (self.__module__, self.__class__.__name__,
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
        except IOError:
            # On some versions of Python, it appears, GzipFile will raise an
            # IOError if you try to read past its end (as opposed to just
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
            raise ValueError('size %d not a multiple of %s' % (size, dtype))

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

            if actualsize < size:
                raise ValueError('size %d is too few bytes for a %s array of '
                                 '%s' % (size, shape, dtype))
            if actualsize < size:
                raise ValueError('size %d is too many bytes for a %s array of '
                                 '%s' % (size, shape, dtype))

        if self.memmap:
            if self._mmap is None:
                # Instantiate Memmap array of the file offset at 0
                # (so we can return slices of it to offset anywhere else into
                # the file)
                memmap = Memmap(self._file,
                                mode=MEMMAP_MODES[self.mode],
                                dtype=np.uint8)

                # Now we immediately discard the memmap array; we are really
                # just using it as a factory function to instantiate the mmap
                # object in a convenient way (may later do away with this
                # usage)
                self._mmap = memmap.base

                # Prevent dorking with self._memmap._mmap by memmap.__del__ in
                # Numpy 1.6 (see
                # https://github.com/numpy/numpy/commit/dcc355a0b179387eeba10c95baf2e1eb21d417c7)
                memmap._mmap = None
                del memmap

            return np.ndarray(shape=shape, dtype=dtype, offset=offset,
                              buffer=self._mmap)
        else:
            count = reduce(lambda x, y: x * y, shape)
            pos = self._file.tell()
            self._file.seek(offset)
            data = _array_from_file(self._file, dtype, count, '')
            data.shape = shape
            self._file.seek(pos)
            return data

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
        # In newer Python versions, GzipFiles support the whence argument, but
        # I don't think it was added until 2.6; instead of assuming it's
        # present, we implement our own support for it here
        if not hasattr(self._file, 'seek'):
            return
        if isinstance(self._file, _GZIP_FILE_TYPES):
            if whence:
                if whence == 1:
                    offset = self._file.offset + offset
                else:
                    raise ValueError('Seek from end not supported')
            self._file.seek(offset)
        else:
            self._file.seek(offset, whence)

        pos = self._file.tell()
        if self.size and pos > self.size:
            warnings.warn('File may have been truncated: actual file length '
                          '(%i) is smaller than the expected size (%i)' %
                          (self.size, pos), AstropyUserWarning)

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

    def _overwrite_existing(self, clobber, fileobj, closed):
        """Overwrite an existing file if ``clobber`` is ``True``, otherwise
        raise an IOError.  The exact behavior of this method depends on the
        _File object state and is only meant for use within the ``_open_*``
        internal methods.
        """

        # The file will be overwritten...
        if ((self.file_like and hasattr(fileobj, 'len') and fileobj.len > 0) or
            (os.path.exists(self.name) and os.path.getsize(self.name) != 0)):
            if clobber:
                if self.file_like and hasattr(fileobj, 'truncate'):
                    fileobj.truncate(0)
                else:
                    if not closed:
                        fileobj.close()
                    os.remove(self.name)
            else:
                raise IOError("File %r already exists." % self.name)

    def _open_fileobj(self, fileobj, mode, clobber):
        """Open a FITS file from a file object or a GzipFile object."""

        closed = fileobj_closed(fileobj)
        fmode = fileobj_mode(fileobj) or PYFITS_MODES[mode]

        if mode == 'ostream':
            self._overwrite_existing(clobber, fileobj, closed)

        if not closed:
            # Although we have a specific mapping in PYFITS_MODES from our
            # custom file modes to raw file object modes, many of the latter
            # can be used appropriately for the former.  So determine whether
            # the modes match up appropriately
            if ((mode in ('readonly', 'denywrite', 'copyonwrite') and
                    not ('r' in fmode or '+' in fmode)) or
                    (mode == 'append' and fmode not in ('ab+', 'rb+')) or
                    (mode == 'ostream' and
                     not ('w' in fmode or 'a' in fmode or '+' in fmode)) or
                    (mode == 'update' and fmode not in ('rb+', 'wb+'))):
                raise ValueError(
                    "Mode argument '%s' does not match mode of the input "
                    "file (%s)." % (mode, fmode))
            self._file = fileobj
        elif isfile(fileobj):
            self._file = fileobj_open(self.name, PYFITS_MODES[mode])
        else:
            self._file = gzip.open(self.name, PYFITS_MODES[mode])

        if fmode == 'ab+':
            # Return to the beginning of the file--in Python 3 when opening in
            # append mode the file pointer is at the end of the file
            self._file.seek(0)

    def _open_filelike(self, fileobj, mode, clobber):
        """Open a FITS file from a file-like object, i.e. one that has
        read and/or write methods.
        """

        self.file_like = True
        self._file = fileobj

        if fileobj_closed(fileobj):
            raise IOError("Cannot read from/write to a closed file-like "
                          "object (%r)." % fileobj)

        if isinstance(fileobj, zipfile.ZipFile):
            self._open_zipfile(fileobj, mode)
            self._file.seek(0)
            # We can bypass any additional checks at this point since now
            # self._file points to the temp file extracted from the zip
            return

        # If there is not seek or tell methods then set the mode to
        # output streaming.
        if (not hasattr(self._file, 'seek') or
            not hasattr(self._file, 'tell')):
            self.mode = mode = 'ostream'

        if mode == 'ostream':
            self._overwrite_existing(clobber, fileobj, False)

        # Any "writeable" mode requires a write() method on the file object
        if (self.mode in ('update', 'append', 'ostream') and
            not hasattr(self._file, 'write')):
            raise IOError("File-like object does not have a 'write' "
                          "method, required for mode '%s'."
                          % self.mode)

        # Any mode except for 'ostream' requires readability
        if self.mode != 'ostream' and not hasattr(self._file, 'read'):
            raise IOError("File-like object does not have a 'read' "
                          "method, required for mode %r."
                          % self.mode)

    def _open_filename(self, filename, mode, clobber):
        """Open a FITS file from a filename string."""

        if mode == 'ostream':
            self._overwrite_existing(clobber, None, True)

        if os.path.exists(self.name):
            with fileobj_open(self.name, 'rb') as f:
                magic = f.read(4)
        else:
            magic = b('')

        ext = os.path.splitext(self.name)[1]

        if ext == '.gz' or magic.startswith(GZIP_MAGIC):
            # Handle gzip files
            self._file = gzip.open(self.name, PYFITS_MODES[mode])
            self.compression = 'gzip'
        elif ext == '.zip' or magic.startswith(PKZIP_MAGIC):
            # Handle zip files
            self._open_zipfile(self.name, mode)
        elif ext == '.bz2' or magic.startswith(BZIP2_MAGIC):
            # Handle bzip2 files
            if mode in ['update', 'append']:
                raise IOError("update and append modes are not supported "
                              "with bzip2 files")
            # bzip2 only supports 'w' and 'r' modes
            bzip2_mode = 'w' if mode == 'ostream' else 'r'
            self._file = bz2.BZ2File(self.name, bzip2_mode)
        else:
            self._file = fileobj_open(self.name, PYFITS_MODES[mode])

        # Make certain we're back at the beginning of the file
        # BZ2File does not support seek when the file is open for writing, but
        # when opening a file for write, bz2.BZ2File always truncates anyway.
        if isinstance(self._file, bz2.BZ2File) and mode == 'ostream':
            pass
        else:
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
            except mmap.error as e:
                warnings.warn('Failed to create mmap: %s; mmap use will be '
                              'disabled' % str(e), AstropyUserWarning)
                del exc
                return False
            try:
                mm.flush()
            except mmap.error:
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
            raise IOError(
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
            raise IOError(
              "Zip files with multiple members are not supported.")
        self._file = tempfile.NamedTemporaryFile(suffix='.fits')
        self._file.write(zfile.read(namelist[0]))

        if close:
            zfile.close()
        self.compression = 'zip'


def _is_random_access_file_backed(fileobj):
    """Returns `True` if fileobj is a `file` or `io.FileIO` object or a
    `gzip.GzipFile` object.

    Although reading from a zip file is supported, this does not include
    support for random access, and we do not yet support reading directly
    from an already opened `zipfile.ZipFile` object.
    """

    return isfile(fileobj) or isinstance(fileobj, _GZIP_FILE_TYPES)
