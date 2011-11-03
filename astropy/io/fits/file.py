from __future__ import division

import gzip
import os
import sys
import tempfile
import urllib
import warnings
import zipfile

import numpy as np
from numpy import memmap as Memmap

from pyfits.util import (Extendable, isreadable, iswritable, isfile,
                         fileobj_name, fileobj_closed, fileobj_mode,
                         _array_from_file, _array_to_file, _write_string,
                         deprecated)


PYTHON_MODES = {'readonly': 'rb', 'copyonwrite': 'rb', 'update': 'rb+',
                'append': 'ab+', 'ostream': 'w'}  # open modes
MEMMAP_MODES = {'readonly': 'r', 'copyonwrite': 'c', 'update': 'r+'}


class _File(object):
    """
    Represents a FITS file on disk (or in some other file-like object).
    """

    __metaclass__ = Extendable

    def __init__(self, fileobj=None, mode='copyonwrite', memmap=False):
        if fileobj is None:
            self.__file = None
            self.closed = False
            self.mode = mode
            self.memmap = memmap
            self.compression = None
            self.readonly = False
            self.writeonly = False
            self.simulateonly = True
            return
        else:
            self.simulateonly = False

        if mode not in PYTHON_MODES:
            raise ValueError("Mode '%s' not recognized" % mode)

        if (isinstance(fileobj, basestring) and mode != 'append' and
            not os.path.exists(fileobj) and
            not os.path.splitdrive(fileobj)[0]):
                #
                # Not writing file and file does not exist on local machine and
                # name does not begin with a drive letter (Windows), try to
                # get it over the web.
                #
            self.name, _ = urllib.urlretrieve(fileobj)
        else:
            self.name = fileobj_name(fileobj)

        self.closed = False
        self.mode = mode
        self.memmap = memmap

        # Underlying fileobj is a file-like object, but an actual file object
        self.file_like = False

        self.compression = None
        if isinstance(fileobj, gzip.GzipFile):
            self.compression = 'gzip'
        elif isinstance(fileobj, zipfile.ZipFile):
            # Reading from zip files is supported but not writing (yet)
            self.compression = 'zip'

        self.readonly = False
        self.writeonly = False
        if (mode in ('readonly', 'copyonwrite') or
                (self.compression and mode == 'update')):
            self.readonly = True
        elif (mode == 'ostream' or
                (self.compression and mode == 'append')):
            self.writeonly = True

        # Initialize the internal self.__file object
        if isfile(fileobj) or isinstance(fileobj, gzip.GzipFile):
            closed = fileobj_closed(fileobj)
            fmode = fileobj_mode(fileobj) or PYTHON_MODES[mode]

            if not closed:
                # In some cases (like on Python 3) a file opened for
                # appending still shows a mode of 'r+', hence the extra
                # check for the append case
                if ((mode == 'append' and fmode not in ('ab+', 'rb+')) or
                    (mode != 'append' and PYTHON_MODES[mode] != fmode)):
                    raise ValueError(
                        "Input mode '%s' (%s) does not match mode of the "
                        "input file (%s)." % (mode, PYTHON_MODES[mode], fmode))
                self.__file = fileobj
            elif isfile(fileobj):
                self.__file = open(self.name, PYTHON_MODES[mode])
                # Return to the beginning of the file--in Python 3 when
                # opening in append mode the file pointer is at the end of
                # the file
                self.__file.seek(0)
            else:
                self.__file = gzip.open(self.name, PYTHON_MODES[mode])
        elif isinstance(fileobj, basestring):
            if os.path.splitext(self.name)[1] == '.gz':
                # Handle gzip files
                if mode in ['update', 'append']:
                    raise IOError(
                          "Writing to gzipped fits files is not currently "
                          "supported")
                self.__file = gzip.open(self.name)
                self.compression = 'gzip'
            elif os.path.splitext(self.name)[1] == '.zip':
                # Handle zip files
                if mode in ['update', 'append']:
                    raise IOError(
                          "Writing to zipped fits files is not currently "
                          "supported")
                zfile = zipfile.ZipFile(self.name)
                namelist = zfile.namelist()
                if len(namelist) != 1:
                    raise IOError(
                      "Zip files with multiple members are not supported.")
                self.__file = tempfile.NamedTemporaryFile(suffix='.fits')
                self.__file.write(zfile.read(namelist[0]))
                zfile.close()
                self.compression = 'zip'
            else:
                self.__file = open(self.name, PYTHON_MODES[mode])
                # Make certain we're back at the beginning of the file
            self.__file.seek(0)
        else:
            # We are dealing with a file like object.
            # Assume it is open.
            self.file_like = True
            self.__file = fileobj

            # If there is not seek or tell methods then set the mode to
            # output streaming.
            if (not hasattr(self.__file, 'seek') or
                not hasattr(self.__file, 'tell')):
                self.mode = mode = 'ostream'

            if (self.mode in ('copyonwrite', 'update', 'append') and
                not hasattr(self.__file, 'write')):
                raise IOError("File-like object does not have a 'write' "
                              "method, required for mode '%s'."
                              % self.mode)

            if self.mode == 'readonly' and not hasattr(self.__file, 'read'):
                raise IOError("File-like object does not have a 'read' "
                              "method, required for mode 'readonly'."
                              % self.mode)

        # For 'ab+' mode, the pointer is at the end after the open in
        # Linux, but is at the beginning in Solaris.
        if (mode == 'ostream' or self.compression or
            not hasattr(self.__file, 'seek')):
            # For output stream start with a truncated file.
            # For compressed files we can't really guess at the size
            self.size = 0
        else:
            pos = self.__file.tell()
            self.__file.seek(0, 2)
            self.size = self.__file.tell()
            self.__file.seek(pos)

        if self.memmap and not isfile(self.__file):
            self.memmap = False
            warnings.warn('Disabling mmap for non-file-backed file-like '
                          'object.')

    def __repr__(self):
        return '<%s.%s %s>' % (self.__module__, self.__class__.__name__,
                               self.__file)

    # Support the 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    @deprecated
    def getfile(self):
        """**Deprecated** Will be going away as soon as I figure out how."""
        return self.__file

    def readable(self):
        if self.writeonly:
            return False
        return isreadable(self.__file)

    def read(self, size=None):
        if not hasattr(self.__file, 'read'):
            raise EOFError
        return self.__file.read(size)

    def readarray(self, size=None, offset=0, dtype=np.uint8, shape=None):
        """
        Similar to file.read(), but returns the contents of the underlying
        file as a numpy array (or mmap'd array if memmap=True) rather than a
        string.

        Usually it's best not to use the `size` argument with this method, but
        it's provided for compatibility.
        """

        if not hasattr(self.__file, 'read'):
            raise EOFError

        if not isinstance(dtype, np.dtype):
            dtype = np.dtype(dtype)

        if size and size % dtype.itemsize != 0:
            raise ValueError('size %d not a multiple of %s' % (size, dtype))

        if isinstance(shape, int):
            shape = (shape,)

        if size and shape:
            actualsize = sum(dim * dtype.itemsize for dim in shape)
            if actualsize < size:
                raise ValueError('size %d is too few bytes for a %s array of '
                                 '%s' % (size, shape, dtype))
            if actualsize < size:
                raise ValueError('size %d is too many bytes for a %s array of '
                                 '%s' % (size, shape, dtype))

        if size and not shape:
            shape = (size // dtype.itemsize,)

        if not (size or shape):
            warnings.warn('No size or shape given to readarray(); assuming a '
                          'shape of (1,)')
            shape = (1,)

        if self.memmap:
            return Memmap(self.__file, offset=offset,
                          mode=MEMMAP_MODES[self.mode], dtype=dtype,
                          shape=shape)
        else:
            count = reduce(lambda x, y: x * y, shape)
            pos = self.__file.tell()
            self.__file.seek(offset)
            data = _array_from_file(self.__file, dtype, count, '')
            data.shape = shape
            self.__file.seek(pos)
            return data

    def writable(self):
        if self.readonly:
            return False
        return iswritable(self.__file)

    def write(self, string):
        if hasattr(self.__file, 'write'):
            _write_string(self.__file, string)

    def writearray(self, array):
        """
        Similar to file.write(), but writes a numpy array instead of a string.

        Also like file.write(), a flush() or close() may be needed before
        the file on disk reflects the data written.
        """

        if hasattr(self.__file, 'write'):
            _array_to_file(array, self.__file)

    def flush(self):
        if hasattr(self.__file, 'flush'):
            self.__file.flush()

    def seek(self, offset, whence=0):
        # In newer Python versions, GzipFiles support the whence argument, but
        # I don't think it was added until 2.6; instead of assuming it's
        # present, we implement our own support for it here
        if not hasattr(self.__file, 'seek'):
            return
        if isinstance(self.__file, gzip.GzipFile):
            if whence:
                if whence == 1:
                    offset = self.__file.offset + offset
                else:
                    raise ValueError('Seek from end not supported')
            self.__file.seek(offset)
        else:
            self.__file.seek(offset, whence)

        pos = self.__file.tell()
        if self.size and pos > self.size:
            warnings.warn('File may have been truncated: actual file length '
                          '(%i) is smaller than the expected size (%i)' %
                          (self.size, pos))

    def tell(self):
        if not hasattr(self.__file, 'tell'):
            raise EOFError
        return self.__file.tell()

    def truncate(self, size=None):
        if hasattr(self.__file, 'truncate'):
            self.__file.truncate(size)

    def close(self):
        """
        Close the 'physical' FITS file.
        """

        if hasattr(self.__file, 'close'):
            self.__file.close()

        self.closed = True

