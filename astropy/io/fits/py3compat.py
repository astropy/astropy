# Licensed under a 3-clause BSD style license - see PYFITS.rst

import sys

PY3 = sys.version_info[0] >= 3

if PY3:  # pragma: py3
    # Stuff to do if Python 3
    import builtins
    import io

    # Bring back the cmp() function
    builtins.cmp = lambda a, b: (a > b) - (a < b)

    # Make the decode_ascii utility function actually work
    from . import util
    import numpy

    def encode_ascii(s):
        if isinstance(s, str):
            return s.encode('ascii')
        elif isinstance(s, numpy.ndarray) and \
                issubclass(s.dtype.type, numpy.str_):
            ns = numpy.char.encode(s, 'ascii').view(type(s))
            if ns.dtype.itemsize != s.dtype.itemsize / 4:
                ns = ns.astype((numpy.bytes_, s.dtype.itemsize / 4))
            return ns
        return s
    util.encode_ascii = encode_ascii

    def decode_ascii(s):
        if isinstance(s, bytes):
            return s.decode('ascii')
        elif (isinstance(s, numpy.ndarray) and
              issubclass(s.dtype.type, numpy.bytes_)):
            # np.char.encode/decode annoyingly don't preserve the type of the
            # array, hence the view() call
            # It also doesn't necessarily preserve widths of the strings,
            # hence the astype()
            ns = numpy.char.decode(s, 'ascii').view(type(s))
            if ns.dtype.itemsize / 4 != s.dtype.itemsize:
                ns = ns.astype((numpy.str_, s.dtype.itemsize))
            return ns
        return s
    util.decode_ascii = decode_ascii

    # Replacements for b and u marks on strings
    def b(s):
        return s.encode('latin-1')

    def u(s):
        return s

    util.b = b
    util.u = u

    # See the docstring for astropy.io.fits.util.fileobj_open for why we need
    # to replace this function
    def fileobj_open(filename, mode):
        return open(filename, mode, buffering=0)
    util.fileobj_open = fileobj_open

    # Support the io.IOBase.readable/writable methods
    from .util import isreadable as _isreadable

    def isreadable(f):
        if hasattr(f, 'readable'):
            return f.readable()
        return _isreadable(f)
    util.isreadable = isreadable

    from .util import iswritable as _iswritable

    def iswritable(f):
        if hasattr(f, 'writable'):
            return f.writable()
        return _iswritable(f)
    util.iswritable = iswritable

    # isfile needs to support the higher-level wrappers around FileIO
    def isfile(f):
        if isinstance(f, io.FileIO):
            return True
        elif hasattr(f, 'buffer'):
            return isfile(f.buffer)
        elif hasattr(f, 'raw'):
            return isfile(f.raw)
        return False
    util.isfile = isfile

    # Here we monkey patch (yes, I know) numpy to fix a few numpy Python 3
    # bugs.  The only behavior that's modified is that bugs are fixed, so that
    # should be OK.

    # Fix chararrays; this is necessary in numpy 1.5.1 and below--hopefully
    # should not be necessary later.  See
    # http://projects.scipy.org/numpy/ticket/1817
    # TODO: Maybe do a version check on numpy for this?  (Note: the fix for
    # this hasn't been accepted in Numpy yet, so a version number check would
    # not be helpful yet...)
    from . import file

    _chararray = numpy.char.chararray

    class chararray(_chararray):
        def __getitem__(self, obj):
                val = numpy.ndarray.__getitem__(self, obj)
                if isinstance(val, numpy.character):
                    temp = val.rstrip()
                    if numpy.char._len(temp) == 0:
                        val = ''
                    else:
                        val = temp
                return val
    for m in [numpy.char, numpy.core.defchararray, numpy.core.records]:
        m.chararray = chararray

    # Fix recarrays with sub-array fields.  See
    # http://projects.scipy.org/numpy/ticket/1766
    # TODO: Same as above, though the fix to this problem hasn't made it into
    # any Numpy release yet either, so we'll have to hold off on a version
    # check
    def _fix_dtype(dtype):
        """
        Numpy has a bug (in Python3 only) that causes a segfault when
        accessing the data of arrays containing nested arrays.  Specifically,
        this happens if the shape of the subarray is not given as a tuple.
        See http://projects.scipy.org/numpy/ticket/1766.
        """

        if not hasattr(dtype, 'fields') or dtype.fields is None:
            return dtype

        formats = []
        offsets = []
        titles = []
        for name in dtype.names:
            field = dtype.fields[name]
            shape = field[0].shape
            if not isinstance(shape, tuple):
                shape = (shape,)
            formats.append((field[0].base, shape))
            offsets.append(field[1])

            # There seems to be no obvious way to extract the titles from
            # a dtype, so this just searches for duplicate fields
            title = None
            for key, dup in dtype.fields.items():
                if key != name and dup == field:
                    title = key
                    break
            titles.append(title)

        return numpy.dtype({'names': dtype.names, 'formats': formats,
                            'offsets': offsets, 'titles': titles})

    _recarray = numpy.recarray

    class recarray(_recarray):
        def __new__(subtype, shape, dtype=None, buf=None, offset=0,
                    strides=None, formats=None, names=None, titles=None,
                    byteorder=None, aligned=False, order='C'):
            if dtype is not None:
                dtype = _fix_dtype(dtype)

            if 'order' in _recarray.__new__.__code__.co_varnames:
                return _recarray.__new__(
                    subtype, shape, dtype, buf, offset, strides, formats,
                    names, titles, byteorder, aligned, order)
            else:
                return _recarray.__new__(
                    subtype, shape, dtype, buf, offset, strides, formats,
                    names, titles, byteorder, aligned)
    numpy.recarray = numpy.core.records.recarray = recarray

    # We also need to patch astropy.io.fits.file._File which can also be
    # affected by the #1766 bug
    old_File = file._File

    class _File(old_File):
        def readarray(self, size=None, offset=0, dtype=numpy.uint8,
                      shape=None):
            if isinstance(dtype, numpy.dtype):
                dtype = _fix_dtype(dtype)
            return old_File.readarray(self, size, offset, dtype, shape)
        readarray.__doc__ = old_File.readarray.__doc__
    file._File = _File

    # Replace astropy.io.fits.util.maketrans and translate with versions that
    # work with Python 3 unicode strings
    util.maketrans = str.maketrans

    def translate(s, table, deletechars):
        if deletechars:
            table = table.copy()
            for c in deletechars:
                table[ord(c)] = None
        return s.translate(table)
    util.translate = translate
else:
    # Stuff to do if not Python 3
    import string
    from . import util
    util.maketrans = string.maketrans

    def b(s):
        return s

    def u(s):
        return unicode(s, 'unicode_escape')

    util.b = b
    util.u = u
