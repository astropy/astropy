import sys

if sys.version_info[0] >= 3:
    # Stuff to do if Python 3
    import builtins
    import io

    # Bring back the cmp() function
    builtins.cmp = lambda a, b: (a > b) - (a < b)


    # Make the decode_ascii utility function actually work
    import pyfits.util
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
    pyfits.util.encode_ascii = encode_ascii

    def decode_ascii(s):
        if isinstance(s, bytes):
            return s.decode('ascii')
        elif isinstance(s, numpy.ndarray) and \
             issubclass(s.dtype.type, numpy.bytes_):
            # np.char.encode/decode annoyingly don't preserve the type of the
            # array, hence the view() call
            # It also doesn't necessarily preserve widths of the strings,
            # hence the astype()
            ns = numpy.char.decode(s, 'ascii').view(type(s))
            if ns.dtype.itemsize / 4 != s.dtype.itemsize:
                ns = ns.astype((numpy.str_, s.dtype.itemsize))
            return ns
        return s
    pyfits.util.decode_ascii = decode_ascii

    # Support the io.IOBase.readable/writable methods
    from pyfits.util import isreadable as _isreadable
    def isreadable(f):
        if hasattr(f, 'readable'):
            return f.readable()
        return _isreadable(f)
    pyfits.util.isreadable = isreadable

    from pyfits.util import iswritable as _iswritable
    def iswritable(f):
        if hasattr(f, 'writable'):
            return f.writable()
        return _iswritable(f)
    pyfits.util.iswritable = iswritable

    # isfile needs to support the higher-level wrappers around FileIO
    def isfile(f):
        if isinstance(f, io.FileIO):
            return True
        elif hasattr(f, 'buffer'):
            return isfile(f.buffer)
        elif hasattr(f, 'raw'):
            return isfile(f.raw)
        return False
    pyfits.util.isfile = isfile

    # Here we monkey patch (yes, I know) numpy to fix a few numpy Python 3
    # bugs.  The only behavior that's modified is that bugs are fixed, so that
    # should be OK.

    # Fix chararrays; this is necessary in numpy 1.5.1 and below--hopefully
    # should not be necessary later.  See
    # http://projects.scipy.org/numpy/ticket/1817
    # TODO: Maybe do a version check on numpy for this?  (Note: the fix for
    # this hasn't been accepted in Numpy yet, so a version number check would
    # not be helpful yet...)
    import pyfits.file

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

        if dtype.fields is None:
            return dtype

        formats = []
        offsets = []
        for field in (dtype.fields[name] for name in dtype.names):
            shape = field[0].shape
            if not isinstance(shape, tuple):
                shape = (shape,)
            formats.append((field[0].base, shape))
            offsets.append(field[1])

        return numpy.dtype({'names': dtype.names, 'formats': formats,
                            'offsets': offsets})

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


    # We also need to patch pyfits.file._File which can also be affected by the
    # #1766 bug
    old_File = pyfits.file._File
    class _File(old_File):
        def readarray(self, size=None, offset=None, dtype=numpy.uint8,
                      shape=None):
            if isinstance(dtype, numpy.dtype):
                dtype = _fix_dtype(dtype)
            return old_File.readarray(self, size, offset, dtype, shape)
        readarray.__doc__ = old_File.readarray.__doc__
    pyfits.file._File = _File


    # Replace pyfits.util.maketrans and translate with versions that work
    # with Python 3 unicode strings
    pyfits.util.maketrans = str.maketrans

    def translate(s, table, deletechars):
        if deletechars:
            table = table.copy()
            for c in deletechars:
                table[ord(c)] = None
        return s.translate(table)
    pyfits.util.translate = translate
else:
    # Stuff to do if not Python 3

    # We want a bytes stand-in so we can do stuff like isinstance(..., bytes)
    # and have it work correctly
    import __builtin__
    if not hasattr(__builtin__, 'bytes'):
        __builtin__.bytes = str

    import string
    import pyfits.util
    pyfits.util.maketrans = string.maketrans
