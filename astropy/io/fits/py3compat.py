# Licensed under a 3-clause BSD style license - see PYFITS.rst

from ...extern import six
from ...utils.compat.numpycompat import NUMPY_LT_1_10

if six.PY3:
    # Stuff to do if Python 3

    # Make the decode_ascii utility function actually work
    from . import util
    import numpy

    def encode_ascii(s):
        if isinstance(s, str):
            return s.encode('ascii')
        elif (isinstance(s, numpy.ndarray) and
              issubclass(s.dtype.type, numpy.str_)):
            ns = numpy.char.encode(s, 'ascii').view(type(s))
            if ns.dtype.itemsize != s.dtype.itemsize / 4:
                ns = ns.astype((numpy.bytes_, s.dtype.itemsize / 4))
            return ns
        elif (isinstance(s, numpy.ndarray) and
              not issubclass(s.dtype.type, numpy.bytes_)):
            raise TypeError('string operation on non-string array')
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
            if s.size == 0:
                # Numpy apparently also has a bug that if a string array is
                # empty calling np.char.decode on it returns an empty float64
                # array wth
                dt = s.dtype.str.replace('S', 'U')
                ns = numpy.array([], dtype=dt).view(type(s))
            else:
                ns = numpy.char.decode(s, 'ascii').view(type(s))
            if ns.dtype.itemsize / 4 != s.dtype.itemsize:
                ns = ns.astype((numpy.str_, s.dtype.itemsize))
            return ns
        elif (isinstance(s, numpy.ndarray) and
              not issubclass(s.dtype.type, numpy.str_)):
            # Don't silently pass through on non-string arrays; we don't want
            # to hide errors where things that are not stringy are attempting
            # to be decoded
            raise TypeError('string operation on non-string array')
        return s
    util.decode_ascii = decode_ascii

    # Here we monkey patch (yes, I know) numpy to fix a few numpy Python 3
    # bugs.  The only behavior that's modified is that bugs are fixed, so that
    # should be OK.

    # Fix chararrays; this is necessary in numpy 1.9.x and below
    # The fix for this is in https://github.com/numpy/numpy/pull/5982 and is
    # available as of Numpy 1.10

    if NUMPY_LT_1_10:
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

        for m in [numpy, numpy.char, numpy.core.defchararray,
                  numpy.core.records]:
            m.chararray = chararray
