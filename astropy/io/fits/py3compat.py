# Licensed under a 3-clause BSD style license - see PYFITS.rst

from ...extern import six

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
