# Licensed under a 3-clause BSD style license - see PYFITS.rst

from ...utils.compat.numpycompat import NUMPY_LT_1_10


import numpy

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
