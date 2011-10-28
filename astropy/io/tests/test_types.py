import re
import sys
import glob
import math

try:
    import StringIO as io
except ImportError:
    import io

try:
    from .. import ascii as asciitable
except ImportError:
    from .. import asciitable

if asciitable.has_numpy:
    import numpy as np

from .common import *

@has_numpy_and_not_has_numpy
def test_types_from_dat(numpy):
    if numpy:
        converters = {'a': [asciitable.convert_numpy(np.float)],
                      'e': [asciitable.convert_numpy(np.str)]}
    else:
        converters = {'a': [asciitable.convert_list(float)],
                      'e': [asciitable.convert_list(str)]}

    dat = asciitable.read(['a b c d e', '1 1 cat 2.1 4.2'], Reader=asciitable.Basic,
                          converters=converters, numpy=numpy)

    reader = asciitable.get_reader(Reader=asciitable.Memory, numpy=numpy)
    reader.read(dat)

    print('numpy=%s' % numpy)
    print('dat=%s' % repr(dat))
    print('reader.table=%s' % repr(reader.table))
    print('types=%s' % repr([x.type for x in reader.cols]))

    assert_true(issubclass(reader.cols[0].type, asciitable.FloatType))
    assert_true(issubclass(reader.cols[1].type, asciitable.IntType))
    assert_true(issubclass(reader.cols[2].type, asciitable.StrType))
    assert_true(issubclass(reader.cols[3].type, asciitable.FloatType))
    assert_true(issubclass(reader.cols[4].type, asciitable.StrType))
    
@has_numpy_and_not_has_numpy
def test_rdb_write_types(numpy):
    dat = asciitable.read(['a b c d', '1 1.0 cat 2.1'], Reader=asciitable.Basic, numpy=numpy)
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.Rdb)
    outs = out.getvalue().splitlines()
    assert_equal(outs[1], 'N\tN\tS\tN')

@has_numpy_and_not_has_numpy
def test_ipac_read_types(numpy):
    table = r"""\
|     ra   |    dec   |   sai   |-----v2---|    sptype        |
|    real  |   float  |   l     |    real  |     char         |
|    unit  |   unit   |   unit  |    unit  |     ergs         |
|    null  |   null   |   null  |    null  |     -999         |
   2.09708   2956        73765    2.06000   B8IVpMnHg
"""
    reader = asciitable.get_reader(Reader=asciitable.Ipac, numpy=numpy)
    dat = reader.read(table)
    types = [asciitable.FloatType,
             asciitable.FloatType,
             asciitable.IntType,
             asciitable.FloatType,
             asciitable.StrType]
    for (col, expected_type) in zip(reader.cols, types):
        assert_equal(col.type, expected_type)


    
