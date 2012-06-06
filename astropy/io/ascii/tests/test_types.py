import re
import sys
import glob
import math

import numpy as np

try:
    import StringIO as io
except ImportError:
    import io

from ... import ascii as asciitable

from .common import (raises,
                     assert_equal, assert_almost_equal, assert_true,
                     setup_function, teardown_function)


def test_types_from_dat():
    converters = {'a': [asciitable.convert_numpy(np.float)],
                  'e': [asciitable.convert_numpy(np.str)]}

    dat = asciitable.read(['a b c d e', '1 1 cat 2.1 4.2'],
                          Reader=asciitable.Basic,
                          converters=converters)

    reader = asciitable.get_reader(Reader=asciitable.Memory)
    reader.read(dat)

    print('dat=%s' % repr(dat))
    print('reader.table=%s' % repr(reader.table))
    print('types=%s' % repr([x.type for x in reader.cols]))

    assert_true(issubclass(reader.cols[0].type, asciitable.FloatType))
    assert_true(issubclass(reader.cols[1].type, asciitable.IntType))
    assert_true(issubclass(reader.cols[2].type, asciitable.StrType))
    assert_true(issubclass(reader.cols[3].type, asciitable.FloatType))
    assert_true(issubclass(reader.cols[4].type, asciitable.StrType))


def test_rdb_write_types():
    dat = asciitable.read(['a b c d', '1 1.0 cat 2.1'],
                          Reader=asciitable.Basic)
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.Rdb)
    outs = out.getvalue().splitlines()
    assert_equal(outs[1], 'N\tN\tS\tN')


def test_ipac_read_types():
    table = r"""\
|     ra   |    dec   |   sai   |-----v2---|    sptype        |
|    real  |   float  |   l     |    real  |     char         |
|    unit  |   unit   |   unit  |    unit  |     ergs         |
|    null  |   null   |   null  |    null  |     -999         |
   2.09708   2956        73765    2.06000   B8IVpMnHg
"""
    reader = asciitable.get_reader(Reader=asciitable.Ipac)
    dat = reader.read(table)
    types = [asciitable.FloatType,
             asciitable.FloatType,
             asciitable.IntType,
             asciitable.FloatType,
             asciitable.StrType]
    for (col, expected_type) in zip(reader.cols, types):
        assert_equal(col.type, expected_type)
