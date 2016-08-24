# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from ....extern.six.moves import cStringIO as StringIO

import numpy as np

from ... import ascii
from ....extern.six.moves import zip

from .common import assert_equal, setup_function, teardown_function


def test_types_from_dat():
    converters = {'a': [ascii.convert_numpy(np.float)],
                  'e': [ascii.convert_numpy(np.str)]}

    dat = ascii.read(['a b c d e', '1 1 cat 2.1 4.2'],
                          Reader=ascii.Basic,
                          converters=converters)

    assert dat['a'].dtype.kind == 'f'
    assert dat['b'].dtype.kind == 'i'
    assert dat['c'].dtype.kind in ('S', 'U')
    assert dat['d'].dtype.kind == 'f'
    assert dat['e'].dtype.kind in ('S', 'U')


def test_rdb_write_types():
    dat = ascii.read(['a b c d', '1 1.0 cat 2.1'],
                          Reader=ascii.Basic)
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.Rdb)
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
    reader = ascii.get_reader(Reader=ascii.Ipac)
    dat = reader.read(table)
    types = [ascii.FloatType,
             ascii.FloatType,
             ascii.IntType,
             ascii.FloatType,
             ascii.StrType]
    for (col, expected_type) in zip(reader.cols, types):
        assert_equal(col.type, expected_type)
