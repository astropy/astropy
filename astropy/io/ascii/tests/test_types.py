# Licensed under a 3-clause BSD style license - see LICENSE.rst


from io import StringIO

import numpy as np

from ... import ascii
from ....tests.helper import pytest

from .common import assert_equal


def test_types_from_dat():
    converters = {'a': [ascii.convert_numpy(float)],
                  'e': [ascii.convert_numpy(str)]}

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
|     ra   |    dec   |   sai   |   sai2  |-----v2---|    sptype        |
|    real  |   float  |   l     |   int   |    real  |     char         |
|    unit  |   unit   |   unit  |   unit  |    unit  |     ergs         |
|    null  |   null   |   null  |   null  |    null  |     -999         |
   2.09708   2956        73765     73765    2.06000   B8IVpMnHg
"""
    reader = ascii.get_reader(Reader=ascii.Ipac)
    dat = reader.read(table)
    types = [ascii.FloatType,
             ascii.FloatType,
             ascii.IntType,
             ascii.IntType,
             ascii.FloatType,
             ascii.StrType]

    for (col, expected_type, dat_col) in zip(reader.cols, types, dat.columns.values()):
        assert_equal(col.type, expected_type)

        # Explicitly check that IPAC numeric types are 64-bits on all platforms (#4684).
        if col.type is ascii.FloatType:
            assert dat_col.dtype.type is np.float64
        elif col.type is ascii.IntType:
            assert dat_col.dtype.type is np.int64
        else:
            assert dat_col.dtype.kind in ('S', 'U')


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_data_type(fast_reader):
    """
    Simple test that data type conversion is yielding int64, float64, and str
    types by default.  On Windows np.int is 32 bit.
    """
    t = ascii.read(['a b c', '1 1.0 x'], format='basic', fast_reader=fast_reader, guess=False)
    assert t['a'].dtype.type is np.int64
    assert t['b'].dtype.type is np.float64
    assert t['c'].dtype.kind in ('S', 'U')
