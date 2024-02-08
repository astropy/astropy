# Licensed under a 3-clause BSD style license - see LICENSE.rst


from io import StringIO

import numpy as np
import pytest

from astropy.io import ascii

from .common import assert_equal


def test_types_from_dat():
    converters = {"a": [ascii.convert_numpy(float)], "e": [ascii.convert_numpy(str)]}

    dat = ascii.read(
        ["a b c d e", "1 1 cat 2.1 4.2"], format="basic", converters=converters
    )

    assert dat["a"].dtype.kind == "f"
    assert dat["b"].dtype.kind == "i"
    assert dat["c"].dtype.kind in ("S", "U")
    assert dat["d"].dtype.kind == "f"
    assert dat["e"].dtype.kind in ("S", "U")


def test_rdb_write_types():
    dat = ascii.read(["a b c d", "1 1.0 cat 2.1"], format="basic")
    out = StringIO()
    ascii.write(dat, out, format="rdb")
    outs = out.getvalue().splitlines()
    assert_equal(outs[1], "N\tN\tS\tN")


def test_ipac_read_types():
    table = r"""\
|     ra   |    dec   |   sai   |-----v2---|    sptype        |
|    real  |   float  |   l     |    real  |     char         |
|    unit  |   unit   |   unit  |    unit  |     ergs         |
|    null  |   null   |   null  |    null  |     -999         |
   2.09708   2956        73765    2.06000   B8IVpMnHg
"""
    reader = ascii.get_reader(reader_cls=ascii.Ipac)
    reader.read(table)
    types = [
        ascii.FloatType,
        ascii.FloatType,
        ascii.IntType,
        ascii.FloatType,
        ascii.StrType,
    ]
    for col, expected_type in zip(reader.cols, types):
        assert_equal(col.type, expected_type)


@pytest.mark.parametrize("fast", [False, "force"])
def test_large_integers(fast):
    # case https://github.com/astropy/astropy/issues/5744
    dat = ascii.read(
        """\
num,cadence
110000000000000001,1
12345,2
""",
        fast_reader=fast,
    )
    # assert `num`` is exactly 64bit int
    assert dat["num"].dtype.kind == "i"
    assert dat["num"].dtype.itemsize == 8
    assert dat["num"][0] == 110000000000000001

    # assert `cadence` is exactly 64bit int, even though the values are small numbers.
    assert dat["cadence"].dtype.kind == "i"
    assert dat["cadence"].dtype.itemsize == 8


def test_col_dtype_in_custom_class():
    """Test code in BaseOutputter._convert_vals to handle Column.dtype
    attribute. See discussion in #11895."""
    dtypes = [np.float32, np.int8, np.int16]

    class TestDtypeHeader(ascii.BasicHeader):
        def get_cols(self, lines):
            super().get_cols(lines)
            for col, dtype in zip(self.cols, dtypes):
                col.dtype = dtype

    class TestDtype(ascii.Basic):
        """
        Basic table Data reader_cls with data type alternating float32, int8
        """

        header_class = TestDtypeHeader

    txt = """
    a b c
    1 2 3
    """
    reader = ascii.get_reader(TestDtype)
    t = reader.read(txt)
    for col, dtype in zip(t.itercols(), dtypes):
        assert col.dtype.type is dtype
