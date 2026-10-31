# Licensed under a 3-clause BSD style license - see LICENSE.rst


from io import StringIO

import numpy as np
import pytest

from astropy.io import ascii
from astropy.io.ascii import ColumnType


class TestColumnType:
    @pytest.mark.parametrize("t", ColumnType)
    def test_top_is_subtype_only_to_itself(self, t: ColumnType) -> None:
        assert ColumnType.ANY.is_subtype(t) == (t is ColumnType.ANY)

    @pytest.mark.parametrize("t", ColumnType)
    def test_all_is_subtype_to_top(self, t: ColumnType) -> None:
        assert t.is_subtype(ColumnType.ANY)

    @pytest.mark.parametrize("t", ColumnType)
    def test_bottom_is_subtype_to_all(self, t: ColumnType) -> None:
        assert ColumnType.UNDEFINED.is_subtype(t)

    @pytest.mark.parametrize("t", ColumnType)
    def test_only_bottom_is_subtype_to_itself(self, t: ColumnType) -> None:
        assert t.is_subtype(ColumnType.UNDEFINED) == (t is ColumnType.UNDEFINED)

    @pytest.mark.parametrize("t", ColumnType)
    def test_reflexivity(self, t: ColumnType) -> None:
        assert t.is_subtype(t)

    @pytest.mark.parametrize(
        "t1, t2, expected",
        [
            # a complete, flat table of all remaining, non-trivial cases
            (ColumnType.BOOL, ColumnType.FLOAT, False),
            (ColumnType.BOOL, ColumnType.INTEGER, False),
            (ColumnType.BOOL, ColumnType.NUMERIC, False),
            (ColumnType.BOOL, ColumnType.STRING, False),
            (ColumnType.FLOAT, ColumnType.BOOL, False),
            (ColumnType.FLOAT, ColumnType.INTEGER, False),
            (ColumnType.FLOAT, ColumnType.NUMERIC, True),
            (ColumnType.FLOAT, ColumnType.STRING, False),
            (ColumnType.INTEGER, ColumnType.BOOL, False),
            (ColumnType.INTEGER, ColumnType.FLOAT, False),
            (ColumnType.INTEGER, ColumnType.NUMERIC, True),
            (ColumnType.INTEGER, ColumnType.STRING, False),
            (ColumnType.NUMERIC, ColumnType.BOOL, False),
            (ColumnType.NUMERIC, ColumnType.FLOAT, False),
            (ColumnType.NUMERIC, ColumnType.INTEGER, False),
            (ColumnType.NUMERIC, ColumnType.STRING, False),
            (ColumnType.STRING, ColumnType.BOOL, False),
            (ColumnType.STRING, ColumnType.FLOAT, False),
            (ColumnType.STRING, ColumnType.INTEGER, False),
            (ColumnType.STRING, ColumnType.NUMERIC, False),
        ],
    )
    def test_is_subtype(
        self,
        t1: ColumnType,
        t2: ColumnType,
        expected: bool,
    ) -> None:
        assert t1.is_subtype(t2) is expected


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
    assert outs[1] == "N\tN\tS\tN"


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
        ascii.ColumnType.FLOAT,
        ascii.ColumnType.FLOAT,
        ascii.ColumnType.INTEGER,
        ascii.ColumnType.FLOAT,
        ascii.ColumnType.STRING,
    ]
    for col, expected_type in zip(reader.cols, types, strict=True):
        assert col.type == expected_type


def test_ipac_read_long_columns():
    """Test for https://github.com/astropy/astropy/issues/15989"""
    test_data = """\
|              oid|expid|cadence|
|             long|    l|    int|
 90000000000000001   123       1
"""
    dat = ascii.read(test_data, format="ipac")

    # assert oid, as a long column, is int64
    oid = dat["oid"]
    assert oid[0] == 90000000000000001
    assert oid.dtype.kind == "i"
    assert oid.dtype.itemsize == 8

    # expid is declared as a long column,
    # the type needs to be int64, even though all
    # the values are within int32 range
    expid = dat["expid"]
    assert expid[0] == 123
    assert expid.dtype.kind == "i"
    assert expid.dtype.itemsize == 8


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
