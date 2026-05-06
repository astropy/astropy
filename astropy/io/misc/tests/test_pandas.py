# Licensed under a 3-clause BSD style license - see LICENSE.rst

from io import BytesIO, StringIO

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io.misc.pandas import connect
from astropy.table import QTable, Table
from astropy.utils.compat.optional_deps import (
    HAS_BS4,
    HAS_HTML5LIB,
    HAS_LXML,
    HAS_OPENPYXL,
)
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

# Check dependencies
pandas = pytest.importorskip("pandas")

HAS_HTML_DEPS = HAS_LXML or (HAS_BS4 and HAS_HTML5LIB)


WRITE_FMTS = [fmt for fmt in connect.PANDAS_FMTS if "write" in connect.PANDAS_FMTS[fmt]]


@pytest.mark.parametrize("fmt", WRITE_FMTS)
def test_read_write_format(fmt):
    """
    Test round-trip through pandas write/read for supported formats.

    :param fmt: format name, e.g. csv, html, json
    :return:
    """
    # Skip the reading tests
    if fmt == "html" and not HAS_HTML_DEPS:
        pytest.skip("Missing lxml or bs4 + html5lib for HTML read/write test")

    if fmt == "excel" and not HAS_OPENPYXL:
        pytest.skip("Missing openpyxl for excel read/write test")

    pandas_fmt = "pandas." + fmt
    # Explicitly provide dtype to avoid casting 'a' to int32.
    # See https://github.com/astropy/astropy/issues/8682
    t = Table(
        [[1, 2, 3], [1.0, 2.5, 5.0], ["a", "b", "c"]], dtype=(np.int64, np.float64, str)
    )
    if fmt == "excel":
        buf = BytesIO()
    else:
        buf = StringIO()
    t.write(buf, format=pandas_fmt)

    buf.seek(0)
    t2 = Table.read(buf, format=pandas_fmt)

    assert t.colnames == t2.colnames

    if fmt == "excel":
        assert np.all(t.as_array() == t2.as_array())
    else:
        assert np.all(t == t2)


@pytest.mark.parametrize("fmt", WRITE_FMTS)
def test_write_overwrite(tmp_path, fmt):
    """Test overwriting."""

    if fmt == "excel":
        tmpfile = tmp_path / "test.xlsx"
        if not HAS_OPENPYXL:
            pytest.skip("Missing openpyxl for excel read/write test")
    else:
        tmpfile = tmp_path / f"test.{fmt}"
    pandas_fmt = f"pandas.{fmt}"

    # Explicitly provide dtype to avoid casting 'a' to int32.
    # See https://github.com/astropy/astropy/issues/8682
    t = Table(
        [[1, 2, 3], [1.0, 2.5, 5.0], ["a", "b", "c"]], dtype=(np.int64, np.float64, str)
    )

    # works when file does not exist
    t.write(tmpfile, format=pandas_fmt)

    # fails when cannot overwrite
    with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
        t.write(tmpfile, format=pandas_fmt, overwrite=False)

    # passes when it can
    t.write(tmpfile, format=pandas_fmt, overwrite=True)


def test_read_fixed_width_format():
    """Test reading with pandas read_fwf()"""
    tbl = """\
    a   b   c
    1  2.0  a
    2  3.0  b"""
    buf = StringIO()
    buf.write(tbl)

    # Explicitly provide converters to avoid casting 'a' to int32.
    # See https://github.com/astropy/astropy/issues/8682
    t = Table.read(
        tbl,
        format="ascii",
        guess=False,
        converters={"a": [ascii.convert_numpy(np.int64)]},
    )

    buf.seek(0)
    t2 = Table.read(buf, format="pandas.fwf")

    assert t.colnames == t2.colnames
    assert np.all(t == t2)


def test_write_with_mixins():
    """Writing a table with mixins just drops them via to_pandas()

    This also tests passing a kwarg to pandas read and write.
    """
    sc = SkyCoord([1, 2], [3, 4], unit="deg")
    q = [5, 6] * u.m
    qt = QTable([[1, 2], q, sc], names=["i", "q", "sc"])

    buf = StringIO()
    qt.write(buf, format="pandas.csv", sep=" ")
    exp = ["i q sc.ra sc.dec", "1 5.0 1.0 3.0", "2 6.0 2.0 4.0"]
    assert buf.getvalue().splitlines() == exp

    # Read it back
    buf.seek(0)
    qt2 = Table.read(buf, format="pandas.csv", sep=" ")
    # Explicitly provide converters to avoid casting 'i' to int32.
    # See https://github.com/astropy/astropy/issues/8682
    exp_t = ascii.read(exp, converters={"i": [ascii.convert_numpy(np.int64)]})
    assert qt2.colnames == exp_t.colnames
    assert np.all(qt2 == exp_t)


def test_excel_multi_sheet_error(tmp_path):
    """Testing reading of excel files with multiple sheets when no sheet_name is specified."""
    import pandas as pd

    if not HAS_OPENPYXL:
        pytest.skip("Missing openpyxl")
    file = tmp_path / "multi.xlsx"

    # create multi-sheet file
    with pd.ExcelWriter(file, engine="openpyxl") as writer:
        pd.DataFrame({"a": [1]}).to_excel(writer, sheet_name="s1")
        pd.DataFrame({"a": [2]}).to_excel(writer, sheet_name="s2")

    with pytest.raises(ValueError, match="Multiple sheets"):
        Table.read(file, format="pandas.excel", sheet_name=None)
