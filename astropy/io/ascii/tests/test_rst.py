# Licensed under a 3-clause BSD style license - see LICENSE.rst

from io import StringIO

import numpy as np
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable

from .common import assert_equal_splitlines


def test_read_normal():
    """Normal SimpleRST Table"""
    table = """
# comment (with blank line above)
======= =========
   Col1      Col2
======= =========
   1.2    "hello"
   2.4  's worlds
======= =========
"""
    reader = ascii.get_reader(reader_cls=ascii.RST)
    dat = reader.read(table)
    assert dat.colnames == ["Col1", "Col2"]
    assert_allclose(dat[1][0], 2.4)
    assert dat[0][1] == '"hello"'
    assert dat[1][1] == "'s worlds"


def test_read_normal_names():
    """Normal SimpleRST Table with provided column names"""
    table = """
# comment (with blank line above)
======= =========
   Col1      Col2
======= =========
   1.2    "hello"
   2.4  's worlds
======= =========
"""
    reader = ascii.get_reader(reader_cls=ascii.RST, names=("name1", "name2"))
    dat = reader.read(table)
    assert dat.colnames == ["name1", "name2"]
    assert_allclose(dat[1][0], 2.4)


def test_read_normal_names_include():
    """Normal SimpleRST Table with provided column names"""
    table = """
# comment (with blank line above)
=======  ========== ======
   Col1     Col2      Col3
=======  ========== ======
   1.2     "hello"       3
   2.4    's worlds      7
=======  ========== ======
"""
    reader = ascii.get_reader(
        reader_cls=ascii.RST,
        names=("name1", "name2", "name3"),
        include_names=("name1", "name3"),
    )
    dat = reader.read(table)
    assert dat.colnames == ["name1", "name3"]
    assert_allclose(dat[1][0], 2.4)
    assert dat[0][1] == 3


def test_read_normal_exclude():
    """Nice, typical SimpleRST table with col name excluded"""
    table = """
======= ==========
  Col1     Col2
======= ==========
  1.2     "hello"
  2.4    's worlds
======= ==========
"""
    reader = ascii.get_reader(reader_cls=ascii.RST, exclude_names=("Col1",))
    dat = reader.read(table)
    assert dat.colnames == ["Col2"]
    assert dat[1][0] == "'s worlds"


def test_read_unbounded_right_column():
    """The right hand column should be allowed to overflow"""
    table = """
# comment (with blank line above)
===== ===== ====
 Col1  Col2 Col3
===== ===== ====
 1.2    2    Hello
 2.4     4   Worlds
===== ===== ====
"""
    reader = ascii.get_reader(reader_cls=ascii.RST)
    dat = reader.read(table)
    assert dat[0][2] == "Hello"
    assert dat[1][2] == "Worlds"


def test_read_unbounded_right_column_header():
    """The right hand column should be allowed to overflow"""
    table = """
# comment (with blank line above)
===== ===== ====
 Col1  Col2 Col3Long
===== ===== ====
 1.2    2    Hello
 2.4     4   Worlds
===== ===== ====
"""
    reader = ascii.get_reader(reader_cls=ascii.RST)
    dat = reader.read(table)
    assert dat.colnames[-1] == "Col3Long"


def test_read_right_indented_table():
    """We should be able to read right indented tables correctly"""
    table = """
# comment (with blank line above)
   ==== ==== ====
   Col1 Col2 Col3
   ==== ==== ====
    3    3.4  foo
    1    4.5  bar
   ==== ==== ====
"""
    reader = ascii.get_reader(reader_cls=ascii.RST)
    dat = reader.read(table)
    assert dat.colnames == ["Col1", "Col2", "Col3"]
    assert dat[0][2] == "foo"
    assert dat[1][0] == 1


def test_trailing_spaces_in_row_definition():
    """Trailing spaces in the row definition column shouldn't matter"""
    table = (
        "\n"
        "# comment (with blank line above)\n"
        "   ==== ==== ====    \n"
        "   Col1 Col2 Col3\n"
        "   ==== ==== ====  \n"
        "    3    3.4  foo\n"
        "    1    4.5  bar\n"
        "   ==== ==== ====  \n"
    )
    # make sure no one accidentally deletes the trailing whitespaces in the
    # table.
    assert len(table) == 151

    reader = ascii.get_reader(reader_cls=ascii.RST)
    dat = reader.read(table)
    assert dat.colnames == ["Col1", "Col2", "Col3"]
    assert dat[0][2] == "foo"
    assert dat[1][0] == 1


table = """\
====== =========== ============ ===========
  Col1    Col2        Col3        Col4
====== =========== ============ ===========
  1.2    "hello"      1           a
  2.4   's worlds          2           2
====== =========== ============ ===========
"""
dat = ascii.read(table, format="rst")


def test_write_normal():
    """Write a table as a normal SimpleRST Table"""
    out = StringIO()
    ascii.write(dat, out, format="rst")
    assert_equal_splitlines(
        out.getvalue(),
        """\
==== ========= ==== ====
Col1      Col2 Col3 Col4
==== ========= ==== ====
 1.2   "hello"    1    a
 2.4 's worlds    2    2
==== ========= ==== ====
""",
    )


def test_rst_with_header_rows():
    """Round-trip a table with header_rows specified"""
    lines = [
        "======= ======== ====",
        "   wave response ints",
        "     nm       ct     ",
        "float64  float32 int8",
        "======= ======== ====",
        "  350.0      1.0    1",
        "  950.0      2.0    2",
        "======= ======== ====",
    ]
    tbl = QTable.read(lines, format="ascii.rst", header_rows=["name", "unit", "dtype"])
    assert tbl["wave"].unit == u.nm
    assert tbl["response"].unit == u.ct
    assert tbl["wave"].dtype == np.float64
    assert tbl["response"].dtype == np.float32
    assert tbl["ints"].dtype == np.int8

    out = StringIO()
    tbl.write(out, format="ascii.rst", header_rows=["name", "unit", "dtype"])
    assert out.getvalue().splitlines() == lines
