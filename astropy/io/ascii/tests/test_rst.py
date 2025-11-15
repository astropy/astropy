# Licensed under a 3-clause BSD style license - see LICENSE.rst

from io import StringIO

from astropy.io import ascii

from .common import assert_almost_equal, assert_equal


def assert_equal_splitlines(arg1, arg2):
    assert_equal(arg1.splitlines(), arg2.splitlines())


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
    reader = ascii.get_reader(Reader=ascii.RST)
    dat = reader.read(table)
    assert_equal(dat.colnames, ["Col1", "Col2"])
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], '"hello"')
    assert_equal(dat[1][1], "'s worlds")


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
    reader = ascii.get_reader(Reader=ascii.RST, names=("name1", "name2"))
    dat = reader.read(table)
    assert_equal(dat.colnames, ["name1", "name2"])
    assert_almost_equal(dat[1][0], 2.4)


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
        Reader=ascii.RST,
        names=("name1", "name2", "name3"),
        include_names=("name1", "name3"),
    )
    dat = reader.read(table)
    assert_equal(dat.colnames, ["name1", "name3"])
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], 3)


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
    reader = ascii.get_reader(Reader=ascii.RST, exclude_names=("Col1",))
    dat = reader.read(table)
    assert_equal(dat.colnames, ["Col2"])
    assert_equal(dat[1][0], "'s worlds")


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
    reader = ascii.get_reader(Reader=ascii.RST)
    dat = reader.read(table)
    assert_equal(dat[0][2], "Hello")
    assert_equal(dat[1][2], "Worlds")


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
    reader = ascii.get_reader(Reader=ascii.RST)
    dat = reader.read(table)
    assert_equal(dat.colnames[-1], "Col3Long")


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
    reader = ascii.get_reader(Reader=ascii.RST)
    dat = reader.read(table)
    assert_equal(dat.colnames, ["Col1", "Col2", "Col3"])
    assert_equal(dat[0][2], "foo")
    assert_equal(dat[1][0], 1)


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

    reader = ascii.get_reader(Reader=ascii.RST)
    dat = reader.read(table)
    assert_equal(dat.colnames, ["Col1", "Col2", "Col3"])
    assert_equal(dat[0][2], "foo")
    assert_equal(dat[1][0], 1)


table = """\
====== =========== ============ ===========
  Col1    Col2        Col3        Col4
====== =========== ============ ===========
  1.2    "hello"      1           a
  2.4   's worlds          2           2
====== =========== ============ ===========
"""
dat = ascii.read(table, Reader=ascii.RST)


def test_write_normal():
    """Write a table as a normal SimpleRST Table"""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.RST)
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


def test_write_with_header_rows_two_lines():
    """Write a table with two header rows (name and unit)"""
    from astropy.table import QTable
    import astropy.units as u

    tbl = QTable({"wave": [350, 950] * u.nm, "response": [0.7, 1.2] * u.count})
    out = StringIO()
    ascii.write(tbl, out, format="rst", header_rows=["name", "unit"])
    result = out.getvalue()

    # Check that we have the expected structure
    lines = result.strip().split("\n")
    assert len(lines) == 7, f"Expected 7 lines, got {len(lines)}"

    # Line 0: separator
    assert "=" in lines[0]
    # Line 1: column names
    assert "wave" in lines[1] and "response" in lines[1]
    # Line 2: units
    assert "nm" in lines[2] and "ct" in lines[2]
    # Line 3: separator
    assert "=" in lines[3]
    # Lines 4-5: data
    assert "350" in lines[4] or "350.0" in lines[4]
    assert "950" in lines[5] or "950.0" in lines[5]
    # Line 6: separator
    assert "=" in lines[6]


def test_write_with_header_rows_three_lines():
    """Write a table with three header rows (name, unit, format)"""
    from astropy.table import Table

    tbl = Table(
        {
            "col1": [1, 2, 3],
            "col2": [4.5, 5.6, 6.7],
            "col3": ["a", "b", "c"],
        }
    )
    tbl["col2"].info.format = ".2f"

    out = StringIO()
    ascii.write(tbl, out, format="rst", header_rows=["name", "format"])
    result = out.getvalue()

    lines = result.strip().split("\n")
    # Expected: separator, name row, format row, separator, 3 data rows, separator = 8 lines
    assert len(lines) == 8, f"Expected 8 lines, got {len(lines)}"

    # Check structure
    assert "=" in lines[0]  # top separator
    assert "col1" in lines[1] and "col2" in lines[1] and "col3" in lines[1]  # names
    assert ".2f" in lines[2]  # format row
    assert "=" in lines[3]  # middle separator
    assert "=" in lines[7]  # bottom separator


def test_read_with_header_rows():
    """Read a table with multiple header rows"""
    table = """\
===== ======== ====
 Col1     Col2 Col3
    m      m/s
===== ======== ====
  1.2      2.3  foo
  2.4      4.5  bar
===== ======== ====
"""
    dat = ascii.read(table, format="rst", header_rows=["name", "unit"])
    assert_equal(dat.colnames, ["Col1", "Col2", "Col3"])
    assert_almost_equal(dat[0][0], 1.2)
    assert_equal(dat[0][2], "foo")
    # Check units were read
    if hasattr(dat["Col1"], "unit"):
        assert str(dat["Col1"].unit) == "m"
        assert str(dat["Col2"].unit) == "m / s"
