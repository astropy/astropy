# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from ... import ascii
from ..core import InconsistentTableError
from .common import (assert_equal, assert_almost_equal,
                     setup_function, teardown_function)
from ....tests.helper import pytest


def assert_equal_splitlines(arg1, arg2):
    assert_equal(arg1.splitlines(), arg2.splitlines())


def test_read_normal():
    """Nice, typical fixed format table"""
    table = """
# comment (with blank line above)
|  Col1  |  Col2   |
|  1.2   | "hello" |
|  2.4   |'s worlds|
"""
    reader = ascii.get_reader(Reader=ascii.FixedWidth)
    dat = reader.read(table)
    assert_equal(dat.colnames, ['Col1', 'Col2'])
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], '"hello"')
    assert_equal(dat[1][1], "'s worlds")


def test_read_normal_names():
    """Nice, typical fixed format table with col names provided"""
    table = """
# comment (with blank line above)
|  Col1  |  Col2   |
|  1.2   | "hello" |
|  2.4   |'s worlds|
"""
    reader = ascii.get_reader(Reader=ascii.FixedWidth,
                                   names=('name1', 'name2'))
    dat = reader.read(table)
    assert_equal(dat.colnames, ['name1', 'name2'])
    assert_almost_equal(dat[1][0], 2.4)


def test_read_normal_names_include():
    """Nice, typical fixed format table with col names provided"""
    table = """
# comment (with blank line above)
|  Col1  |  Col2   |  Col3 |
|  1.2   | "hello" |     3 |
|  2.4   |'s worlds|     7 |
"""
    reader = ascii.get_reader(Reader=ascii.FixedWidth,
                                   names=('name1', 'name2', 'name3'),
                                   include_names=('name1', 'name3'))
    dat = reader.read(table)
    assert_equal(dat.colnames, ['name1', 'name3'])
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], 3)


def test_read_normal_exclude():
    """Nice, typical fixed format table with col name excluded"""
    table = """
# comment (with blank line above)
|  Col1  |  Col2   |
|  1.2   | "hello" |
|  2.4   |'s worlds|
"""
    reader = ascii.get_reader(Reader=ascii.FixedWidth,
                                   exclude_names=('Col1',))
    dat = reader.read(table)
    assert_equal(dat.colnames, ['Col2'])
    assert_equal(dat[1][0], "'s worlds")


def test_read_weird():
    """Weird input table with data values chopped by col extent """
    table = """
  Col1  |  Col2 |
  1.2       "hello"
  2.4   sdf's worlds
"""
    reader = ascii.get_reader(Reader=ascii.FixedWidth)
    dat = reader.read(table)
    assert_equal(dat.colnames, ['Col1', 'Col2'])
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], '"hel')
    assert_equal(dat[1][1], "df's wo")


def test_read_double():
    """Table with double delimiters"""
    table = """
|| Name ||   Phone ||         TCP||
|  John  | 555-1234 |192.168.1.10X|
|  Mary  | 555-2134 |192.168.1.12X|
|   Bob  | 555-4527 | 192.168.1.9X|
"""
    dat = ascii.read(table, Reader=ascii.FixedWidth, guess=False)
    assert_equal(tuple(dat.dtype.names), ('Name', 'Phone', 'TCP'))
    assert_equal(dat[1][0], "Mary")
    assert_equal(dat[0][1], "555-1234")
    assert_equal(dat[2][2], "192.168.1.9")


def test_read_space_delimiter():
    """Table with space delimiter"""
    table = """
 Name  --Phone-    ----TCP-----
 John  555-1234    192.168.1.10
 Mary  555-2134    192.168.1.12
  Bob  555-4527     192.168.1.9
"""
    dat = ascii.read(table, Reader=ascii.FixedWidth, guess=False,
                          delimiter=' ')
    assert_equal(tuple(dat.dtype.names), ('Name', '--Phone-', '----TCP-----'))
    assert_equal(dat[1][0], "Mary")
    assert_equal(dat[0][1], "555-1234")
    assert_equal(dat[2][2], "192.168.1.9")


def test_read_no_header_autocolumn():
    """Table with no header row and auto-column naming"""
    table = """
|  John  | 555-1234 |192.168.1.10|
|  Mary  | 555-2134 |192.168.1.12|
|   Bob  | 555-4527 | 192.168.1.9|
"""
    dat = ascii.read(table, Reader=ascii.FixedWidth, guess=False,
                          header_start=None, data_start=0)
    assert_equal(tuple(dat.dtype.names), ('col1', 'col2', 'col3'))
    assert_equal(dat[1][0], "Mary")
    assert_equal(dat[0][1], "555-1234")
    assert_equal(dat[2][2], "192.168.1.9")


def test_read_no_header_names():
    """Table with no header row and with col names provided.  Second
    and third rows also have hanging spaces after final |."""
    table = """
|  John  | 555-1234 |192.168.1.10|
|  Mary  | 555-2134 |192.168.1.12|
|   Bob  | 555-4527 | 192.168.1.9|
"""
    dat = ascii.read(table, Reader=ascii.FixedWidth, guess=False,
                          header_start=None, data_start=0,
                          names=('Name', 'Phone', 'TCP'))
    assert_equal(tuple(dat.dtype.names), ('Name', 'Phone', 'TCP'))
    assert_equal(dat[1][0], "Mary")
    assert_equal(dat[0][1], "555-1234")
    assert_equal(dat[2][2], "192.168.1.9")


def test_read_no_header_autocolumn_NoHeader():
    """Table with no header row and auto-column naming"""
    table = """
|  John  | 555-1234 |192.168.1.10|
|  Mary  | 555-2134 |192.168.1.12|
|   Bob  | 555-4527 | 192.168.1.9|
"""
    dat = ascii.read(table, Reader=ascii.FixedWidthNoHeader)
    assert_equal(tuple(dat.dtype.names), ('col1', 'col2', 'col3'))
    assert_equal(dat[1][0], "Mary")
    assert_equal(dat[0][1], "555-1234")
    assert_equal(dat[2][2], "192.168.1.9")


def test_read_no_header_names_NoHeader():
    """Table with no header row and with col names provided.  Second
    and third rows also have hanging spaces after final |."""
    table = """
|  John  | 555-1234 |192.168.1.10|
|  Mary  | 555-2134 |192.168.1.12|
|   Bob  | 555-4527 | 192.168.1.9|
"""
    dat = ascii.read(table, Reader=ascii.FixedWidthNoHeader,
                          names=('Name', 'Phone', 'TCP'))
    assert_equal(tuple(dat.dtype.names), ('Name', 'Phone', 'TCP'))
    assert_equal(dat[1][0], "Mary")
    assert_equal(dat[0][1], "555-1234")
    assert_equal(dat[2][2], "192.168.1.9")


def test_read_col_starts():
    """Table with no delimiter with column start and end values specified."""
    table = """
#    5   9     17  18      28
#    |   |       ||         |
  John   555- 1234 192.168.1.10
  Mary   555- 2134 192.168.1.12
   Bob   555- 4527  192.168.1.9
"""
    dat = ascii.read(table, Reader=ascii.FixedWidthNoHeader,
                          names=('Name', 'Phone', 'TCP'),
                          col_starts=(0, 9, 18),
                          col_ends=(5, 17, 28),
                          )
    assert_equal(tuple(dat.dtype.names), ('Name', 'Phone', 'TCP'))
    assert_equal(dat[0][1], "555- 1234")
    assert_equal(dat[1][0], "Mary")
    assert_equal(dat[1][2], "192.168.1.")
    assert_equal(dat[2][2], "192.168.1")  # col_end=28 cuts this column off


table = """\
| Col1 |  Col2     |  Col3     |  Col4     |
| 1.2  | "hello"   |  1        |  a        |
| 2.4  | 's worlds |         2 |         2 |
"""
dat = ascii.read(table, Reader=ascii.FixedWidth)


def test_write_normal():
    """Write a table as a normal fixed width table."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidth)
    assert_equal_splitlines(out.getvalue(), """\
| Col1 |      Col2 | Col3 | Col4 |
|  1.2 |   "hello" |    1 |    a |
|  2.4 | 's worlds |    2 |    2 |
""")


def test_write_fill_values():
    """Write a table as a normal fixed width table."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidth,
                     fill_values=('a', 'N/A'))
    assert_equal_splitlines(out.getvalue(), """\
| Col1 |      Col2 | Col3 | Col4 |
|  1.2 |   "hello" |    1 |  N/A |
|  2.4 | 's worlds |    2 |    2 |
""")


def test_write_no_pad():
    """Write a table as a fixed width table with no padding."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidth,
                     delimiter_pad=None)
    assert_equal_splitlines(out.getvalue(), """\
|Col1|     Col2|Col3|Col4|
| 1.2|  "hello"|   1|   a|
| 2.4|'s worlds|   2|   2|
""")


def test_write_no_bookend():
    """Write a table as a fixed width table with no bookend."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidth, bookend=False)
    assert_equal_splitlines(out.getvalue(), """\
Col1 |      Col2 | Col3 | Col4
 1.2 |   "hello" |    1 |    a
 2.4 | 's worlds |    2 |    2
""")


def test_write_no_delimiter():
    """Write a table as a fixed width table with no delimiter."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidth, bookend=False,
                     delimiter=None)
    assert_equal_splitlines(out.getvalue(), """\
Col1       Col2  Col3  Col4
 1.2    "hello"     1     a
 2.4  's worlds     2     2
""")


def test_write_noheader_normal():
    """Write a table as a normal fixed width table."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidthNoHeader)
    assert_equal_splitlines(out.getvalue(), """\
| 1.2 |   "hello" | 1 | a |
| 2.4 | 's worlds | 2 | 2 |
""")


def test_write_noheader_no_pad():
    """Write a table as a fixed width table with no padding."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidthNoHeader,
                     delimiter_pad=None)
    assert_equal_splitlines(out.getvalue(), """\
|1.2|  "hello"|1|a|
|2.4|'s worlds|2|2|
""")


def test_write_noheader_no_bookend():
    """Write a table as a fixed width table with no bookend."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidthNoHeader,
                     bookend=False)
    assert_equal_splitlines(out.getvalue(), """\
1.2 |   "hello" | 1 | a
2.4 | 's worlds | 2 | 2
""")


def test_write_noheader_no_delimiter():
    """Write a table as a fixed width table with no delimiter."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidthNoHeader, bookend=False,
                     delimiter=None)
    assert_equal_splitlines(out.getvalue(), """\
1.2    "hello"  1  a
2.4  's worlds  2  2
""")


def test_write_formats():
    """Write a table as a fixed width table with no delimiter."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidth,
                     formats={'Col1': '%-8.3f', 'Col2': '%-15s'})
    assert_equal_splitlines(out.getvalue(), """\
|     Col1 |            Col2 | Col3 | Col4 |
| 1.200    | "hello"         |    1 |    a |
| 2.400    | 's worlds       |    2 |    2 |
""")


def test_read_twoline_normal():
    """Typical fixed format table with two header lines (with some cruft
    thrown in to test column positioning"""
    table = """
  Col1    Col2
  ----  ---------
   1.2xx"hello"
  2.4   's worlds
"""
    dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine)
    assert_equal(dat.dtype.names, ('Col1', 'Col2'))
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], '"hello"')
    assert_equal(dat[1][1], "'s worlds")


def test_read_twoline_ReST():
    """Read restructured text table"""
    table = """
======= ===========
  Col1    Col2
======= ===========
  1.2   "hello"
  2.4   's worlds
======= ===========
"""
    dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine,
                          header_start=1, position_line=2, data_end=-1)
    assert_equal(dat.dtype.names, ('Col1', 'Col2'))
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], '"hello"')
    assert_equal(dat[1][1], "'s worlds")


def test_read_twoline_human():
    """Read text table designed for humans and test having position line
    before the header line"""
    table = """
+------+----------+
| Col1 |   Col2   |
+------|----------+
|  1.2 | "hello"  |
|  2.4 | 's worlds|
+------+----------+
"""
    dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine,
                          delimiter='+',
                          header_start=1, position_line=0,
                          data_start=3, data_end=-1)
    assert_equal(dat.dtype.names, ('Col1', 'Col2'))
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], '"hello"')
    assert_equal(dat[1][1], "'s worlds")


def test_read_twoline_fail():
    """Test failure if too many different character are on position line.

    The position line shall consist of only one character in addition to
    the delimiter.
    """
    table = """
| Col1 |   Col2   |
|------|==========|
|  1.2 | "hello"  |
|  2.4 | 's worlds|
"""
    with pytest.raises(InconsistentTableError) as excinfo:
        dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine,
                         delimiter='|', guess=False)
    assert 'Position line should only contain delimiters and one other character' in str(excinfo.value)


def test_read_twoline_wrong_marker():
    '''Test failure when position line uses characters prone to ambiguity

    Characters in position line must be part an allowed set because
    normal letters or numbers will lead to ambiguous tables.
    '''
    table = """
| Col1 |   Col2   |
|aaaaaa|aaaaaaaaaa|
|  1.2 | "hello"  |
|  2.4 | 's worlds|
"""
    with pytest.raises(InconsistentTableError) as excinfo:
        dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine,
                         delimiter='|', guess=False)
    assert 'Characters in position line must be part' in str(excinfo.value)


def test_write_twoline_normal():
    """Write a table as a normal fixed width table."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidthTwoLine)
    assert_equal_splitlines(out.getvalue(), """\
Col1      Col2 Col3 Col4
---- --------- ---- ----
 1.2   "hello"    1    a
 2.4 's worlds    2    2
""")


def test_write_twoline_no_pad():
    """Write a table as a fixed width table with no padding."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidthTwoLine,
                     delimiter_pad=' ', position_char='=')
    assert_equal_splitlines(out.getvalue(), """\
Col1        Col2   Col3   Col4
====   =========   ====   ====
 1.2     "hello"      1      a
 2.4   's worlds      2      2
""")


def test_write_twoline_no_bookend():
    """Write a table as a fixed width table with no bookend."""
    out = StringIO()
    ascii.write(dat, out, Writer=ascii.FixedWidthTwoLine,
                     bookend=True, delimiter='|')
    assert_equal_splitlines(out.getvalue(), """\
|Col1|     Col2|Col3|Col4|
|----|---------|----|----|
| 1.2|  "hello"|   1|   a|
| 2.4|'s worlds|   2|   2|
""")
