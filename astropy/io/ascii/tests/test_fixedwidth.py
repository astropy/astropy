# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import glob
import numpy as np
from ... import ascii as asciitable

io = asciitable.core.io

from .common import (raises,
                     assert_equal, assert_almost_equal, assert_true,
                     setup_function, teardown_function)


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
    reader = asciitable.get_reader(Reader=asciitable.FixedWidth)
    dat = reader.read(table)
    assert_equal(reader.header.colnames, ('Col1', 'Col2'))
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
    reader = asciitable.get_reader(Reader=asciitable.FixedWidth,
                                   names=('name1', 'name2'))
    dat = reader.read(table)
    assert_equal(reader.header.colnames, ('name1', 'name2'))
    assert_almost_equal(dat[1][0], 2.4)


def test_read_normal_names_include():
    """Nice, typical fixed format table with col names provided"""
    table = """
# comment (with blank line above)
|  Col1  |  Col2   |  Col3 |
|  1.2   | "hello" |     3 |
|  2.4   |'s worlds|     7 |
"""
    reader = asciitable.get_reader(Reader=asciitable.FixedWidth,
                                   names=('name1', 'name2', 'name3'),
                                   include_names=('name1', 'name3'))
    dat = reader.read(table)
    assert_equal(reader.header.colnames, ('name1', 'name3'))
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
    reader = asciitable.get_reader(Reader=asciitable.FixedWidth,
                                   exclude_names=('Col1',))
    dat = reader.read(table)
    assert_equal(reader.header.colnames, ('Col2',))
    assert_equal(dat[1][0], "'s worlds")


def test_read_weird():
    """Weird input table with data values chopped by col extent """
    table = """
  Col1  |  Col2 |
  1.2       "hello"
  2.4   sdf's worlds
"""
    reader = asciitable.get_reader(Reader=asciitable.FixedWidth)
    dat = reader.read(table)
    assert_equal(reader.header.colnames, ('Col1', 'Col2'))
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidth, guess=False)
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidth, guess=False,
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidth, guess=False,
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidth, guess=False,
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidthNoHeader)
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidthNoHeader,
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidthNoHeader,
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
dat = asciitable.read(table, Reader=asciitable.FixedWidth)


def test_write_normal():
    """Write a table as a normal fixed width table."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidth)
    assert_equal_splitlines(out.getvalue(), """\
| Col1 |      Col2 | Col3 | Col4 |
|  1.2 |   "hello" |    1 |    a |
|  2.4 | 's worlds |    2 |    2 |
""")


def test_write_no_pad():
    """Write a table as a fixed width table with no padding."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidth,
                     delimiter_pad=None)
    assert_equal_splitlines(out.getvalue(), """\
|Col1|     Col2|Col3|Col4|
| 1.2|  "hello"|   1|   a|
| 2.4|'s worlds|   2|   2|
""")


def test_write_no_bookend():
    """Write a table as a fixed width table with no bookend."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidth, bookend=False)
    assert_equal_splitlines(out.getvalue(), """\
Col1 |      Col2 | Col3 | Col4
 1.2 |   "hello" |    1 |    a
 2.4 | 's worlds |    2 |    2
""")


def test_write_no_delimiter():
    """Write a table as a fixed width table with no delimiter."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidth, bookend=False,
                     delimiter=None)
    assert_equal_splitlines(out.getvalue(), """\
Col1       Col2  Col3  Col4
 1.2    "hello"     1     a
 2.4  's worlds     2     2
""")


def test_write_noheader_normal():
    """Write a table as a normal fixed width table."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidthNoHeader)
    assert_equal_splitlines(out.getvalue(), """\
| 1.2 |   "hello" | 1 | a |
| 2.4 | 's worlds | 2 | 2 |
""")


def test_write_noheader_no_pad():
    """Write a table as a fixed width table with no padding."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidthNoHeader,
                     delimiter_pad=None)
    assert_equal_splitlines(out.getvalue(), """\
|1.2|  "hello"|1|a|
|2.4|'s worlds|2|2|
""")


def test_write_noheader_no_bookend():
    """Write a table as a fixed width table with no bookend."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidthNoHeader,
                     bookend=False)
    assert_equal_splitlines(out.getvalue(), """\
1.2 |   "hello" | 1 | a
2.4 | 's worlds | 2 | 2
""")


def test_write_noheader_no_delimiter():
    """Write a table as a fixed width table with no delimiter."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidthNoHeader, bookend=False,
                     delimiter=None)
    assert_equal_splitlines(out.getvalue(), """\
1.2    "hello"  1  a
2.4  's worlds  2  2
""")


def test_write_formats():
    """Write a table as a fixed width table with no delimiter."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidth,
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidthTwoLine)
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidthTwoLine,
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
    dat = asciitable.read(table, Reader=asciitable.FixedWidthTwoLine,
                          delimiter='+',
                          header_start=1, position_line=0,
                          data_start=3, data_end=-1)
    assert_equal(dat.dtype.names, ('Col1', 'Col2'))
    assert_almost_equal(dat[1][0], 2.4)
    assert_equal(dat[0][1], '"hello"')
    assert_equal(dat[1][1], "'s worlds")


def test_write_twoline_normal():
    """Write a table as a normal fixed width table."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidthTwoLine)
    assert_equal_splitlines(out.getvalue(), """\
Col1      Col2 Col3 Col4
---- --------- ---- ----
 1.2   "hello"    1    a
 2.4 's worlds    2    2
""")


def test_write_twoline_no_pad():
    """Write a table as a fixed width table with no padding."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidthTwoLine,
                     delimiter_pad=' ', position_char='=')
    assert_equal_splitlines(out.getvalue(), """\
Col1        Col2   Col3   Col4
====   =========   ====   ====
 1.2     "hello"      1      a
 2.4   's worlds      2      2
""")


def test_write_twoline_no_bookend():
    """Write a table as a fixed width table with no bookend."""
    out = io.StringIO()
    asciitable.write(dat, out, Writer=asciitable.FixedWidthTwoLine,
                     bookend=True, delimiter='|')
    assert_equal_splitlines(out.getvalue(), """\
|Col1|     Col2|Col3|Col4|
|----|---------|----|----|
| 1.2|  "hello"|   1|   a|
| 2.4|'s worlds|   2|   2|
""")
