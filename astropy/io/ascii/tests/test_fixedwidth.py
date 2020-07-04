# Licensed under a 3-clause BSD style license - see LICENSE.rst


from io import StringIO

import numpy as np
import pytest

from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError
from .common import (assert_equal, assert_almost_equal)


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


def test_read_detect_col_starts_or_ends():
    """Table with no delimiter with only column start or end values specified"""
    table = """
#1       9        19                <== Column start indexes
#|       |         |                <== Column start positions
#<------><--------><------------->  <== Inferred column positions
  John   555- 1234 192.168.1.10
  Mary   555- 2134 192.168.1.123
   Bob   555- 4527  192.168.1.9
   Bill  555-9875  192.255.255.255
"""
    for kwargs in ({'col_starts': (1, 9, 19)},
                   {'col_ends': (8, 18, 33)}):
        dat = ascii.read(table,
                         Reader=ascii.FixedWidthNoHeader,
                         names=('Name', 'Phone', 'TCP'),
                         **kwargs)
        assert_equal(tuple(dat.dtype.names), ('Name', 'Phone', 'TCP'))
        assert_equal(dat[0][1], "555- 1234")
        assert_equal(dat[1][0], "Mary")
        assert_equal(dat[1][2], "192.168.1.123")
        assert_equal(dat[3][2], "192.255.255.255")


table = """\
| Col1 |  Col2     |  Col3     |  Col4     |
| 1.2  | "hello"   |  1        |  a        |
| 2.4  | 's worlds |         2 |         2 |
"""
dat = ascii.read(table, Reader=ascii.FixedWidth)


def test_read_header_col_onesided():
    """Table with only start or end positions of columns demarked by header"""
    table = "\n".join([' Col1   Col2         Col3        Col4      ',
                       '  1.2    "hello"      1           a         ',
                       "  2.4    's worlds           2           2  "])

    tab = ascii.read(table, Reader=ascii.FixedWidth, col_starts='from_header')
    assert_equal(tab.colnames, dat.colnames)
    for c in dat.colnames:
        if dat[c].dtype.kind == 'f':
            assert np.allclose(dat[c], tab[c])
        else:
            assert np.array_equal(dat[c], tab[c])

    tab = ascii.read(table, Reader=ascii.FixedWidth, col_ends='from_header')
    assert_equal(tab.colnames, dat.colnames)
    assert np.allclose(dat['Col1'], tab['Col1'])

    # Compare col_starts set to col_ends to both to none.
    table = "\n".join(["Col01  Col02   Col3 ",
                       "    100  1.2345XA",
                       "    200123.4e+12A03b"])

    tab = ascii.read(table, format='fixed_width', col_starts='from_header')
    assert np.array_equal(tab['Col01'], np.arange(1, 3) * 100)
    assert np.allclose(tab['Col02'], np.array([1.2345, 1234]))
    assert np.array_equal(tab['Col3'], np.array(['XA', '2A03b']))

    tab = ascii.read(table, format='fixed_width', col_ends='from_header')
    assert np.array_equal(tab['Col01'], np.arange(1, 3))
    assert np.array_equal(tab['Col02'], np.array(['00  1.2', '00123.4']))
    assert np.array_equal(tab['Col3'], np.array(['345XA', 'e+12A03']))

    tab = ascii.read(table, format='fixed_width', col_starts='from_header', col_ends='from_header')
    assert np.array_equal(tab['Col01'], np.arange(1, 3))
    assert np.allclose(tab['Col02'], np.array([1.2, 123.4]))
    assert np.array_equal(tab['Col3'], np.array(['XA', '2A03']))

    tab = ascii.read(table, format='fixed_width', delimiter=' ')
    assert np.array_equal(tab['Col01'], np.arange(1, 3))
    assert np.allclose(tab['Col02'], np.array([1.2, 123.4]))
    assert np.array_equal(tab['Col3'], np.array(['XA', '2A03']))

    # Note that delimiter=None is split differently than blank!
    tab = ascii.read(table, format='fixed_width', delimiter=None)
    assert np.array_equal(tab['Col01'], np.arange(1, 3))
    assert np.array_equal(tab['Col02'], np.array(['0  1.', '0123.']))
    assert np.array_equal(tab['Col3'], np.array(['345X', 'e+12']))


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
        ascii.read(table, Reader=ascii.FixedWidthTwoLine,
                   delimiter='|', guess=False)
    assert 'Position line should only contain delimiters and one other character' in str(
        excinfo.value)


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
        ascii.read(table, Reader=ascii.FixedWidthTwoLine,
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


def test_fixedwidthnoheader_splitting():
    """Test fix in #8511 where data_start is being ignored"""
    tbl = """\
AAA y z
1 2 3
4 5 6
7 8 9
"""
    names = ['a', 'b', 'c']
    dat = ascii.read(tbl, data_start=1, data_end=3,
                     delimiter=' ', names=names,
                     format='fixed_width_no_header')
    assert dat.colnames == names
    assert np.all(dat['a'] == [1, 4])
    assert np.all(dat['b'] == [2, 5])
    assert np.all(dat['c'] == [3, 6])
