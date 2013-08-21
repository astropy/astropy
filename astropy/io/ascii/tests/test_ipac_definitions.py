# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

from ..ui import read
from ..ipac import Ipac, IpacFormatError, IpacFormatErrorDBMS
from ....tests.helper import pytest
from ... import ascii
from ....table import Table

try:
    import StringIO as io
except ImportError:
    import io


DATA = '''
|   a  |   b   |
| char | char  |
ABBBBBBABBBBBBBA
'''


def test_ipac_default():
    # default should be ignore
    table = read(DATA, Reader=Ipac)
    assert table['a'][0] == 'BBBBBB'
    assert table['b'][0] == 'BBBBBBB'


def test_ipac_ignore():
    table = read(DATA, Reader=Ipac, definition='ignore')
    assert table['a'][0] == 'BBBBBB'
    assert table['b'][0] == 'BBBBBBB'


def test_ipac_left():
    table = read(DATA, Reader=Ipac, definition='left')
    assert table['a'][0] == 'BBBBBBA'
    assert table['b'][0] == 'BBBBBBBA'


def test_ipac_right():
    table = read(DATA, Reader=Ipac, definition='right')
    assert table['a'][0] == 'ABBBBBB'
    assert table['b'][0] == 'ABBBBBBB'


def test_too_long_colname_default():
    table = Table([[3]], names=['a1234567890123456789012345678901234567890'])
    out = io.StringIO()
    with pytest.raises(IpacFormatError):
        ascii.write(table, out, Writer=Ipac)


def test_too_long_colname_strict():
    table = Table([[3]], names=['a1234567890123456'])
    out = io.StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


def test_too_long_colname_notstrict():
    table = Table([[3]], names=['a1234567890123456789012345678901234567890'])
    out = io.StringIO()
    with pytest.raises(IpacFormatError):
        ascii.write(table, out, Writer=Ipac, DBMS=False)


@pytest.mark.parametrize(("strict_", "Err"), [(True, IpacFormatErrorDBMS), (False, IpacFormatError)])
def test_non_alfnum_colname(strict_, Err):
    table = Table([[3]], names=['a123456789 01234'])
    out = io.StringIO()
    with pytest.raises(Err):
        ascii.write(table, out, Writer=Ipac, DBMS=strict_)


def test_colname_starswithnumber_strict():
    table = Table([[3]], names=['a123456789 01234'])
    out = io.StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


def test_double_colname_strict():
    table = Table([[3], [1]], names=['DEC', 'dec'])
    out = io.StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


@pytest.mark.parametrize('colname', ['x', 'y', 'z', 'X', 'Y', 'Z'])
def test_reserved_colname_strict(colname):
    table = Table([['reg']], names=[colname])
    out = io.StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


def test_too_long_comment(recwarn):
    table = Table([[3]])
    table.meta['comments'] = ['a' * 79]
    out = io.StringIO()
    ascii.write(table, out, Writer=Ipac)
    w = recwarn.pop(UserWarning)
    assert 'Comment string > 78 characters was automatically wrapped.' == str(w.message)
    expected_out = """\
\\ aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
\\ a
|col0|
|long|
|    |
|null|
    3 
"""
    assert out.getvalue().splitlines() == expected_out.splitlines()

def test_out_with_nonstring_null():
    dirname = os.path.dirname(os.path.abspath(__file__))
    tab = Table.read(os.path.join(dirname, 't', 'withnullvalues.vot'))
    out = io.StringIO()
    ascii.write(tab, out, Writer=Ipac)
    expected_out = """\
|Prog|    Name|       Obs|   tos|Type|   Vel|n_Vel|  Flow|n_Flow| Fhigh|n_Fhigh|Conf|    RAJ2000|    DEJ2000|
|char|    char|      char|  long|char|double| char|  long|  char|  long|   char|char|       char|       char|
|    |        |   "Y:M:D"|     s|    |km / s|     |   MHz|      |   MHz|       |    |    "h:m:s"|    "d:m:s"|
| N/A|     N/A|       N/A|999999| N/A| 1e+20|    N|999999|     N|999999|      N| N/A|        N/A|        N/A|
 N032 AFGL2591 2003-12-06  16800  MAP   -5.5     L  80578      L 203407       L  6Cp 20:29:24.87 +40:11:19.8 
 N032 AFGL2591 2004-05-15   1140  MAP   -5.5     L  80578      L 203407       L  6Dp 20:29:24.87 +40:11:19.8 
 N032 AFGL2591 2004-05-15   7740  MAP   -5.5     L  80578      L 203407       L  6Dp 20:29:24.87 +40:11:19.8 
 N032 AFGL2591 2004-05-16  24060  MAP   -5.5     L  80578      L 203407       L  6Dp 20:29:24.87 +40:11:19.8 
 N032 AFGL2591 2004-05-17   3060  MAP   -5.5     L  80578      L 203407       L  6Dp 20:29:24.87 +40:11:19.8 
 PB3F AFGL2591 2005-12-25   4050  MAP    0.0     L  86610      L 230538       L  6Cq 20:29:24.80 +40:11:19.0 
 PB3F AFGL2591 2005-12-26   4050  MAP    0.0     L  86610      L 230538       L  6Cq 20:29:24.80 +40:11:19.0 
 P04A AFGL2591 2006-02-03   9000  MAP   -5.5     L  80578      L 203407       L  6Aq 20:29:24.87 +40:11:19.5 
 PB3F AFGL2591 2006-06-03    675  MAP    0.0     L  86610      L 230538       L  5Dq 20:29:24.80 +40:11:19.0 
 PB3F AFGL2591 2006-06-04   7425  MAP    0.0     L  86610      L 230538       L  5Dq 20:29:24.80 +40:11:19.0 
 PB4A AFGL2591 2007-02-19   9450  MAP   -5.5     L  80126      L  81030       L  6Aq 20:29:24.87 +40:11:19.5 
 PA4A AFGL2591 2007-03-04  10800  MAP   -5.5     L 203055      L 203739       L  6Aq 20:29:24.87 +40:11:19.5 
 PA4A AFGL2591 2007-03-13  10935  MAP   -5.5     L 203055      L 203739       L  6Bq 20:29:24.87 +40:11:19.5 
 PA4A AFGL2591 2007-03-15   9450  MAP   -5.5     L 203055      L 203739       L  6Bq 20:29:24.87 +40:11:19.5 """
    assert out.getvalue().splitlines() == expected_out.splitlines()
