# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import os

from ..ui import read
from ..ipac import Ipac, IpacFormatError, IpacFormatErrorDBMS
from ....tests.helper import pytest, catch_warnings
from ... import ascii
from ....table import Table

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO


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
    out = StringIO()
    with pytest.raises(IpacFormatError):
        ascii.write(table, out, Writer=Ipac)


def test_too_long_colname_strict():
    table = Table([[3]], names=['a1234567890123456'])
    out = StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


def test_too_long_colname_notstrict():
    table = Table([[3]], names=['a1234567890123456789012345678901234567890'])
    out = StringIO()
    with pytest.raises(IpacFormatError):
        ascii.write(table, out, Writer=Ipac, DBMS=False)


@pytest.mark.parametrize(("strict_", "Err"), [(True, IpacFormatErrorDBMS), (False, IpacFormatError)])
def test_non_alfnum_colname(strict_, Err):
    table = Table([[3]], names=['a123456789 01234'])
    out = StringIO()
    with pytest.raises(Err):
        ascii.write(table, out, Writer=Ipac, DBMS=strict_)


def test_colname_starswithnumber_strict():
    table = Table([[3]], names=['a123456789 01234'])
    out = StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


def test_double_colname_strict():
    table = Table([[3], [1]], names=['DEC', 'dec'])
    out = StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


@pytest.mark.parametrize('colname', ['x', 'y', 'z', 'X', 'Y', 'Z'])
def test_reserved_colname_strict(colname):
    table = Table([['reg']], names=[colname])
    out = StringIO()
    with pytest.raises(IpacFormatErrorDBMS):
        ascii.write(table, out, Writer=Ipac, DBMS=True)


def test_too_long_comment():
    with catch_warnings(UserWarning) as w:
        table = Table([[3]])
        table.meta['comments'] = ['a' * 79]
        out = StringIO()
        ascii.write(table, out, Writer=Ipac)
    w = w[0]
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
    assert out.getvalue().strip().splitlines() == expected_out.splitlines()

def test_out_with_nonstring_null():
    # unmasked tables don't have fill_values
    table = Table([[3]], masked=True)
    table['col0'].fill_value = -99999
    out = StringIO()
    ascii.write(table, out, Writer=Ipac)
    expected_out = """\
|  col0|
|  long|
|      |
|-99999|
      3
"""
    assert out.getvalue().strip().splitlines() == expected_out.splitlines()
