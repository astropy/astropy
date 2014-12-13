try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from ...tests.helper import pytest
from ...table import QTable, col_setattr, col_getattr
from ... import units as u

# ISSUES
# - Overlap with SkyCoord.name => SkyCoord.frame.name => 'icrs' by default.
#   This gets clobbered by directly setting `frame` attribute.

def test_attributes(mixin_cols):
    """
    Required attributes for a column can be set.
    """
    m = mixin_cols['m']
    col_setattr(m, 'name', 'a')
    assert col_getattr(m, 'name') == 'a'

    col_setattr(m, 'description', 'a')
    assert col_getattr(m, 'description') == 'a'

    if not isinstance(m, u.Quantity):
        col_setattr(m, 'unit', u.m)
    assert col_getattr(m, 'unit') is u.m

    col_setattr(m, 'format', 'a')
    assert col_getattr(m, 'format') == 'a'

    col_setattr(m, 'meta', {'a': 1})
    assert col_getattr(m, 'meta') == {'a': 1}

    with pytest.raises(AttributeError):
        col_setattr(m, 'bad_attr', 1)

    with pytest.raises(AttributeError):
        col_getattr(m, 'bad_attr')


def check_mixin_type(table, table_col, in_col):
    if isinstance(in_col, u.Quantity) and type(table) is not QTable:
        assert type(table_col) is table.ColumnClass
    else:
        assert type(table_col) is type(in_col)


def test_make_table(table_types, mixin_cols):
    """
    Make a table with the columns in mixin_cols, which is an ordered dict of
    three cols: 'a' and 'b' are table_types.Column type, and 'm' is a mixin.
    """
    t = table_types.Table(mixin_cols)
    check_mixin_type(t, t['m'], mixin_cols['m'])

    cols = list(mixin_cols.values())
    t = table_types.Table(cols, names=('a', 'b', 'c', 'm'))
    check_mixin_type(t, t['m'], mixin_cols['m'])

    t = table_types.Table(cols)
    check_mixin_type(t, t['col3'], mixin_cols['m'])


def test_io_ascii_write(mixin_cols):
    """
    Test that table with mixin column can be written by io.ascii for
    every pure Python writer.  No validation of the output is done,
    this just confirms no exceptions.
    """
    from ...io.ascii.connect import _get_connectors_table
    t = QTable(mixin_cols)
    for fmt in _get_connectors_table():
        if fmt['Write'] and '.fast_' not in fmt['Format']:
            out = StringIO()
            t.write(out, format=fmt['Format'])
