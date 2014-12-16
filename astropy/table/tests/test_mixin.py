try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

try:
    import h5py
except ImportError:
    HAS_H5PY = False
else:
    HAS_H5PY = True

import numpy as np

from ...tests.helper import pytest
from ...table import QTable, col_setattr, col_getattr, join
from ... import units as u
from .conftest import MIXIN_COLS

# ISSUES / TODO
# - Test that slicing / indexing table gives right values and col attrs inherit
# - Test hstack, vstack, groups
# - Add column to table => makes copy and copies col attrs if existent

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

    # Make sure in_col got copied and creating table did not touch it
    assert col_getattr(in_col, 'name') is None


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


def test_io_ascii_write():
    """
    Test that table with mixin column can be written by io.ascii for
    every pure Python writer.  No validation of the output is done,
    this just confirms no exceptions.
    """
    from ...io.ascii.connect import _get_connectors_table
    t = QTable(MIXIN_COLS)
    for fmt in _get_connectors_table():
        if fmt['Write'] and '.fast_' not in fmt['Format']:
            out = StringIO()
            t.write(out, format=fmt['Format'])


def test_io_write_fail(mixin_cols):
    """
    Test that table with mixin column cannot be written by io.votable,
    io.fits, and io.misc.hdf5
    every pure Python writer.  No validation of the output is done,
    this just confirms no exceptions.
    """
    t = QTable(mixin_cols)
    for fmt in ('fits', 'votable', 'hdf5'):
        if fmt == 'hdf5' and not HAS_H5PY:
            continue
        out = StringIO()
        with pytest.raises(ValueError) as err:
            t.write(out, format=fmt)
        assert 'cannot write table with mixin column(s)' in str(err.value)


def test_join(table_types):
    """
    Join tables with mixin cols.  Use column "i" as proxy for what the
    result should be for each mixin.
    """
    t1 = table_types.Table()
    t1['a'] = table_types.Column(['a', 'b', 'b', 'c'])
    t1['i'] = table_types.Column([0, 1, 2, 3])
    for name, col in MIXIN_COLS.items():
        t1[name] = col

    t2 = table_types.Table(t1)
    t2['a'] = ['b', 'c', 'a', 'd']

    for join_type in ('inner', 'left'):
        t12 = join(t1, t2, keys='a', join_type=join_type)
        idx1 = t12['i_1']
        idx2 = t12['i_2']
        for name, col in MIXIN_COLS.items():
            name1 = name + '_1'
            name2 = name + '_2'
            if name == 'skycoord':
                assert np.all(t12[name1].ra == col[idx1].ra)
                assert np.all(t12[name1].dec == col[idx1].dec)
                assert np.all(t12[name2].ra == col[idx2].ra)
                assert np.all(t12[name2].dec == col[idx2].dec)
            elif name == 'quantity':
                if table_types.Table is QTable:
                    assert np.all(t12[name1].value == col[idx1].value)
                    assert np.all(t12[name2].value == col[idx2].value)
            else:
                assert np.all(t12[name1] == col[idx1])
                assert np.all(t12[name2] == col[idx2])

    for join_type in ('outer', 'right'):
        with pytest.raises(ValueError) as exc:
            t12 = join(t1, t2, keys='a', join_type=join_type)
        assert 'join requires masking column' in str(exc.value)

    with pytest.raises(ValueError) as exc:
        t12 = join(t1, t2, keys=['a', 'skycoord'])
    assert 'not allowed as a key column' in str(exc.value)
