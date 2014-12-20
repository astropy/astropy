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
from ...table import Table, QTable, join, hstack, vstack
from ..column import col_setattr, col_getattr
from ... import units as u
from ... import coordinates
from .. import table_helpers
from .conftest import MIXIN_COLS


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

    for name, col in MIXIN_COLS.items():
        col_setattr(t1[name], 'description', name)
        col_setattr(t2[name], 'description', name + '2')

    for join_type in ('inner', 'left'):
        t12 = join(t1, t2, keys='a', join_type=join_type)
        idx1 = t12['i_1']
        idx2 = t12['i_2']
        for name, col in MIXIN_COLS.items():
            name1 = name + '_1'
            name2 = name + '_2'
            assert_table_name_col_equal(t12, name1, col[idx1])
            assert_table_name_col_equal(t12, name2, col[idx2])
            assert col_getattr(t12[name1], 'description') == name
            assert col_getattr(t12[name2], 'description') == name + '2'

    for join_type in ('outer', 'right'):
        with pytest.raises(ValueError) as exc:
            t12 = join(t1, t2, keys='a', join_type=join_type)
        assert 'join requires masking column' in str(exc.value)

    with pytest.raises(ValueError) as exc:
        t12 = join(t1, t2, keys=['a', 'skycoord'])
    assert 'not allowed as a key column' in str(exc.value)

def test_hstack(table_types):
    """
    Hstack tables with mixin cols.  Use column "i" as proxy for what the
    result should be for each mixin.
    """
    t1 = table_types.Table()
    t1['i'] = table_types.Column([0, 1, 2, 3])
    for name, col in MIXIN_COLS.items():
        t1[name] = col
        col_setattr(t1[name], 'description', name)
        col_setattr(t1[name], 'meta', {'a': 1})

    for join_type in ('inner', 'outer'):
        for chop in (True, False):
            t2 = table_types.Table(t1)
            if chop:
                t2 = t2[:-1]
                if join_type == 'outer':
                    with pytest.raises(ValueError) as exc:
                        t12 = hstack([t1, t2], join_type=join_type)
                    assert 'hstack requires masking column' in str(exc.value)
                    continue

            t12 = hstack([t1, t2], join_type=join_type)
            idx1 = t12['i_1']
            idx2 = t12['i_2']
            for name, col in MIXIN_COLS.items():
                name1 = name + '_1'
                name2 = name + '_2'
                assert_table_name_col_equal(t12, name1, col[idx1])
                assert_table_name_col_equal(t12, name2, col[idx2])
                for attr in ('description', 'meta'):
                    assert col_getattr(t1[name], attr) == col_getattr(t12[name1], attr)
                    assert col_getattr(t2[name], attr) == col_getattr(t12[name2], attr)


def assert_table_name_col_equal(t, name, col):
    """
    Assert all(t[name] == col), with special handling for known mixin cols.
    """
    if isinstance(col, coordinates.SkyCoord):
        assert np.all(t[name].ra == col.ra)
        assert np.all(t[name].dec == col.dec)
    elif isinstance(col, u.Quantity):
        if type(t) is QTable:
            assert np.all(t[name].value == col.value)
    elif isinstance(col, table_helpers.ArrayWrapper):
        assert np.all(t[name].data == col.data)
    else:
        assert np.all(t[name] == col)

def test_get_items(mixin_cols):
    """
    Test that slicing / indexing table gives right values and col attrs inherit
    """
    attrs = ('name', 'unit', 'dtype', 'format', 'description', 'meta')
    m = mixin_cols['m']
    col_setattr(m, 'name', 'm')
    col_setattr(m, 'format', '{0}')
    col_setattr(m, 'description', 'd')
    col_setattr(m, 'meta', {'a': 1})
    t = QTable([m])
    for item in ([1, 3], np.array([0, 2]), slice(1, 3)):
        t2 = t[item]
        assert_table_name_col_equal(t2, 'm', m[item])
        for attr in attrs:
            assert col_getattr(t2['m'], attr) == col_getattr(m, attr)

def test_add_column(mixin_cols):
    """
    Test that adding a column preserves values and attributes
    """
    attrs = ('name', 'unit', 'dtype', 'format', 'description', 'meta')
    m = mixin_cols['m']
    assert col_getattr(m, 'name') is None

    # Make sure adding column in various ways doesn't touch
    t = QTable([m], names=['a'])
    assert col_getattr(m, 'name') is None

    t['new'] = m
    assert col_getattr(m, 'name') is None

    col_setattr(m, 'name', 'm')
    col_setattr(m, 'format', '{0}')
    col_setattr(m, 'description', 'd')
    col_setattr(m, 'meta', {'a': 1})
    t = QTable([m])

    # Add columns m2 and m3 by two different methods and test expected equality
    t['m2'] = m
    col_setattr(m, 'name', 'm3')
    t.add_columns([m], copy=True)
    col_setattr(m, 'name', 'm4')
    t.add_columns([m], copy=False)
    for name in ('m2', 'm3', 'm4'):
        assert_table_name_col_equal(t, 'm', t[name])
        for attr in attrs:
            if attr != 'name':
                assert col_getattr(t['m'], attr) == col_getattr(t[name], attr)

def test_vstack():
    """
    Vstack tables with mixin cols.
    """
    t1 = QTable(MIXIN_COLS)
    t2 = QTable(MIXIN_COLS)
    with pytest.raises(NotImplementedError):
        vstack([t1, t2])

def test_insert_row(mixin_cols):
    """
    Test inserting a row, which only works for BaseColumn and Quantity
    """
    t = QTable(mixin_cols)
    col_setattr(t['m'], 'description', 'd')
    if isinstance(t['m'], u.Quantity):
        t.insert_row(1, t[-1])
        assert t[1] == t[-1]
        assert col_getattr(t['m'], 'description') == 'd'
    else:
        with pytest.raises(ValueError) as exc:
            t.insert_row(1, t[-1])
        assert "Unable to insert row" in str(exc.value)

def test_convert_np_array(mixin_cols):
    """
    Test that converting to numpy array creates an object dtype and that
    each instance in the array has the expected type.
    """
    t = QTable(mixin_cols)
    ta = t.as_array()
    m = mixin_cols['m']
    dtype_kind = m.dtype.kind if hasattr(m, 'dtype') else 'O'
    assert ta['m'].dtype.kind == dtype_kind

def test_assignment_and_copy():
    """
    Test that assignment of an int, slice, and fancy index works.
    Along the way test that copying table works.
    """
    for name in ('quantity', 'arraywrap'):
        m = MIXIN_COLS[name]
        t0 = QTable([m], names=['m'])
        for i0, i1 in ((1, 2),
                       (slice(0, 2), slice(1, 3)),
                       (np.array([1, 2]), np.array([2, 3]))):
            t = t0.copy()
            t['m'][i0] = m[i1]
            if name == 'arraywrap':
                assert np.all(t['m'].data[i0] == m.data[i1])
                assert np.all(t0['m'].data[i0] == m.data[i0])
                assert np.all(t0['m'].data[i0] != t['m'].data[i0])
            else:
                assert np.all(t['m'][i0] == m[i1])
                assert np.all(t0['m'][i0] == m[i0])
                assert np.all(t0['m'][i0] != t['m'][i0])

def test_grouping():
    """
    Test grouping with mixin columns.  Raises not yet implemented error.
    """
    t = QTable(MIXIN_COLS)
    t['index'] = ['a', 'b', 'b', 'c']
    with pytest.raises(NotImplementedError):
        t.group_by('index')

def test_conversion_qtable_table():
    """
    Test that a table round trips from QTable => Table => QTable
    """
    qt = QTable(MIXIN_COLS)
    names = qt.colnames
    for name in names:
        col_setattr(qt[name], 'description', name)

    t = Table(qt)
    for name in names:
        assert col_getattr(t[name], 'description') == name
        if name == 'quantity':
            assert np.all(t['quantity'] == qt['quantity'].value)
            assert np.all(t['quantity'].unit is qt['quantity'].unit)
            assert isinstance(t['quantity'], t.ColumnClass)
        else:
            assert_table_name_col_equal(t, name, qt[name])

    qt2 = QTable(qt)
    for name in names:
        assert col_getattr(qt2[name], 'description') == name
        assert_table_name_col_equal(qt2, name, qt[name])

@pytest.mark.xfail
def test_column_rename():
    qt = QTable(MIXIN_COLS)
    names = qt.colnames
    for name in names:
        qt.rename_column(name, name + '2')
    assert qt.colnames == [name + '2' for name in names]
