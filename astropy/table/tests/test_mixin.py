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
from ...table import QTable, col_setattr, col_getattr, join, hstack, vstack
from ... import units as u
from ... import coordinates
from .. import table_helpers
from .conftest import MIXIN_COLS

# ISSUES / TODO
# - Test groups
# - Check attributes in join outputs
# - Test convert QTable <=> Table
# - Assignment
# - Copy
# - Array subsetting

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
            assert_table_name_col_equal(t12, name1, col[idx1])
            assert_table_name_col_equal(t12, name2, col[idx2])

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

def test_assignment():
    for name in ('quantity', 'arraywrap'):
        m = MIXIN_COLS[name]
        t = QTable([m], names=['m'])
        t['m'][0:2] = m[1:3]
        if name == 'arraywrap':
            assert np.all(t['m'].data[0:2] == m.data[1:3])
        else:
            assert np.all(t['m'][0:2] == m[1:3])
