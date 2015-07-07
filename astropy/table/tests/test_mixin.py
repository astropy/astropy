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

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

import numpy as np
import copy

from ...extern.six.moves import cPickle as pickle
from ...tests.helper import pytest
from ...table import Table, QTable, join, hstack, vstack, Column, NdarrayMixin
from ... import time
from ... import coordinates
from ... import units as u
from .. import table_helpers
from .conftest import MIXIN_COLS


def test_attributes(mixin_cols):
    """
    Required attributes for a column can be set.
    """
    m = mixin_cols['m']
    m.info.name = 'a'
    assert m.info.name == 'a'

    m.info.description = 'a'
    assert m.info.description == 'a'

    # Cannot set unit for these classes
    if isinstance(m, (u.Quantity, coordinates.SkyCoord, time.Time)):
        with pytest.raises(AttributeError):
            m.info.unit = u.m
    else:
        m.info.unit = u.m
        assert m.info.unit is u.m

    m.info.format = 'a'
    assert m.info.format == 'a'

    m.info.meta = {'a': 1}
    assert m.info.meta == {'a': 1}

    with pytest.raises(AttributeError):
        m.info.bad_attr = 1

    with pytest.raises(AttributeError):
        m.info.bad_attr


def check_mixin_type(table, table_col, in_col):
    if ((isinstance(in_col, u.Quantity) and type(table) is not QTable)
            or isinstance(in_col, Column)):
        assert type(table_col) is table.ColumnClass
    else:
        assert type(table_col) is type(in_col)

    # Make sure in_col got copied and creating table did not touch it
    assert in_col.info.name is None

def test_make_table(table_types, mixin_cols):
    """
    Make a table with the columns in mixin_cols, which is an ordered dict of
    three cols: 'a' and 'b' are table_types.Column type, and 'm' is a mixin.
    """
    t = table_types.Table(mixin_cols)
    check_mixin_type(t, t['m'], mixin_cols['m'])

    cols = list(mixin_cols.values())
    t = table_types.Table(cols, names=('i', 'a', 'b', 'm'))
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
        if fmt['Format'] == 'ascii.ecsv' and not HAS_YAML:
            continue
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
        t1[name].info.description = name
        t2[name].info.description = name + '2'

    for join_type in ('inner', 'left'):
        t12 = join(t1, t2, keys='a', join_type=join_type)
        idx1 = t12['i_1']
        idx2 = t12['i_2']
        for name, col in MIXIN_COLS.items():
            name1 = name + '_1'
            name2 = name + '_2'
            assert_table_name_col_equal(t12, name1, col[idx1])
            assert_table_name_col_equal(t12, name2, col[idx2])
            assert t12[name1].info.description == name
            assert t12[name2].info.description == name + '2'

    for join_type in ('outer', 'right'):
        with pytest.raises(ValueError) as exc:
            t12 = join(t1, t2, keys='a', join_type=join_type)
        assert 'join requires masking column' in str(exc.value)

    with pytest.raises(ValueError) as exc:
        t12 = join(t1, t2, keys=['a', 'skycoord'])
    assert 'not allowed as a key column' in str(exc.value)

    # Join does work for a mixin which is a subclass of np.ndarray
    t12 = join(t1, t2, keys=['quantity'])
    assert np.all(t12['a_1'] == t1['a'])

def test_hstack(table_types):
    """
    Hstack tables with mixin cols.  Use column "i" as proxy for what the
    result should be for each mixin.
    """
    t1 = table_types.Table()
    t1['i'] = table_types.Column([0, 1, 2, 3])
    for name, col in MIXIN_COLS.items():
        t1[name] = col
        t1[name].info.description = name
        t1[name].info.meta = {'a': 1}

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
                    assert getattr(t1[name].info, attr) == getattr(t12[name1].info, attr)
                    assert getattr(t2[name].info, attr) == getattr(t12[name2].info, attr)


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
    m.info.name = 'm'
    m.info.format = '{0}'
    m.info.description = 'd'
    m.info.meta = {'a': 1}
    t = QTable([m])
    for item in ([1, 3], np.array([0, 2]), slice(1, 3)):
        t2 = t[item]
        m2 = m[item]
        assert_table_name_col_equal(t2, 'm', m[item])
        for attr in attrs:
            assert getattr(t2['m'].info, attr) == getattr(m.info, attr)
            assert getattr(m2.info, attr) == getattr(m.info, attr)

def test_info_preserved_pickle_copy_init(mixin_cols):
    """
    Test copy, pickle, and init from class roundtrip preserve info.  This
    tests not only the mixin classes but a regular column as well.
    """
    def pickle_roundtrip(c):
        return pickle.loads(pickle.dumps(c))

    def init_from_class(c):
        return c.__class__(c)

    attrs = ('name', 'unit', 'dtype', 'format', 'description', 'meta')
    for colname in ('i', 'm'):
        m = mixin_cols[colname]
        m.info.name = colname
        m.info.format = '{0}'
        m.info.description = 'd'
        m.info.meta = {'a': 1}
        for func in (copy.copy, copy.deepcopy, pickle_roundtrip, init_from_class):
            m2 = func(m)
            for attr in attrs:
                assert getattr(m2.info, attr) == getattr(m.info, attr)


def test_add_column(mixin_cols):
    """
    Test that adding a column preserves values and attributes
    """
    attrs = ('name', 'unit', 'dtype', 'format', 'description', 'meta')
    m = mixin_cols['m']
    assert m.info.name is None

    # Make sure adding column in various ways doesn't touch
    t = QTable([m], names=['a'])
    assert m.info.name is None

    t['new'] = m
    assert m.info.name is None

    m.info.name = 'm'
    m.info.format = '{0}'
    m.info.description = 'd'
    m.info.meta = {'a': 1}
    t = QTable([m])

    # Add columns m2 and m3 by two different methods and test expected equality
    t['m2'] = m
    m.info.name = 'm3'
    t.add_columns([m], copy=True)
    m.info.name = 'm4'
    t.add_columns([m], copy=False)
    for name in ('m2', 'm3', 'm4'):
        assert_table_name_col_equal(t, 'm', t[name])
        for attr in attrs:
            if attr != 'name':
                assert getattr(t['m'].info, attr) == getattr(t[name].info, attr)

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
    t['m'].info.description = 'd'
    if isinstance(t['m'], (u.Quantity, Column)):
        t.insert_row(1, t[-1])
        assert t[1] == t[-1]
        assert t['m'].info.description == 'd'
    else:
        with pytest.raises(ValueError) as exc:
            t.insert_row(1, t[-1])
        assert "Unable to insert row" in str(exc.value)

def test_insert_row_bad_unit():
    """
    Insert a row into a QTable with the wrong unit
    """
    t = QTable([[1] * u.m])
    with pytest.raises(ValueError) as exc:
        t.insert_row(0, (2 * u.m / u.s,))
    assert "'m / s' (speed) and 'm' (length) are not convertible" in str(exc.value)

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
        qt[name].info.description = name

    t = Table(qt)
    for name in names:
        assert t[name].info.description == name
        if name == 'quantity':
            assert np.all(t['quantity'] == qt['quantity'].value)
            assert np.all(t['quantity'].unit is qt['quantity'].unit)
            assert isinstance(t['quantity'], t.ColumnClass)
        else:
            assert_table_name_col_equal(t, name, qt[name])

    qt2 = QTable(qt)
    for name in names:
        assert qt2[name].info.description == name
        assert_table_name_col_equal(qt2, name, qt[name])

@pytest.mark.xfail
def test_column_rename():
    qt = QTable(MIXIN_COLS)
    names = qt.colnames
    for name in names:
        qt.rename_column(name, name + '2')
    assert qt.colnames == [name + '2' for name in names]


def test_setitem_as_column_name():
    """
    Test for mixin-related regression described in #3321.
    """
    t = Table()
    t['a'] = ['x', 'y']
    t['b'] = 'b'  # Previously was failing with KeyError
    assert np.all(t['a'] == ['x', 'y'])
    assert np.all(t['b'] == ['b', 'b'])


def test_quantity_representation():
    """
    Test that table representation of quantities does not have unit
    """
    t = QTable([[1, 2] * u.m])
    assert t.pformat() == ['col0',
                           ' m  ',
                           '----',
                           ' 1.0',
                           ' 2.0']


def test_skycoord_representation():
    """
    Test that skycoord representation works, both in the way that the
    values are output and in changing the frame representation.
    """
    # With no unit we get "None" in the unit row
    c = coordinates.SkyCoord([0], [1], [0], representation='cartesian')
    t = Table([c])
    assert t.pformat() == ['     col0     ',
                           'None,None,None',
                           '--------------',
                           '   0.0,1.0,0.0']

    # Test that info works with a dynamically changed representation
    c = coordinates.SkyCoord([0], [1], [0], unit='m', representation='cartesian')
    t = Table([c])
    assert t.pformat() == ['    col0   ',
                           '   m,m,m   ',
                           '-----------',
                           '0.0,1.0,0.0']

    t['col0'].representation = 'unitspherical'
    assert t.pformat() == ['  col0  ',
                           'deg,deg ',
                           '--------',
                           '90.0,0.0']

    t['col0'].representation = 'cylindrical'
    assert t.pformat() == ['    col0    ',
                           '  m,deg,m   ',
                           '------------',
                           '1.0,90.0,0.0']


def test_ndarray_mixin():
    """
    Test directly adding a plain structured array into a table instead of the
    view as an NdarrayMixin.  Once added as an NdarrayMixin then all the previous
    tests apply.
    """
    a = np.array([(1, 'a'), (2, 'b'), (3, 'c'), (4, 'd')],
                 dtype='<i4,|S1')
    b = np.array([(10, 'aa'), (20, 'bb'), (30, 'cc'), (40, 'dd')],
                 dtype=[('x', 'i4'), ('y', 'S2')])
    c = np.rec.fromrecords([(100, 'raa'), (200, 'rbb'), (300, 'rcc'), (400, 'rdd')],
                           names=['rx', 'ry'])
    d = np.arange(8).reshape(4, 2).view(NdarrayMixin)

    # Add one during initialization and the next as a new column.
    t = Table([a], names=['a'])
    t['b'] = b
    t['c'] = c
    t['d'] = d

    assert isinstance(t['a'], NdarrayMixin)

    assert t['a'][1][1] == 'b'
    assert t['a'][2][0] == 3

    assert t[1]['a'][1] == 'b'
    assert t[2]['a'][0] == 3

    assert isinstance(t['b'], NdarrayMixin)

    assert t['b'][1]['x'] == 20
    assert t['b'][1]['y'] == 'bb'

    assert t[1]['b']['x'] == 20
    assert t[1]['b']['y'] == 'bb'

    assert isinstance(t['c'], NdarrayMixin)

    assert t['c'][1]['rx'] == 200
    assert t['c'][1]['ry'] == 'rbb'

    assert t[1]['c']['rx'] == 200
    assert t[1]['c']['ry'] == 'rbb'

    assert isinstance(t['d'], NdarrayMixin)

    assert t['d'][1][0] == 2
    assert t['d'][1][1] == 3

    assert t[1]['d'][0] == 2
    assert t[1]['d'][1] == 3

    assert t.pformat() == ['   a         b           c       d [2] ',
                           '-------- ---------- ------------ ------',
                           "(1, 'a') (10, 'aa') (100, 'raa') 0 .. 1",
                           "(2, 'b') (20, 'bb') (200, 'rbb') 2 .. 3",
                           "(3, 'c') (30, 'cc') (300, 'rcc') 4 .. 5",
                           "(4, 'd') (40, 'dd') (400, 'rdd') 6 .. 7"]
