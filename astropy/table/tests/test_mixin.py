# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

try:
    import h5py  # pylint: disable=W0611 # noqa
except ImportError:
    HAS_H5PY = False
else:
    HAS_H5PY = True

try:
    import yaml  # pylint: disable=W0611 # noqa
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

import copy
import pickle
from io import StringIO

import pytest
import numpy as np

from astropy.coordinates import EarthLocation
from astropy.table import Table, QTable, join, hstack, vstack, Column, NdarrayMixin
from astropy.table import serialize
from astropy import time
from astropy import coordinates
from astropy import units as u
from astropy.table.column import BaseColumn
from astropy.table import table_helpers
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.metadata import MergeConflictWarning

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
    # We check for QuantityInfo rather than just isinstance(col, u.Quantity)
    # since we want to treat EarthLocation as a mixin, even though it is
    # a Quantity subclass.
    if ((isinstance(in_col.info, u.QuantityInfo) and type(table) is not QTable)
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
    from astropy.io.ascii.connect import _get_connectors_table
    t = QTable(MIXIN_COLS)
    for fmt in _get_connectors_table():
        if fmt['Format'] == 'ascii.ecsv' and not HAS_YAML:
            continue
        if fmt['Write'] and '.fast_' not in fmt['Format']:
            out = StringIO()
            t.write(out, format=fmt['Format'])


def test_votable_quantity_write(tmpdir):
    """
    Test that table with Quantity mixin column can be round-tripped by
    io.votable.  Note that FITS and HDF5 mixin support are tested (much more
    thoroughly) in their respective subpackage tests
    (io/fits/tests/test_connect.py and io/misc/tests/test_hdf5.py).
    """
    t = QTable()
    t['a'] = u.Quantity([1, 2, 4], unit='Angstrom')

    filename = str(tmpdir.join('table-tmp'))
    t.write(filename, format='votable', overwrite=True)
    qt = QTable.read(filename, format='votable')
    assert isinstance(qt['a'], u.Quantity)
    assert qt['a'].unit == 'Angstrom'


@pytest.mark.remote_data
@pytest.mark.parametrize('table_types', (Table, QTable))
def test_io_time_write_fits_standard(tmpdir, table_types):
    """
    Test that table with Time mixin columns can be written by io.fits.
    Validation of the output is done. Test that io.fits writes a table
    containing Time mixin columns that can be partially round-tripped
    (metadata scale, location).

    Note that we postpone checking the "local" scale, since that cannot
    be done with format 'cxcsec', as it requires an epoch.
    """
    t = table_types([[1, 2], ['string', 'column']])
    for scale in time.STANDARD_TIME_SCALES:
        t['a'+scale] = time.Time([[1, 2], [3, 4]], format='cxcsec',
                                 scale=scale, location=EarthLocation(
                                     -2446354, 4237210, 4077985, unit='m'))
        t['b'+scale] = time.Time(['1999-01-01T00:00:00.123456789',
                                  '2010-01-01T00:00:00'], scale=scale)
    t['c'] = [3., 4.]

    filename = str(tmpdir.join('table-tmp'))

    # Show that FITS format succeeds
    with pytest.warns(
            AstropyUserWarning,
            match='Time Column "btai" has no specified location, '
            'but global Time Position is present'):
        t.write(filename, format='fits', overwrite=True)
    with pytest.warns(
            AstropyUserWarning,
            match='Time column reference position "TRPOSn" is not specified'):
        tm = table_types.read(filename, format='fits', astropy_native=True)

    for scale in time.STANDARD_TIME_SCALES:
        for ab in ('a', 'b'):
            name = ab + scale

            # Assert that the time columns are read as Time
            assert isinstance(tm[name], time.Time)

            # Assert that the scales round-trip
            assert tm[name].scale == t[name].scale

            # Assert that the format is jd
            assert tm[name].format == 'jd'

            # Assert that the location round-trips
            assert tm[name].location == t[name].location

            # Finally assert that the column data round-trips
            assert (tm[name] == t[name]).all()

    for name in ('col0', 'col1', 'c'):
        # Assert that the non-time columns are read as Column
        assert isinstance(tm[name], Column)

        # Assert that the non-time columns' data round-trips
        assert (tm[name] == t[name]).all()

    # Test for conversion of time data to its value, as defined by its format
    for scale in time.STANDARD_TIME_SCALES:
        for ab in ('a', 'b'):
            name = ab + scale
            t[name].info.serialize_method['fits'] = 'formatted_value'

    t.write(filename, format='fits', overwrite=True)
    tm = table_types.read(filename, format='fits')

    for scale in time.STANDARD_TIME_SCALES:
        for ab in ('a', 'b'):
            name = ab + scale

            assert not isinstance(tm[name], time.Time)
            assert (tm[name] == t[name].value).all()


@pytest.mark.parametrize('table_types', (Table, QTable))
def test_io_time_write_fits_local(tmpdir, table_types):
    """
    Test that table with a Time mixin with scale local can also be written
    by io.fits. Like ``test_io_time_write_fits_standard`` above, but avoiding
    ``cxcsec`` format, which requires an epoch and thus cannot be used for a
    local time scale.
    """
    t = table_types([[1, 2], ['string', 'column']])
    t['a_local'] = time.Time([[50001, 50002], [50003, 50004]],
                             format='mjd', scale='local',
                             location=EarthLocation(-2446354, 4237210, 4077985,
                                                    unit='m'))
    t['b_local'] = time.Time(['1999-01-01T00:00:00.123456789',
                              '2010-01-01T00:00:00'], scale='local')
    t['c'] = [3., 4.]

    filename = str(tmpdir.join('table-tmp'))

    # Show that FITS format succeeds

    with pytest.warns(AstropyUserWarning, match='Time Column "b_local" has no specified location'):
        t.write(filename, format='fits', overwrite=True)

    with pytest.warns(AstropyUserWarning, match='Time column reference position "TRPOSn" is not specified.'):
        tm = table_types.read(filename, format='fits', astropy_native=True)

    for ab in ('a', 'b'):
        name = ab + '_local'

        # Assert that the time columns are read as Time
        assert isinstance(tm[name], time.Time)

        # Assert that the scales round-trip
        assert tm[name].scale == t[name].scale

        # Assert that the format is jd
        assert tm[name].format == 'jd'

        # Assert that the location round-trips
        assert tm[name].location == t[name].location

        # Finally assert that the column data round-trips
        assert (tm[name] == t[name]).all()

    for name in ('col0', 'col1', 'c'):
        # Assert that the non-time columns are read as Column
        assert isinstance(tm[name], Column)

        # Assert that the non-time columns' data round-trips
        assert (tm[name] == t[name]).all()

    # Test for conversion of time data to its value, as defined by its format.
    for ab in ('a', 'b'):
        name = ab + '_local'
        t[name].info.serialize_method['fits'] = 'formatted_value'

    t.write(filename, format='fits', overwrite=True)
    tm = table_types.read(filename, format='fits')

    for ab in ('a', 'b'):
        name = ab + '_local'

        assert not isinstance(tm[name], time.Time)
        assert (tm[name] == t[name].value).all()


def test_votable_mixin_write_fail(mixin_cols):
    """
    Test that table with mixin columns (excluding Quantity) cannot be written by
    io.votable.
    """
    t = QTable(mixin_cols)
    # Only do this test if there are unsupported column types (i.e. anything besides
    # BaseColumn and Quantity class instances).
    unsupported_cols = t.columns.not_isinstance((BaseColumn, u.Quantity))

    if not unsupported_cols:
        pytest.skip("no unsupported column types")

    out = StringIO()
    with pytest.raises(ValueError) as err:
        t.write(out, format='votable')
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
        with pytest.raises(NotImplementedError) as exc:
            t12 = join(t1, t2, keys='a', join_type=join_type)
        assert 'join requires masking column' in str(exc.value)

    with pytest.raises(TypeError) as exc:
        t12 = join(t1, t2, keys=['a', 'skycoord'])
    assert 'one or more key columns are not sortable' in str(exc.value)

    # Join does work for a mixin which is a subclass of np.ndarray
    with pytest.warns(MergeConflictWarning,
                      match="In merged column 'quantity' the 'description' attribute does not match"):
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
                    with pytest.raises(NotImplementedError) as exc:
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
            assert np.all(t[name] == col)
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

    # Add columns m2, m3, m4 by two different methods and test expected equality
    t['m2'] = m
    m.info.name = 'm3'
    t.add_columns([m], copy=True)
    m.info.name = 'm4'
    t.add_columns([m], copy=False)
    for name in ('m2', 'm3', 'm4'):
        assert_table_name_col_equal(t, name, m)
        for attr in attrs:
            if attr != 'name':
                assert getattr(t['m'].info, attr) == getattr(t[name].info, attr)
    # Also check that one can set using a scalar.
    s = m[0]
    if type(s) is type(m):
        # We're not going to worry about testing classes for which scalars
        # are a different class than the real array (and thus loose info, etc.)
        t['s'] = m[0]
        assert_table_name_col_equal(t, 's', m[0])
        for attr in attrs:
            if attr != 'name':
                assert getattr(t['m'].info, attr) == getattr(t['s'].info, attr)

    # While we're add it, also check a length-1 table.
    t = QTable([m[1:2]], names=['m'])
    if type(s) is type(m):
        t['s'] = m[0]
        assert_table_name_col_equal(t, 's', m[0])
        for attr in attrs:
            if attr != 'name':
                assert getattr(t['m'].info, attr) == getattr(t['s'].info, attr)


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
    if isinstance(t['m'], (u.Quantity, Column, time.Time)):
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
    c = coordinates.SkyCoord([0], [1], [0], representation_type='cartesian')
    t = Table([c])
    assert t.pformat() == ['     col0     ',
                           'None,None,None',
                           '--------------',
                           '   0.0,1.0,0.0']

    # Test that info works with a dynamically changed representation
    c = coordinates.SkyCoord([0], [1], [0], unit='m', representation_type='cartesian')
    t = Table([c])
    assert t.pformat() == ['    col0   ',
                           '   m,m,m   ',
                           '-----------',
                           '0.0,1.0,0.0']

    t['col0'].representation_type = 'unitspherical'
    assert t.pformat() == ['  col0  ',
                           'deg,deg ',
                           '--------',
                           '90.0,0.0']

    t['col0'].representation_type = 'cylindrical'
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
                 dtype='<i4,' + ('|U1'))
    b = np.array([(10, 'aa'), (20, 'bb'), (30, 'cc'), (40, 'dd')],
                 dtype=[('x', 'i4'), ('y', ('U2'))])
    c = np.rec.fromrecords([(100, 'raa'), (200, 'rbb'), (300, 'rcc'), (400, 'rdd')],
                           names=['rx', 'ry'])
    d = np.arange(8).reshape(4, 2).view(NdarrayMixin)

    # Add one during initialization and the next as a new column.
    t = Table([a], names=['a'])
    t['b'] = b
    t['c'] = c
    t['d'] = d

    assert isinstance(t['a'], NdarrayMixin)

    assert t['a'][1][1] == a[1][1]
    assert t['a'][2][0] == a[2][0]

    assert t[1]['a'][1] == a[1][1]
    assert t[2]['a'][0] == a[2][0]

    assert isinstance(t['b'], NdarrayMixin)

    assert t['b'][1]['x'] == b[1]['x']
    assert t['b'][1]['y'] == b[1]['y']

    assert t[1]['b']['x'] == b[1]['x']
    assert t[1]['b']['y'] == b[1]['y']

    assert isinstance(t['c'], NdarrayMixin)

    assert t['c'][1]['rx'] == c[1]['rx']
    assert t['c'][1]['ry'] == c[1]['ry']

    assert t[1]['c']['rx'] == c[1]['rx']
    assert t[1]['c']['ry'] == c[1]['ry']

    assert isinstance(t['d'], NdarrayMixin)

    assert t['d'][1][0] == d[1][0]
    assert t['d'][1][1] == d[1][1]

    assert t[1]['d'][0] == d[1][0]
    assert t[1]['d'][1] == d[1][1]

    assert t.pformat() == ['   a         b           c       d [2] ',
                           '-------- ---------- ------------ ------',
                           "(1, 'a') (10, 'aa') (100, 'raa') 0 .. 1",
                           "(2, 'b') (20, 'bb') (200, 'rbb') 2 .. 3",
                           "(3, 'c') (30, 'cc') (300, 'rcc') 4 .. 5",
                           "(4, 'd') (40, 'dd') (400, 'rdd') 6 .. 7"]


def test_possible_string_format_functions():
    """
    The QuantityInfo info class for Quantity implements a
    possible_string_format_functions() method that overrides the
    standard pprint._possible_string_format_functions() function.
    Test this.
    """
    t = QTable([[1, 2] * u.m])
    t['col0'].info.format = '%.3f'
    assert t.pformat() == [' col0',
                           '  m  ',
                           '-----',
                           '1.000',
                           '2.000']

    t['col0'].info.format = 'hi {:.3f}'
    assert t.pformat() == ['  col0  ',
                           '   m    ',
                           '--------',
                           'hi 1.000',
                           'hi 2.000']

    t['col0'].info.format = '.4f'
    assert t.pformat() == [' col0 ',
                           '  m   ',
                           '------',
                           '1.0000',
                           '2.0000']


def test_rename_mixin_columns(mixin_cols):
    """
    Rename a mixin column.
    """
    t = QTable(mixin_cols)
    tc = t.copy()
    t.rename_column('m', 'mm')
    assert t.colnames == ['i', 'a', 'b', 'mm']
    if isinstance(t['mm'], table_helpers.ArrayWrapper):
        assert np.all(t['mm'].data == tc['m'].data)
    elif isinstance(t['mm'], coordinates.SkyCoord):
        assert np.all(t['mm'].ra == tc['m'].ra)
        assert np.all(t['mm'].dec == tc['m'].dec)
    else:
        assert np.all(t['mm'] == tc['m'])


def test_represent_mixins_as_columns_unit_fix():
    """
    If the unit is invalid for a column that gets serialized this would
    cause an exception.  Fixed in #7481.
    """
    t = Table({'a': [1, 2]}, masked=True)
    t['a'].unit = 'not a valid unit'
    t['a'].mask[1] = True
    serialize.represent_mixins_as_columns(t)
