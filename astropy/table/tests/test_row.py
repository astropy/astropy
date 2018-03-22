# Licensed under a 3-clause BSD style license - see LICENSE.rst


import sys

import pytest
import numpy as np

from ... import table
from ...table import Row
from ... import units as u
from .conftest import MaskedTable


def test_masked_row_with_object_col():
    """
    Numpy < 1.8 has a bug in masked array that prevents access a row if there is
    a column with object type.
    """
    t = table.Table([[1]], dtype=['O'], masked=True)
    t['col0'].mask = False
    assert t[0]['col0'] == 1
    t['col0'].mask = True
    assert t[0]['col0'] is np.ma.masked


@pytest.mark.usefixtures('table_types')
class TestRow():
    def _setup(self, table_types):
        self._table_type = table_types.Table
        self._column_type = table_types.Column

    @property
    def t(self):
        # py.test wants to run this method once before table_types is run
        # to set Table and Column.  In this case just return None, which would
        # cause any downstream test to fail if this happened in any other context.
        if self._column_type is None:
            return None
        if not hasattr(self, '_t'):
            a = self._column_type(name='a', data=[1, 2, 3], dtype='i8')
            b = self._column_type(name='b', data=[4, 5, 6], dtype='i8')
            self._t = self._table_type([a, b])
        return self._t

    def test_subclass(self, table_types):
        """Row is subclass of ndarray and Row"""
        self._setup(table_types)
        c = Row(self.t, 2)
        assert isinstance(c, Row)

    def test_values(self, table_types):
        """Row accurately reflects table values and attributes"""
        self._setup(table_types)
        table = self.t
        row = table[1]
        assert row['a'] == 2
        assert row['b'] == 5
        assert row[0] == 2
        assert row[1] == 5
        assert row.meta is table.meta
        assert row.colnames == table.colnames
        assert row.columns is table.columns
        with pytest.raises(IndexError):
            row[2]
        if sys.byteorder == 'little':
            assert str(row.dtype) == "[('a', '<i8'), ('b', '<i8')]"
        else:
            assert str(row.dtype) == "[('a', '>i8'), ('b', '>i8')]"

    def test_ref(self, table_types):
        """Row is a reference into original table data"""
        self._setup(table_types)
        table = self.t
        row = table[1]
        row['a'] = 10
        if table_types.Table is not MaskedTable:
            assert table['a'][1] == 10

    def test_left_equal(self, table_types):
        """Compare a table row to the corresponding structured array row"""
        self._setup(table_types)
        np_t = self.t.as_array()
        if table_types.Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(row == np_row)

    def test_left_not_equal(self, table_types):
        """Compare a table row to the corresponding structured array row"""
        self._setup(table_types)
        np_t = self.t.as_array()
        np_t['a'] = [0, 0, 0]
        if table_types.Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(row != np_row)

    def test_right_equal(self, table_types):
        """Test right equal"""
        self._setup(table_types)
        np_t = self.t.as_array()
        if table_types.Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(np_row == row)

    def test_convert_numpy_array(self, table_types):
        self._setup(table_types)
        d = self.t[1]

        np_data = np.array(d)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d.as_void())
        assert np_data is not d.as_void()
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d.as_void())
        assert np_data is not d.as_void()
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[(str('c'), 'i8'), (str('d'), 'i8')])

    def test_format_row(self, table_types):
        """Test formatting row"""
        self._setup(table_types)
        table = self.t
        row = table[0]
        assert repr(row).splitlines() == ['<{0} {1}{2}>'
                                          .format(row.__class__.__name__,
                                                  'index=0',
                                                  ' masked=True' if table.masked else ''),
                                          '  a     b  ',
                                          'int64 int64',
                                          '----- -----',
                                          '    1     4']
        assert str(row).splitlines() == [' a   b ',
                                         '--- ---',
                                         '  1   4']

        assert row._repr_html_().splitlines() == ['<i>{0} {1}{2}</i>'
                                                  .format(row.__class__.__name__,
                                                          'index=0',
                                                          ' masked=True' if table.masked else ''),
                                                  '<table id="table{0}">'.format(id(table)),
                                                  '<thead><tr><th>a</th><th>b</th></tr></thead>',
                                                  '<thead><tr><th>int64</th><th>int64</th></tr></thead>',
                                                  '<tr><td>1</td><td>4</td></tr>',
                                                  '</table>']

    def test_as_void(self, table_types):
        """Test the as_void() method"""
        self._setup(table_types)
        table = self.t
        row = table[0]

        # If masked then with no masks, issue numpy/numpy#483 should come
        # into play.  Make sure as_void() code is working.
        row_void = row.as_void()
        if table.masked:
            assert isinstance(row_void, np.ma.mvoid)
        else:
            assert isinstance(row_void, np.void)
        assert row_void['a'] == 1
        assert row_void['b'] == 4

        # Confirm row is a view of table but row_void is not.
        table['a'][0] = -100
        assert row['a'] == -100
        assert row_void['a'] == 1

        # Make sure it works for a table that has masked elements
        if table.masked:
            table['a'].mask = True

            # row_void is not a view, need to re-make
            assert row_void['a'] == 1
            row_void = row.as_void()  # but row is a view
            assert row['a'] is np.ma.masked

    def test_row_and_as_void_with_objects(self, table_types):
        """Test the deprecated data property and as_void() method"""
        t = table_types.Table([[{'a': 1}, {'b': 2}]], names=('a',))
        assert t[0][0] == {'a': 1}
        assert t[0]['a'] == {'a': 1}
        assert t[0].as_void()[0] == {'a': 1}
        assert t[0].as_void()['a'] == {'a': 1}

    def test_bounds_checking(self, table_types):
        """Row gives index error upon creation for out-of-bounds index"""
        self._setup(table_types)
        for ibad in (-5, -4, 3, 4):
            with pytest.raises(IndexError):
                self.t[ibad]


def test_row_tuple_column_slice():
    """
    Test getting and setting a row using a tuple or list of column names
    """
    t = table.QTable([[1, 2, 3] * u.m,
                      [10., 20., 30.],
                      [100., 200., 300.],
                      ['x', 'y', 'z']], names=['a', 'b', 'c', 'd'])
    # Get a row for index=1
    r1 = t[1]
    # Column slice with tuple of col names
    r1_abc = r1['a', 'b', 'c']  # Row object for these cols
    r1_abc_repr = ['<Row index=1>',
                   '   a       b       c   ',
                   '   m                   ',
                   'float64 float64 float64',
                   '------- ------- -------',
                   '    2.0    20.0   200.0']
    assert repr(r1_abc).splitlines() == r1_abc_repr

    # Column slice with list of col names
    r1_abc = r1[['a', 'b', 'c']]
    assert repr(r1_abc).splitlines() == r1_abc_repr

    # Make sure setting on a tuple or slice updates parent table and row
    r1['c'] = 1000
    r1['a', 'b'] = 1000 * u.cm, 100.
    assert r1['a'] == 10 * u.m
    assert r1['b'] == 100
    assert t['a'][1] == 10 * u.m
    assert t['b'][1] == 100.
    assert t['c'][1] == 1000

    # Same but using a list of column names instead of tuple
    r1[['a', 'b']] = 2000 * u.cm, 200.
    assert r1['a'] == 20 * u.m
    assert r1['b'] == 200
    assert t['a'][1] == 20 * u.m
    assert t['b'][1] == 200.

    # Set column slice of column slice
    r1_abc['a', 'c'] = -1 * u.m, -10
    assert t['a'][1] == -1 * u.m
    assert t['b'][1] == 200.
    assert t['c'][1] == -10.

    # Bad column name
    with pytest.raises(KeyError) as err:
        t[1]['a', 'not_there']
    assert "KeyError: 'not_there'" in str(err)

    # Too many values
    with pytest.raises(ValueError) as err:
        t[1]['a', 'b'] = 1 * u.m, 2, 3
    assert 'right hand side must be a sequence' in str(err)

    # Something without a length
    with pytest.raises(ValueError) as err:
        t[1]['a', 'b'] = 1
    assert 'right hand side must be a sequence' in str(err)


def test_row_tuple_column_slice_transaction():
    """
    Test that setting a row that fails part way through does not
    change the table at all.
    """
    t = table.QTable([[10., 20., 30.],
                      [1, 2, 3] * u.m], names=['a', 'b'])
    tc = t.copy()

    # First one succeeds but second fails.
    with pytest.raises(ValueError) as err:
        t[1]['a', 'b'] = (-1, -1 * u.s)  # Bad unit
    assert "'s' (time) and 'm' (length) are not convertible" in str(err)
    assert t[1] == tc[1]
