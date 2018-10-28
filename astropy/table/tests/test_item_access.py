# Licensed under a 3-clause BSD style license - see LICENSE.rst


""" Verify item access API in:
https://github.com/astropy/astropy/wiki/Table-item-access-definition
"""

import pytest
import numpy as np


@pytest.mark.usefixtures('table_data')
class BaseTestItems():
    pass


@pytest.mark.usefixtures('table_data')
class TestTableColumnsItems(BaseTestItems):

    def test_by_name(self, table_data):
        """Access TableColumns by name and show that item access returns
        a Column that refers to underlying table data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        assert self.tc['a'].name == 'a'
        assert self.tc['a'][1] == 2
        assert self.tc['a'].description == 'da'
        assert self.tc['a'].format == '%i'
        assert self.tc['a'].meta == {'ma': 1}
        assert self.tc['a'].unit == 'ua'
        assert self.tc['a'].attrs_equal(table_data.COLS[0])
        assert isinstance(self.tc['a'], table_data.Column)

        self.tc['b'][1] = 0
        assert self.t['b'][1] == 0

    def test_by_position(self, table_data):
        """Access TableColumns by position and show that item access returns
        a Column that refers to underlying table data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        assert self.tc[1].name == 'b'
        assert np.all(self.tc[1].data == table_data.COLS[1].data)
        assert self.tc[1].description == 'db'
        assert self.tc[1].format == '%d'
        assert self.tc[1].meta == {'mb': 1}
        assert self.tc[1].unit == 'ub'
        assert self.tc[1].attrs_equal(table_data.COLS[1])
        assert isinstance(self.tc[1], table_data.Column)

        assert self.tc[2].unit == 'ub'

        self.tc[1][1] = 0
        assert self.t['b'][1] == 0

    def test_mult_columns(self, table_data):
        """Access TableColumns with "fancy indexing" and showed that returned
        TableColumns object still references original data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        tc2 = self.tc['b', 'c']
        assert tc2[1].name == 'c'
        assert tc2[1][1] == 8
        assert tc2[0].name == 'b'
        assert tc2[0][1] == 5

        tc2['c'][1] = 0
        assert self.tc['c'][1] == 0
        assert self.t['c'][1] == 0

    def test_column_slice(self, table_data):
        """Access TableColumns with slice and showed that returned
        TableColumns object still references original data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        tc2 = self.tc[1:3]
        assert tc2[1].name == 'c'
        assert tc2[1][1] == 8
        assert tc2[0].name == 'b'
        assert tc2[0][1] == 5

        tc2['c'][1] = 0
        assert self.tc['c'][1] == 0
        assert self.t['c'][1] == 0


@pytest.mark.usefixtures('table_data')
class TestTableItems(BaseTestItems):

    @pytest.mark.parametrize("idx", [1, np.int64(1), np.array(1)])
    def test_column(self, table_data, idx):
        """Column access returns REFERENCE to data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        a = self.t['a']
        assert a[idx] == 2
        a[idx] = 0
        assert self.t['a'][idx] == 0

    @pytest.mark.parametrize("idx", [1, np.int64(1), np.array(1)])
    def test_row(self, table_data, idx):
        """Row  access returns REFERENCE to data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        row = self.t[idx]
        assert row['a'] == 2
        assert row[idx] == 5
        assert row.columns['a'].attrs_equal(table_data.COLS[0])
        assert row.columns['b'].attrs_equal(table_data.COLS[1])
        assert row.columns['c'].attrs_equal(table_data.COLS[2])

        # Check that setting by col index sets the table and row value
        row[idx] = 0
        assert row[idx] == 0
        assert row['b'] == 0
        assert self.t['b'][idx] == 0
        assert self.t[idx]['b'] == 0

        # Check that setting by col name sets the table and row value
        row['a'] = 0
        assert row[0] == 0
        assert row['a'] == 0
        assert self.t['a'][1] == 0
        assert self.t[1]['a'] == 0

    def test_empty_iterable_item(self, table_data):
        """
        Table item access with [], (), or np.array([]) returns the same table
        with no rows.
        """
        self.t = table_data.Table(table_data.COLS)
        for item in [], (), np.array([]):
            t2 = self.t[item]
            assert not t2
            assert len(t2) == 0
            assert t2['a'].attrs_equal(table_data.COLS[0])
            assert t2['b'].attrs_equal(table_data.COLS[1])
            assert t2['c'].attrs_equal(table_data.COLS[2])

    def test_table_slice(self, table_data):
        """Table slice returns REFERENCE to data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        t2 = self.t[1:3]
        assert np.all(t2['a'] == table_data.DATA['a'][1:3])
        assert t2['a'].attrs_equal(table_data.COLS[0])
        assert t2['b'].attrs_equal(table_data.COLS[1])
        assert t2['c'].attrs_equal(table_data.COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t['a'] == np.array([1, 0, 3]))
        assert t2.masked == self.t.masked
        assert t2._column_class == self.t._column_class
        assert isinstance(t2, table_data.Table)

    def test_fancy_index_slice(self, table_data):
        """Table fancy slice returns COPY of data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        slice = np.array([0, 2])
        t2 = self.t[slice]
        assert np.all(t2['a'] == table_data.DATA['a'][slice])
        assert t2['a'].attrs_equal(table_data.COLS[0])
        assert t2['b'].attrs_equal(table_data.COLS[1])
        assert t2['c'].attrs_equal(table_data.COLS[2])
        t2['a'][0] = 0

        assert np.all(self.t.as_array() == table_data.DATA)
        assert np.any(t2['a'] != table_data.DATA['a'][slice])
        assert t2.masked == self.t.masked
        assert t2._column_class == self.t._column_class
        assert isinstance(t2, table_data.Table)

    def test_list_index_slice(self, table_data):
        """Table list index slice returns COPY of data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        slice = [0, 2]
        t2 = self.t[slice]
        assert np.all(t2['a'] == table_data.DATA['a'][slice])
        assert t2['a'].attrs_equal(table_data.COLS[0])
        assert t2['b'].attrs_equal(table_data.COLS[1])
        assert t2['c'].attrs_equal(table_data.COLS[2])
        t2['a'][0] = 0

        assert np.all(self.t.as_array() == table_data.DATA)
        assert np.any(t2['a'] != table_data.DATA['a'][slice])
        assert t2.masked == self.t.masked
        assert t2._column_class == self.t._column_class
        assert isinstance(t2, table_data.Table)

    def test_select_columns(self, table_data):
        """Select columns returns COPY of data and all column
        attributes"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        # try both lists and tuples
        for columns in (('a', 'c'), ['a', 'c']):
            t2 = self.t[columns]
            assert np.all(t2['a'] == table_data.DATA['a'])
            assert np.all(t2['c'] == table_data.DATA['c'])
            assert t2['a'].attrs_equal(table_data.COLS[0])
            assert t2['c'].attrs_equal(table_data.COLS[2])
            t2['a'][0] = 0
            assert np.all(self.t.as_array() == table_data.DATA)
            assert np.any(t2['a'] != table_data.DATA['a'])
            assert t2.masked == self.t.masked
            assert t2._column_class == self.t._column_class

    def test_select_columns_fail(self, table_data):
        """Selecting a column that doesn't exist fails"""
        self.t = table_data.Table(table_data.COLS)

        with pytest.raises(KeyError) as err:
            self.t[['xxxx']]
        assert "KeyError: 'xxxx'" in str(err)

        with pytest.raises(KeyError) as err:
            self.t[['xxxx', 'yyyy']]
        assert "KeyError: 'xxxx'" in str(err)

    def test_np_where(self, table_data):
        """Select rows using output of np.where"""
        t = table_data.Table(table_data.COLS)
        # Select last two rows
        rows = np.where(t['a'] > 1.5)
        t2 = t[rows]
        assert np.all(t2['a'] == [2, 3])
        assert np.all(t2['b'] == [5, 6])
        assert isinstance(t2, table_data.Table)

        # Select no rows
        rows = np.where(t['a'] > 100)
        t2 = t[rows]
        assert len(t2) == 0
        assert isinstance(t2, table_data.Table)

    def test_np_integers(self, table_data):
        """
        Select rows using numpy integers.  This is a regression test for a
        py 3.3 failure mode
        """
        t = table_data.Table(table_data.COLS)
        idxs = np.random.randint(len(t), size=2)
        item = t[idxs[1]]

    def test_select_bad_column(self, table_data):
        """Select column name that does not exist"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        with pytest.raises(ValueError):
            self.t['a', 1]
