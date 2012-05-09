""" Verify item access API in:
https://github.com/astropy/astropy/wiki/Table-item-access-definition
"""
from .. import Table, Column
import numpy as np
import pytest

COLS = [Column('a', [1, 2, 3], description='da',
               format='fa', meta={'ma': 1}, units='ua'),
        Column('b', [4, 5, 6], description='db',
               format='fb', meta={'mb': 1}, units='ub'),
        Column('c', [7, 8, 9], description='dc',
               format='fc', meta={'mc': 1}, units='ub')]
DATA = Table(COLS)


class BaseTestItems():

    def setup_method(self, method):
        self.t = Table(COLS)
        self.tc = self.t.columns


class TestTableColumnsItems(BaseTestItems):

    def test_by_name(self):
        """Access TableColumns by name and show that item access returns
        a Column that refers to underlying table data"""
        assert self.tc['a'].name == 'a'
        assert self.tc['a'][1] == 2
        assert self.tc['a'].description == 'da'
        assert self.tc['a'].format == 'fa'
        assert self.tc['a'].meta == {'ma': 1}
        assert self.tc['a'].units == 'ua'
        assert self.tc['a'].attrs_equal(COLS[0])
        assert isinstance(self.tc['a'], Column)

        self.tc['b'][1] = 0
        assert self.t['b'][1] == 0

    def test_by_position(self):
        """Access TableColumns by position and show that item access returns
        a Column that refers to underlying table data"""
        assert self.tc[1].name == 'b'
        assert np.all(self.tc[1].data == COLS[1].data)
        assert self.tc[1].description == 'db'
        assert self.tc[1].format == 'fb'
        assert self.tc[1].meta == {'mb': 1}
        assert self.tc[1].units == 'ub'
        assert self.tc[1].attrs_equal(COLS[1])
        assert isinstance(self.tc[1], Column)

        self.tc[1][1] = 0
        assert self.t['b'][1] == 0

    def test_mult_columns(self):
        """Access TableColumns with "fancy indexing" and showed that returned
        TableColumns object still references original data"""
        tc2 = self.tc['b', 'c']
        assert tc2[1].name == 'c'
        assert tc2[1][1] == 8
        assert tc2[0].name == 'b'
        assert tc2[0][1] == 5

        tc2['c'][1] = 0
        assert self.tc['c'][1] == 0
        assert self.t['c'][1] == 0

    def test_column_slice(self):
        """Access TableColumns with slice and showed that returned
        TableColumns object still references original data"""
        tc2 = self.tc[1:3]
        assert tc2[1].name == 'c'
        assert tc2[1][1] == 8
        assert tc2[0].name == 'b'
        assert tc2[0][1] == 5

        tc2['c'][1] = 0
        assert self.tc['c'][1] == 0
        assert self.t['c'][1] == 0


class TestTableItems(BaseTestItems):

    def test_column(self):
        """Column access returns REFERENCE to data"""
        a = self.t['a']
        assert a[1] == 2
        a[1] = 0
        assert self.t['a'][1] == 0

    def test_row(self):
        """Row  access returns REFERENCE to data"""
        row = self.t[1]
        assert row['a'] == 2
        assert row[1] == 5
        assert row.columns['a'].attrs_equal(COLS[0])
        assert row.columns['b'].attrs_equal(COLS[1])
        assert row.columns['c'].attrs_equal(COLS[2])
        row[1] = 0
        assert row[1] == 0
        assert self.t['b'][1] == 0

    def test_table_slice(self):
        """Table slice returns REFERENCE to data"""
        t2 = self.t[1:3]
        assert np.all(t2['a'] == DATA['a'][1:3])
        assert t2['a'].attrs_equal(COLS[0])
        assert t2['b'].attrs_equal(COLS[1])
        assert t2['c'].attrs_equal(COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t['a'] == np.array([1, 0, 3]))

    def test_fancy_index_slice(self):
        """Table fancy slice returns COPY of data"""
        slice = np.array([0, 2])
        t2 = self.t[slice]
        assert np.all(t2['a'] == DATA['a'][slice])
        assert t2['a'].attrs_equal(COLS[0])
        assert t2['b'].attrs_equal(COLS[1])
        assert t2['c'].attrs_equal(COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t._data == DATA)
        assert np.any(t2['a'] != DATA['a'])

    def test_list_index_slice(self):
        """Table list index slice returns COPY of data"""
        slice = [0, 2]
        t2 = self.t[slice]
        assert np.all(t2['a'] == DATA['a'][slice])
        assert t2['a'].attrs_equal(COLS[0])
        assert t2['b'].attrs_equal(COLS[1])
        assert t2['c'].attrs_equal(COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t._data == DATA)
        assert np.any(t2['a'] != DATA['a'])

    def test_select_columns(self):
        """Select columns returns COPY of data and all column
        attributes"""
        t2 = self.t['a', 'c']
        assert np.all(t2['a'] == DATA['a'])
        assert np.all(t2['c'] == DATA['c'])
        assert t2['a'].attrs_equal(COLS[0])
        assert t2['c'].attrs_equal(COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t._data == DATA)
        assert np.any(t2['a'] != DATA['a'])

    def test_select_bad_column(self):
        """Select column name that does not exist"""
        with pytest.raises(ValueError):
            self.t['a', 1]
