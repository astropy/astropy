# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" Verify item access API in:
https://github.com/astropy/astropy/wiki/Table-item-access-definition
"""
from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False] if numpy_lt_1p5 else [False, True])
def table_data(request):
    class TableData:
        def __init__(self, request):
            self.Table = MaskedTable if request.param else table.Table
            self.Column = table.MaskedColumn if request.param else table.Column
            self.COLS = [
                self.Column(name='a', data=[1, 2, 3], description='da',
                            format='fa', meta={'ma': 1}, units='ua'),
                self.Column(name='b', data=[4, 5, 6], description='db',
                            format='fb', meta={'mb': 1}, units='ub'),
                self.Column(name='c', data=[7, 8, 9], description='dc',
                            format='fc', meta={'mc': 1}, units='ub')]
            self.DATA = self.Table(self.COLS)
    return TableData(request)


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
        assert self.tc['a'].format == 'fa'
        assert self.tc['a'].meta == {'ma': 1}
        assert self.tc['a'].units == 'ua'
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
        assert self.tc[1].format == 'fb'
        assert self.tc[1].meta == {'mb': 1}
        assert self.tc[1].units == 'ub'
        assert self.tc[1].attrs_equal(table_data.COLS[1])
        assert isinstance(self.tc[1], table_data.Column)

        assert self.tc[2].units == 'ub'

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

    def test_column(self, table_data):
        """Column access returns REFERENCE to data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        a = self.t['a']
        assert a[1] == 2
        a[1] = 0
        assert self.t['a'][1] == 0

    def test_row(self, table_data):
        """Row  access returns REFERENCE to data"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        row = self.t[1]
        assert row['a'] == 2
        assert row[1] == 5
        assert row.columns['a'].attrs_equal(table_data.COLS[0])
        assert row.columns['b'].attrs_equal(table_data.COLS[1])
        assert row.columns['c'].attrs_equal(table_data.COLS[2])
        row[1] = 0
        assert row[1] == 0
        if table_data.Table is not MaskedTable:
            # numpy.core.ma.mvoid makes a copy so this test is skipped for masked table
            assert self.t['b'][1] == 0

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
        assert np.all(self.t._data == table_data.DATA)
        assert np.any(t2['a'] != table_data.DATA['a'][slice])
        assert t2.masked == self.t.masked
        assert t2._column_class == self.t._column_class

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
        assert np.all(self.t._data == table_data.DATA)
        assert np.any(t2['a'] != table_data.DATA['a'][slice])
        assert t2.masked == self.t.masked
        assert t2._column_class == self.t._column_class


    def test_select_columns(self, table_data):
        """Select columns returns COPY of data and all column
        attributes"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        t2 = self.t['a', 'c']
        assert np.all(t2['a'] == table_data.DATA['a'])
        assert np.all(t2['c'] == table_data.DATA['c'])
        assert t2['a'].attrs_equal(table_data.COLS[0])
        assert t2['c'].attrs_equal(table_data.COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t._data == table_data.DATA)
        assert np.any(t2['a'] != table_data.DATA['a'])
        assert t2.masked == self.t.masked
        assert t2._column_class == self.t._column_class


    def test_np_where(self, table_data):
        """Select rows using output of np.where"""
        t = table_data.Table(table_data.COLS)
        # Select last two rows
        rows = np.where(t['a'] > 1.5)
        t2 = t[rows]
        assert np.all(t2['a'] == [2, 3])
        assert np.all(t2['b'] == [5, 6])

        # Select no rows
        rows = np.where(t['a'] > 100)
        t2 = t[rows]
        assert len(t2) == 0

    def test_select_bad_column(self, table_data):
        """Select column name that does not exist"""
        self.t = table_data.Table(table_data.COLS)
        self.tc = self.t.columns

        with pytest.raises(ValueError):
            self.t['a', 1]
