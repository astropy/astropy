# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" Verify item access API in:
https://github.com/astropy/astropy/wiki/Table-item-access-definition
"""
from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')

# Dummy init of Table, DATA for pyflakes and to be sure test fixture is working
Table = None
Column = None
DATA = None
COLS = None


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False] if numpy_lt_1p5 else [False, True])
def set_global_Table_DATA(request):
    global Table, Column, DATA, COLS

    Table = MaskedTable if request.param else table.Table
    Column = table.MaskedColumn if request.param else table.Column
    COLS = [Column(name='a', data=[1, 2, 3], description='da',
                   format='fa', meta={'ma': 1}, units='ua'),
            Column(name='b', data=[4, 5, 6], description='db',
                   format='fb', meta={'mb': 1}, units='ub'),
            Column(name='c', data=[7, 8, 9], description='dc',
                   format='fc', meta={'mc': 1}, units='ub')]
    DATA = Table(COLS)


@pytest.mark.usefixtures('set_global_Table_DATA')
class BaseTestItems():

    def setup_method(self, method):
        pass


@pytest.mark.usefixtures('set_global_Table_DATA')
class TestTableColumnsItems(BaseTestItems):

    def test_by_name(self):
        """Access TableColumns by name and show that item access returns
        a Column that refers to underlying table data"""
        self.t = Table(COLS)
        self.tc = self.t.columns

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
        self.t = Table(COLS)
        self.tc = self.t.columns

        assert self.tc[1].name == 'b'
        assert np.all(self.tc[1].data == COLS[1].data)
        assert self.tc[1].description == 'db'
        assert self.tc[1].format == 'fb'
        assert self.tc[1].meta == {'mb': 1}
        assert self.tc[1].units == 'ub'
        assert self.tc[1].attrs_equal(COLS[1])
        assert isinstance(self.tc[1], Column)

        assert self.tc[2].units == 'ub'

        self.tc[1][1] = 0
        assert self.t['b'][1] == 0

    def test_mult_columns(self):
        """Access TableColumns with "fancy indexing" and showed that returned
        TableColumns object still references original data"""
        self.t = Table(COLS)
        self.tc = self.t.columns

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
        self.t = Table(COLS)
        self.tc = self.t.columns

        tc2 = self.tc[1:3]
        assert tc2[1].name == 'c'
        assert tc2[1][1] == 8
        assert tc2[0].name == 'b'
        assert tc2[0][1] == 5

        tc2['c'][1] = 0
        assert self.tc['c'][1] == 0
        assert self.t['c'][1] == 0


@pytest.mark.usefixtures('set_global_Table_DATA')
class TestTableItems(BaseTestItems):

    def test_column(self):
        """Column access returns REFERENCE to data"""
        self.t = Table(COLS)
        self.tc = self.t.columns

        a = self.t['a']
        assert a[1] == 2
        a[1] = 0
        assert self.t['a'][1] == 0

    def test_row(self):
        """Row  access returns REFERENCE to data"""
        self.t = Table(COLS)
        self.tc = self.t.columns

        row = self.t[1]
        assert row['a'] == 2
        assert row[1] == 5
        assert row.columns['a'].attrs_equal(COLS[0])
        assert row.columns['b'].attrs_equal(COLS[1])
        assert row.columns['c'].attrs_equal(COLS[2])
        row[1] = 0
        assert row[1] == 0
        if Table is not MaskedTable:
            # numpy.core.ma.mvoid makes a copy so this test is skipped for masked table
            assert self.t['b'][1] == 0

    def test_table_slice(self):
        """Table slice returns REFERENCE to data"""
        self.t = Table(COLS)
        self.tc = self.t.columns

        t2 = self.t[1:3]
        assert np.all(t2['a'] == DATA['a'][1:3])
        assert t2['a'].attrs_equal(COLS[0])
        assert t2['b'].attrs_equal(COLS[1])
        assert t2['c'].attrs_equal(COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t['a'] == np.array([1, 0, 3]))

    def test_fancy_index_slice(self):
        """Table fancy slice returns COPY of data"""
        self.t = Table(COLS)
        self.tc = self.t.columns

        slice = np.array([0, 2])
        t2 = self.t[slice]
        assert np.all(t2['a'] == DATA['a'][slice])
        assert t2['a'].attrs_equal(COLS[0])
        assert t2['b'].attrs_equal(COLS[1])
        assert t2['c'].attrs_equal(COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t._data == DATA)
        assert np.any(t2['a'] != DATA['a'][slice])

    def test_list_index_slice(self):
        """Table list index slice returns COPY of data"""
        self.t = Table(COLS)
        self.tc = self.t.columns

        slice = [0, 2]
        t2 = self.t[slice]
        assert np.all(t2['a'] == DATA['a'][slice])
        assert t2['a'].attrs_equal(COLS[0])
        assert t2['b'].attrs_equal(COLS[1])
        assert t2['c'].attrs_equal(COLS[2])
        t2['a'][0] = 0
        assert np.all(self.t._data == DATA)
        assert np.any(t2['a'] != DATA['a'][slice])

    def test_select_columns(self):
        """Select columns returns COPY of data and all column
        attributes"""
        self.t = Table(COLS)
        self.tc = self.t.columns

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
        self.t = Table(COLS)
        self.tc = self.t.columns

        with pytest.raises(ValueError):
            self.t['a', 1]
