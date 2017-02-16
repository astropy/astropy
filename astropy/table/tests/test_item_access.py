# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

""" Verify item access API in:
https://github.com/astropy/astropy/wiki/Table-item-access-definition
"""
import numpy as np

from ...units import Quantity
from ...tests.helper import pytest
from ..table import QTable


def attrs_equal(col1, col2):
    try:
        return col1.attrs_equal(col2)
    except AttributeError:
        pass

    # copied from Column,attrs_equal
    attrs = ('name', 'unit', 'dtype', 'format', 'description', 'meta')
    return all(getattr(col1.info, x) == getattr(col2.info, x) for x in attrs)


class BaseTestItems():
    def _setup(self, table_types):
        # use reals to ensure no dtype conversion is done with Quantity.
        self.cols = [
            table_types.Column(name='a', data=[1., 2., 3.], description='da',
                               format='d', meta={'ma': 1}, unit='ua'),
            table_types.Column(name='b', data=[4., 5., 6.], description='db',
                               format='f', meta={'mb': 1}, unit='ub'),
            table_types.Column(name='c', data=[7., 8., 9.], description='dc',
                               format='g', meta={'mc': 1}, unit='ub')]
        self.table = table_types.Table(self.cols)
        if table_types.Table is not QTable:
            self.expected_column_cls = table_types.Column
        else:
            self.expected_column_cls = Quantity

    def q_if_needed(self, item, name='a'):
        """Make an item into a quantity if we're comparing to a quantity."""
        col = self.table[name]
        return item * col.unit if isinstance(col, Quantity) else item


@pytest.mark.usefixtures('table_types')
class TestTableColumnsItems(BaseTestItems):

    def test_by_name(self, table_types):
        """Access TableColumns by name and show that item access returns
        a Column that refers to underlying table data"""
        self._setup(table_types)
        tc = self.table.columns

        assert tc['a'].info.name == 'a'
        assert tc['a'][1] == self.q_if_needed(2)
        assert tc['a'].info.description == 'da'
        assert tc['a'].info.format == 'd'
        assert tc['a'].info.meta == {'ma': 1}
        assert tc['a'].info.unit == 'ua'
        assert attrs_equal(tc['a'], self.cols[0])
        assert isinstance(tc['a'], self.expected_column_cls)

        tc['b'][1] = 0
        assert self.table['b'][1] == 0

    def test_by_position(self, table_types):
        """Access TableColumns by position and show that item access returns
        a Column that refers to underlying table data"""
        self._setup(table_types)
        tc = self.table.columns

        assert tc[1].info.name == 'b'
        assert np.all(tc[1].data == self.cols[1].data)
        assert tc[1].info.description == 'db'
        assert tc[1].info.format == 'f'
        assert tc[1].info.meta == {'mb': 1}
        assert tc[1].info.unit == 'ub'
        assert attrs_equal(tc[1], self.cols[1])
        assert isinstance(tc[1], self.expected_column_cls)

        assert tc[2].unit == 'ub'

        tc[1][1] = 0
        assert self.table['b'][1] == 0

    def test_mult_columns(self, table_types):
        """Access TableColumns with "fancy indexing" and showed that returned
        TableColumns object still references original data"""
        self._setup(table_types)
        tc = self.table.columns

        tc2 = tc['b', 'c']
        assert tc2[1].info.name == 'c'
        assert tc2[1][1] == self.q_if_needed(8, 'c')
        assert tc2[0].info.name == 'b'
        assert tc2[0][1] == self.q_if_needed(5, 'b')

        tc2['c'][1] = 0
        assert tc['c'][1] == 0
        assert self.table['c'][1] == 0

    def test_column_slice(self, table_types):
        """Access TableColumns with slice and showed that returned
        TableColumns object still references original data"""
        self._setup(table_types)
        tc = self.table.columns

        tc2 = tc[1:3]
        assert tc2[1].info.name == 'c'
        assert tc2[1][1] == self.q_if_needed(8, 'c')
        assert tc2[0].info.name == 'b'
        assert tc2[0][1] == self.q_if_needed(5, 'b')

        tc2['c'][1] = 0
        assert tc['c'][1] == 0
        assert self.table['c'][1] == 0


@pytest.mark.usefixtures('table_types')
class TestTableItems(BaseTestItems):

    def test_column(self, table_types):
        """Column access returns REFERENCE to data"""
        self._setup(table_types)

        a = self.table['a']
        assert a[1] == self.q_if_needed(2)
        a[1] = 0
        assert self.table['a'][1] == 0

    def test_row(self, table_types):
        """Row  access returns REFERENCE to data"""
        self._setup(table_types)

        row = self.table[1]
        assert row['a'] == self.q_if_needed(2)
        assert row[1] == self.q_if_needed(5, 'b')
        assert attrs_equal(row.columns['a'], self.cols[0])
        assert attrs_equal(row.columns['b'], self.cols[1])
        assert attrs_equal(row.columns['c'], self.cols[2])

        # Check that setting by col index sets the table and row value
        row[1] = 0
        assert row[1] == 0
        assert row['b'] == 0
        assert self.table['b'][1] == 0
        assert self.table[1]['b'] == 0

        # Check that setting by col name sets the table and row value
        row['a'] = 0
        assert row[0] == 0
        assert row['a'] == 0
        assert self.table['a'][1] == 0
        assert self.table[1]['a'] == 0

    def test_empty_iterable_item(self, table_types):
        """
        Table item access with [], (), or np.array([]) returns the same table
        with no rows.
        """
        self._setup(table_types)
        for item in [], (), np.array([]):
            t2 = self.table[item]
            assert not t2
            assert len(t2) == 0
            assert attrs_equal(t2['a'], self.cols[0])
            assert attrs_equal(t2['b'], self.cols[1])
            assert attrs_equal(t2['c'], self.cols[2])

    def test_table_slice(self, table_types):
        """Table slice returns REFERENCE to data"""
        self._setup(table_types)

        t2 = self.table[1:3]
        assert np.all(t2['a'] == self.table['a'][1:3])
        assert attrs_equal(t2['a'], self.cols[0])
        assert attrs_equal(t2['b'], self.cols[1])
        assert attrs_equal(t2['c'], self.cols[2])
        t2['a'][0] = 0
        assert np.all(self.table['a'] == self.q_if_needed(np.array([1, 0, 3])))
        assert t2.masked == self.table.masked
        assert t2._column_class == self.table._column_class
        assert isinstance(t2, table_types.Table)

    def test_fancy_index_slice(self, table_types):
        """Table fancy slice returns COPY of data"""
        self._setup(table_types)

        slc = np.array([0, 2])
        t2 = self.table[slc]
        assert np.all(t2['a'] == self.table['a'][slc])
        assert attrs_equal(t2['a'], self.cols[0])
        assert attrs_equal(t2['b'], self.cols[1])
        assert attrs_equal(t2['c'], self.cols[2])
        t2['a'][0] = 0

        assert np.all(self.table.as_array() == self.table)
        assert np.any(t2['a'] != self.table['a'][slc])
        assert t2.masked == self.table.masked
        assert t2._column_class == self.table._column_class
        assert isinstance(t2, table_types.Table)

    def test_list_index_slice(self, table_types):
        """Table list index slice returns COPY of data"""
        self._setup(table_types)

        slc = [0, 2]
        t2 = self.table[slc]
        assert np.all(t2['a'] == self.table['a'][slc])
        assert attrs_equal(t2['a'], self.cols[0])
        assert attrs_equal(t2['b'], self.cols[1])
        assert attrs_equal(t2['c'], self.cols[2])
        t2['a'][0] = 0

        assert np.all(self.table.as_array() == self.table)
        assert np.any(t2['a'] != self.table['a'][slc])
        assert t2.masked == self.table.masked
        assert t2._column_class == self.table._column_class
        assert isinstance(t2, table_types.Table)

    def test_select_columns(self, table_types):
        """Select columns returns COPY of data and all column
        attributes"""
        self._setup(table_types)

        # try both lists and tuples
        for columns in (('a', 'c'), ['a', 'c']):
            t2 = self.table[columns]
            assert np.all(t2['a'] == self.table['a'])
            assert np.all(t2['c'] == self.table['c'])
            assert attrs_equal(t2['a'], self.cols[0])
            assert attrs_equal(t2['c'], self.cols[2])
            t2['a'][0] = 0
            assert np.all(self.table.as_array() == self.table)
            assert np.any(t2['a'] != self.table['a'])
            assert t2.masked == self.table.masked
            assert t2._column_class == self.table._column_class

    def test_select_columns_fail(self, table_types):
        """Selecting a column that doesn't exist fails"""
        self._setup(table_types)

        with pytest.raises(ValueError) as err:
            self.table[['xxxx']]
        assert 'Slice name(s) xxxx not valid column name(s)' in str(err)

        with pytest.raises(ValueError) as err:
            self.table[['xxxx', 'yyyy']]
        assert 'Slice name(s) xxxx, yyyy not valid column name(s)' in str(err)

    def test_np_where(self, table_types):
        """Select rows using output of np.where"""
        self._setup(table_types)

        # Select last two rows
        rows = np.where(self.table['a'] > self.q_if_needed(1.5))
        t2 = self.table[rows]
        assert np.all(t2['a'] == self.q_if_needed([2, 3]))
        assert np.all(t2['b'] == self.q_if_needed([5, 6], 'b'))
        assert isinstance(t2, table_types.Table)

        # Select no rows
        rows = np.where(self.table['a'] > self.q_if_needed(100))
        t2 = self.table[rows]
        assert len(t2) == 0
        assert isinstance(t2, table_types.Table)

    def test_np_integers(self, table_types):
        """
        Select rows using numpy integers.  This is a regression test for a
        py 3.3 failure mode
        """
        self._setup(table_types)
        idxs = np.random.randint(len(self.table), size=2)
        self.table[idxs[1]]

    def test_select_bad_column(self, table_types):
        """Select column name that does not exist"""
        self._setup(table_types)

        with pytest.raises(ValueError):
            self.table['a', 1]
