# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import copy
import gc
import sys

import numpy as np
from numpy.testing import assert_allclose

from ...extern import six
from ...io import fits
from ...tests.helper import pytest, assert_follows_unicode_guidelines
from ...utils.data import get_pkg_data_filename
from ... import table
from ... import units as u
from .conftest import MaskedTable

try:
    import pandas
except ImportError:
    HAS_PANDAS = False
else:
    HAS_PANDAS = True


class SetupData(object):
    def _setup(self, table_types):
        self._table_type = table_types.Table
        self._column_type = table_types.Column

    @property
    def a(self):
        if self._column_type is not None:
            if not hasattr(self, '_a'):
                self._a = self._column_type(
                    [1, 2, 3], name='a', format='%d',
                    meta={'aa': [0, 1, 2, 3, 4]})
            return self._a

    @property
    def b(self):
        if self._column_type is not None:
            if not hasattr(self, '_b'):
                self._b = self._column_type(
                    [4, 5, 6], name='b', format='%d', meta={'aa': 1})
            return self._b

    @property
    def c(self):
        if self._column_type is not None:
            if not hasattr(self, '_c'):
                self._c = self._column_type([7, 8, 9], 'c')
            return self._c

    @property
    def d(self):
        if self._column_type is not None:
            if not hasattr(self, '_d'):
                self._d = self._column_type([7, 8, 7], 'd')
            return self._d

    @property
    def obj(self):
        if self._column_type is not None:
            if not hasattr(self, '_obj'):
                self._obj = self._column_type([1, 'string', 3], 'obj', dtype='O')
            return self._obj

    @property
    def t(self):
        if self._table_type is not None:
            if not hasattr(self, '_t'):
                self._t = self._table_type([self.a, self.b])
            return self._t


@pytest.mark.usefixtures('table_types')
class TestSetTableColumn(SetupData):

    def test_set_row(self, table_types):
        """Set a row from a tuple of values"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t[1] = (20, 21)
        assert t['a'][0] == 1
        assert t['a'][1] == 20
        assert t['a'][2] == 3
        assert t['b'][0] == 4
        assert t['b'][1] == 21
        assert t['b'][2] == 6

    def test_set_row_existing(self, table_types):
        """Set a row from another existing row"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t[0] = t[1]
        assert t[0][0] == 2
        assert t[0][1] == 5

    def test_set_row_fail_1(self, table_types):
        """Set a row from an incorrectly-sized or typed set of values"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        with pytest.raises(ValueError):
            t[1] = (20, 21, 22)
        with pytest.raises(TypeError):
            t[1] = 0

    def test_set_row_fail_2(self, table_types):
        """Set a row from an incorrectly-typed tuple of values"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        with pytest.raises(ValueError):
            t[1] = ('abc', 'def')

    def test_set_new_col_new_table(self, table_types):
        """Create a new column in empty table using the item access syntax"""
        self._setup(table_types)
        t = table_types.Table()
        t['aa'] = self.a
        # Test that the new column name is 'aa' and that the values match
        assert np.all(t['aa'] == self.a)
        assert t.colnames == ['aa']

    def test_set_new_col_new_table_quantity(self, table_types):
        """Create a new column (from a quantity) in empty table using the item access syntax"""
        self._setup(table_types)
        t = table_types.Table()

        t['aa'] = np.array([1,2,3]) * u.m
        assert np.all(t['aa'] == np.array([1,2,3]))
        assert t['aa'].unit == u.m

        t['bb'] = 3 * u.m
        assert np.all(t['bb'] == 3)
        assert t['bb'].unit == u.m

    def test_set_new_col_existing_table(self, table_types):
        """Create a new column in an existing table using the item access syntax"""
        self._setup(table_types)
        t = table_types.Table([self.a])

        # Add a column
        t['bb'] = self.b
        assert np.all(t['bb'] == self.b)
        assert t.colnames == ['a', 'bb']
        assert t['bb'].meta == self.b.meta
        assert t['bb'].format == self.b.format

        # Add another column
        t['c'] = t['a']
        assert np.all(t['c'] == t['a'])
        assert t.colnames == ['a', 'bb', 'c']
        assert t['c'].meta == t['a'].meta
        assert t['c'].format == t['a'].format

        # Add a multi-dimensional column
        t['d'] = table_types.Column(np.arange(12).reshape(3, 2, 2))
        assert t['d'].shape == (3, 2, 2)
        assert t['d'][0, 0, 1] == 1

        # Add column from a list
        t['e'] = ['hello', 'the', 'world']
        assert np.all(t['e'] == np.array(['hello', 'the', 'world']))

        # Make sure setting existing column still works
        t['e'] = ['world', 'hello', 'the']
        assert np.all(t['e'] == np.array(['world', 'hello', 'the']))

        # Add a column via broadcasting
        t['f'] = 10
        assert np.all(t['f'] == 10)

        # Add a column from a Quantity
        t['g'] = np.array([1,2,3]) * u.m
        assert np.all(t['g'].data == np.array([1,2,3]))
        assert t['g'].unit == u.m

        # Add a column from a (scalar) Quantity
        t['g'] = 3 * u.m
        assert np.all(t['g'].data == 3)
        assert t['g'].unit == u.m

    def test_set_new_unmasked_col_existing_table(self, table_types):
        """Create a new column in an existing table using the item access syntax"""
        self._setup(table_types)
        t = table_types.Table([self.a])  # masked or unmasked
        b = table.Column(name='b', data=[1, 2, 3])  # unmasked
        t['b'] = b
        assert np.all(t['b'] == b)

    def test_set_new_masked_col_existing_table(self, table_types):
        """Create a new column in an existing table using the item access syntax"""
        self._setup(table_types)
        t = table_types.Table([self.a])  # masked or unmasked
        b = table.MaskedColumn(name='b', data=[1, 2, 3])  # masked
        t['b'] = b
        assert np.all(t['b'] == b)

    def test_set_new_col_existing_table_fail(self, table_types):
        """Generate failure when creating a new column using the item access syntax"""
        self._setup(table_types)
        t = table_types.Table([self.a])
        # Wrong size
        with pytest.raises(ValueError):
            t['b'] = [1, 2]


@pytest.mark.usefixtures('table_types')
class TestEmptyData():

    def test_1(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', dtype=int, length=100))
        assert len(t['a']) == 100

    def test_2(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', dtype=int, shape=(3, ), length=100))
        assert len(t['a']) == 100

    def test_3(self, table_types):
        t = table_types.Table()  # length is not given
        t.add_column(table_types.Column(name='a', dtype=int))
        assert len(t['a']) == 0

    def test_4(self, table_types):
        t = table_types.Table()  # length is not given
        t.add_column(table_types.Column(name='a', dtype=int, shape=(3, 4)))
        assert len(t['a']) == 0

    def test_5(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a'))  # dtype is not specified
        assert len(t['a']) == 0

    def test_add_via_setitem_and_slice(self, table_types):
        """Test related to #3023 where a MaskedColumn is created with name=None
        and then gets changed to name='a'.  After PR #2790 this test fails
        without the #3023 fix."""
        t = table_types.Table()
        t['a'] = table_types.Column([1, 2, 3])
        t2 = t[:]
        assert t2.colnames == t.colnames


@pytest.mark.usefixtures('table_types')
class TestNewFromColumns():

    def test_simple(self, table_types):
        cols = [table_types.Column(name='a', data=[1, 2, 3]),
                table_types.Column(name='b', data=[4, 5, 6], dtype=np.float32)]
        t = table_types.Table(cols)
        assert np.all(t['a'].data == np.array([1, 2, 3]))
        assert np.all(t['b'].data == np.array([4, 5, 6], dtype=np.float32))
        assert type(t['b'][1]) == np.float32

    def test_from_np_array(self, table_types):
        cols = [table_types.Column(name='a', data=np.array([1, 2, 3], dtype=np.int64),
                       dtype=np.float64),
                table_types.Column(name='b', data=np.array([4, 5, 6], dtype=np.float32))]
        t = table_types.Table(cols)
        assert np.all(t['a'] == np.array([1, 2, 3], dtype=np.float64))
        assert np.all(t['b'] == np.array([4, 5, 6], dtype=np.float32))
        assert type(t['a'][1]) == np.float64
        assert type(t['b'][1]) == np.float32

    def test_size_mismatch(self, table_types):
        cols = [table_types.Column(name='a', data=[1, 2, 3]),
                table_types.Column(name='b', data=[4, 5, 6, 7])]
        with pytest.raises(ValueError):
            table_types.Table(cols)

    def test_name_none(self, table_types):
        """Column with name=None can init a table whether or not names are supplied"""
        c = table_types.Column(data=[1, 2], name='c')
        d = table_types.Column(data=[3, 4])
        t = table_types.Table([c, d], names=(None, 'd'))
        assert t.colnames == ['c', 'd']
        t = table_types.Table([c, d])
        assert t.colnames == ['c', 'col1']

@pytest.mark.usefixtures('table_types')
class TestReverse():

    def test_reverse(self, table_types):
        t = table_types.Table([[1, 2, 3],
                   ['a', 'b', 'cc']])
        t.reverse()
        assert np.all(t['col0'] == np.array([3, 2, 1]))
        assert np.all(t['col1'] == np.array(['cc', 'b', 'a']))

        t2 = table_types.Table(t, copy=False)
        assert np.all(t2['col0'] == np.array([3, 2, 1]))
        assert np.all(t2['col1'] == np.array(['cc', 'b', 'a']))

        t2 = table_types.Table(t, copy=True)
        assert np.all(t2['col0'] == np.array([3, 2, 1]))
        assert np.all(t2['col1'] == np.array(['cc', 'b', 'a']))

        t2.sort('col0')
        assert np.all(t2['col0'] == np.array([1, 2, 3]))
        assert np.all(t2['col1'] == np.array(['a', 'b', 'cc']))

    def test_reverse_big(self, table_types):
        x = np.arange(10000)
        y = x + 1
        t = table_types.Table([x, y], names=('x', 'y'))
        t.reverse()
        assert np.all(t['x'] == x[::-1])
        assert np.all(t['y'] == y[::-1])


@pytest.mark.usefixtures('table_types')
class TestColumnAccess():

    def test_1(self, table_types):
        t = table_types.Table()
        with pytest.raises(KeyError):
            t['a']

    def test_2(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[1, 2, 3]))
        assert np.all(t['a'] == np.array([1, 2, 3]))
        with pytest.raises(KeyError):
            t['b']  # column does not exist


@pytest.mark.usefixtures('table_types')
class TestAddLength(SetupData):

    def test_right_length(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.add_column(self.b)

    def test_too_long(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        with pytest.raises(ValueError):
            t.add_column(table_types.Column(name='b', data=[4, 5, 6, 7]))  # data too long

    def test_too_short(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        with pytest.raises(ValueError):
            t.add_column(table_types.Column(name='b', data=[4, 5]))  # data too short


@pytest.mark.usefixtures('table_types')
class TestAddPosition(SetupData):

    def test_1(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a, 0)

    def test_2(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a, 1)

    def test_3(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a, -1)

    def test_5(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        with pytest.raises(ValueError):
            t.index_column('b')

    def test_6(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a)
        t.add_column(self.b)
        assert t.columns.keys() == ['a', 'b']

    def test_7(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.add_column(self.b, t.index_column('a'))
        assert t.columns.keys() == ['b', 'a']

    def test_8(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.add_column(self.b, t.index_column('a') + 1)
        assert t.columns.keys() == ['a', 'b']

    def test_9(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a)
        t.add_column(self.b, t.index_column('a') + 1)
        t.add_column(self.c, t.index_column('b'))
        assert t.columns.keys() == ['a', 'c', 'b']

    def test_10(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a)
        ia = t.index_column('a')
        t.add_column(self.b, ia + 1)
        t.add_column(self.c, ia)
        assert t.columns.keys() == ['c', 'a', 'b']


@pytest.mark.usefixtures('table_types')
class TestInitFromTable(SetupData):

    def test_from_table_cols(self, table_types):
        """Ensure that using cols from an existing table gives
        a clean copy.
        """
        self._setup(table_types)
        t = self.t
        cols = t.columns
        # Construct Table with cols via Table._new_from_cols
        t2a = table_types.Table([cols['a'], cols['b'], self.c])

        # Construct with add_column
        t2b = table_types.Table()
        t2b.add_column(cols['a'])
        t2b.add_column(cols['b'])
        t2b.add_column(self.c)

        t['a'][1] = 20
        t['b'][1] = 21
        for t2 in [t2a, t2b]:
            t2['a'][2] = 10
            t2['b'][2] = 11
            t2['c'][2] = 12
            t2.columns['a'].meta['aa'][3] = 10
            assert np.all(t['a'] == np.array([1, 20, 3]))
            assert np.all(t['b'] == np.array([4, 21, 6]))
            assert np.all(t2['a'] == np.array([1, 2, 10]))
            assert np.all(t2['b'] == np.array([4, 5, 11]))
            assert np.all(t2['c'] == np.array([7, 8, 12]))
            assert t2['a'].name == 'a'
            assert t2.columns['a'].meta['aa'][3] == 10
            assert t.columns['a'].meta['aa'][3] == 3


@pytest.mark.usefixtures('table_types')
class TestAddColumns(SetupData):

    def test_add_columns1(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_columns([self.a, self.b, self.c])
        assert t.colnames == ['a', 'b', 'c']

    def test_add_columns2(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.add_columns([self.c, self.d])
        assert t.colnames == ['a', 'b', 'c', 'd']
        assert np.all(t['c'] == np.array([7, 8, 9]))

    def test_add_columns3(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.add_columns([self.c, self.d], indexes=[1, 0])
        assert t.colnames == ['d', 'a', 'c', 'b']

    def test_add_columns4(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.add_columns([self.c, self.d], indexes=[0, 0])
        assert t.colnames == ['c', 'd', 'a', 'b']

    def test_add_columns5(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.add_columns([self.c, self.d], indexes=[2, 2])
        assert t.colnames == ['a', 'b', 'c', 'd']

    def test_add_duplicate_column(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a)
        with pytest.raises(ValueError):
            t.add_column(table_types.Column(name='a', data=[0, 1, 2]))
        t.add_column(table_types.Column(name='a', data=[0, 1, 2]),
                     rename_duplicate=True)
        t.add_column(self.b)
        t.add_column(self.c)
        assert t.colnames == ['a', 'a_1', 'b', 'c']
        t.add_column(table_types.Column(name='a', data=[0, 1, 2]),
                     rename_duplicate=True)
        assert t.colnames == ['a', 'a_1', 'b', 'c', 'a_2']

        # test adding column from a separate Table
        t1 = table_types.Table()
        t1.add_column(self.a)
        with pytest.raises(ValueError):
            t.add_column(t1['a'])
        t.add_column(t1['a'], rename_duplicate=True)

        t1['a'][0] = 100  # Change original column
        assert t.colnames == ['a', 'a_1', 'b', 'c', 'a_2', 'a_3']
        assert t1.colnames == ['a']

        # Check new column didn't change (since name conflict forced a copy)
        assert t['a_3'][0] == self.a[0]

    def test_add_duplicate_columns(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b, self.c])
        with pytest.raises(ValueError):
            t.add_columns([table_types.Column(name='a', data=[0, 1, 2]), table_types.Column(name='b', data=[0, 1, 2])])
        t.add_columns([table_types.Column(name='a', data=[0, 1, 2]),
                       table_types.Column(name='b', data=[0, 1, 2])],
                      rename_duplicate=True)
        t.add_column(self.d)
        assert t.colnames == ['a', 'b', 'c', 'a_1', 'b_1', 'd']


@pytest.mark.usefixtures('table_types')
class TestAddRow(SetupData):

    @property
    def b(self):
        if self._column_type is not None:
            if not hasattr(self, '_b'):
                self._b = self._column_type(name='b', data=[4.0, 5.1, 6.2])
            return self._b

    @property
    def c(self):
        if self._column_type is not None:
            if not hasattr(self, '_c'):
                self._c = self._column_type(name='c', data=['7', '8', '9'])
            return self._c

    @property
    def d(self):
        if self._column_type is not None:
            if not hasattr(self, '_d'):
                self._d = self._column_type(name='d', data=[[1, 2], [3, 4], [5, 6]])
            return self._d

    @property
    def t(self):
        if self._table_type is not None:
            if not hasattr(self, '_t'):
                self._t = self._table_type([self.a, self.b, self.c])
            return self._t

    def test_add_none_to_empty_table(self, table_types):
        self._setup(table_types)
        t = table_types.Table(names=('a', 'b', 'c'), dtype=('(2,)i', 'S4', 'O'))
        t.add_row()
        assert np.all(t['a'][0] == [0, 0])
        assert t['b'][0] == b''
        assert t['c'][0] == 0
        t.add_row()
        assert np.all(t['a'][1] == [0, 0])
        assert t['b'][1] == b''
        assert t['c'][1] == 0

    def test_add_stuff_to_empty_table(self, table_types):
        self._setup(table_types)
        t = table_types.Table(names=('a', 'b', 'obj'), dtype=('(2,)i', 'S8', 'O'))
        t.add_row([[1, 2], 'hello', 'world'])
        assert np.all(t['a'][0] == [1, 2])
        assert t['b'][0] == b'hello'
        assert t['obj'][0] == 'world'
        # Make sure it is not repeating last row but instead
        # adding zeros (as documented)
        t.add_row()
        assert np.all(t['a'][1] == [0, 0])
        assert t['b'][1] == b''
        assert t['obj'][1] == 0

    def test_add_table_row(self, table_types):
        self._setup(table_types)
        t = self.t
        t['d'] = self.d
        t2 = table_types.Table([self.a, self.b, self.c, self.d])
        t.add_row(t2[0])
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 1]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 4.0]))
        assert np.all(t['c'] == np.array(['7', '8', '9', '7']))
        assert np.all(t['d'] == np.array([[1, 2], [3, 4], [5, 6], [1, 2]]))

    def test_add_table_row_obj(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b, self.obj])
        t.add_row([1, 4.0, [10]])
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 1]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 4.0]))
        assert np.all(t['obj'] == np.array([1, 'string', 3, [10]], dtype='O'))

    def test_add_with_tuple(self, table_types):
        self._setup(table_types)
        t = self.t
        t.add_row((4, 7.2, '1'))
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 4]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 7.2]))
        assert np.all(t['c'] == np.array(['7', '8', '9', '1']))

    def test_add_with_list(self, table_types):
        self._setup(table_types)
        t = self.t
        t.add_row([4, 7.2, '10'])
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 4]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 7.2]))
        assert np.all(t['c'] == np.array(['7', '8', '9', '1']))

    def test_add_with_dict(self, table_types):
        self._setup(table_types)
        t = self.t
        t.add_row({'a': 4, 'b': 7.2})
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 4]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 7.2]))
        if t.masked:
            assert np.all(t['c'] == np.array(['7', '8', '9', '7']))
        else:
            assert np.all(t['c'] == np.array(['7', '8', '9', '']))

    def test_add_with_none(self, table_types):
        self._setup(table_types)
        t = self.t
        t.add_row()
        assert len(t) == 4
        assert np.all(t['a'].data == np.array([1, 2, 3, 0]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 0.0]))
        assert np.all(t['c'].data == np.array(['7', '8', '9', '']))

    def test_add_missing_column(self, table_types):
        self._setup(table_types)
        t = self.t
        with pytest.raises(ValueError):
            t.add_row({'bad_column': 1})

    def test_wrong_size_tuple(self, table_types):
        self._setup(table_types)
        t = self.t
        with pytest.raises(ValueError):
            t.add_row((1, 2))

    def test_wrong_vals_type(self, table_types):
        self._setup(table_types)
        t = self.t
        with pytest.raises(TypeError):
            t.add_row(1)

    def test_add_row_failures(self, table_types):
        self._setup(table_types)
        t = self.t
        t_copy = table_types.Table(t, copy=True)
        # Wrong number of columns
        try:
            t.add_row([1,2,3,4])
        except ValueError:
            pass
        assert len(t) == 3
        assert np.all(t.as_array() == t_copy.as_array())
        # Wrong data type
        try:
            t.add_row(['one',2,3])
        except ValueError:
            pass
        assert len(t) == 3
        assert np.all(t.as_array() == t_copy.as_array())

    def test_insert_table_row(self, table_types):
        """
        Light testing of Table.insert_row() method.  The deep testing is done via
        the add_row() tests which calls insert_row(index=len(self), ...), so
        here just test that the added index parameter is handled correctly.
        """
        self._setup(table_types)
        row = (10, 40.0, 'x', [10, 20])
        for index in range(-3, 4):
            indices = np.insert(np.arange(3), index, 3)
            t = table_types.Table([self.a, self.b, self.c, self.d])
            t2 = t.copy()
            t.add_row(row)  # By now we know this works
            t2.insert_row(index, row)
            for name in t.colnames:
                if t[name].dtype.kind == 'f':
                    assert np.allclose(t[name][indices], t2[name])
                else:
                    assert np.all(t[name][indices] == t2[name])

        for index in (-4, 4):
            t = table_types.Table([self.a, self.b, self.c, self.d])
            with pytest.raises(IndexError):
                t.insert_row(index, row)


@pytest.mark.usefixtures('table_types')
class TestTableColumn(SetupData):

    def test_column_view(self, table_types):
        self._setup(table_types)
        t = self.t
        a = t.columns['a']
        a[2] = 10
        assert t['a'][2] == 10


@pytest.mark.usefixtures('table_types')
class TestArrayColumns(SetupData):

    def test_1d(self, table_types):
        self._setup(table_types)
        b = table_types.Column(name='b', dtype=int, shape=(2, ), length=3)
        t = table_types.Table([self.a])
        t.add_column(b)
        assert t['b'].shape == (3, 2)
        assert t['b'][0].shape == (2, )

    def test_2d(self, table_types):
        self._setup(table_types)
        b = table_types.Column(name='b', dtype=int, shape=(2, 4), length=3)
        t = table_types.Table([self.a])
        t.add_column(b)
        assert t['b'].shape == (3, 2, 4)
        assert t['b'][0].shape == (2, 4)

    def test_3d(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        b = table_types.Column(name='b', dtype=int, shape=(2, 4, 6), length=3)
        t.add_column(b)
        assert t['b'].shape == (3, 2, 4, 6)
        assert t['b'][0].shape == (2, 4, 6)


@pytest.mark.usefixtures('table_types')
class TestRemove(SetupData):

    @property
    def t(self):
        if self._table_type is not None:
            if not hasattr(self, '_t'):
                self._t = self._table_type([self.a])
            return self._t

    @property
    def t2(self):
        if self._table_type is not None:
            if not hasattr(self, '_t2'):
                self._t2 = self._table_type([self.a, self.b, self.c])
            return self._t2

    def test_1(self, table_types):
        self._setup(table_types)
        self.t.remove_columns('a')
        assert self.t.columns.keys() == []
        assert self.t.as_array() is None

    def test_2(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.remove_columns('a')
        assert self.t.columns.keys() == ['b']
        assert self.t.dtype.names == ('b',)
        assert np.all(self.t['b'] == np.array([4, 5, 6]))

    def test_3(self, table_types):
        """Check remove_columns works for a single column with a name of
        more than one character.  Regression test against #2699"""
        self._setup(table_types)
        self.t['new_column'] = self.t['a']
        assert 'new_column' in self.t.columns.keys()
        self.t.remove_columns('new_column')
        assert 'new_column' not in self.t.columns.keys()

    def test_remove_nonexistent_row(self, table_types):
        self._setup(table_types)
        with pytest.raises(IndexError):
            self.t.remove_row(4)

    def test_remove_row_0(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        self.t.remove_row(0)
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['b'] == np.array([5, 6]))

    def test_remove_row_1(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        self.t.remove_row(1)
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['a'] == np.array([1, 3]))

    def test_remove_row_2(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        self.t.remove_row(2)
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['c'] == np.array([7, 8]))

    def test_remove_row_slice(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        self.t.remove_rows(slice(0, 2, 1))
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['c'] == np.array([9]))

    def test_remove_row_list(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        self.t.remove_rows([0, 2])
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['c'] == np.array([8]))

    def test_remove_row_preserves_meta(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.remove_rows([0, 2])
        assert self.t['a'].meta == {'aa': [0, 1, 2, 3, 4]}
        assert self.t.dtype == np.dtype([(str('a'), 'int'),
                                         (str('b'), 'int')])

    def test_delitem1(self, table_types):
        self._setup(table_types)
        del self.t['a']
        assert self.t.columns.keys() == []
        assert self.t.as_array() is None

    def test_delitem2(self, table_types):
        self._setup(table_types)
        del self.t2['b']
        assert self.t2.colnames == ['a', 'c']

    def test_delitems(self, table_types):
        self._setup(table_types)
        del self.t2['a', 'b']
        assert self.t2.colnames == ['c']

    def test_delitem_fail(self, table_types):
        self._setup(table_types)
        with pytest.raises(KeyError):
            del self.t['d']


@pytest.mark.usefixtures('table_types')
class TestKeep(SetupData):

    def test_1(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.keep_columns([])
        assert t.columns.keys() == []
        assert t.as_array() is None

    def test_2(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.keep_columns('b')
        assert t.columns.keys() == ['b']
        assert t.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([4, 5, 6]))


@pytest.mark.usefixtures('table_types')
class TestRename(SetupData):

    def test_1(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.rename_column('a', 'b')
        assert t.columns.keys() == ['b']
        assert t.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([1, 2, 3]))

    def test_2(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.rename_column('a', 'c')
        t.rename_column('b', 'a')
        assert t.columns.keys() == ['c', 'a']
        assert t.dtype.names == ('c', 'a')
        if t.masked:
            assert t.mask.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))

    def test_rename_by_attr(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t['a'].name = 'c'
        t['b'].name = 'a'
        assert t.columns.keys() == ['c', 'a']
        assert t.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))


@pytest.mark.usefixtures('table_types')
class TestSort():

    def test_single(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[2, 1, 3]))
        t.add_column(table_types.Column(name='b', data=[6, 5, 4]))
        t.add_column(table_types.Column(name='c', data=[(1, 2), (3, 4), (4, 5)]))
        assert np.all(t['a'] == np.array([2, 1, 3]))
        assert np.all(t['b'] == np.array([6, 5, 4]))
        t.sort('a')
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['b'] == np.array([5, 6, 4]))
        assert np.all(t['c'] == np.array([[3, 4],
                                          [1, 2],
                                          [4, 5]]))
        t.sort('b')
        assert np.all(t['a'] == np.array([3, 1, 2]))
        assert np.all(t['b'] == np.array([4, 5, 6]))
        assert np.all(t['c'] == np.array([[4, 5],
                                          [3, 4],
                                          [1, 2]]))

    def test_single_big(self, table_types):
        """Sort a big-ish table with a non-trivial sort order"""
        x = np.arange(10000)
        y = np.sin(x)
        t = table_types.Table([x, y], names=('x', 'y'))
        t.sort('y')
        idx = np.argsort(y)
        assert np.all(t['x'] == x[idx])
        assert np.all(t['y'] == y[idx])

    def test_multiple(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[2, 1, 3, 2, 3, 1]))
        t.add_column(table_types.Column(name='b', data=[6, 5, 4, 3, 5, 4]))
        assert np.all(t['a'] == np.array([2, 1, 3, 2, 3, 1]))
        assert np.all(t['b'] == np.array([6, 5, 4, 3, 5, 4]))
        t.sort(['a', 'b'])
        assert np.all(t['a'] == np.array([1, 1, 2, 2, 3, 3]))
        assert np.all(t['b'] == np.array([4, 5, 3, 6, 4, 5]))
        t.sort(['b', 'a'])
        assert np.all(t['a'] == np.array([2, 1, 3, 1, 3, 2]))
        assert np.all(t['b'] == np.array([3, 4, 4, 5, 5, 6]))

    def test_multiple_with_bytes(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='firstname', data=[b"Max", b"Jo", b"John"]))
        t.add_column(table_types.Column(name='name', data=[b"Miller", b"Miller", b"Jackson"]))
        t.add_column(table_types.Column(name='tel', data=[12, 15, 19]))
        t.sort(['name','firstname'])
        assert np.all([t['firstname'] == np.array([b"John", b"Jo", b"Max"])])
        assert np.all([t['name'] == np.array([b"Jackson", b"Miller", b"Miller"])])
        assert np.all([t['tel'] == np.array([19, 15, 12])])

    def test_multiple_with_unicode(self, table_types):
        # Before Numpy 1.6.2, sorting with multiple column names
        # failed when a unicode column was present.
        t = table_types.Table()
        t.add_column(table_types.Column(
            name='firstname',
            data=[six.text_type(x) for x in ["Max", "Jo", "John"]]))
        t.add_column(table_types.Column(
            name='name',
            data=[six.text_type(x) for x in ["Miller", "Miller", "Jackson"]]))
        t.add_column(table_types.Column(name='tel', data=[12, 15, 19]))
        t.sort(['name','firstname'])
        assert np.all([t['firstname'] == np.array(
            [six.text_type(x) for x in ["John", "Jo", "Max"]])])
        assert np.all([t['name'] == np.array(
            [six.text_type(x) for x in ["Jackson", "Miller", "Miller"]])])
        assert np.all([t['tel'] == np.array([19, 15, 12])])

    def test_argsort(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[2, 1, 3, 2, 3, 1]))
        t.add_column(table_types.Column(name='b', data=[6, 5, 4, 3, 5, 4]))
        assert np.all(t.argsort() == t.as_array().argsort())
        i0 = t.argsort('a')
        i1 = t.as_array().argsort(order=['a'])
        assert np.all(t['a'][i0] == t['a'][i1])
        i0 = t.argsort(['a', 'b'])
        i1 = t.as_array().argsort(order=['a', 'b'])
        assert np.all(t['a'][i0] == t['a'][i1])
        assert np.all(t['b'][i0] == t['b'][i1])

    def test_argsort_bytes(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='firstname', data=[b"Max", b"Jo", b"John"]))
        t.add_column(table_types.Column(name='name', data=[b"Miller", b"Miller", b"Jackson"]))
        t.add_column(table_types.Column(name='tel', data=[12, 15, 19]))
        assert np.all(t.argsort(['name', 'firstname']) == np.array([2, 1, 0]))

    def test_argsort_unicode(self, table_types):
        # Before Numpy 1.6.2, sorting with multiple column names
        # failed when a unicode column was present.
        t = table_types.Table()
        t.add_column(table_types.Column(
            name='firstname',
            data=[six.text_type(x) for x in ["Max", "Jo", "John"]]))
        t.add_column(table_types.Column(
            name='name',
            data=[six.text_type(x) for x in ["Miller", "Miller", "Jackson"]]))
        t.add_column(table_types.Column(name='tel', data=[12, 15, 19]))
        assert np.all(t.argsort(['name', 'firstname']) == np.array([2, 1, 0]))

    def test_rebuild_column_view_then_rename(self, table_types):
        """
        Issue #2039 where renaming fails after any method that calls
        _rebuild_table_column_view (this includes sort and add_row).
        """
        t = table_types.Table([[1]], names=('a',))
        assert t.colnames == ['a']
        assert t.dtype.names == ('a',)

        t.add_row((2,))
        assert t.colnames == ['a']
        assert t.dtype.names == ('a',)

        t.rename_column('a', 'b')
        assert t.colnames == ['b']
        assert t.dtype.names == ('b',)

        t.sort('b')
        assert t.colnames == ['b']
        assert t.dtype.names == ('b',)

        t.rename_column('b', 'c')
        assert t.colnames == ['c']
        assert t.dtype.names == ('c',)


@pytest.mark.usefixtures('table_types')
class TestIterator():

    def test_iterator(self, table_types):
        d = np.array([(2, 1),
                      (3, 6),
                      (4, 5)], dtype=[(str('a'), 'i4'), (str('b'), 'i4')])
        t = table_types.Table(d)
        if t.masked:
            with pytest.raises(ValueError):
                t[0] == d[0]
        else:
            for row, np_row in zip(t, d):
                assert np.all(row == np_row)


@pytest.mark.usefixtures('table_types')
class TestSetMeta():

    def test_set_meta(self, table_types):
        d = table_types.Table(names=('a', 'b'))
        d.meta['a'] = 1
        d.meta['b'] = 1
        d.meta['c'] = 1
        d.meta['d'] = 1
        assert list(d.meta.keys()) == ['a', 'b', 'c', 'd']


@pytest.mark.usefixtures('table_types')
class TestConvertNumpyArray():

    def test_convert_numpy_array(self, table_types):
        d = table_types.Table([[1, 2], [3, 4]], names=('a', 'b'))

        np_data = np.array(d)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d.as_array())
        assert not np_data is d.as_array()
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d.as_array())
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[(str('c'), 'i8'), (str('d'), 'i8')])

    def test_as_array_byteswap(self, table_types):
        """Test for https://github.com/astropy/astropy/pull/4080"""

        byte_orders = ('>', '<')
        native_order = byte_orders[sys.byteorder == 'little']

        for order in byte_orders:
            col = table_types.Column([1.0, 2.0], name='a', dtype=order + 'f8')
            t = table_types.Table([col])
            arr = t.as_array()
            assert arr['a'].dtype.byteorder in (native_order, '=')
            arr = t.as_array(keep_byteorder=True)
            if order == native_order:
                assert arr['a'].dtype.byteorder in (order, '=')
            else:
                assert arr['a'].dtype.byteorder == order

    def test_byteswap_fits_array(self, table_types):
        """
        Test for https://github.com/astropy/astropy/pull/4080, demonstrating
        that FITS tables are converted to native byte order.
        """

        non_native_order = ('>', '<')[sys.byteorder != 'little']

        filename = get_pkg_data_filename('data/tb.fits',
                                         'astropy.io.fits.tests')
        t = table_types.Table.read(filename)
        arr = t.as_array()

        for idx in range(len(arr.dtype)):
            assert arr.dtype[idx].byteorder != non_native_order

        with fits.open(filename) as hdul:
            data = hdul[1].data
            for colname in data.columns.names:
                assert np.all(data[colname] == arr[colname])

            arr2 = t.as_array(keep_byteorder=True)
            for colname in data.columns.names:
                assert (data[colname].dtype.byteorder ==
                        arr2[colname].dtype.byteorder)


def _assert_copies(t, t2, deep=True):
    assert t.colnames == t2.colnames
    np.testing.assert_array_equal(t.as_array(), t2.as_array())
    assert t.meta == t2.meta

    for col, col2 in zip(t.columns.values(), t2.columns.values()):
        if deep:
            assert not np.may_share_memory(col, col2)
        else:
            assert np.may_share_memory(col, col2)


def test_copy():
    t = table.Table([[1, 2, 3], [2, 3, 4]], names=['x', 'y'])
    t2 = t.copy()
    _assert_copies(t, t2)


def test_copy_masked():
    t = table.Table([[1, 2, 3], [2, 3, 4]], names=['x', 'y'], masked=True,
                    meta={'name': 'test'})
    t['x'].mask == [True, False, True]
    t2 = t.copy()
    _assert_copies(t, t2)


def test_copy_protocol():
    t = table.Table([[1, 2, 3], [2, 3, 4]], names=['x', 'y'])

    t2 = copy.copy(t)
    t3 = copy.deepcopy(t)

    _assert_copies(t, t2, deep=False)
    _assert_copies(t, t3)

def test_disallow_inequality_comparisons():
    """
    Regression test for #828 - disallow comparison operators on whole Table
    """

    t = table.Table()

    with pytest.raises(TypeError):
        t > 2

    with pytest.raises(TypeError):
        t < 1.1

    with pytest.raises(TypeError):
        t >= 5.5

    with pytest.raises(TypeError):
        t <= -1.1

def test_equality():

    t = table.Table.read([' a b  c  d',
                          ' 2 c 7.0 0',
                          ' 2 b 5.0 1',
                          ' 2 b 6.0 2',
                          ' 2 a 4.0 3',
                          ' 0 a 0.0 4',
                          ' 1 b 3.0 5',
                          ' 1 a 2.0 6',
                          ' 1 a 1.0 7',
                         ], format='ascii')

    # All rows are equal
    assert np.all(t==t)

    # Assert no rows are different
    assert not np.any(t!=t)

    # Check equality result for a given row
    assert np.all((t == t[3]) == np.array([0,0,0,1,0,0,0,0], dtype=bool))

    # Check inequality result for a given row
    assert np.all((t != t[3]) == np.array([1,1,1,0,1,1,1,1], dtype=bool))

    t2 = table.Table.read([' a b  c  d',
                           ' 2 c 7.0 0',
                           ' 2 b 5.0 1',
                           ' 3 b 6.0 2',
                           ' 2 a 4.0 3',
                           ' 0 a 1.0 4',
                           ' 1 b 3.0 5',
                           ' 1 c 2.0 6',
                           ' 1 a 1.0 7',
                          ], format='ascii')

    # In the above cases, Row.__eq__ gets called, but now need to make sure
    # Table.__eq__ also gets called.
    assert np.all((t == t2) == np.array([1,1,0,1,0,1,0,1], dtype=bool))
    assert np.all((t != t2) == np.array([0,0,1,0,1,0,1,0], dtype=bool))

    # Check that comparing to a structured array works
    assert np.all((t == t2.as_array()) == np.array([1,1,0,1,0,1,0,1], dtype=bool))
    assert np.all((t.as_array() == t2) == np.array([1,1,0,1,0,1,0,1], dtype=bool))


def test_equality_masked():

    t = table.Table.read([' a b  c  d',
                          ' 2 c 7.0 0',
                          ' 2 b 5.0 1',
                          ' 2 b 6.0 2',
                          ' 2 a 4.0 3',
                          ' 0 a 0.0 4',
                          ' 1 b 3.0 5',
                          ' 1 a 2.0 6',
                          ' 1 a 1.0 7',
                         ], format='ascii')

    # Make into masked table
    t = table.Table(t, masked=True)

    # All rows are equal
    assert np.all(t==t)

    # Assert no rows are different
    assert not np.any(t!=t)

    # Check equality result for a given row
    assert np.all((t == t[3]) == np.array([0,0,0,1,0,0,0,0], dtype=bool))

    # Check inequality result for a given row
    assert np.all((t != t[3]) == np.array([1,1,1,0,1,1,1,1], dtype=bool))

    t2 = table.Table.read([' a b  c  d',
                           ' 2 c 7.0 0',
                           ' 2 b 5.0 1',
                           ' 3 b 6.0 2',
                           ' 2 a 4.0 3',
                           ' 0 a 1.0 4',
                           ' 1 b 3.0 5',
                           ' 1 c 2.0 6',
                           ' 1 a 1.0 7',
                          ], format='ascii')

    # In the above cases, Row.__eq__ gets called, but now need to make sure
    # Table.__eq__ also gets called.
    assert np.all((t == t2) == np.array([1,1,0,1,0,1,0,1], dtype=bool))
    assert np.all((t != t2) == np.array([0,0,1,0,1,0,1,0], dtype=bool))

    # Check that masking a value causes the row to differ
    t.mask['a'][0] = True
    assert np.all((t == t2) == np.array([0,1,0,1,0,1,0,1], dtype=bool))
    assert np.all((t != t2) == np.array([1,0,1,0,1,0,1,0], dtype=bool))

    # Check that comparing to a structured array works
    assert np.all((t == t2.as_array()) == np.array([0,1,0,1,0,1,0,1], dtype=bool))


@pytest.mark.xfail
def test_equality_masked_bug():
    """
    This highlights a Numpy bug. Once it works, it can be moved into the
    test_equality_masked test. Related Numpy bug report:

      https://github.com/numpy/numpy/issues/3840
    """

    t = table.Table.read([' a b  c  d',
                          ' 2 c 7.0 0',
                          ' 2 b 5.0 1',
                          ' 2 b 6.0 2',
                          ' 2 a 4.0 3',
                          ' 0 a 0.0 4',
                          ' 1 b 3.0 5',
                          ' 1 a 2.0 6',
                          ' 1 a 1.0 7',
                         ], format='ascii')

    t = table.Table(t, masked=True)

    t2 = table.Table.read([' a b  c  d',
                           ' 2 c 7.0 0',
                           ' 2 b 5.0 1',
                           ' 3 b 6.0 2',
                           ' 2 a 4.0 3',
                           ' 0 a 1.0 4',
                           ' 1 b 3.0 5',
                           ' 1 c 2.0 6',
                           ' 1 a 1.0 7',
                          ], format='ascii')

    assert np.all((t.as_array() == t2) == np.array([0,1,0,1,0,1,0,1], dtype=bool))


# Check that the meta descriptor is working as expected. The MetaBaseTest class
# takes care of defining all the tests, and we simply have to define the class
# and any minimal set of args to pass.

from ...utils.tests.test_metadata import MetaBaseTest


class TestMetaTable(MetaBaseTest):
    test_class = table.Table
    args = ()


def test_unicode_column_names(table_types):
    """
    Test that unicode column names are accepted.  Only do this for
    Python 2 since strings are unicode already in Python 3.
    """
    if six.PY2:
        t = table_types.Table([[1]], names=(six.text_type('a'),))
        assert t.colnames == ['a']
        t[six.text_type('b')] = 0.0
        assert t.colnames == ['a', 'b']


def test_unicode_content():
    # If we don't have unicode literals then return
    if isinstance('', bytes):
        return

    # Define unicode literals
    string_a = 'астрономическая питона'
    string_b = 'миллиарды световых лет'

    a = table.Table(
        [[string_a, 2],
         [string_b, 3]],
        names=('a', 'b'))

    assert string_a in six.text_type(a)
    # This only works because the coding of this file is utf-8, which
    # matches the default encoding of Table.__str__
    assert string_a.encode('utf-8') in bytes(a)


def test_unicode_policy():
    t = table.Table.read([' a b  c  d',
                          ' 2 c 7.0 0',
                          ' 2 b 5.0 1',
                          ' 2 b 6.0 2',
                          ' 2 a 4.0 3',
                          ' 0 a 0.0 4',
                          ' 1 b 3.0 5',
                          ' 1 a 2.0 6',
                          ' 1 a 1.0 7',
                         ], format='ascii')
    assert_follows_unicode_guidelines(t)


def test_unicode_bytestring_conversion(table_types):
    t = table_types.Table([['abc'], ['def'], [1]], dtype=('S', 'U', 'i'))
    assert t['col0'].dtype.kind == 'S'
    assert t['col1'].dtype.kind == 'U'
    assert t['col2'].dtype.kind == 'i'

    t1 = t.copy()
    t1.convert_unicode_to_bytestring()
    assert t1['col0'].dtype.kind == 'S'
    assert t1['col1'].dtype.kind == 'S'
    assert t1['col2'].dtype.kind == 'i'
    assert t1['col0'][0] == 'abc'.encode('ascii')
    assert t1['col1'][0] == 'def'.encode('ascii')
    assert t1['col2'][0] == 1

    t1 = t.copy()
    t1.convert_bytestring_to_unicode()
    assert t1['col0'].dtype.kind == 'U'
    assert t1['col1'].dtype.kind == 'U'
    assert t1['col2'].dtype.kind == 'i'
    assert t1['col0'][0] == six.text_type('abc')
    assert t1['col1'][0] == six.text_type('def')
    assert t1['col2'][0] == 1


def test_table_deletion():
    """
    Regression test for the reference cycle discussed in
    https://github.com/astropy/astropy/issues/2877
    """

    deleted = set()

    # A special table subclass which leaves a record when it is finalized
    class TestTable(table.Table):
        def __del__(self):
            deleted.add(id(self))

    t = TestTable({'a': [1, 2, 3]})
    the_id = id(t)
    assert t['a'].parent_table is t

    del t

    # Cleanup
    gc.collect()

    assert the_id in deleted

def test_nested_iteration():
    """
    Regression test for issue 3358 where nested iteration over a single table fails.
    """
    t = table.Table([[0, 1]], names=['a'])
    out = []
    for r1 in t:
        for r2 in t:
            out.append((r1['a'], r2['a']))
    assert out == [(0, 0), (0, 1), (1, 0), (1, 1)]


def test_table_init_from_degenerate_arrays(table_types):
    t = table_types.Table(np.array([]))
    assert len(t.columns) == 0

    with pytest.raises(ValueError):
        t = table_types.Table(np.array(0))

    t = table_types.Table(np.array([1, 2, 3]))
    assert len(t.columns) == 3


@pytest.mark.skipif('not HAS_PANDAS')
class TestPandas(object):

    def test_simple(self):

        t = table.Table()

        for endian in ['<', '>']:
            for kind in ['f', 'i']:
                for byte in ['2','4','8']:
                    dtype = np.dtype(endian + kind + byte)
                    x = np.array([1,2,3], dtype=dtype)
                    t[endian + kind + byte] = x

        t['u'] = ['a','b','c']
        t['s'] = [b'a', b'b', b'c']

        d = t.to_pandas()

        for column in t.columns:
            if column == 'u':
                assert np.all(t['u'] == np.array(['a','b','c']))
                assert d[column].dtype == np.dtype("O")  # upstream feature of pandas
            elif column == 's':
                assert np.all(t['s'] == np.array([b'a',b'b',b'c']))
                assert d[column].dtype == np.dtype("O")  # upstream feature of pandas
            else:
                # We should be able to compare exact values here
                assert np.all(t[column] == d[column])
                if t[column].dtype.byteorder in ('=', '|'):
                    assert d[column].dtype == t[column].dtype
                else:
                    assert d[column].dtype == t[column].byteswap().newbyteorder().dtype


        # Regression test for astropy/astropy#1156 - the following code gave a
        # ValueError: Big-endian buffer not supported on little-endian
        # compiler. We now automatically swap the endian-ness to native order
        # upon adding the arrays to the data frame.
        d[['<i4','>i4']]
        d[['<f4','>f4']]

        t2 = table.Table.from_pandas(d)

        for column in t.columns:
            if column in ('u', 's'):
                assert np.all(t[column] == t2[column])
            else:
                assert_allclose(t[column], t2[column])
            if t[column].dtype.byteorder in ('=', '|'):
                assert t[column].dtype == t2[column].dtype
            else:
                assert t[column].byteswap().newbyteorder().dtype == t2[column].dtype

    def test_2d(self):

        t = table.Table()
        t['a'] = [1,2,3]
        t['b'] = np.ones((3,2))

        with pytest.raises(ValueError) as exc:
            t.to_pandas()
        assert exc.value.args[0] == "Cannot convert a table with multi-dimensional columns to a pandas DataFrame"

    def test_mixin(self):

        from ...coordinates import SkyCoord

        t = table.Table()
        t['c'] = SkyCoord([1,2,3], [4,5,6], unit='deg')

        with pytest.raises(ValueError) as exc:
            t.to_pandas()
        assert exc.value.args[0] == "Cannot convert a table with mixin columns to a pandas DataFrame"

    def test_masking(self):

        t = table.Table(masked=True)

        t['a'] = [1, 2, 3]
        t['a'].mask = [True, False, True]

        t['b'] = [1., 2., 3.]
        t['b'].mask = [False, False, True]

        t['u'] = ['a','b','c']
        t['u'].mask = [False, True, False]

        t['s'] = [b'a', b'b', b'c']
        t['s'].mask = [False, True, False]

        d = t.to_pandas()

        t2 = table.Table.from_pandas(d)

        for name, column in t.columns.items():
            assert np.all(column.data == t2[name].data)
            assert np.all(column.mask == t2[name].mask)
            # Masked integer type comes back as float.  Nothing we can do about this.
            if column.dtype.kind == 'i':
                assert t2[name].dtype.kind == 'f'
            else:
                if column.dtype.byteorder in ('=', '|'):
                    assert column.dtype == t2[name].dtype
                else:
                    assert column.byteswap().newbyteorder().dtype == t2[name].dtype


@pytest.mark.usefixtures('table_types')
class TestReplaceColumn(SetupData):
    def test_fail_replace_column(self, table_types):
        """Raise exception when trying to replace column via table.columns object"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])

        with pytest.raises(ValueError):
            t.columns['a'] = [1, 2, 3]

        with pytest.raises(ValueError):
            t.replace_column('not there', [1, 2, 3])

    def test_replace_column(self, table_types):
        """Replace existing column with a new column"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        ta = t['a']
        tb = t['b']

        vals = [1.2, 3.4, 5.6]
        for col in (vals,
                    table_types.Column(vals),
                    table_types.Column(vals, name='a'),
                    table_types.Column(vals, name='b')):
            t.replace_column('a', col)
            assert np.all(t['a'] == vals)
            assert t['a'] is not ta  # New a column
            assert t['b'] is tb  # Original b column unchanged
            assert t.colnames == ['a', 'b']
            assert t['a'].meta == {}
            assert t['a'].format is None
