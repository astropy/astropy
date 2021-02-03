# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.tests.test_metadata import MetaBaseTest
import gc
import sys
import copy
from io import StringIO
from collections import OrderedDict
import pickle

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from astropy.io import fits
from astropy.table import (Table, QTable, MaskedColumn, TableReplaceWarning,
                           TableAttribute)
from astropy.tests.helper import assert_follows_unicode_guidelines
from astropy.coordinates import SkyCoord

from astropy.utils.data import get_pkg_data_filename
from astropy import table
from astropy import units as u
from astropy.time import Time, TimeDelta
from .conftest import MaskedTable, MIXIN_COLS

try:
    import pandas  # noqa
except ImportError:
    HAS_PANDAS = False
else:
    HAS_PANDAS = True

try:
    import yaml  # noqa
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


class SetupData:
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
        with pytest.raises(ValueError):
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

        t['aa'] = np.array([1, 2, 3]) * u.m
        assert np.all(t['aa'] == np.array([1, 2, 3]))
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
        t['g'] = np.array([1, 2, 3]) * u.m
        assert np.all(t['g'].data == np.array([1, 2, 3]))
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

    def test_scalar(self, table_types):
        """Test related to #3811 where setting empty tables to scalar values
        should raise an error instead of having an error raised when accessing
        the table."""
        t = table_types.Table()
        with pytest.raises(TypeError, match='Empty table cannot have column set to scalar value'):
            t.add_column(0)

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
        assert type(t['b'][1]) is np.float32

    def test_from_np_array(self, table_types):
        cols = [table_types.Column(name='a', data=np.array([1, 2, 3], dtype=np.int64),
                                   dtype=np.float64),
                table_types.Column(name='b', data=np.array([4, 5, 6], dtype=np.float32))]
        t = table_types.Table(cols)
        assert np.all(t['a'] == np.array([1, 2, 3], dtype=np.float64))
        assert np.all(t['b'] == np.array([4, 5, 6], dtype=np.float32))
        assert type(t['a'][1]) is np.float64
        assert type(t['b'][1]) is np.float32

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

    def test_reverse_mixin(self):
        """Test reverse for a mixin with no item assignment, fix for #9836"""
        sc = SkyCoord([1, 2], [3, 4], unit='deg')
        t = Table([[2, 1], sc], names=['a', 'sc'])
        t.reverse()
        assert np.all(t['a'] == [1, 2])
        assert np.allclose(t['sc'].ra.to_value('deg'), [2, 1])


@pytest.mark.usefixtures('table_types')
class TestRound():

    def test_round_int(self, table_types):
        t = table_types.Table([['a', 'b', 'c'],
                               [1.11, 2.3, 3.0],
                               [1.123456, 2.9876, 3.901]])
        t.round()
        assert np.all(t['col0'] == ['a', 'b', 'c'])
        assert np.all(t['col1'] == [1., 2., 3.])
        assert np.all(t['col2'] == [1., 3., 4.])

    def test_round_dict(self, table_types):
        t = table_types.Table([['a', 'b', 'c'],
                               [1.5, 2.5, 3.0111],
                               [1.123456, 2.9876, 3.901]])

        t.round({'col1': 0, 'col2': 3})
        assert np.all(t['col0'] == ['a', 'b', 'c'])
        assert np.all(t['col1'] == [2.0, 2.0, 3.0])
        assert np.all(t['col2'] == [1.123, 2.988, 3.901])

    def test_round_invalid(self, table_types):
        t = table_types.Table([[1, 2, 3]])
        with pytest.raises(ValueError, match="'decimals' argument must be an int or a dict"):
            t.round(0.5)

    def test_round_kind(self, table_types):
        for typecode in 'bBhHiIlLqQpPefdgFDG':  # AllInteger, AllFloat
            arr = np.array([4, 16], dtype=typecode)
            t = Table([arr])
            col0 = t['col0']
            t.round(decimals=-1)  # Round to nearest 10
            assert np.all(t['col0'] == [0, 20])
            assert t['col0'] is col0


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

    def test_itercols(self, table_types):
        names = ['a', 'b', 'c']
        t = table_types.Table([[1], [2], [3]], names=names)
        for name, col in zip(names, t.itercols()):
            assert name == col.name
            assert isinstance(col, table_types.Column)


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
        assert t.colnames == ['a', 'b']

    def test_7(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.add_column(self.b, t.index_column('a'))
        assert t.colnames == ['b', 'a']

    def test_8(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.add_column(self.b, t.index_column('a') + 1)
        assert t.colnames == ['a', 'b']

    def test_9(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a)
        t.add_column(self.b, t.index_column('a') + 1)
        t.add_column(self.c, t.index_column('b'))
        assert t.colnames == ['a', 'c', 'b']

    def test_10(self, table_types):
        self._setup(table_types)
        t = table_types.Table()
        t.add_column(self.a)
        ia = t.index_column('a')
        t.add_column(self.b, ia + 1)
        t.add_column(self.c, ia)
        assert t.colnames == ['c', 'a', 'b']


@pytest.mark.usefixtures('table_types')
class TestAddName(SetupData):

    def test_override_name(self, table_types):
        self._setup(table_types)
        t = table_types.Table()

        # Check that we can override the name of the input column in the Table
        t.add_column(self.a, name='b')
        t.add_column(self.b, name='a')
        assert t.colnames == ['b', 'a']
        # Check that we did not change the name of the input column
        assert self.a.info.name == 'a'
        assert self.b.info.name == 'b'

        # Now test with an input column from another table
        t2 = table_types.Table()
        t2.add_column(t['a'], name='c')
        assert t2.colnames == ['c']
        # Check that we did not change the name of the input column
        assert t.colnames == ['b', 'a']

        # Check that we can give a name if none was present
        col = table_types.Column([1, 2, 3])
        t.add_column(col, name='c')
        assert t.colnames == ['b', 'a', 'c']

    def test_default_name(self, table_types):
        t = table_types.Table()
        col = table_types.Column([1, 2, 3])
        t.add_column(col)
        assert t.colnames == ['col0']


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

    def test_add_columns6(self, table_types):
        """Check that we can override column names."""
        self._setup(table_types)
        t = table_types.Table()
        t.add_columns([self.a, self.b, self.c], names=['b', 'c', 'a'])
        assert t.colnames == ['b', 'c', 'a']

    def test_add_columns7(self, table_types):
        """Check that default names are used when appropriate."""
        t = table_types.Table()
        col0 = table_types.Column([1, 2, 3])
        col1 = table_types.Column([4, 5, 3])
        t.add_columns([col0, col1])
        assert t.colnames == ['col0', 'col1']

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

        # Check that rename_duplicate=True is ok if there are no duplicates
        t.add_column(table_types.Column(name='q', data=[0, 1, 2]),
                     rename_duplicate=True)
        assert t.colnames == ['a', 'a_1', 'b', 'c', 'a_2', 'a_3', 'q']

    def test_add_duplicate_columns(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b, self.c])
        with pytest.raises(ValueError):
            t.add_columns([table_types.Column(name='a', data=[0, 1, 2]),
                           table_types.Column(name='b', data=[0, 1, 2])])
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
        assert t['b'][0] == ''
        assert t['c'][0] == 0
        t.add_row()
        assert np.all(t['a'][1] == [0, 0])
        assert t['b'][1] == ''
        assert t['c'][1] == 0

    def test_add_stuff_to_empty_table(self, table_types):
        self._setup(table_types)
        t = table_types.Table(names=('a', 'b', 'obj'), dtype=('(2,)i', 'S8', 'O'))
        t.add_row([[1, 2], 'hello', 'world'])
        assert np.all(t['a'][0] == [1, 2])
        assert t['b'][0] == 'hello'
        assert t['obj'][0] == 'world'
        # Make sure it is not repeating last row but instead
        # adding zeros (as documented)
        t.add_row()
        assert np.all(t['a'][1] == [0, 0])
        assert t['b'][1] == ''
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

    def test_add_qtable_row_multidimensional(self):
        q = [[1, 2], [3, 4]] * u.m
        qt = table.QTable([q])
        qt.add_row(([5, 6] * u.km,))
        assert np.all(qt['col0'] == [[1, 2], [3, 4], [5000, 6000]] * u.m)

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
        assert np.all(t['c'] == np.array(['7', '8', '9', '10']))

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
            t.add_row([1, 2, 3, 4])
        except ValueError:
            pass
        assert len(t) == 3
        assert np.all(t.as_array() == t_copy.as_array())
        # Wrong data type
        try:
            t.add_row(['one', 2, 3])
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
        assert self.t.colnames == []
        assert self.t.as_array().size == 0
        # Regression test for gh-8640
        assert not self.t
        assert isinstance(self.t == None, np.ndarray)  # noqa
        assert (self.t == None).size == 0  # noqa

    def test_2(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.remove_columns('a')
        assert self.t.colnames == ['b']
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
        assert self.t.dtype == np.dtype([('a', 'int'),
                                         ('b', 'int')])

    def test_delitem_row(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        del self.t[1]
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['a'] == np.array([1, 3]))

    @pytest.mark.parametrize("idx", [[0, 2], np.array([0, 2])])
    def test_delitem_row_list(self, table_types, idx):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        del self.t[idx]
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['c'] == np.array([8]))

    def test_delitem_row_slice(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.add_column(self.c)
        del self.t[0:2]
        assert self.t.colnames == ['a', 'b', 'c']
        assert np.all(self.t['c'] == np.array([9]))

    def test_delitem_row_fail(self, table_types):
        self._setup(table_types)
        with pytest.raises(IndexError):
            del self.t[4]

    def test_delitem_row_float(self, table_types):
        self._setup(table_types)
        with pytest.raises(IndexError):
            del self.t[1.]

    def test_delitem1(self, table_types):
        self._setup(table_types)
        del self.t['a']
        assert self.t.colnames == []
        assert self.t.as_array().size == 0
        # Regression test for gh-8640
        assert not self.t
        assert isinstance(self.t == None, np.ndarray)  # noqa
        assert (self.t == None).size == 0  # noqa

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
        assert t.colnames == []
        assert t.as_array().size == 0
        # Regression test for gh-8640
        assert not t
        assert isinstance(t == None, np.ndarray)  # noqa
        assert (t == None).size == 0  # noqa

    def test_2(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.keep_columns('b')
        assert t.colnames == ['b']
        assert t.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([4, 5, 6]))


@pytest.mark.usefixtures('table_types')
class TestRename(SetupData):

    def test_1(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.rename_column('a', 'b')
        assert t.colnames == ['b']
        assert t.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([1, 2, 3]))

    def test_2(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.rename_column('a', 'c')
        t.rename_column('b', 'a')
        assert t.colnames == ['c', 'a']
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
        assert t.colnames == ['c', 'a']
        assert t.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))

    def test_rename_columns(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b, self.c])
        t.rename_columns(('a', 'b', 'c'), ('aa', 'bb', 'cc'))
        assert t.colnames == ['aa', 'bb', 'cc']
        t.rename_columns(['bb', 'cc'], ['b', 'c'])
        assert t.colnames == ['aa', 'b', 'c']
        with pytest.raises(TypeError):
            t.rename_columns(('aa'), ['a'])
        with pytest.raises(ValueError):
            t.rename_columns(['a'], ['b', 'c'])


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

    @pytest.mark.parametrize('create_index', [False, True])
    def test_single_reverse(self, table_types, create_index):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[2, 1, 3]))
        t.add_column(table_types.Column(name='b', data=[6, 5, 4]))
        t.add_column(table_types.Column(name='c', data=[(1, 2), (3, 4), (4, 5)]))
        assert np.all(t['a'] == np.array([2, 1, 3]))
        assert np.all(t['b'] == np.array([6, 5, 4]))
        t.sort('a', reverse=True)
        assert np.all(t['a'] == np.array([3, 2, 1]))
        assert np.all(t['b'] == np.array([4, 6, 5]))
        assert np.all(t['c'] == np.array([[4, 5],
                                          [1, 2],
                                          [3, 4]]))
        t.sort('b', reverse=True)
        assert np.all(t['a'] == np.array([2, 1, 3]))
        assert np.all(t['b'] == np.array([6, 5, 4]))
        assert np.all(t['c'] == np.array([[1, 2],
                                          [3, 4],
                                          [4, 5]]))

    def test_single_big(self, table_types):
        """Sort a big-ish table with a non-trivial sort order"""
        x = np.arange(10000)
        y = np.sin(x)
        t = table_types.Table([x, y], names=('x', 'y'))
        t.sort('y')
        idx = np.argsort(y)
        assert np.all(t['x'] == x[idx])
        assert np.all(t['y'] == y[idx])

    @pytest.mark.parametrize('reverse', [True, False])
    def test_empty_reverse(self, table_types, reverse):
        t = table_types.Table([[], []], dtype=['f4', 'U1'])
        t.sort('col1', reverse=reverse)

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
        t.sort(('a', 'b'))
        assert np.all(t['a'] == np.array([1, 1, 2, 2, 3, 3]))
        assert np.all(t['b'] == np.array([4, 5, 3, 6, 4, 5]))

    def test_multiple_reverse(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[2, 1, 3, 2, 3, 1]))
        t.add_column(table_types.Column(name='b', data=[6, 5, 4, 3, 5, 4]))
        assert np.all(t['a'] == np.array([2, 1, 3, 2, 3, 1]))
        assert np.all(t['b'] == np.array([6, 5, 4, 3, 5, 4]))
        t.sort(['a', 'b'], reverse=True)
        assert np.all(t['a'] == np.array([3, 3, 2, 2, 1, 1]))
        assert np.all(t['b'] == np.array([5, 4, 6, 3, 5, 4]))
        t.sort(['b', 'a'], reverse=True)
        assert np.all(t['a'] == np.array([2, 3, 1, 3, 1, 2]))
        assert np.all(t['b'] == np.array([6, 5, 5, 4, 4, 3]))
        t.sort(('a', 'b'), reverse=True)
        assert np.all(t['a'] == np.array([3, 3, 2, 2, 1, 1]))
        assert np.all(t['b'] == np.array([5, 4, 6, 3, 5, 4]))

    def test_multiple_with_bytes(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='firstname', data=[b"Max", b"Jo", b"John"]))
        t.add_column(table_types.Column(name='name', data=[b"Miller", b"Miller", b"Jackson"]))
        t.add_column(table_types.Column(name='tel', data=[12, 15, 19]))
        t.sort(['name', 'firstname'])
        assert np.all([t['firstname'] == np.array([b"John", b"Jo", b"Max"])])
        assert np.all([t['name'] == np.array([b"Jackson", b"Miller", b"Miller"])])
        assert np.all([t['tel'] == np.array([19, 15, 12])])

    def test_multiple_with_unicode(self, table_types):
        # Before Numpy 1.6.2, sorting with multiple column names
        # failed when a unicode column was present.
        t = table_types.Table()
        t.add_column(table_types.Column(
            name='firstname',
            data=[str(x) for x in ["Max", "Jo", "John"]]))
        t.add_column(table_types.Column(
            name='name',
            data=[str(x) for x in ["Miller", "Miller", "Jackson"]]))
        t.add_column(table_types.Column(name='tel', data=[12, 15, 19]))
        t.sort(['name', 'firstname'])
        assert np.all([t['firstname'] == np.array(
            [str(x) for x in ["John", "Jo", "Max"]])])
        assert np.all([t['name'] == np.array(
            [str(x) for x in ["Jackson", "Miller", "Miller"]])])
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

    @pytest.mark.parametrize('add_index', [False, True])
    def test_argsort_reverse(self, table_types, add_index):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[2, 1, 3, 2, 3, 1]))
        t.add_column(table_types.Column(name='b', data=[6, 5, 4, 3, 5, 4]))
        if add_index:
            t.add_index('a')
        assert np.all(t.argsort(reverse=True) == np.array([4, 2, 0, 3, 1, 5]))
        i0 = t.argsort('a', reverse=True)
        i1 = np.array([4, 2, 3, 0, 5, 1])
        assert np.all(t['a'][i0] == t['a'][i1])
        i0 = t.argsort(['a', 'b'], reverse=True)
        i1 = np.array([4, 2, 0, 3, 1, 5])
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
            data=[str(x) for x in ["Max", "Jo", "John"]]))
        t.add_column(table_types.Column(
            name='name',
            data=[str(x) for x in ["Miller", "Miller", "Jackson"]]))
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
                      (4, 5)], dtype=[('a', 'i4'), ('b', 'i4')])
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
        assert np_data is not d.as_array()
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d.as_array())
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[('c', 'i8'), ('d', 'i8')])

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

        with fits.open(filename, character_as_bytes=True) as hdul:
            data = hdul[1].data
            for colname in data.columns.names:
                assert np.all(data[colname] == arr[colname])

            arr2 = t.as_array(keep_byteorder=True)
            for colname in data.columns.names:
                assert (data[colname].dtype.byteorder
                        == arr2[colname].dtype.byteorder)


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


def test_values_equal_part1():

    col1 = [1, 2]
    col2 = [1.0, 2.0]
    col3 = ['a', 'b']
    t1 = table.Table([col1, col2, col3], names=['a', 'b', 'c'])
    t2 = table.Table([col1, col2], names=['a', 'b'])
    t3 = table.table_helpers.simple_table()
    tm = t1.copy()
    tm['time'] = Time([1, 2], format='cxcsec')
    tm1 = tm.copy()
    tm1['time'][0] = np.ma.masked

    tq = table.table_helpers.simple_table()
    tq['quantity'] = [1., 2., 3.] * u.m

    tsk = table.table_helpers.simple_table()
    tsk['sk'] = SkyCoord(1, 2, unit='deg')
    eqsk = tsk.values_equal(tsk)
    for col in eqsk.itercols():
        assert np.all(col)

    with pytest.raises(ValueError, match='cannot compare tables with different column names'):
        t2.values_equal(t1)

    with pytest.raises(ValueError, match='unable to compare column a'):
        # Shape mismatch
        t3.values_equal(t1)

    with pytest.raises(ValueError, match='unable to compare column c'):
        # Type mismatch in column c causes FutureWarning
        t1.values_equal(2)

    with pytest.raises(ValueError, match='unable to compare column c'):
        t1.values_equal([1, 2])

    eq = t2.values_equal(t2)
    for col in eq.colnames:
        assert np.all(eq[col] == [True, True])

    eq1 = tm1.values_equal(tm)
    for col in eq1.colnames:
        assert np.all(eq1[col] == [True, True])

    eq2 = tq.values_equal(tq)
    for col in eq2.colnames:
        assert np.all(eq2[col] == [True, True, True])

    eq3 = t2.values_equal(2)
    for col in eq3.colnames:
        assert np.all(eq3[col] == [False, True])

    eq4 = t2.values_equal([1, 2])
    for col in eq4.colnames:
        assert np.all(eq4[col] == [True, True])

    # Compare table to its first row
    t = table.Table(rows=[(1, 'a'),
                          (1, 'b')])
    eq = t.values_equal(t[0])
    assert np.all(eq['col0'] == [True, True])
    assert np.all(eq['col1'] == [True, False])


def test_rows_equal():

    t = table.Table.read([' a b  c  d',
                          ' 2 c 7.0 0',
                          ' 2 b 5.0 1',
                          ' 2 b 6.0 2',
                          ' 2 a 4.0 3',
                          ' 0 a 0.0 4',
                          ' 1 b 3.0 5',
                          ' 1 a 2.0 6',
                          ' 1 a 1.0 7'],
                         format='ascii')

    # All rows are equal
    assert np.all(t == t)

    # Assert no rows are different
    assert not np.any(t != t)

    # Check equality result for a given row
    assert np.all((t == t[3]) == np.array([0, 0, 0, 1, 0, 0, 0, 0], dtype=bool))

    # Check inequality result for a given row
    assert np.all((t != t[3]) == np.array([1, 1, 1, 0, 1, 1, 1, 1], dtype=bool))

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
    assert np.all((t == t2) == np.array([1, 1, 0, 1, 0, 1, 0, 1], dtype=bool))
    assert np.all((t != t2) == np.array([0, 0, 1, 0, 1, 0, 1, 0], dtype=bool))

    # Check that comparing to a structured array works
    assert np.all((t == t2.as_array()) == np.array([1, 1, 0, 1, 0, 1, 0, 1], dtype=bool))
    assert np.all((t.as_array() == t2) == np.array([1, 1, 0, 1, 0, 1, 0, 1], dtype=bool))


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
    assert np.all(t == t)

    # Assert no rows are different
    assert not np.any(t != t)

    # Check equality result for a given row
    assert np.all((t == t[3]) == np.array([0, 0, 0, 1, 0, 0, 0, 0], dtype=bool))

    # Check inequality result for a given row
    assert np.all((t != t[3]) == np.array([1, 1, 1, 0, 1, 1, 1, 1], dtype=bool))

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
    assert np.all((t == t2) == np.array([1, 1, 0, 1, 0, 1, 0, 1], dtype=bool))
    assert np.all((t != t2) == np.array([0, 0, 1, 0, 1, 0, 1, 0], dtype=bool))

    # Check that masking a value causes the row to differ
    t.mask['a'][0] = True
    assert np.all((t == t2) == np.array([0, 1, 0, 1, 0, 1, 0, 1], dtype=bool))
    assert np.all((t != t2) == np.array([1, 0, 1, 0, 1, 0, 1, 0], dtype=bool))

    # Check that comparing to a structured array works
    assert np.all((t == t2.as_array()) == np.array([0, 1, 0, 1, 0, 1, 0, 1], dtype=bool))


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

    assert np.all((t.as_array() == t2) == np.array([0, 1, 0, 1, 0, 1, 0, 1], dtype=bool))


# Check that the meta descriptor is working as expected. The MetaBaseTest class
# takes care of defining all the tests, and we simply have to define the class
# and any minimal set of args to pass.


class TestMetaTable(MetaBaseTest):
    test_class = table.Table
    args = ()


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

    assert string_a in str(a)
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


@pytest.mark.parametrize('uni', ['питона', 'ascii'])
def test_unicode_bytestring_conversion(table_types, uni):
    """
    Test converting columns to all unicode or all bytestring.  Thi
    makes two columns, one which is unicode (str in Py3) and one which
    is bytes (UTF-8 encoded).  There are two code paths in the conversions,
    a faster one where the data are actually ASCII and a slower one where
    UTF-8 conversion is required.  This tests both via the ``uni`` param.
    """
    byt = uni.encode('utf-8')
    t = table_types.Table([[byt], [uni], [1]], dtype=('S', 'U', 'i'))
    assert t['col0'].dtype.kind == 'S'
    assert t['col1'].dtype.kind == 'U'
    assert t['col2'].dtype.kind == 'i'
    t['col0'].description = 'col0'
    t['col1'].description = 'col1'
    t['col0'].meta['val'] = 'val0'
    t['col1'].meta['val'] = 'val1'

    # Unicode to bytestring
    t1 = t.copy()
    t1.convert_unicode_to_bytestring()
    assert t1['col0'].dtype.kind == 'S'
    assert t1['col1'].dtype.kind == 'S'
    assert t1['col2'].dtype.kind == 'i'

    # Meta made it through
    assert t1['col0'].description == 'col0'
    assert t1['col1'].description == 'col1'
    assert t1['col0'].meta['val'] == 'val0'
    assert t1['col1'].meta['val'] == 'val1'

    # Need to de-fang the automatic unicode sandwiching of Table
    assert np.array(t1['col0'])[0] == byt
    assert np.array(t1['col1'])[0] == byt
    assert np.array(t1['col2'])[0] == 1

    # Bytestring to unicode
    t1 = t.copy()
    t1.convert_bytestring_to_unicode()
    assert t1['col0'].dtype.kind == 'U'
    assert t1['col1'].dtype.kind == 'U'
    assert t1['col2'].dtype.kind == 'i'

    # Meta made it through
    assert t1['col0'].description == 'col0'
    assert t1['col1'].description == 'col1'
    assert t1['col0'].meta['val'] == 'val0'
    assert t1['col1'].meta['val'] == 'val1'

    # No need to de-fang the automatic unicode sandwiching of Table here, but
    # do just for consistency to prove things are working.
    assert np.array(t1['col0'])[0] == uni
    assert np.array(t1['col1'])[0] == uni
    assert np.array(t1['col2'])[0] == 1


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
class TestPandas:

    def test_simple(self):

        t = table.Table()

        for endian in ['<', '>', '=']:
            for kind in ['f', 'i']:
                for byte in ['2', '4', '8']:
                    dtype = np.dtype(endian + kind + byte)
                    x = np.array([1, 2, 3], dtype=dtype)
                    t[endian + kind + byte] = x.newbyteorder(endian)

        t['u'] = ['a', 'b', 'c']
        t['s'] = ['a', 'b', 'c']

        d = t.to_pandas()

        for column in t.columns:
            if column == 'u':
                assert np.all(t['u'] == np.array(['a', 'b', 'c']))
                assert d[column].dtype == np.dtype("O")  # upstream feature of pandas
            elif column == 's':
                assert np.all(t['s'] == np.array(['a', 'b', 'c']))
                assert d[column].dtype == np.dtype("O")  # upstream feature of pandas
            else:
                # We should be able to compare exact values here
                assert np.all(t[column] == d[column])
                if t[column].dtype.isnative:
                    assert d[column].dtype == t[column].dtype
                else:
                    assert d[column].dtype == t[column].byteswap().newbyteorder().dtype

        # Regression test for astropy/astropy#1156 - the following code gave a
        # ValueError: Big-endian buffer not supported on little-endian
        # compiler. We now automatically swap the endian-ness to native order
        # upon adding the arrays to the data frame.
        # Explicitly testing little/big/native endian separately -
        # regression for a case in astropy/astropy#11286 not caught by #3729.
        d[['<i4', '>i4']]
        d[['<f4', '>f4']]

        t2 = table.Table.from_pandas(d)

        for column in t.columns:
            if column in ('u', 's'):
                assert np.all(t[column] == t2[column])
            else:
                assert_allclose(t[column], t2[column])
            if t[column].dtype.isnative:
                assert t[column].dtype == t2[column].dtype
            else:
                assert t[column].byteswap().newbyteorder().dtype == t2[column].dtype

    @pytest.mark.parametrize('unsigned', ['u', ''])
    @pytest.mark.parametrize('bits', [8, 16, 32, 64])
    def test_nullable_int(self, unsigned, bits):
        np_dtype = f'{unsigned}int{bits}'
        c = MaskedColumn([1, 2], mask=[False, True], dtype=np_dtype)
        t = Table([c])
        df = t.to_pandas()
        pd_dtype = np_dtype.replace('i', 'I').replace('u', 'U')
        assert str(df['col0'].dtype) == pd_dtype
        t2 = Table.from_pandas(df)
        assert str(t2['col0'].dtype) == np_dtype
        assert np.all(t2['col0'].mask == [False, True])
        assert np.all(t2['col0'] == c)

    def test_2d(self):

        t = table.Table()
        t['a'] = [1, 2, 3]
        t['b'] = np.ones((3, 2))

        with pytest.raises(ValueError,
                           match='Cannot convert a table with multidimensional columns'):
            t.to_pandas()

    def test_mixin_pandas(self):
        t = table.QTable()
        for name in sorted(MIXIN_COLS):
            if name != 'ndarray':
                t[name] = MIXIN_COLS[name]

        t['dt'] = TimeDelta([0, 2, 4, 6], format='sec')

        tp = t.to_pandas()
        t2 = table.Table.from_pandas(tp)

        assert np.allclose(t2['quantity'], [0, 1, 2, 3])
        assert np.allclose(t2['longitude'], [0., 1., 5., 6.])
        assert np.allclose(t2['latitude'], [5., 6., 10., 11.])
        assert np.allclose(t2['skycoord.ra'], [0, 1, 2, 3])
        assert np.allclose(t2['skycoord.dec'], [0, 1, 2, 3])
        assert np.allclose(t2['arraywrap'], [0, 1, 2, 3])
        assert np.allclose(t2['arrayswap'], [0, 1, 2, 3])
        assert np.allclose(t2['earthlocation.y'], [0, 110708, 547501, 654527], rtol=0, atol=1)

        # For pandas, Time, TimeDelta are the mixins that round-trip the class
        assert isinstance(t2['time'], Time)
        assert np.allclose(t2['time'].jyear, [2000, 2001, 2002, 2003])
        assert np.all(t2['time'].isot == ['2000-01-01T12:00:00.000',
                                          '2000-12-31T18:00:00.000',
                                          '2002-01-01T00:00:00.000',
                                          '2003-01-01T06:00:00.000'])
        assert t2['time'].format == 'isot'

        # TimeDelta
        assert isinstance(t2['dt'], TimeDelta)
        assert np.allclose(t2['dt'].value, [0, 2, 4, 6])
        assert t2['dt'].format == 'sec'

    def test_to_pandas_index(self):
        import pandas as pd
        row_index = pd.RangeIndex(0, 2, 1)
        tm_index = pd.DatetimeIndex(['1998-01-01', '2002-01-01'],
                                    dtype='datetime64[ns]',
                                    name='tm', freq=None)

        tm = Time([1998, 2002], format='jyear')
        x = [1, 2]
        t = table.QTable([tm, x], names=['tm', 'x'])
        tp = t.to_pandas()
        assert np.all(tp.index == row_index)

        tp = t.to_pandas(index='tm')
        assert np.all(tp.index == tm_index)

        t.add_index('tm')
        tp = t.to_pandas()
        assert np.all(tp.index == tm_index)
        # Make sure writing to pandas didn't hack the original table
        assert t['tm'].info.indices

        tp = t.to_pandas(index=True)
        assert np.all(tp.index == tm_index)

        tp = t.to_pandas(index=False)
        assert np.all(tp.index == row_index)

        with pytest.raises(ValueError) as err:
            t.to_pandas(index='not a column')
        assert 'index must be None, False' in str(err.value)

    def test_mixin_pandas_masked(self):
        tm = Time([1, 2, 3], format='cxcsec')
        dt = TimeDelta([1, 2, 3], format='sec')
        tm[1] = np.ma.masked
        dt[1] = np.ma.masked
        t = table.QTable([tm, dt], names=['tm', 'dt'])

        tp = t.to_pandas()
        assert np.all(tp['tm'].isnull() == [False, True, False])
        assert np.all(tp['dt'].isnull() == [False, True, False])

        t2 = table.Table.from_pandas(tp)

        assert np.all(t2['tm'].mask == tm.mask)
        assert np.ma.allclose(t2['tm'].jd, tm.jd, rtol=1e-14, atol=1e-14)

        assert np.all(t2['dt'].mask == dt.mask)
        assert np.ma.allclose(t2['dt'].jd, dt.jd, rtol=1e-14, atol=1e-14)

    def test_from_pandas_index(self):
        tm = Time([1998, 2002], format='jyear')
        x = [1, 2]
        t = table.Table([tm, x], names=['tm', 'x'])
        tp = t.to_pandas(index='tm')

        t2 = table.Table.from_pandas(tp)
        assert t2.colnames == ['x']

        t2 = table.Table.from_pandas(tp, index=True)
        assert t2.colnames == ['tm', 'x']
        assert np.allclose(t2['tm'].jyear, tm.jyear)

    @pytest.mark.parametrize('use_nullable_int', [True, False])
    def test_masking(self, use_nullable_int):

        t = table.Table(masked=True)

        t['a'] = [1, 2, 3]
        t['a'].mask = [True, False, True]

        t['b'] = [1., 2., 3.]
        t['b'].mask = [False, False, True]

        t['u'] = ['a', 'b', 'c']
        t['u'].mask = [False, True, False]

        t['s'] = ['a', 'b', 'c']
        t['s'].mask = [False, True, False]

        # https://github.com/astropy/astropy/issues/7741
        t['Source'] = [2584290278794471936, 2584290038276303744,
                       2584288728310999296]
        t['Source'].mask = [False, False, False]

        if use_nullable_int:  # Default
            # No warning with the default use_nullable_int=True
            d = t.to_pandas(use_nullable_int=use_nullable_int)
        else:
            with pytest.warns(TableReplaceWarning,
                              match=r"converted column 'a' from int(32|64) to float64"):
                d = t.to_pandas(use_nullable_int=use_nullable_int)

        t2 = table.Table.from_pandas(d)

        for name, column in t.columns.items():
            assert np.all(column.data == t2[name].data)
            if hasattr(t2[name], 'mask'):
                assert np.all(column.mask == t2[name].mask)

            if column.dtype.kind == 'i':
                if np.any(column.mask) and not use_nullable_int:
                    assert t2[name].dtype.kind == 'f'
                else:
                    assert t2[name].dtype.kind == 'i'

                assert_array_equal(column.data,
                                   t2[name].data.astype(column.dtype))
            else:
                if column.dtype.byteorder in ('=', '|'):
                    assert column.dtype == t2[name].dtype
                else:
                    assert column.byteswap().newbyteorder().dtype == t2[name].dtype

    def test_units(self):
        import pandas as pd
        import astropy.units as u

        df = pd.DataFrame({'x': [1, 2, 3], 't': [1.3, 1.2, 1.8]})
        t = table.Table.from_pandas(df, units={'x': u.m, 't': u.s})

        assert t['x'].unit == u.m
        assert t['t'].unit == u.s

        # test error if not a mapping
        with pytest.raises(TypeError):
            table.Table.from_pandas(df, units=[u.m, u.s])

        # test warning is raised if additional columns in units dict
        with pytest.warns(UserWarning) as record:
            table.Table.from_pandas(df, units={'x': u.m, 't': u.s, 'y': u.m})
        assert len(record) == 1
        assert "{'y'}" in record[0].message.args[0]


@pytest.mark.usefixtures('table_types')
class TestReplaceColumn(SetupData):
    def test_fail_replace_column(self, table_types):
        """Raise exception when trying to replace column via table.columns object"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])

        with pytest.raises(ValueError,
                           match=r"Cannot replace column 'a'.  Use "
                           "Table.replace_column.. instead."):
            t.columns['a'] = [1, 2, 3]

        with pytest.raises(ValueError, match=r"column name not there is not in the table"):
            t.replace_column('not there', [1, 2, 3])

        with pytest.raises(ValueError, match=r"length of new column must match table length"):
            t.replace_column('a', [1, 2])

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

        # Special case: replacing the only column can resize table
        del t['b']
        assert len(t) == 3
        t['a'] = [1, 2]
        assert len(t) == 2

    def test_replace_index_column(self, table_types):
        """Replace index column and generate expected exception"""
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.add_index('a')

        with pytest.raises(ValueError) as err:
            t.replace_column('a', [1, 2, 3])
        assert err.value.args[0] == 'cannot replace a table index column'

    def test_replace_column_no_copy(self):
        t = Table([[1, 2], [3, 4]], names=['a', 'b'])
        a = np.array([1.5, 2.5])
        t.replace_column('a', a, copy=False)
        assert t['a'][0] == a[0]
        t['a'][0] = 10
        assert t['a'][0] == a[0]

    def test_replace_with_masked_col_with_units_in_qtable(self):
        """This is a small regression from #8902"""
        t = QTable([[1, 2], [3, 4]], names=['a', 'b'])
        t['a'] = MaskedColumn([5, 6], unit='m')
        assert isinstance(t['a'], u.Quantity)


class Test__Astropy_Table__():
    """
    Test initializing a Table subclass from a table-like object that
    implements the __astropy_table__ interface method.
    """

    class SimpleTable:
        def __init__(self):
            self.columns = [[1, 2, 3],
                            [4, 5, 6],
                            [7, 8, 9] * u.m]
            self.names = ['a', 'b', 'c']
            self.meta = OrderedDict([('a', 1), ('b', 2)])

        def __astropy_table__(self, cls, copy, **kwargs):
            a, b, c = self.columns
            c.info.name = 'c'
            cols = [table.Column(a, name='a'),
                    table.MaskedColumn(b, name='b'),
                    c]
            names = [col.info.name for col in cols]
            return cls(cols, names=names, copy=copy, meta=kwargs or self.meta)

    def test_simple_1(self):
        """Make a SimpleTable and convert to Table, QTable with copy=False, True"""
        for table_cls in (table.Table, table.QTable):
            col_c_class = u.Quantity if table_cls is table.QTable else table.Column
            for cpy in (False, True):
                st = self.SimpleTable()
                # Test putting in a non-native kwarg `extra_meta` to Table initializer
                t = table_cls(st, copy=cpy, extra_meta='extra!')
                assert t.colnames == ['a', 'b', 'c']
                assert t.meta == {'extra_meta': 'extra!'}
                assert np.all(t['a'] == st.columns[0])
                assert np.all(t['b'] == st.columns[1])
                vals = t['c'].value if table_cls is table.QTable else t['c']
                assert np.all(st.columns[2].value == vals)

                assert isinstance(t['a'], table.Column)
                assert isinstance(t['b'], table.MaskedColumn)
                assert isinstance(t['c'], col_c_class)
                assert t['c'].unit is u.m
                assert type(t) is table_cls

                # Copy being respected?
                t['a'][0] = 10
                assert st.columns[0][0] == 1 if cpy else 10

    def test_simple_2(self):
        """Test converting a SimpleTable and changing column names and types"""
        st = self.SimpleTable()
        dtypes = [np.int32, np.float32, np.float16]
        names = ['a', 'b', 'c']
        meta = OrderedDict([('c', 3)])
        t = table.Table(st, dtype=dtypes, names=names, meta=meta)
        assert t.colnames == names
        assert all(col.dtype.type is dtype
                   for col, dtype in zip(t.columns.values(), dtypes))

        # The supplied meta is overrides the existing meta.  Changed in astropy 3.2.
        assert t.meta != st.meta
        assert t.meta == meta

    def test_kwargs_exception(self):
        """If extra kwargs provided but without initializing with a table-like
        object, exception is raised"""
        with pytest.raises(TypeError) as err:
            table.Table([[1]], extra_meta='extra!')
        assert '__init__() got unexpected keyword argument' in str(err.value)


def test_table_meta_copy():
    """
    Test no copy vs light (key) copy vs deep copy of table meta for different
    situations.  #8404.
    """
    t = table.Table([[1]])
    meta = {1: [1, 2]}

    # Assigning meta directly implies using direct object reference
    t.meta = meta
    assert t.meta is meta

    # Table slice implies key copy, so values are unchanged
    t2 = t[:]
    assert t2.meta is not t.meta  # NOT the same OrderedDict object but equal
    assert t2.meta == t.meta
    assert t2.meta[1] is t.meta[1]  # Value IS the list same object

    # Table init with copy=False implies key copy
    t2 = table.Table(t, copy=False)
    assert t2.meta is not t.meta  # NOT the same OrderedDict object but equal
    assert t2.meta == t.meta
    assert t2.meta[1] is t.meta[1]  # Value IS the same list object

    # Table init with copy=True implies deep copy
    t2 = table.Table(t, copy=True)
    assert t2.meta is not t.meta  # NOT the same OrderedDict object but equal
    assert t2.meta == t.meta
    assert t2.meta[1] is not t.meta[1]  # Value is NOT the same list object


def test_table_meta_copy_with_meta_arg():
    """
    Test no copy vs light (key) copy vs deep copy of table meta when meta is
    supplied as a table init argument.  #8404.
    """
    meta = {1: [1, 2]}
    meta2 = {2: [3, 4]}
    t = table.Table([[1]], meta=meta, copy=False)
    assert t.meta is meta

    t = table.Table([[1]], meta=meta)  # default copy=True
    assert t.meta is not meta
    assert t.meta == meta

    # Test initializing from existing table with meta with copy=False
    t2 = table.Table(t, meta=meta2, copy=False)
    assert t2.meta is meta2
    assert t2.meta != t.meta  # Change behavior in #8404

    # Test initializing from existing table with meta with default copy=True
    t2 = table.Table(t, meta=meta2)
    assert t2.meta is not meta2
    assert t2.meta != t.meta  # Change behavior in #8404

    # Table init with copy=True and empty dict meta gets that empty dict
    t2 = table.Table(t, copy=True, meta={})
    assert t2.meta == {}

    # Table init with copy=True and kwarg meta=None gets the original table dict.
    # This is a somewhat ambiguous case because it could be interpreted as the
    # user wanting NO meta set on the output.  This could be implemented by inspecting
    # call args.
    t2 = table.Table(t, copy=True, meta=None)
    assert t2.meta == t.meta

    # Test initializing empty table with meta with copy=False
    t = table.Table(meta=meta, copy=False)
    assert t.meta is meta
    assert t.meta[1] is meta[1]

    # Test initializing empty table with meta with default copy=True (deepcopy meta)
    t = table.Table(meta=meta)
    assert t.meta is not meta
    assert t.meta == meta
    assert t.meta[1] is not meta[1]


def test_replace_column_qtable():
    """Replace existing Quantity column with a new column in a QTable"""
    a = [1, 2, 3] * u.m
    b = [4, 5, 6]
    t = table.QTable([a, b], names=['a', 'b'])

    ta = t['a']
    tb = t['b']
    ta.info.meta = {'aa': [0, 1, 2, 3, 4]}
    ta.info.format = '%f'

    t.replace_column('a', a.to('cm'))
    assert np.all(t['a'] == ta)
    assert t['a'] is not ta  # New a column
    assert t['b'] is tb  # Original b column unchanged
    assert t.colnames == ['a', 'b']
    assert t['a'].info.meta is None
    assert t['a'].info.format is None


def test_replace_update_column_via_setitem():
    """
    Test table update like ``t['a'] = value``.  This leverages off the
    already well-tested ``replace_column`` and in-place update
    ``t['a'][:] = value``, so this testing is fairly light.
    """
    a = [1, 2] * u.m
    b = [3, 4]
    t = table.QTable([a, b], names=['a', 'b'])
    assert isinstance(t['a'], u.Quantity)

    # Inplace update
    ta = t['a']
    t['a'] = 5 * u.m
    assert np.all(t['a'] == [5, 5] * u.m)
    assert t['a'] is ta

    # Replace
    t['a'] = [5, 6]
    assert np.all(t['a'] == [5, 6])
    assert isinstance(t['a'], table.Column)
    assert t['a'] is not ta


def test_replace_update_column_via_setitem_warnings_normal():
    """
    Test warnings related to table replace change in #5556:
    Normal warning-free replace
    """
    t = table.Table([[1, 2, 3], [4, 5, 6]], names=['a', 'b'])
    with table.conf.set_temp('replace_warnings',
                             ['refcount', 'attributes', 'slice']):
        t['a'] = 0  # in-place update
        t['a'] = [10, 20, 30]  # replace column


def test_replace_update_column_via_setitem_warnings_slice():
    """
    Test warnings related to table replace change in #5556:
    Replace a slice, one warning.
    """
    t = table.Table([[1, 2, 3], [4, 5, 6]], names=['a', 'b'])
    with table.conf.set_temp('replace_warnings',
                             ['refcount', 'attributes', 'slice']):
        t2 = t[:2]

        t2['a'] = 0  # in-place slice update
        assert np.all(t['a'] == [0, 0, 3])

        with pytest.warns(TableReplaceWarning, match="replaced column 'a' "
                          "which looks like an array slice") as w:
            t2['a'] = [10, 20]  # replace slice
        assert len(w) == 1


def test_replace_update_column_via_setitem_warnings_attributes():
    """
    Test warnings related to table replace change in #5556:
    Lost attributes.
    """
    t = table.Table([[1, 2, 3], [4, 5, 6]], names=['a', 'b'])
    t['a'].unit = 'm'

    with pytest.warns(TableReplaceWarning, match=r"replaced column 'a' "
                      r"and column attributes \['unit'\]") as w:
        with table.conf.set_temp('replace_warnings',
                                 ['refcount', 'attributes', 'slice']):
            t['a'] = [10, 20, 30]
    assert len(w) == 1


def test_replace_update_column_via_setitem_warnings_refcount():
    """
    Test warnings related to table replace change in #5556:
    Reference count changes.
    """
    t = table.Table([[1, 2, 3], [4, 5, 6]], names=['a', 'b'])
    ta = t['a']  # noqa : Generate an extra reference to original column

    with pytest.warns(TableReplaceWarning, match="replaced column 'a' and the "
                      "number of references") as w:
        with table.conf.set_temp('replace_warnings',
                                 ['refcount', 'attributes', 'slice']):
            t['a'] = [10, 20, 30]
    assert len(w) == 1


def test_replace_update_column_via_setitem_warnings_always():
    """
    Test warnings related to table replace change in #5556:
    Test 'always' setting that raises warning for any replace.
    """
    from inspect import currentframe, getframeinfo

    t = table.Table([[1, 2, 3], [4, 5, 6]], names=['a', 'b'])

    with table.conf.set_temp('replace_warnings', ['always']):
        t['a'] = 0  # in-place slice update

        with pytest.warns(TableReplaceWarning, match="replaced column 'a'") as w:
            frameinfo = getframeinfo(currentframe())
            t['a'] = [10, 20, 30]  # replace column
        assert len(w) == 1

        # Make sure the warning points back to the user code line
        assert w[0].lineno == frameinfo.lineno + 1
        assert 'test_table' in w[0].filename


def test_replace_update_column_via_setitem_replace_inplace():
    """
    Test the replace_inplace config option related to #5556.  In this
    case no replace is done.
    """
    t = table.Table([[1, 2, 3], [4, 5, 6]], names=['a', 'b'])
    ta = t['a']
    t['a'].unit = 'm'

    with table.conf.set_temp('replace_inplace', True):
        with table.conf.set_temp('replace_warnings',
                                 ['always', 'refcount', 'attributes', 'slice']):
            t['a'] = 0  # in-place update
            assert ta is t['a']

            t['a'] = [10, 20, 30]  # normally replaces column, but not now
            assert ta is t['a']
            assert np.all(t['a'] == [10, 20, 30])


def test_primary_key_is_inherited():
    """Test whether a new Table inherits the primary_key attribute from
    its parent Table. Issue #4672"""

    t = table.Table([(2, 3, 2, 1), (8, 7, 6, 5)], names=('a', 'b'))
    t.add_index('a')
    original_key = t.primary_key

    # can't test if tuples are equal, so just check content
    assert original_key[0] == 'a'

    t2 = t[:]
    t3 = t.copy()
    t4 = table.Table(t)

    # test whether the reference is the same in the following
    assert original_key == t2.primary_key
    assert original_key == t3.primary_key
    assert original_key == t4.primary_key

    # just test one element, assume rest are equal if assert passes
    assert t.loc[1] == t2.loc[1]
    assert t.loc[1] == t3.loc[1]
    assert t.loc[1] == t4.loc[1]


def test_qtable_read_for_ipac_table_with_char_columns():
    '''Test that a char column of a QTable is assigned no unit and not
    a dimensionless unit, otherwise conversion of reader output to
    QTable fails.'''
    t1 = table.QTable([["A"]], names="B")
    out = StringIO()
    t1.write(out, format="ascii.ipac")
    t2 = table.QTable.read(out.getvalue(), format="ascii.ipac", guess=False)
    assert t2["B"].unit is None


def test_create_table_from_final_row():
    """Regression test for issue #8422: passing the last row of a table into
    Table should return a new table containing that row."""
    t1 = table.Table([(1, 2)], names=['col'])
    row = t1[-1]
    t2 = table.Table(row)['col']
    assert t2[0] == 2


def test_key_values_in_as_array():
    # Test for cheking column slicing using key_values in Table.as_array()
    data_rows = [(1, 2.0, 'x'),
                 (4, 5.0, 'y'),
                 (5, 8.2, 'z')]
    # Creating a table with three columns
    t1 = table.Table(rows=data_rows, names=('a', 'b', 'c'),
                     meta={'name': 'first table'},
                     dtype=('i4', 'f8', 'S1'))
    # Values of sliced column a,b is stored in a numpy array
    a = np.array([(1, 2.), (4, 5.), (5, 8.2)],
                 dtype=[('a', '<i4'), ('b', '<f8')])
    # Values fo sliced column c is stored in a numpy array
    b = np.array([(b'x',), (b'y',), (b'z',)], dtype=[('c', 'S1')])
    # Comparing initialised array with sliced array using Table.as_array()
    assert np.array_equal(a, t1.as_array(names=['a', 'b']))
    assert np.array_equal(b, t1.as_array(names=['c']))


def test_tolist():
    t = table.Table([[1, 2, 3], [1.1, 2.2, 3.3], [b'foo', b'bar', b'hello']],
                    names=('a', 'b', 'c'))
    assert t['a'].tolist() == [1, 2, 3]
    assert_array_equal(t['b'].tolist(), [1.1, 2.2, 3.3])
    assert t['c'].tolist() == ['foo', 'bar', 'hello']

    assert isinstance(t['a'].tolist()[0], int)
    assert isinstance(t['b'].tolist()[0], float)
    assert isinstance(t['c'].tolist()[0], str)

    t = table.Table([[[1, 2], [3, 4]],
                     [[b'foo', b'bar'], [b'hello', b'world']]],
                    names=('a', 'c'))

    assert t['a'].tolist() == [[1, 2], [3, 4]]
    assert t['c'].tolist() == [['foo', 'bar'], ['hello', 'world']]
    assert isinstance(t['a'].tolist()[0][0], int)
    assert isinstance(t['c'].tolist()[0][0], str)


class MyTable(Table):
    foo = TableAttribute()
    bar = TableAttribute(default=[])
    baz = TableAttribute(default=1)


def test_table_attribute():
    assert repr(MyTable.baz) == '<TableAttribute name=baz default=1>'

    t = MyTable([[1, 2]])
    # __attributes__ created on the fly on the first access of an attribute
    assert '__attributes__' not in t.meta
    assert t.foo is None
    assert '__attributes__' in t.meta
    t.bar.append(2.0)
    assert t.bar == [2.0]
    assert t.baz == 1

    t.baz = 'baz'
    assert t.baz == 'baz'

    # Table attributes round-trip through pickle
    tp = pickle.loads(pickle.dumps(t))
    assert tp.foo is None
    assert tp.baz == 'baz'
    assert tp.bar == [2.0]

    # Allow initialization of attributes in table creation, with / without data
    for data in None, [[1, 2]]:
        t2 = MyTable(data, foo=3, bar='bar', baz='baz')
        assert t2.foo == 3
        assert t2.bar == 'bar'
        assert t2.baz == 'baz'

    # Initializing from an existing MyTable works, with and without kwarg attrs
    t3 = MyTable(t2)
    assert t3.foo == 3
    assert t3.bar == 'bar'
    assert t3.baz == 'baz'

    t3 = MyTable(t2, foo=5, bar='fubar')
    assert t3.foo == 5
    assert t3.bar == 'fubar'
    assert t3.baz == 'baz'


@pytest.mark.skipif('not HAS_YAML')
def test_table_attribute_ecsv():
    # Table attribute round-trip through ECSV
    t = MyTable([[1, 2]], bar=[2.0], baz='baz')
    out = StringIO()
    t.write(out, format='ascii.ecsv')
    t2 = MyTable.read(out.getvalue(), format='ascii.ecsv')
    assert t2.foo is None
    assert t2.bar == [2.0]
    assert t2.baz == 'baz'


def test_table_attribute_fail():
    # Code raises ValueError(f'{attr} not allowed as TableAttribute') but in this
    # context it gets re-raised as a RuntimeError during class definition.
    with pytest.raises(RuntimeError, match='Error calling __set_name__'):
        class MyTable2(Table):
            descriptions = TableAttribute()  # Conflicts with init arg

    with pytest.raises(RuntimeError, match='Error calling __set_name__'):
        class MyTable3(Table):
            colnames = TableAttribute()  # Conflicts with built-in property


def test_set_units_fail():
    dat = [[1.0, 2.0], ['aa', 'bb']]
    with pytest.raises(ValueError, match='sequence of unit values must match number of columns'):
        Table(dat, units=[u.m])
    with pytest.raises(ValueError, match='invalid column name c for setting unit attribute'):
        Table(dat, units={'c': u.m})


def test_set_units():
    dat = [[1.0, 2.0], ['aa', 'bb'], [3, 4]]
    exp_units = (u.m, None, None)
    for cls in Table, QTable:
        for units in ({'a': u.m, 'c': ''}, exp_units):
            qt = cls(dat, units=units, names=['a', 'b', 'c'])
            if cls is QTable:
                assert isinstance(qt['a'], u.Quantity)
                assert isinstance(qt['b'], table.Column)
                assert isinstance(qt['c'], table.Column)
            for col, unit in zip(qt.itercols(), exp_units):
                assert col.info.unit is unit


def test_set_descriptions():
    dat = [[1.0, 2.0], ['aa', 'bb']]
    exp_descriptions = ('my description', None)
    for cls in Table, QTable:
        for descriptions in ({'a': 'my description'}, exp_descriptions):
            qt = cls(dat, descriptions=descriptions, names=['a', 'b'])
            for col, description in zip(qt.itercols(), exp_descriptions):
                assert col.info.description == description


def test_set_units_from_row():
    text = ['a,b',
            ',s',
            '1,2',
            '3,4']
    units = Table.read(text, format='ascii', data_start=1, data_end=2)[0]
    t = Table.read(text, format='ascii', data_start=2, units=units)
    assert isinstance(units, table.Row)
    assert t['a'].info.unit is None
    assert t['b'].info.unit is u.s


def test_set_units_descriptions_read():
    """Test setting units and descriptions via Table.read.  The test here
    is less comprehensive because the implementation is exactly the same
    as for Table.__init__ (calling Table._set_column_attribute) """
    for cls in Table, QTable:
        t = cls.read(['a b', '1 2'],
                     format='ascii',
                     units=[u.m, u.s],
                     descriptions=['hi', 'there'])
        assert t['a'].info.unit is u.m
        assert t['b'].info.unit is u.s
        assert t['a'].info.description == 'hi'
        assert t['b'].info.description == 'there'


def test_broadcasting_8933():
    """Explicitly check re-work of code related to broadcasting in #8933"""
    t = table.Table([[1, 2]])  # Length=2 table
    t['a'] = [[3, 4]]  # Can broadcast if ndim > 1 and shape[0] == 1
    t['b'] = 5
    t['c'] = [1]  # Treat as broadcastable scalar, not length=1 array (which would fail)
    assert np.all(t['a'] == [[3, 4], [3, 4]])
    assert np.all(t['b'] == [5, 5])
    assert np.all(t['c'] == [1, 1])

    # Test that broadcasted column is writeable
    t['c'][1] = 10
    assert np.all(t['c'] == [1, 10])


def test_custom_masked_column_in_nonmasked_table():
    """Test the refactor and change in column upgrades introduced
    in 95902650f.  This fixes a regression introduced by #8789
    (Change behavior of Table regarding masked columns)."""
    class MyMaskedColumn(table.MaskedColumn):
        pass

    class MySubMaskedColumn(MyMaskedColumn):
        pass

    class MyColumn(table.Column):
        pass

    class MySubColumn(MyColumn):
        pass

    class MyTable(table.Table):
        Column = MyColumn
        MaskedColumn = MyMaskedColumn

    a = table.Column([1])
    b = table.MaskedColumn([2], mask=[True])
    c = MyMaskedColumn([3], mask=[True])
    d = MySubColumn([4])
    e = MySubMaskedColumn([5], mask=[True])

    # Two different pathways for making table
    t1 = MyTable([a, b, c, d, e], names=['a', 'b', 'c', 'd', 'e'])
    t2 = MyTable()
    t2['a'] = a
    t2['b'] = b
    t2['c'] = c
    t2['d'] = d
    t2['e'] = e

    for t in (t1, t2):
        assert type(t['a']) is MyColumn
        assert type(t['b']) is MyMaskedColumn  # upgrade
        assert type(t['c']) is MyMaskedColumn
        assert type(t['d']) is MySubColumn
        assert type(t['e']) is MySubMaskedColumn  # sub-class not downgraded


def test_sort_with_mutable_skycoord():
    """Test sorting a table that has a mutable column such as SkyCoord.

    In this case the sort is done in-place
    """
    t = Table([[2, 1], SkyCoord([4, 3], [6, 5], unit='deg,deg')], names=['a', 'sc'])
    meta = {'a': [1, 2]}
    ta = t['a']
    tsc = t['sc']
    t['sc'].info.meta = meta
    t.sort('a')
    assert np.all(t['a'] == [1, 2])
    assert np.allclose(t['sc'].ra.to_value(u.deg), [3, 4])
    assert np.allclose(t['sc'].dec.to_value(u.deg), [5, 6])
    assert t['a'] is ta
    assert t['sc'] is tsc

    # Prior to astropy 4.1 this was a deep copy of SkyCoord column; after 4.1
    # it is a reference.
    t['sc'].info.meta['a'][0] = 100
    assert meta['a'][0] == 100


def test_sort_with_non_mutable():
    """Test sorting a table that has a non-mutable column.
    """
    t = Table([[2, 1], [3, 4]], names=['a', 'b'])
    ta = t['a']
    tb = t['b']
    t['b'].setflags(write=False)
    meta = {'a': [1, 2]}
    t['b'].info.meta = meta
    t.sort('a')
    assert np.all(t['a'] == [1, 2])
    assert np.all(t['b'] == [4, 3])
    assert ta is t['a']
    assert tb is not t['b']

    # Prior to astropy 4.1 this was a deep copy of SkyCoord column; after 4.1
    # it is a reference.
    t['b'].info.meta['a'][0] = 100
    assert meta['a'][0] == 1


def test_init_with_list_of_masked_arrays():
    """Test the fix for #8977"""
    m0 = np.ma.array([0, 1, 2], mask=[True, False, True])
    m1 = np.ma.array([3, 4, 5], mask=[False, True, False])
    mc = [m0, m1]

    # Test _init_from_list
    t = table.Table([mc], names=['a'])

    # Test add_column
    t['b'] = [m1, m0]

    assert t['a'].shape == (2, 3)
    assert np.all(t['a'][0] == m0)
    assert np.all(t['a'][1] == m1)
    assert np.all(t['a'][0].mask == m0.mask)
    assert np.all(t['a'][1].mask == m1.mask)

    assert t['b'].shape == (2, 3)
    assert np.all(t['b'][0] == m1)
    assert np.all(t['b'][1] == m0)
    assert np.all(t['b'][0].mask == m1.mask)
    assert np.all(t['b'][1].mask == m0.mask)


def test_data_to_col_convert_strategy():
    """Test the update to how data_to_col works (#8972), using the regression
    example from #8971.
    """
    t = table.Table([[0, 1]])
    t['a'] = 1
    t['b'] = np.int64(2)  # Failed previously
    assert np.all(t['a'] == [1, 1])
    assert np.all(t['b'] == [2, 2])


def test_rows_with_mixins():
    """Test for #9165 to allow adding a list of mixin objects.
    Also test for fix to #9357 where group_by() failed due to
    mixin object not having info.indices set to [].
    """
    tm = Time([1, 2], format='cxcsec')
    q = [1, 2] * u.m
    mixed1 = [1 * u.m, 2]  # Mixed input, fails to convert to Quantity
    mixed2 = [2, 1 * u.m]  # Mixed input, not detected as potential mixin
    rows = [(1, q[0], tm[0]),
            (2, q[1], tm[1])]
    t = table.QTable(rows=rows)
    t['a'] = [q[0], q[1]]
    t['b'] = [tm[0], tm[1]]
    t['m1'] = mixed1
    t['m2'] = mixed2

    assert np.all(t['col1'] == q)
    assert np.all(t['col2'] == tm)
    assert np.all(t['a'] == q)
    assert np.all(t['b'] == tm)
    assert np.all(t['m1'][ii] == mixed1[ii] for ii in range(2))
    assert np.all(t['m2'][ii] == mixed2[ii] for ii in range(2))
    assert type(t['m1']) is table.Column
    assert t['m1'].dtype is np.dtype(object)
    assert type(t['m2']) is table.Column
    assert t['m2'].dtype is np.dtype(object)

    # Ensure group_by() runs without failing for sortable columns.
    # The columns 'm1', and 'm2' are object dtype and not sortable.
    for name in ['col0', 'col1', 'col2', 'a', 'b']:
        t.group_by(name)

    # For good measure include exactly the failure in #9357 in which the
    # list of Time() objects is in the Table initializer.
    mjds = [Time(58000, format="mjd")]
    t = Table([mjds, ["gbt"]], names=("mjd", "obs"))
    t.group_by("obs")


def test_iterrows():
    dat = [(1, 2, 3),
           (4, 5, 6),
           (7, 8, 6)]
    t = table.Table(rows=dat, names=('a', 'b', 'c'))
    c_s = []
    a_s = []
    for c, a in t.iterrows('c', 'a'):
        a_s.append(a)
        c_s.append(c)
    assert np.all(t['a'] == a_s)
    assert np.all(t['c'] == c_s)

    rows = [row for row in t.iterrows()]
    assert rows == dat

    with pytest.raises(ValueError, match='d is not a valid column name'):
        t.iterrows('d')


def test_values_and_types():
    dat = [(1, 2, 3),
           (4, 5, 6),
           (7, 8, 6)]
    t = table.Table(rows=dat, names=('a', 'b', 'c'))
    assert isinstance(t.values(), type(OrderedDict().values()))
    assert isinstance(t.columns.values(), type(OrderedDict().values()))
    assert isinstance(t.columns.keys(), type(OrderedDict().keys()))
    for i in t.values():
        assert isinstance(i, table.column.Column)


def test_items():
    dat = [(1, 2, 3),
           (4, 5, 6),
           (7, 8, 9)]
    t = table.Table(rows=dat, names=('a', 'b', 'c'))

    assert isinstance(t.items(), type(OrderedDict({}).items()))

    for i in list(t.items()):
        assert isinstance(i, tuple)
