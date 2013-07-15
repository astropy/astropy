# Licensed under a 3-clause BSD style license - see LICENSE.rst
from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table

NUMPY_LT_1P5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False] if NUMPY_LT_1P5 else [False, True])
def table_types(request):
    class TableTypes:
        def __init__(self, request):
            self.Table = MaskedTable if request.param else table.Table
            self.Column = table.MaskedColumn if request.param else table.Column
    return TableTypes(request)


class SetupData(object):
    def _setup(self, table_types):
        self._table_type = table_types.Table
        self._column_type = table_types.Column

    @property
    def a(self):
        if self._column_type is not None:
            if not hasattr(self, '_a'):
                self._a = self._column_type(
                    name='a', data=[1, 2, 3], format='%d',
                    meta={'aa': [0, 1, 2, 3, 4]})
            return self._a

    @property
    def b(self):
        if self._column_type is not None:
            if not hasattr(self, '_b'):
                self._b = self._column_type(
                    name='b', data=[4, 5, 6], format='%d', meta={'aa': 1})
            return self._b

    @property
    def c(self):
        if self._column_type is not None:
            if not hasattr(self, '_c'):
                self._c = self._column_type(name='c', data=[7, 8, 9])
            return self._c

    @property
    def d(self):
        if self._column_type is not None:
            if not hasattr(self, '_d'):
                self._d = self._column_type(name='d', data=[7, 8, 7])
            return self._d

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
        """Set a row from an incorrectly-sized set of values"""
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
        """Column with name=None can init a table IFF names are supplied"""
        c = table_types.Column(data=[1, 2])
        table_types.Table([c], names=('c',))
        with pytest.raises(TypeError):
            table_types.Table([c])


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
        t.add_column(self.b)
        t.add_column(self.c)
        assert t.colnames == ['a', 'b', 'c']

    def test_add_duplicate_columns(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b, self.c])
        with pytest.raises(ValueError):
            t.add_columns([table_types.Column(name='a', data=[0, 1, 2]), table_types.Column(name='b', data=[0, 1, 2])])
        t.add_column(self.d)
        assert t.colnames == ['a', 'b', 'c', 'd']


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
    def t(self):
        if self._table_type is not None:
            if not hasattr(self, '_t'):
                self._t = self._table_type([self.a, self.b, self.c])
            return self._t

    def test_add_none_to_empty_table(self, table_types):
        self._setup(table_types)
        t = table_types.Table(names=('a', 'b'), dtype=('i', 'S4'))
        t.add_row()
        assert t['a'][0] == 0
        assert t['b'][0] == b''
        t.add_row()
        assert t['a'][1] == 0
        assert t['b'][1] == b''

    def test_add_stuff_to_empty_table(self, table_types):
        self._setup(table_types)
        t = table_types.Table(names=('a', 'b'), dtype=('i', 'S8'))
        t.add_row([1, 'hello'])
        assert t['a'][0] == 1
        assert t['b'][0] == b'hello'
        # Make sure it is not repeating last row but instead
        # adding zeros (as documented)
        t.add_row()
        assert t['a'][1] == 0
        assert t['b'][1] == b''

    def test_add_table_row(self, table_types):
        self._setup(table_types)
        t = self.t
        t2 = table_types.Table([self.a, self.b, self.c])
        t.add_row(t2[0])
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 1]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 4.0]))
        assert np.all(t['c'] == np.array(['7', '8', '9', '7']))

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

    def test_add_without_own_fails(self, table_types):
        """Add row to a table that doesn't own the data"""
        self._setup(table_types)
        data = np.array([(1, 2, 3),
                         (3, 4, 5)],
                        dtype='i4')
        t = table_types.Table(data, copy=False)
        if not t.masked:
            with pytest.raises(ValueError):
                t.add_row([6, 7, 8])


@pytest.mark.usefixtures('table_types')
class TestTableColumn(SetupData):

    def test_column_view(self, table_types):
        self._setup(table_types)
        t = self.t
        a = t.columns['a']
        a[2] = 10
        assert t._data['a'][2] == 10


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
        assert self.t._data is None

    def test_2(self, table_types):
        self._setup(table_types)
        self.t.add_column(self.b)
        self.t.remove_columns('a')
        assert self.t.columns.keys() == ['b']
        assert self.t._data.dtype.names == ('b',)
        assert np.all(self.t['b'] == np.array([4, 5, 6]))

    def test_remove_nonexistent_row(self, table_types):
        self._setup(table_types)
        with pytest.raises(ValueError):
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
        assert self.t.dtype == np.dtype([('a', '<i8'),
                                         ('b', '<i8')])

    def test_delitem1(self, table_types):
        self._setup(table_types)
        del self.t['a']
        assert self.t.columns.keys() == []
        assert self.t._data is None

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
        assert t._data is None

    def test_2(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.keep_columns('b')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([4, 5, 6]))


@pytest.mark.usefixtures('table_types')
class TestRename(SetupData):

    def test_1(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a])
        t.rename_column('a', 'b')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([1, 2, 3]))

    def test_2(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t.rename_column('a', 'c')
        t.rename_column('b', 'a')
        assert t.columns.keys() == ['c', 'a']
        assert t._data.dtype.names == ('c', 'a')
        if t.masked:
            assert t._data.mask.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))

    def test_rename_by_attr(self, table_types):
        self._setup(table_types)
        t = table_types.Table([self.a, self.b])
        t['a'].name = 'c'
        t['b'].name = 'a'
        assert t.columns.keys() == ['c', 'a']
        assert t._data.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))


@pytest.mark.usefixtures('table_types')
class TestSort():

    def test_single(self, table_types):
        t = table_types.Table()
        t.add_column(table_types.Column(name='a', data=[2, 1, 3]))
        t.add_column(table_types.Column(name='b', data=[6, 5, 4]))
        assert np.all(t['a'] == np.array([2, 1, 3]))
        assert np.all(t['b'] == np.array([6, 5, 4]))
        t.sort('a')
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['b'] == np.array([5, 6, 4]))
        t.sort('b')
        assert np.all(t['a'] == np.array([3, 1, 2]))
        assert np.all(t['b'] == np.array([4, 5, 6]))

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
            assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d._data)
            assert np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[('c', 'i8'), ('d', 'i8')])


def _assert_copies(t, t2, deep=True):
    assert t.colnames == t2.colnames
    np.testing.assert_array_equal(t._data, t2._data)
    assert t.meta == t2.meta

    if deep:
        assert not np.may_share_memory(t._data, t2._data)
    else:
        assert np.may_share_memory(t._data, t2._data)


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
    from copy import copy, deepcopy

    t = table.Table([[1, 2, 3], [2, 3, 4]], names=['x', 'y'])

    t2 = copy(t)
    t3 = deepcopy(t)

    _assert_copies(t, t2, deep=False)
    _assert_copies(t, t3)
