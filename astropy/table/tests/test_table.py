from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')

# Dummy init of Table, DATA for pyflakes and to be sure test fixture is working
Table = None
Column = None


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False] if numpy_lt_1p5 else [False, True])
def set_global_Table(request):
    global Table, Column

    Table = MaskedTable if request.param else table.Table
    Column = table.MaskedColumn if request.param else table.Column


class SetupData(object):
    @property
    def a(self):
        if Column is not None:
            if not hasattr(self, '_a'):
                self._a = Column('a', [1, 2, 3], meta={'aa': np.arange(5)})
            return self._a

    @property
    def b(self):
        if Column is not None:
            if not hasattr(self, '_b'):
                self._b = Column('b', [4, 5, 6])
            return self._b

    @property
    def c(self):
        if Column is not None:
            if not hasattr(self, '_c'):
                self._c = Column('c', [7, 8, 9])
            return self._c

    @property
    def d(self):
        if Column is not None:
            if not hasattr(self, '_d'):
                self._d = Column('d', [7, 8, 7])
            return self._d

    @property
    def t(self):
        if Table is not None:
            if not hasattr(self, '_t'):
                self._t = Table([self.a, self.b])
            return self._t


@pytest.mark.usefixtures('set_global_Table')
class TestEmptyData():

    def test_1(self):
        t = Table()
        t.add_column(Column('a', dtype=int, length=100))
        assert len(t['a']) == 100

    def test_2(self):
        t = Table()
        t.add_column(Column('a', dtype=int, shape=(3, ), length=100))
        assert len(t['a']) == 100

    def test_3(self):
        t = Table()  # length is not given
        t.add_column(Column('a', dtype=int))
        assert len(t['a']) == 0

    def test_4(self):
        t = Table()  # length is not given
        t.add_column(Column('a', dtype=int, shape=(3, 4)))
        assert len(t['a']) == 0

    def test_5(self):
        t = Table()
        t.add_column(Column('a'))  # dtype is not specified
        assert len(t['a']) == 0


@pytest.mark.usefixtures('set_global_Table')
class TestNewFromColumns():

    def test_simple(self):
        cols = [Column('a', [1, 2, 3]),
                Column('b', [4, 5, 6], dtype=np.float32)]
        t = Table(cols)
        assert np.all(t['a'].data == np.array([1, 2, 3]))
        assert np.all(t['b'].data == np.array([4, 5, 6], dtype=np.float32))
        assert type(t['b'][1]) == np.float32

    def test_from_np_array(self):
        cols = [Column('a', np.array([1, 2, 3], dtype=np.int64),
                       dtype=np.float64),
                Column('b', np.array([4, 5, 6], dtype=np.float32))]
        t = Table(cols)
        assert np.all(t['a'] == np.array([1, 2, 3], dtype=np.float64))
        assert np.all(t['b'] == np.array([4, 5, 6], dtype=np.float32))
        assert type(t['a'][1]) == np.float64
        assert type(t['b'][1]) == np.float32

    def test_size_mismatch(self):
        cols = [Column('a', [1, 2, 3]),
                Column('b', [4, 5, 6, 7])]
        with pytest.raises(ValueError):
            Table(cols)


@pytest.mark.usefixtures('set_global_Table')
class TestReverse():

    def test_reverse(self):
        t = Table([[1, 2, 3],
                   ['a', 'b', 'cc']])
        t.reverse()
        assert np.all(t['col0'] == np.array([3, 2, 1]))
        assert np.all(t['col1'] == np.array(['cc', 'b', 'a']))

        t2 = Table(t, copy=False)
        assert np.all(t2['col0'] == np.array([3, 2, 1]))
        assert np.all(t2['col1'] == np.array(['cc', 'b', 'a']))

        t2 = Table(t, copy=True)
        assert np.all(t2['col0'] == np.array([3, 2, 1]))
        assert np.all(t2['col1'] == np.array(['cc', 'b', 'a']))

        t2.sort('col0')
        assert np.all(t2['col0'] == np.array([1, 2, 3]))
        assert np.all(t2['col1'] == np.array(['a', 'b', 'cc']))


@pytest.mark.usefixtures('set_global_Table')
class TestColumnAccess():

    def test_1(self):
        t = Table()
        with pytest.raises(KeyError):
            t['a']

    def test_2(self):
        t = Table()
        t.add_column(Column('a', [1, 2, 3]))
        assert np.all(t['a'] == np.array([1, 2, 3]))
        with pytest.raises(KeyError):
            t['b']  # column does not exist


@pytest.mark.usefixtures('set_global_Table')
class TestAddLength(SetupData):

    def test_right_length(self):
        t = Table([self.a])
        t.add_column(self.b)

    def test_too_long(self):
        t = Table([self.a])
        with pytest.raises(ValueError):
            t.add_column(Column('b', [4, 5, 6, 7]))  # data too long

    def test_too_short(self):
        t = Table([self.a])
        with pytest.raises(ValueError):
            t.add_column(Column('b', [4, 5]))  # data too short


@pytest.mark.usefixtures('set_global_Table')
class TestAddPosition(SetupData):

    def test_1(self):
        t = Table()
        t.add_column(self.a, 0)

    def test_2(self):
        t = Table()
        t.add_column(self.a, 1)

    def test_3(self):
        t = Table()
        t.add_column(self.a, -1)

    def test_5(self):
        t = Table()
        with pytest.raises(ValueError):
            t.index_column('b')

    def test_6(self):
        t = Table()
        t.add_column(self.a)
        t.add_column(self.b)
        assert t.columns.keys() == ['a', 'b']

    def test_7(self):
        t = Table([self.a])
        t.add_column(self.b, t.index_column('a'))
        assert t.columns.keys() == ['b', 'a']

    def test_8(self):
        t = Table([self.a])
        t.add_column(self.b, t.index_column('a') + 1)
        assert t.columns.keys() == ['a', 'b']

    def test_9(self):
        t = Table()
        t.add_column(self.a)
        t.add_column(self.b, t.index_column('a') + 1)
        t.add_column(self.c, t.index_column('b'))
        assert t.columns.keys() == ['a', 'c', 'b']

    def test_10(self):
        t = Table()
        t.add_column(self.a)
        ia = t.index_column('a')
        t.add_column(self.b, ia + 1)
        t.add_column(self.c, ia)
        assert t.columns.keys() == ['c', 'a', 'b']


@pytest.mark.usefixtures('set_global_Table')
class TestInitFromTable(SetupData):

    def test_from_table_cols(self):
        """Ensure that using cols from an existing table gives
        a clean copy.
        """
        t = self.t
        cols = t.columns
        # Construct Table with cols via Table._new_from_cols
        t2a = Table([cols['a'], cols['b'], self.c])

        # Construct with add_column
        t2b = Table()
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


@pytest.mark.usefixtures('set_global_Table')
class TestAddColumns(SetupData):

    def test_add_columns1(self):
        t = Table()
        t.add_columns([self.a, self.b, self.c])
        assert t.colnames == ['a', 'b', 'c']

    def test_add_columns2(self):
        t = Table([self.a, self.b])
        t.add_columns([self.c, self.d])
        assert t.colnames == ['a', 'b', 'c', 'd']
        assert np.all(t['c'] == np.array([7, 8, 9]))

    def test_add_columns3(self):
        t = Table([self.a, self.b])
        t.add_columns([self.c, self.d], indexes=[1, 0])
        assert t.colnames == ['d', 'a', 'c', 'b']

    def test_add_columns4(self):
        t = Table([self.a, self.b])
        t.add_columns([self.c, self.d], indexes=[0, 0])
        assert t.colnames == ['c', 'd', 'a', 'b']

    def test_add_columns5(self):
        t = Table([self.a, self.b])
        t.add_columns([self.c, self.d], indexes=[2, 2])
        assert t.colnames == ['a', 'b', 'c', 'd']

    def test_add_duplicate_column(self):
        t = Table()
        t.add_column(self.a)
        with pytest.raises(ValueError):
            t.add_column(Column('a', [0, 1, 2]))
        t.add_column(self.b)
        t.add_column(self.c)
        assert t.colnames == ['a', 'b', 'c']

    def test_add_duplicate_columns(self):
        t = Table([self.a, self.b, self.c])
        with pytest.raises(ValueError):
            t.add_columns([Column('a', [0, 1, 2]), Column('b', [0, 1, 2])])
        t.add_column(self.d)
        assert t.colnames == ['a', 'b', 'c', 'd']


@pytest.mark.usefixtures('set_global_Table')
class TestAddRow(SetupData):

    @property
    def b(self):
        if Column is not None:
            if not hasattr(self, '_b'):
                self._b = Column('b', [4.0, 5.1, 6.2])
            return self._b

    @property
    def c(self):
        if Column is not None:
            if not hasattr(self, '_c'):
                self._c = Column('c', ['7', '8', '9'])
            return self._c

    @property
    def t(self):
        if Table is not None:
            if not hasattr(self, '_t'):
                self._t = Table([self.a, self.b, self.c])
            return self._t

    def test_add_table_row(self):
        t = self.t
        t2 = Table([self.a, self.b, self.c])
        t.add_row(t2[0])
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 1]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 4.0]))
        assert np.all(t['c'] == np.array(['7', '8', '9', '7']))

    def test_add_with_tuple(self):
        t = self.t
        t.add_row((4, 7.2, '1'))
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 4]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 7.2]))
        assert np.all(t['c'] == np.array(['7', '8', '9', '1']))

    def test_add_with_list(self):
        t = self.t
        t.add_row([4, 7.2, '10'])
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 4]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 7.2]))
        assert np.all(t['c'] == np.array(['7', '8', '9', '1']))

    def test_add_with_dict(self):
        t = self.t
        t.add_row({'a': 4, 'b': 7.2})
        assert len(t) == 4
        assert np.all(t['a'] == np.array([1, 2, 3, 4]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 7.2]))
        if t.masked:
            assert np.all(t['c'] == np.array(['7', '8', '9', '7']))
        else:
            assert np.all(t['c'] == np.array(['7', '8', '9', '']))

    def test_add_with_none(self):
        t = self.t
        t.add_row()
        assert len(t) == 4
        if t.masked:
            assert np.all(t['a'].data == np.array([1, 2, 3, 1]))
            assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 4.0]))
            assert np.all(t['c'].data == np.array(['7', '8', '9', '7']))
        else:
            assert np.all(t['a'].data == np.array([1, 2, 3, 0]))
            assert np.allclose(t['b'], np.array([4.0, 5.1, 6.2, 0.0]))
            assert np.all(t['c'].data == np.array(['7', '8', '9', '']))

    def test_add_missing_column(self):
        t = self.t
        with pytest.raises(ValueError):
            t.add_row({'bad_column': 1})

    def test_wrong_size_tuple(self):
        t = self.t
        with pytest.raises(ValueError):
            t.add_row((1, 2))

    def test_wrong_vals_type(self):
        t = self.t
        with pytest.raises(TypeError):
            t.add_row(1)

    def test_add_without_own_fails(self):
        """Add row to a table that doesn't own the data"""
        data = np.array([(1, 2, 3),
                         (3, 4, 5)],
                        dtype='i4')
        t = Table(data, copy=False)
        if not t.masked:
            with pytest.raises(ValueError):
                t.add_row([6, 7, 8])


@pytest.mark.usefixtures('set_global_Table')
class TestTableColumn(SetupData):

    def test_column_view(self):
        t = self.t
        a = t.columns['a']
        a[2] = 10
        assert t._data['a'][2] == 10


@pytest.mark.usefixtures('set_global_Table')
class TestArrayColumns(SetupData):

    def test_1d(self):
        b = Column('b', dtype=int, shape=(2, ), length=3)
        t = Table([self.a])
        t.add_column(b)
        assert t['b'].shape == (3, 2)
        assert t['b'][0].shape == (2, )

    def test_2d(self):
        b = Column('b', dtype=int, shape=(2, 4), length=3)
        t = Table([self.a])
        t.add_column(b)
        assert t['b'].shape == (3, 2, 4)
        assert t['b'][0].shape == (2, 4)

    def test_3d(self):
        t = Table([self.a])
        b = Column('b', dtype=int, shape=(2, 4, 6), length=3)
        t.add_column(b)
        assert t['b'].shape == (3, 2, 4, 6)
        assert t['b'][0].shape == (2, 4, 6)


@pytest.mark.usefixtures('set_global_Table')
class TestRemove(SetupData):

    @property
    def t(self):
        if Table is not None:
            if not hasattr(self, '_t'):
                self._t = Table([self.a])
            return self._t

    @property
    def t2(self):
        if Table is not None:
            if not hasattr(self, '_t2'):
                self._t2 = Table([self.a, self.b, self.c])
            return self._t2

    def test_1(self):
        self.t.remove_columns('a')
        assert self.t.columns.keys() == []
        assert self.t._data is None

    def test_2(self):
        self.t.add_column(self.b)
        self.t.remove_columns('a')
        assert self.t.columns.keys() == ['b']
        assert self.t._data.dtype.names == ('b',)
        assert np.all(self.t['b'] == np.array([4, 5, 6]))

    def test_delitem1(self):
        del self.t['a']
        assert self.t.columns.keys() == []
        assert self.t._data is None

    def test_delitem2(self):
        del self.t2['b']
        assert self.t2.colnames == ['a', 'c']

    def test_delitems(self):
        del self.t2['a', 'b']
        assert self.t2.colnames == ['c']

    def test_delitem_fail(self):
        with pytest.raises(KeyError):
            del self.t['d']


@pytest.mark.usefixtures('set_global_Table')
class TestKeep(SetupData):

    def test_1(self):
        t = Table([self.a, self.b])
        t.keep_columns([])
        assert t.columns.keys() == []
        assert t._data is None

    def test_2(self):
        t = Table([self.a, self.b])
        t.keep_columns('b')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([4, 5, 6]))


@pytest.mark.usefixtures('set_global_Table')
class TestRename(SetupData):

    def test_1(self):
        t = Table([self.a])
        t.rename_column('a', 'b')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([1, 2, 3]))

    def test_2(self):
        t = Table([self.a, self.b])
        t.rename_column('a', 'c')
        t.rename_column('b', 'a')
        assert t.columns.keys() == ['c', 'a']
        assert t._data.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))

    def test_rename_by_attr(self):
        t = Table([self.a, self.b])
        t['a'].name = 'c'
        t['b'].name = 'a'
        assert t.columns.keys() == ['c', 'a']
        assert t._data.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))


@pytest.mark.usefixtures('set_global_Table')
class TestSort():

    def test_single(self):
        t = Table()
        t.add_column(Column('a', [2, 1, 3]))
        t.add_column(Column('b', [6, 5, 4]))
        assert np.all(t['a'] == np.array([2, 1, 3]))
        assert np.all(t['b'] == np.array([6, 5, 4]))
        t.sort('a')
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['b'] == np.array([5, 6, 4]))
        t.sort('b')
        assert np.all(t['a'] == np.array([3, 1, 2]))
        assert np.all(t['b'] == np.array([4, 5, 6]))

    def test_multiple(self):
        t = Table()
        t.add_column(Column('a', [2, 1, 3, 2, 3, 1]))
        t.add_column(Column('b', [6, 5, 4, 3, 5, 4]))
        assert np.all(t['a'] == np.array([2, 1, 3, 2, 3, 1]))
        assert np.all(t['b'] == np.array([6, 5, 4, 3, 5, 4]))
        t.sort(['a', 'b'])
        assert np.all(t['a'] == np.array([1, 1, 2, 2, 3, 3]))
        assert np.all(t['b'] == np.array([4, 5, 3, 6, 4, 5]))
        t.sort(['b', 'a'])
        assert np.all(t['a'] == np.array([2, 1, 3, 1, 3, 2]))
        assert np.all(t['b'] == np.array([3, 4, 4, 5, 5, 6]))


@pytest.mark.usefixtures('set_global_Table')
class TestIterator():

    def test_iterator(self):
        d = np.array([(2, 1),
                      (3, 6),
                      (4, 5)], dtype=[('a', 'i4'), ('b', 'i4')])
        t = Table(d)
        if t.masked:
            with pytest.raises(ValueError):
                t[0] == d[0]
        else:
            for row, np_row in zip(t, d):
                assert np.all(row == np_row)


@pytest.mark.usefixtures('set_global_Table')
class TestSetMeta():

    def test_set_meta(self):
        d = Table(names=('a', 'b'))
        d.meta['a'] = 1
        d.meta['b'] = 1
        d.meta['c'] = 1
        d.meta['d'] = 1
        assert list(d.meta.keys()) == ['a', 'b', 'c', 'd']


@pytest.mark.usefixtures('set_global_Table')
class TestConvertNumpyArray():

    def test_convert_numpy_array(self):
        d = Table([[1, 2], [3, 4]], names=('a', 'b'))

        np_data = np.array(d)
        if Table is not MaskedTable:
            assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        if Table is not MaskedTable:
            assert np.all(np_data == d._data)
            assert np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[('c', 'i8'), ('d', 'i8')])
