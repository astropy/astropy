import pytest
import numpy as np

from .. import Table, Column


class TestEmptyData():

    def test_1(self):
        t = Table()
        t.append_column(Column('a', datatype=int, length=100))
        assert len(t['a']) == 100

    def test_2(self):
        t = Table()
        t.append_column(Column('a', datatype=int, shape=(3, ), length=100))
        assert len(t['a']) == 100

    def test_3(self):
        t = Table()  # length is not given
        t.append_column(Column('a', datatype=int))
        assert len(t['a']) == 0

    def test_4(self):
        t = Table()  # length is not given
        t.append_column(Column('a', datatype=int, shape=(3, 4)))
        assert len(t['a']) == 0

    def test_5(self):
        t = Table()
        t.append_column(Column('a'))  # dtype is not specified
        assert len(t['a']) == 0


class TestNewFromColumns():

    def test_simple(self):
        cols = [Column('a', [1, 2, 3]),
                Column('b', [4, 5, 6], datatype=np.float32)]
        t = Table(cols)
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['b'] == np.array([4, 5, 6], dtype=np.float32))
        assert type(t['b'][1]) == np.float32

    def test_size_mismatch(self):
        cols = [Column('a', [1, 2, 3]),
                Column('b', [4, 5, 6, 7])]
        with pytest.raises(ValueError):
            Table(cols)


class TestColumnAccess():

    def test_1(self):
        t = Table()
        with pytest.raises(KeyError):
            t['a']

    def test_2(self):
        t = Table()
        t.append_column(Column('a', [1, 2, 3]))
        assert np.all(t['a'] == np.array([1, 2, 3]))
        with pytest.raises(KeyError):
            t['b']  # column does not exist


class TestAddLength():
    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])
        self.b = Column('b', [4, 5, 6])

    def test_right_length(self):
        t = Table([self.a])
        t.append_column(self.b)

    def test_too_long(self):
        t = Table([self.a])
        with pytest.raises(ValueError):
            t.append_column(Column('b', [4, 5, 6, 7]))  # data too long

    def test_too_short(self):
        t = Table([self.a])
        with pytest.raises(ValueError):
            t.append_column(Column('b', [4, 5]))  # data too short


class TestAddPosition():
    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])
        self.b = Column('b', [4, 5, 6])
        self.c = Column('c', [7, 8, 9])

    def test_1(self):
        t = Table()
        t.insert_column(0, self.a)

    def test_2(self):
        t = Table()
        t.insert_column(1, self.a)

    def test_3(self):
        t = Table()
        t.insert_column(-1, self.a)

    def test_5(self):
        t = Table()
        with pytest.raises(ValueError):
            t.index_column('b')

    def test_6(self):
        t = Table()
        t.append_column(self.a)
        t.append_column(self.b)
        assert t.columns.keys() == ['a', 'b']

    def test_7(self):
        t = Table([self.a])
        t.insert_column(t.index_column('a'), self.b)
        assert t.columns.keys() == ['b', 'a']

    def test_8(self):
        t = Table([self.a])
        t.insert_column(t.index_column('a') + 1, self.b)
        assert t.columns.keys() == ['a', 'b']

    def test_9(self):
        t = Table()
        t.append_column(self.a)
        t.insert_column(t.index_column('a') + 1, self.b)
        t.insert_column(t.index_column('b'), self.c)
        assert t.columns.keys() == ['a', 'c', 'b']

    def test_10(self):
        t = Table()
        t.append_column(self.a)
        ia = t.index_column('a')
        t.insert_column(ia + 1, self.b)
        t.insert_column(ia, self.c)
        assert t.columns.keys() == ['c', 'a', 'b']


class TestArrayColumns():

    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])

    def test_1d(self):
        b = Column('b', datatype=int, shape=(2, ), length=3)
        t = Table([self.a])
        t.append_column(b)
        assert t['b'].shape == (3, 2)
        assert t['b'][0].shape == (2, )

    def test_2d(self):
        b = Column('b', datatype=int, shape=(2, 4), length=3)
        t = Table([self.a])
        t.append_column(b)
        assert t['b'].shape == (3, 2, 4)
        assert t['b'][0].shape == (2, 4)

    def test_3d(self):
        t = Table([self.a])
        b = Column('b', datatype=int, shape=(2, 4, 6), length=3)
        t.append_column(b)
        assert t['b'].shape == (3, 2, 4, 6)
        assert t['b'][0].shape == (2, 4, 6)


class TestRemove():

    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])
        self.b = Column('b', [4, 5, 6])

    def test_1(self):
        t = Table([self.a])
        t.remove_columns('a')
        assert t.columns.keys() == []
        assert t._data is None

    def test_2(self):
        t = Table([self.a])
        t.append_column(self.b)
        t.remove_columns('a')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([4, 5, 6]))


class TestKeep():

    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])
        self.b = Column('b', [4, 5, 6])

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


class TestRename():

    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])
        self.b = Column('b', [4, 5, 6])

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


class TestSort():

    def test_single(self):
        t = Table()
        t.append_column(Column('a', [2, 1, 3]))
        t.append_column(Column('b', [6, 5, 4]))
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
        t.append_column(Column('a', [2, 1, 3, 2, 3, 1]))
        t.append_column(Column('b', [6, 5, 4, 3, 5, 4]))
        assert np.all(t['a'] == np.array([2, 1, 3, 2, 3, 1]))
        assert np.all(t['b'] == np.array([6, 5, 4, 3, 5, 4]))
        t.sort(['a', 'b'])
        assert np.all(t['a'] == np.array([1, 1, 2, 2, 3, 3]))
        assert np.all(t['b'] == np.array([4, 5, 3, 6, 4, 5]))
        t.sort(['b', 'a'])
        assert np.all(t['a'] == np.array([2, 1, 3, 1, 3, 2]))
        assert np.all(t['b'] == np.array([3, 4, 4, 5, 5, 6]))
