import pytest
import numpy as np

from .. import Table, ArgumentError


class TestEmptyData():

    def test_1(self):
        t = Table(length=100)
        t.add_column('a', dtype=int)

    def test_2(self):
        t = Table(length=100)
        t.add_column('a', dtype=int, shape=(3, ))

    def test_3(self):
        t = Table()  # length is not given
        with pytest.raises(ArgumentError):
            t.add_column('a', dtype=int)

    def test_4(self):
        t = Table()  # length is not given
        with pytest.raises(ArgumentError):
            t.add_column('a', dtype=int, shape=(3, 4))

    def test_5(self):
        t = Table()
        with pytest.raises(ArgumentError):
            t.add_column('a')  # dtype is not specified


class TestColumnAccess():

    def test_1(self):
        t = Table()
        with pytest.raises(KeyError):
            t['a']

    def test_2(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t['a']
        with pytest.raises(KeyError):
            t['b']  # column does not exist


class TestAddLength():

    def test_right_length(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6])

    def test_too_long(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        with pytest.raises(ValueError):
            t.add_column('b', [4, 5, 6, 7])  # data is too long

    def test_too_short(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        with pytest.raises(ValueError):
            t.add_column('b', [4, 5])  # data is too short


class TestAddPosition():

    def test_1(self):
        t = Table()
        t.add_column('a', [1, 2, 3], position=0)

    def test_2(self):
        t = Table()
        with pytest.raises(ValueError):
            t.add_column('a', [1, 2, 3], position=1)  # invalid position

    def test_3(self):
        t = Table()
        with pytest.raises(ValueError):
            t.add_column('a', [1, 2, 3], position=-1)  # invalid position

    def test_4(self):
        t = Table()
        with pytest.raises(KeyError):
            t.add_column('a', [1, 2, 3], before='b')  # 'b' does not exist

    def test_5(self):
        t = Table()
        with pytest.raises(KeyError):
            t.add_column('a', [1, 2, 3], after='b')  # 'b' does not exist

    def test_6(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6])
        assert t.columns.keys() == ['a', 'b']

    def test_7(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6], before='a')
        assert t.columns.keys() == ['b', 'a']

    def test_8(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6], after='a')
        assert t.columns.keys() == ['a', 'b']

    def test_9(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6], after='a')
        t.add_column('c', [7, 8, 9], before='b')
        assert t.columns.keys() == ['a', 'c', 'b']

    def test_10(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6], after='a')
        t.add_column('c', [7, 8, 9], before='a')
        assert t.columns.keys() == ['c', 'a', 'b']


class TestArrayColumns():

    def test_1d(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', dtype=int, shape=(2, ))
        assert t['b'].shape == (3, 2)
        assert t['b'][0].shape == (2, )

    def test_2d(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', dtype=int, shape=(2, 4))
        assert t['b'].shape == (3, 2, 4)
        assert t['b'][0].shape == (2, 4)

    def test_3d(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', dtype=int, shape=(2, 4, 6))
        assert t['b'].shape == (3, 2, 4, 6)
        assert t['b'][0].shape == (2, 4, 6)


class TestRemove():

    def test_1(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.remove_columns('a')
        assert t.columns.keys() == []
        assert t._data is None

    def test_2(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6])
        t.remove_columns('a')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([4, 5, 6]))


class TestKeep():

    def test_1(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6])
        t.keep_columns([])
        assert t.columns.keys() == []
        assert t._data is None

    def test_2(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6])
        t.keep_columns('b')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([4, 5, 6]))


class TestRename():

    def test_1(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.rename_column('a', 'b')
        assert t.columns.keys() == ['b']
        assert t._data.dtype.names == ('b',)
        assert np.all(t['b'] == np.array([1, 2, 3]))

    def test_2(self):
        t = Table()
        t.add_column('a', [1, 2, 3])
        t.add_column('b', [4, 5, 6])
        t.rename_column('a', 'c')
        t.rename_column('b', 'a')
        assert t.columns.keys() == ['c', 'a']
        assert t._data.dtype.names == ('c', 'a')
        assert np.all(t['c'] == np.array([1, 2, 3]))
        assert np.all(t['a'] == np.array([4, 5, 6]))
