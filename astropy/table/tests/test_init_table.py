from __future__ import print_function  # For print debugging with python 2 or 3

import pytest
import numpy as np

from .. import Table, Column


class BaseInitFrom():

    def test_basic_init(self):
        t = Table(self.data, names=('a', 'b', 'c'))
        assert t.colnames == ['a', 'b', 'c']
        assert np.all(t['a'] == np.array([1, 3]))
        assert np.all(t['b'] == np.array([2, 4]))
        assert np.all(t['c'] == np.array([3, 5]))

    def test_set_dtypes(self):
        t = Table(self.data, names=('a', 'b', 'c'), dtypes=('i4', 'f4', 'f8'))
        assert t.colnames == ['a', 'b', 'c']
        assert np.all(t['a'] == np.array([1, 3], dtype='i4'))
        assert np.all(t['b'] == np.array([2, 4], dtype='f4'))
        assert np.all(t['c'] == np.array([3, 5], dtype='f8'))
        assert t['a'].dtype.type == np.int32
        assert t['b'].dtype.type == np.float32
        assert t['c'].dtype.type == np.float64

    def test_names_dtypes_mismatch(self):
        with pytest.raises(ValueError):
            Table(self.data, names=('a',), dtypes=('i4', 'f4'))


class BaseInitFromListLike(BaseInitFrom):

    @pytest.mark.xfail
    def test_names_cols_mismatch(self):
        with pytest.raises(ValueError):
            Table(self.data, names=['a'], dtypes=[int])

    def test_names_copy_false(self):
        with pytest.raises(ValueError):
            Table(self.data, names=['a'], dtypes=[int], copy=False)


class BaseInitFromDictLike(BaseInitFrom):

    @pytest.mark.xfail
    def test_select_names(self):
        t = Table(self.data, names=('b', 'a'))
        assert t.colnames == ['b', 'a']
        assert np.all(t['a'] == np.array([1, 3]))
        assert np.all(t['b'] == np.array([2, 4]))


class TestInitFromNdarrayHomo(BaseInitFromListLike):

    def setup_method(self, method):
        self.data = np.array([(1, 2, 3),
                             (3, 4, 5)],
                            dtype='i4')

    def test_default_names(self):
        t = Table(self.data)
        assert t.colnames == ['col0', 'col1', 'col2']

    def test_ndarray_ref(self):
        """Init with ndarray and copy=False and show that table uses reference
        to input ndarray"""
        t = Table(self.data, copy=False)
        t['col1'][1] = 0
        assert t._data['col1'][1] == 0
        assert self.data[1][1] == 0
        # NOTE: assert np.all(t._data == self.data) fails because when
        # homogenous array is viewcast to structured then the == is False

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['a', None, 'c'], dtypes=[None, None, 'f8'])
        assert t.colnames == ['a', 'col1', 'c']
        assert t['a'].dtype.type == np.int32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.float64

    def test_partial_names_ref(self):
        t = Table(self.data, names=['a', None, 'c'])
        assert t.colnames == ['a', 'col1', 'c']
        assert t['a'].dtype.type == np.int32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.int32


class TestInitFromList(BaseInitFromListLike):

    def setup_method(self, method):
        self.data = [Column('a', [1, 3]),
                     [2, 4],
                     np.array([3, 5], dtype='i8')]

    def test_default_names(self):
        t = Table(self.data)
        assert t.colnames == ['a', 'col1', 'col2']

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['b', None, 'c'], dtypes=['f4', None, 'f8'])
        assert t.colnames == ['b', 'col1', 'c']
        assert t['b'].dtype.type == np.float32
        assert t['col1'].dtype.type == np.int64
        assert t['c'].dtype.type == np.float64

    @pytest.mark.xfail
    def test_no_copy(self):
        with pytest.raises(ValueError):
            Table(self.data, copy=False)


class TestInitFromNdarrayStruct(BaseInitFromDictLike):

    def setup_method(self, method):
        self.data = np.array([(1, 2, 3),
                              (3, 4, 5)],
                             dtype=[('a', 'i8'), ('b', 'i4'), ('c', 'i8')])

    def test_ndarray_ref(self):
        """Init with ndarray and copy=False and show that table uses reference
        to input ndarray"""
        t = Table(self.data, copy=False)
        assert np.all(t._data == self.data)
        t['a'][1] = 0
        assert t._data['a'][1] == 0
        assert self.data['a'][1] == 0
        assert np.all(t._data == self.data)

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['e', None, 'd'], dtypes=['f4', None, 'f8'])
        assert t.colnames == ['e', 'b', 'd']
        assert t['e'].dtype.type == np.float32
        assert t['b'].dtype.type == np.int32
        assert t['d'].dtype.type == np.float64

    def test_partial_names_ref(self):
        t = Table(self.data, names=['e', None, 'd'], copy=False)
        assert t.colnames == ['e', 'b', 'd']
        assert t['e'].dtype.type == np.int64
        assert t['b'].dtype.type == np.int32
        assert t['d'].dtype.type == np.int64


class TestInitFromDict(BaseInitFromDictLike):

    def setup_method(self, method):
        self.data = dict([('a', Column('a', [1, 3])),
                          ('b', [2, 4]),
                          ('c', np.array([3, 5], dtype='i8'))])


class TestInitFromTable(BaseInitFromDictLike):

    def setup_method(self, method):
        arr = np.array([(1, 2, 3),
                        (3, 4, 5)],
                       dtype=[('a', 'i8'), ('b', 'i8'), ('c', 'f8')])
        self.data = Table(arr, meta={'comments': ['comment1', 'comment2']})

    def test_data_meta_copy(self):
        t = Table(self.data)
        assert t.meta['comments'][0] == 'comment1'
        t['a'][1] = 8
        t.meta['comments'][1] = 'new comment2'
        assert self.data.meta['comments'][1] == 'comment2'
        assert np.all(t['a'] == np.array([1, 8]))
        assert np.all(self.data['a'] == np.array([1, 3]))
        assert t['c'].name == 'c'

    def test_table_ref(self):
        t = Table(self.data, copy=False)
        assert np.all(t._data == self.data._data)
        t['a'][1] = 0
        assert t._data['a'][1] == 0
        assert self.data._data['a'][1] == 0
        assert np.all(t._data == self.data._data)
