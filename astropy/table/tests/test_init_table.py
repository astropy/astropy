from __future__ import print_function  # For print debugging with python 2 or 3

import pytest
import numpy as np

from .. import Table, Column
from astropy.utils import OrderedDict


def names_match_cols(t):
    return 


class BaseInitFrom():

    def test_basic_init(self):
        t = Table(self.data, names=('a', 'b', 'c'))
        assert t.colnames == ['a', 'b', 'c']
        assert np.all(t['a'] == np.array([1, 3]))
        assert np.all(t['b'] == np.array([2, 4]))
        assert np.all(t['c'] == np.array([3, 5]))
        assert all(t[name].name == name for name in t.colnames)

    def test_set_dtypes(self):
        t = Table(self.data, names=('a', 'b', 'c'), dtypes=('i4', 'f4', 'f8'))
        assert t.colnames == ['a', 'b', 'c']
        assert np.all(t['a'] == np.array([1, 3], dtype='i4'))
        assert np.all(t['b'] == np.array([2, 4], dtype='f4'))
        assert np.all(t['c'] == np.array([3, 5], dtype='f8'))
        assert t['a'].dtype.type == np.int32
        assert t['b'].dtype.type == np.float32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_names_dtypes_mismatch(self):
        with pytest.raises(ValueError):
            Table(self.data, names=('a',), dtypes=('i4', 'f4', 'i4'))

    def test_names_cols_mismatch(self):
        with pytest.raises(ValueError):
            Table(self.data, names=('a',), dtypes=('i4'))


class BaseInitFromListLike(BaseInitFrom):

    def test_names_cols_mismatch(self):
        with pytest.raises(ValueError):
            Table(self.data, names=['a'], dtypes=[int])

    def test_names_copy_false(self):
        with pytest.raises(ValueError):
            Table(self.data, names=['a'], dtypes=[int], copy=False)


class BaseInitFromDictLike(BaseInitFrom):
    pass


class TestInitFromNdarrayHomo(BaseInitFromListLike):

    def setup_method(self, method):
        self.data = np.array([(1, 2, 3),
                              (3, 4, 5)],
                             dtype='i4')

    def test_default_names(self):
        t = Table(self.data)
        assert t.colnames == ['col0', 'col1', 'col2']

    def test_ndarray_ref(self):
        """Init with ndarray and copy=False and show that ValueError is raised
        to input ndarray"""
        t = Table(self.data, copy=False)
        t['col1'][1] = 0
        assert t._data['col1'][1] == 0
        assert t['col1'][1] == 0
        assert self.data[1][1] == 0
        # NOTE: assert np.all(t._data == self.data) fails because when
        # homogenous array is viewcast to structured then the == is False

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['a', None, 'c'], dtypes=[None, None, 'f8'])
        assert t.colnames == ['a', 'col1', 'c']
        assert t['a'].dtype.type == np.int32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_ref(self):
        t = Table(self.data, names=['a', None, 'c'])
        assert t.colnames == ['a', 'col1', 'c']
        assert t['a'].dtype.type == np.int32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.int32
        assert all(t[name].name == name for name in t.colnames)


class TestInitFromListOfLists(BaseInitFromListLike):

    def setup_method(self, method):
        self.data = [(1, 3),
                     Column('col1', [2, 4]),
                     np.array([3, 5])]

    def test_default_names(self):
        t = Table(self.data)
        assert t.colnames == ['col0', 'col1', 'col2']
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['b', None, 'c'],
                  dtypes=['f4', None, 'f8'])
        assert t.colnames == ['b', 'col1', 'c']
        assert t['b'].dtype.type == np.float32
        assert t['col1'].dtype.type == np.int64
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_bad_data(self):
        with pytest.raises(ValueError):
            Table([[1, 2],
                   [3, 4, 5]])


class TestInitFromColsList(BaseInitFromListLike):

    def setup_method(self, method):
        self.data = [Column('x', [1, 3]),
                     np.array([2, 4]),
                     np.array([3, 5], dtype='i8')]

    def test_default_names(self):
        t = Table(self.data)
        assert t.colnames == ['x', 'col1', 'col2']
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['b', None, 'c'], dtypes=['f4', None, 'f8'])
        assert t.colnames == ['b', 'col1', 'c']
        assert t['b'].dtype.type == np.float32
        assert t['col1'].dtype.type == np.int64
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_ref(self):
        with pytest.raises(ValueError):
            Table(self.data, copy=False)


class TestInitFromNdarrayStruct(BaseInitFromDictLike):

    def setup_method(self, method):
        self.data = np.array([(1, 2, 3),
                              (3, 4, 5)],
                             dtype=[('x', 'i8'), ('y', 'i4'), ('z', 'i8')])

    def test_ndarray_ref(self):
        """Init with ndarray and copy=False and show that table uses reference
        to input ndarray"""
        t = Table(self.data, copy=False)
        assert np.all(t._data == self.data)
        t['x'][1] = 0
        assert t._data['x'][1] == 0
        assert self.data['x'][1] == 0
        assert np.all(t._data == self.data)
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['e', None, 'd'], dtypes=['f4', None, 'f8'])
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.float32
        assert t['y'].dtype.type == np.int32
        assert t['d'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_ref(self):
        t = Table(self.data, names=['e', None, 'd'], copy=False)
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.int64
        assert t['y'].dtype.type == np.int32
        assert t['d'].dtype.type == np.int64
        assert all(t[name].name == name for name in t.colnames)


class TestInitFromDict(BaseInitFromDictLike):

    def setup_method(self, method):
        self.data = dict([('a', Column('x', [1, 3])),
                          ('b', [2, 4]),
                          ('c', np.array([3, 5], dtype='i8'))])


class TestInitFromOrderedDict(BaseInitFromDictLike):

    def setup_method(self, method):
        self.data = OrderedDict([('a', Column('x', [1, 3])),
                                 ('b', [2, 4]),
                                 ('c', np.array([3, 5], dtype='i8'))])

    def test_col_order(self):
        t = Table(self.data)
        assert t.colnames == ['a', 'b', 'c']


class TestInitFromTable(BaseInitFromDictLike):

    def setup_method(self, method):
        arr = np.array([(1, 2, 3),
                        (3, 4, 5)],
                       dtype=[('x', 'i8'), ('y', 'i8'), ('z', 'f8')])
        self.data = Table(arr, meta={'comments': ['comment1', 'comment2']})

    def test_data_meta_copy(self):
        t = Table(self.data)
        assert t.meta['comments'][0] == 'comment1'
        t['x'][1] = 8
        t.meta['comments'][1] = 'new comment2'
        assert self.data.meta['comments'][1] == 'comment2'
        assert np.all(t['x'] == np.array([1, 8]))
        assert np.all(self.data['x'] == np.array([1, 3]))
        assert t['z'].name == 'z'
        assert all(t[name].name == name for name in t.colnames)

    def test_table_ref(self):
        t = Table(self.data, copy=False)
        assert np.all(t._data == self.data._data)
        t['x'][1] = 0
        assert t._data['x'][1] == 0
        assert self.data._data['x'][1] == 0
        assert np.all(t._data == self.data._data)
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self):
        t = Table(self.data, names=['e', None, 'd'], dtypes=['f4', None, 'i8'])
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.float32
        assert t['y'].dtype.type == np.int64
        assert t['d'].dtype.type == np.int64
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_ref(self):
        t = Table(self.data, names=['e', None, 'd'], copy=False)
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.int64
        assert t['y'].dtype.type == np.int64
        assert t['d'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)
