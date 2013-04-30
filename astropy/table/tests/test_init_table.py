# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function  # For print debugging with python 2 or 3

from distutils import version

import numpy as np

from ...tests.helper import pytest
from ... import table
from ...table import Column
from ...utils import OrderedDict

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False] if numpy_lt_1p5 else [False, True])
def table_type(request):
    return MaskedTable if request.param else table.Table


@pytest.mark.usefixtures('table_type')
class BaseInitFrom():
    def _setup(self, table_type):
        pass

    def test_basic_init(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=('a', 'b', 'c'))
        assert t.colnames == ['a', 'b', 'c']
        assert np.all(t['a'] == np.array([1, 3]))
        assert np.all(t['b'] == np.array([2, 4]))
        assert np.all(t['c'] == np.array([3, 5]))
        assert all(t[name].name == name for name in t.colnames)

    def test_set_dtypes(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=('a', 'b', 'c'), dtypes=('i4', 'f4', 'f8'))
        assert t.colnames == ['a', 'b', 'c']
        assert np.all(t['a'] == np.array([1, 3], dtype='i4'))
        assert np.all(t['b'] == np.array([2, 4], dtype='f4'))
        assert np.all(t['c'] == np.array([3, 5], dtype='f8'))
        assert t['a'].dtype.type == np.int32
        assert t['b'].dtype.type == np.float32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_names_dtypes_mismatch(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=('a',), dtypes=('i4', 'f4', 'i4'))

    def test_names_cols_mismatch(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=('a',), dtypes=('i4'))


@pytest.mark.usefixtures('table_type')
class BaseInitFromListLike(BaseInitFrom):

    def test_names_cols_mismatch(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=['a'], dtypes=[int])

    def test_names_copy_false(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=['a'], dtypes=[int], copy=False)


@pytest.mark.usefixtures('table_type')
class BaseInitFromDictLike(BaseInitFrom):
    pass


@pytest.mark.usefixtures('table_type')
class TestInitFromNdarrayHomo(BaseInitFromListLike):

    def setup_method(self, method):
        self.data = np.array([(1, 2, 3),
                              (3, 4, 5)],
                             dtype='i4')

    def test_default_names(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        assert t.colnames == ['col0', 'col1', 'col2']

    def test_ndarray_ref(self, table_type):
        """Init with ndarray and copy=False and show that ValueError is raised
        to input ndarray"""
        self._setup(table_type)
        t = table_type(self.data, copy=False)
        t['col1'][1] = 0
        assert t._data['col1'][1] == 0
        assert t['col1'][1] == 0
        assert self.data[1][1] == 0
        # NOTE: assert np.all(t._data == self.data) fails because when
        # homogenous array is viewcast to structured then the == is False

    def test_partial_names_dtypes(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['a', None, 'c'], dtypes=[None, None, 'f8'])
        assert t.colnames == ['a', 'col1', 'c']
        assert t['a'].dtype.type == np.int32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_ref(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['a', None, 'c'])
        assert t.colnames == ['a', 'col1', 'c']
        assert t['a'].dtype.type == np.int32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.int32
        assert all(t[name].name == name for name in t.colnames)


@pytest.mark.usefixtures('table_type')
class TestInitFromListOfLists(BaseInitFromListLike):

    def setup_method(self, method):
        self._setup(table_type)
        self.data = [(np.int32(1), np.int32(3)),
                     Column(name='col1', data=[2, 4], dtype=np.int32),
                     np.array([3, 5], dtype=np.int32)]

    def test_default_names(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        assert t.colnames == ['col0', 'col1', 'col2']
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['b', None, 'c'],
                  dtypes=['f4', None, 'f8'])
        assert t.colnames == ['b', 'col1', 'c']
        assert t['b'].dtype.type == np.float32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_bad_data(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type([[1, 2],
                   [3, 4, 5]])


@pytest.mark.usefixtures('table_type')
class TestInitFromListOfDicts(BaseInitFromListLike):

    def _setup(self, table_type):
        self.data = [{'a': 1, 'b': 2, 'c': 3},
                     {'a': 3, 'b': 4, 'c': 5}]

    def test_names(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        assert all(colname in set(['a', 'b', 'c']) for colname in t.colnames)

    def test_names_ordered(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=('c', 'b', 'a'))
        assert t.colnames == ['c', 'b', 'a']

    def test_bad_data(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type([{'a': 1, 'b': 2, 'c': 3},
                   {'a': 2, 'b': 4}])


@pytest.mark.usefixtures('table_type')
class TestInitFromColsList(BaseInitFromListLike):

    def _setup(self, table_type):
        self.data = [Column([1, 3], name='x', dtype=np.int32),
                     np.array([2, 4], dtype=np.int32),
                     np.array([3, 5], dtype='i8')]

    def test_default_names(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        assert t.colnames == ['x', 'col1', 'col2']
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['b', None, 'c'], dtypes=['f4', None, 'f8'])
        assert t.colnames == ['b', 'col1', 'c']
        assert t['b'].dtype.type == np.float32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_ref(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, copy=False)


@pytest.mark.usefixtures('table_type')
class TestInitFromNdarrayStruct(BaseInitFromDictLike):

    def _setup(self, table_type):
        self.data = np.array([(1, 2, 3),
                              (3, 4, 5)],
                             dtype=[('x', 'i8'), ('y', 'i4'), ('z', 'i8')])

    def test_ndarray_ref(self, table_type):
        """Init with ndarray and copy=False and show that table uses reference
        to input ndarray"""
        self._setup(table_type)
        t = table_type(self.data, copy=False)
        assert np.all(t._data == self.data)
        t['x'][1] = 0
        assert t._data['x'][1] == 0
        assert self.data['x'][1] == 0
        assert np.all(t._data == self.data)
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['e', None, 'd'], dtypes=['f4', None, 'f8'])
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.float32
        assert t['y'].dtype.type == np.int32
        assert t['d'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_ref(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['e', None, 'd'], copy=False)
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.int64
        assert t['y'].dtype.type == np.int32
        assert t['d'].dtype.type == np.int64
        assert all(t[name].name == name for name in t.colnames)


@pytest.mark.usefixtures('table_type')
class TestInitFromDict(BaseInitFromDictLike):

    def _setup(self, table_type):
        self.data = dict([('a', Column([1, 3], name='x')),
                          ('b', [2, 4]),
                          ('c', np.array([3, 5], dtype='i8'))])


@pytest.mark.usefixtures('table_type')
class TestInitFromOrderedDict(BaseInitFromDictLike):

    def _setup(self, table_type):
        self.data = OrderedDict([('a', Column(name='x', data=[1, 3])),
                                 ('b', [2, 4]),
                                 ('c', np.array([3, 5], dtype='i8'))])

    def test_col_order(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        assert t.colnames == ['a', 'b', 'c']


@pytest.mark.usefixtures('table_type')
class TestInitFromTable(BaseInitFromDictLike):

    def _setup(self, table_type):
        arr = np.array([(1, 2, 3),
                        (3, 4, 5)],
                       dtype=[('x', 'i8'), ('y', 'i8'), ('z', 'f8')])
        self.data = table_type(arr, meta={'comments': ['comment1', 'comment2']})

    def test_data_meta_copy(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        assert t.meta['comments'][0] == 'comment1'
        t['x'][1] = 8
        t.meta['comments'][1] = 'new comment2'
        assert self.data.meta['comments'][1] == 'comment2'
        assert np.all(t['x'] == np.array([1, 8]))
        assert np.all(self.data['x'] == np.array([1, 3]))
        assert t['z'].name == 'z'
        assert all(t[name].name == name for name in t.colnames)

    def test_table_ref(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, copy=False)
        assert np.all(t._data == self.data._data)
        t['x'][1] = 0
        assert t._data['x'][1] == 0
        assert self.data._data['x'][1] == 0
        assert np.all(t._data == self.data._data)
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtypes(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['e', None, 'd'], dtypes=['f4', None, 'i8'])
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.float32
        assert t['y'].dtype.type == np.int64
        assert t['d'].dtype.type == np.int64
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_ref(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['e', None, 'd'], copy=False)
        assert t.colnames == ['e', 'y', 'd']
        assert t['e'].dtype.type == np.int64
        assert t['y'].dtype.type == np.int64
        assert t['d'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_init_from_columns(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        t2 = table_type(t.columns['z', 'x', 'y'])
        assert t2.colnames == ['z', 'x', 'y']
        assert t2._data.dtype.names == ('z', 'x', 'y')

    def test_init_from_columns_slice(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        t2 = table_type(t.columns[0:2])
        assert t2.colnames == ['x', 'y']
        assert t2._data.dtype.names == ('x', 'y')

    def test_init_from_columns_mix(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        t2 = table_type([t.columns[0], t.columns['z']])
        assert t2.colnames == ['x', 'z']
        assert t2._data.dtype.names == ('x', 'z')


@pytest.mark.usefixtures('table_type')
class TestInitFromNone():
    # Note table_table.TestEmptyData tests initializing a completely empty
    # table and adding data.

    def test_data_none_with_cols(self, table_type):
        t = table_type(names=('a', 'b'))
        assert len(t['a']) == 0
        assert len(t['b']) == 0
        assert t.colnames == ['a', 'b']
        t = table_type(names=('a', 'b'), dtypes=('f4', 'i4'))
        assert t['a'].dtype.type == np.float32
        assert t['b'].dtype.type == np.int32
        assert t.colnames == ['a', 'b']
