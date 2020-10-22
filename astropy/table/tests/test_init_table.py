# Licensed under a 3-clause BSD style license - see LICENSE.rst

from collections import OrderedDict, UserDict
from collections.abc import Mapping

import pytest
import numpy as np

from astropy.table import Column, TableColumns, Table, MaskedColumn


class TestTableColumnsInit():
    def test_init(self):
        """Test initialisation with lists, tuples, dicts of arrays
        rather than Columns [regression test for #2647]"""
        x1 = np.arange(10.)
        x2 = np.arange(5.)
        x3 = np.arange(7.)
        col_list = [('x1', x1), ('x2', x2), ('x3', x3)]
        tc_list = TableColumns(col_list)
        for col in col_list:
            assert col[0] in tc_list
            assert tc_list[col[0]] is col[1]

        col_tuple = (('x1', x1), ('x2', x2), ('x3', x3))
        tc_tuple = TableColumns(col_tuple)
        for col in col_tuple:
            assert col[0] in tc_tuple
            assert tc_tuple[col[0]] is col[1]

        col_dict = dict([('x1', x1), ('x2', x2), ('x3', x3)])
        tc_dict = TableColumns(col_dict)
        for col in tc_dict.keys():
            assert col in tc_dict
            assert tc_dict[col] is col_dict[col]

        columns = [Column(col[1], name=col[0]) for col in col_list]
        tc = TableColumns(columns)
        for col in columns:
            assert col.name in tc
            assert tc[col.name] is col


# pytest.mark.usefixtures('table_type')
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

    def test_set_dtype(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=('a', 'b', 'c'), dtype=('i4', 'f4', 'f8'))
        assert t.colnames == ['a', 'b', 'c']
        assert np.all(t['a'] == np.array([1, 3], dtype='i4'))
        assert np.all(t['b'] == np.array([2, 4], dtype='f4'))
        assert np.all(t['c'] == np.array([3, 5], dtype='f8'))
        assert t['a'].dtype.type == np.int32
        assert t['b'].dtype.type == np.float32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_names_dtype_mismatch(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=('a',), dtype=('i4', 'f4', 'i4'))

    def test_names_cols_mismatch(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=('a',), dtype=('i4'))


@pytest.mark.usefixtures('table_type')
class BaseInitFromListLike(BaseInitFrom):

    def test_names_cols_mismatch(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=['a'], dtype=[int])

    def test_names_copy_false(self, table_type):
        self._setup(table_type)
        with pytest.raises(ValueError):
            table_type(self.data, names=['a'], dtype=[int], copy=False)


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
        """Init with ndarray and copy=False and show that this is a reference
        to input ndarray"""
        self._setup(table_type)
        t = table_type(self.data, copy=False)
        t['col1'][1] = 0
        assert t.as_array()['col1'][1] == 0
        assert t['col1'][1] == 0
        assert self.data[1][1] == 0

    def test_partial_names_dtype(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['a', None, 'c'], dtype=[None, None, 'f8'])
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

    def setup_method(self, table_type):
        self._setup(table_type)
        self.data = [(np.int32(1), np.int32(3)),
                     Column(name='col1', data=[2, 4], dtype=np.int32),
                     np.array([3, 5], dtype=np.int32)]

    def test_default_names(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        assert t.colnames == ['col0', 'col1', 'col2']
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtype(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['b', None, 'c'],
                       dtype=['f4', None, 'f8'])
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

    def test_missing_data_init_from_dict(self, table_type):
        dat = [{'a': 1, 'b': 2},
               {'a': 2, 'c': 4}]
        for rows in [False, True]:
            t = table_type(rows=dat) if rows else table_type(dat)

            assert np.all(t['a'] == [1, 2])
            assert np.all(t['b'].mask == [False, True])
            assert np.all(t['b'].data == [2, 2])
            assert np.all(t['c'].mask == [True, False])
            assert np.all(t['c'].data == [4, 4])

            assert type(t['a']) is (MaskedColumn if t.masked else Column)
            assert type(t['b']) is MaskedColumn
            assert type(t['c']) is MaskedColumn


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

    def test_partial_names_dtype(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['b', None, 'c'], dtype=['f4', None, 'f8'])
        assert t.colnames == ['b', 'col1', 'c']
        assert t['b'].dtype.type == np.float32
        assert t['col1'].dtype.type == np.int32
        assert t['c'].dtype.type == np.float64
        assert all(t[name].name == name for name in t.colnames)

    def test_ref(self, table_type):
        """Test that initializing from a list of columns can be done by reference"""
        self._setup(table_type)
        t = table_type(self.data, copy=False)
        t['x'][0] = 100
        assert self.data[0][0] == 100


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

        t['x'][1] = 0  # Column-wise assignment
        t[0]['y'] = 0  # Row-wise assignment
        assert self.data['x'][1] == 0
        assert self.data['y'][0] == 0
        assert np.all(np.array(t) == self.data)
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtype(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['e', None, 'd'], dtype=['f4', None, 'f8'])
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
class TestInitFromMapping(BaseInitFromDictLike):

    def _setup(self, table_type):
        self.data = UserDict([('a', Column([1, 3], name='x')),
                              ('b', [2, 4]),
                              ('c', np.array([3, 5], dtype='i8'))])
        assert isinstance(self.data, Mapping)
        assert not isinstance(self.data, dict)


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
class TestInitFromRow(BaseInitFromDictLike):

    def _setup(self, table_type):
        arr = np.array([(1, 2, 3),
                        (3, 4, 5)],
                       dtype=[('x', 'i8'), ('y', 'i8'), ('z', 'f8')])
        self.data = table_type(arr, meta={'comments': ['comment1', 'comment2']})

    def test_init_from_row(self, table_type):
        self._setup(table_type)
        t = table_type(self.data[0])

        # Values and meta match original
        assert t.meta['comments'][0] == 'comment1'
        for name in t.colnames:
            assert np.all(t[name] == self.data[name][0:1])
        assert all(t[name].name == name for name in t.colnames)

        # Change value in new instance and check that original is the same
        t['x'][0] = 8
        t.meta['comments'][1] = 'new comment2'
        assert np.all(t['x'] == np.array([8]))
        assert np.all(self.data['x'] == np.array([1, 3]))
        assert self.data.meta['comments'][1] == 'comment2'


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
        t['x'][1] = 0
        assert t['x'][1] == 0
        assert self.data['x'][1] == 0
        assert np.all(t.as_array() == self.data.as_array())
        assert all(t[name].name == name for name in t.colnames)

    def test_partial_names_dtype(self, table_type):
        self._setup(table_type)
        t = table_type(self.data, names=['e', None, 'd'], dtype=['f4', None, 'i8'])
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
        assert t2.dtype.names == ('z', 'x', 'y')

    def test_init_from_columns_slice(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        t2 = table_type(t.columns[0:2])
        assert t2.colnames == ['x', 'y']
        assert t2.dtype.names == ('x', 'y')

    def test_init_from_columns_mix(self, table_type):
        self._setup(table_type)
        t = table_type(self.data)
        t2 = table_type([t.columns[0], t.columns['z']])
        assert t2.colnames == ['x', 'z']
        assert t2.dtype.names == ('x', 'z')


@pytest.mark.usefixtures('table_type')
class TestInitFromNone():
    # Note table_table.TestEmptyData tests initializing a completely empty
    # table and adding data.

    def test_data_none_with_cols(self, table_type):
        """
        Test different ways of initing an empty table
        """
        np_t = np.empty(0, dtype=[('a', 'f4', (2,)),
                                  ('b', 'i4')])
        for kwargs in ({'names': ('a', 'b')},
                       {'names': ('a', 'b'), 'dtype': (('f4', (2,)), 'i4')},
                       {'dtype': [('a', 'f4', (2,)), ('b', 'i4')]},
                       {'dtype': np_t.dtype}):
            t = table_type(**kwargs)
            assert t.colnames == ['a', 'b']
            assert len(t['a']) == 0
            assert len(t['b']) == 0
            if 'dtype' in kwargs:
                assert t['a'].dtype.type == np.float32
                assert t['b'].dtype.type == np.int32
                assert t['a'].shape[1:] == (2,)


@pytest.mark.usefixtures('table_types')
class TestInitFromRows():

    def test_init_with_rows(self, table_type):
        for rows in ([[1, 'a'], [2, 'b']],
                     [(1, 'a'), (2, 'b')],
                     ((1, 'a'), (2, 'b'))):
            t = table_type(rows=rows, names=('a', 'b'))
            assert np.all(t['a'] == [1, 2])
            assert np.all(t['b'] == ['a', 'b'])
            assert t.colnames == ['a', 'b']
            assert t['a'].dtype.kind == 'i'
            assert t['b'].dtype.kind in ('S', 'U')
            # Regression test for
            # https://github.com/astropy/astropy/issues/3052
            assert t['b'].dtype.str.endswith('1')

        rows = np.arange(6).reshape(2, 3)
        t = table_type(rows=rows, names=('a', 'b', 'c'), dtype=['f8', 'f4', 'i8'])
        assert np.all(t['a'] == [0, 3])
        assert np.all(t['b'] == [1, 4])
        assert np.all(t['c'] == [2, 5])
        assert t.colnames == ['a', 'b', 'c']
        assert t['a'].dtype.str.endswith('f8')
        assert t['b'].dtype.str.endswith('f4')
        assert t['c'].dtype.str.endswith('i8')

    def test_init_with_rows_and_data(self, table_type):
        with pytest.raises(ValueError) as err:
            table_type(data=[[1]], rows=[[1]])
        assert "Cannot supply both `data` and `rows` values" in str(err.value)


@pytest.mark.parametrize('has_data', [True, False])
def test_init_table_with_names_and_structured_dtype(has_data):
    """Test fix for #10393"""
    arr = np.ones(2, dtype=np.dtype([('a', 'i4'), ('b', 'f4')]))
    data_args = [arr] if has_data else []
    t = Table(*data_args, names=['x', 'y'], dtype=arr.dtype)
    assert t.colnames == ['x', 'y']
    assert str(t['x'].dtype) == 'int32'
    assert str(t['y'].dtype) == 'float32'
    assert len(t) == (2 if has_data else 0)


@pytest.mark.usefixtures('table_type')
def test_init_and_ref_from_multidim_ndarray(table_type):
    """
    Test that initializing from an ndarray structured array with
    a multi-dim column works for both copy=False and True and that
    the referencing is as expected.
    """
    for copy in (False, True):
        nd = np.array([(1, [10, 20]),
                       (3, [30, 40])],
                      dtype=[('a', 'i8'), ('b', 'i8', (2,))])
        t = table_type(nd, copy=copy)
        assert t.colnames == ['a', 'b']
        assert t['a'].shape == (2,)
        assert t['b'].shape == (2, 2)
        t['a'][0] = -200
        t['b'][1][1] = -100
        if copy:
            assert nd['a'][0] == 1
            assert nd['b'][1][1] == 40
        else:
            assert nd['a'][0] == -200
            assert nd['b'][1][1] == -100


@pytest.mark.usefixtures('table_type')
@pytest.mark.parametrize('copy', [False, True])
def test_init_and_ref_from_dict(table_type, copy):
    """
    Test that initializing from a dict works for both copy=False and True and that
    the referencing is as expected.
    """
    x1 = np.arange(10.)
    x2 = np.zeros(10)
    col_dict = dict([('x1', x1), ('x2', x2)])
    t = table_type(col_dict, copy=copy)
    assert set(t.colnames) == set(['x1', 'x2'])
    assert t['x1'].shape == (10,)
    assert t['x2'].shape == (10,)
    t['x1'][0] = -200
    t['x2'][1] = -100
    if copy:
        assert x1[0] == 0.
        assert x2[1] == 0.
    else:
        assert x1[0] == -200
        assert x2[1] == -100


def test_add_none_object_column():
    """Test fix for a problem introduced in #10636 (see
    https://github.com/astropy/astropy/pull/10636#issuecomment-676847515)
    """
    t = Table(data={'a': [1, 2, 3]})
    t['b'] = None
    assert all(val is None for val in t['b'])
    assert t['b'].dtype.kind == 'O'


@pytest.mark.usefixtures('table_type')
def test_init_from_row_OrderedDict(table_type):
    row1 = OrderedDict([('b', 1), ('a', 0)])
    row2 = {'a': 10, 'b': 20}
    rows12 = [row1, row2]
    row3 = dict([('b', 1), ('a', 0)])
    row4 = dict([('b', 11), ('a', 10)])
    rows34 = [row3, row4]
    t1 = table_type(rows=rows12)
    t2 = table_type(rows=rows34)
    t3 = t2[sorted(t2.colnames)]
    assert t1.colnames == ['b', 'a']
    assert t2.colnames == ['b', 'a']
    assert t3.colnames == ['a', 'b']


def test_init_from_rows_as_generator():
    rows = ((1 + ii, 2 + ii) for ii in range(2))
    t = Table(rows=rows)
    assert np.all(t['col0'] == [1, 2])
    assert np.all(t['col1'] == [2, 3])
