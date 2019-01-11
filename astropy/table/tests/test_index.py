# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from .test_table import SetupData
from astropy.table.bst import BST, FastRBT, FastBST
from astropy.table.sorted_array import SortedArray
from astropy.table.soco import SCEngine, HAS_SOCO
from astropy.table.table import QTable, Row, Table
from astropy import units as u
from astropy.time import Time
from astropy.table.column import BaseColumn

try:
    import bintrees
except ImportError:
    HAS_BINTREES = False
else:
    HAS_BINTREES = True


if HAS_BINTREES:
    available_engines = [BST, FastBST, FastRBT, SortedArray]
else:
    available_engines = [BST, SortedArray]

if HAS_SOCO:
    available_engines.append(SCEngine)


@pytest.fixture(params=available_engines)
def engine(request):
    return request.param


_col = [1, 2, 3, 4, 5]


@pytest.fixture(params=[
    _col,
    u.Quantity(_col),
    Time(_col, format='jyear'),
])
def main_col(request):
    return request.param


def assert_col_equal(col, array):
    if isinstance(col, Time):
        assert np.all(col == Time(array, format='jyear'))
    else:
        assert np.all(col == col.__class__(array))


@pytest.mark.usefixtures('table_types')
class TestIndex(SetupData):
    def _setup(self, main_col, table_types):
        super()._setup(table_types)
        self.main_col = main_col
        if isinstance(main_col, u.Quantity):
            self._table_type = QTable
        if not isinstance(main_col, list):
            self._column_type = lambda x: x  # don't change mixin type
        self.mutable = isinstance(main_col, (list, u.Quantity))

    def make_col(self, name, lst):
        return self._column_type(lst, name=name)

    def make_val(self, val):
        if isinstance(self.main_col, Time):
            return Time(val, format='jyear')
        return val

    @property
    def t(self):
        if not hasattr(self, '_t'):
            self._t = self._table_type()
            self._t['a'] = self._column_type(self.main_col)
            self._t['b'] = self._column_type([4.0, 5.1, 6.2, 7.0, 1.1])
            self._t['c'] = self._column_type(['7', '8', '9', '10', '11'])
        return self._t

    @pytest.mark.parametrize("composite", [False, True])
    def test_table_index(self, main_col, table_types, composite, engine):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index(('a', 'b') if composite else 'a', engine=engine)
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

        if not self.mutable:
            return

        # test altering table columns
        t['a'][0] = 4
        t.add_row((6, 6.0, '7'))
        t['a'][3] = 10
        t.remove_row(2)
        t.add_row((4, 5.0, '9'))

        assert_col_equal(t['a'], np.array([4, 2, 10, 5, 6, 4]))
        assert np.allclose(t['b'], np.array([4.0, 5.1, 7.0, 1.1, 6.0, 5.0]))
        assert np.all(t['c'].data == np.array(['7', '8', '10', '11', '7', '9']))
        index = t.indices[0]
        l = list(index.data.items())

        if composite:
            assert np.all(l == [((2, 5.1), [1]),
                                ((4, 4.0), [0]),
                                ((4, 5.0), [5]),
                                ((5, 1.1), [3]),
                                ((6, 6.0), [4]),
                                ((10, 7.0), [2])])
        else:
            assert np.all(l == [((2,), [1]),
                                ((4,), [0, 5]),
                                ((5,), [3]),
                                ((6,), [4]),
                                ((10,), [2])])
        t.remove_indices('a')
        assert len(t.indices) == 0

    def test_table_slicing(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', engine=engine)
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

        for slice_ in ([0, 2], np.array([0, 2])):
            t2 = t[slice_]
            # t2 should retain an index on column 'a'
            assert len(t2.indices) == 1
            assert_col_equal(t2['a'], [1, 3])

            # the index in t2 should reorder row numbers after slicing
            assert np.all(t2.indices[0].sorted_data() == [0, 1])
            # however, this index should be a deep copy of t1's index
            assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

    def test_remove_rows(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        if not self.mutable:
            return
        t = self.t
        t.add_index('a', engine=engine)

        # remove individual row
        t2 = t.copy()
        t2.remove_rows(2)
        assert_col_equal(t2['a'], [1, 2, 4, 5])
        assert np.all(t2.indices[0].sorted_data() == [0, 1, 2, 3])

        # remove by list, ndarray, or slice
        for cut in ([0, 2, 4], np.array([0, 2, 4]), slice(0, 5, 2)):
            t2 = t.copy()
            t2.remove_rows(cut)
            assert_col_equal(t2['a'], [2, 4])
            assert np.all(t2.indices[0].sorted_data() == [0, 1])

        with pytest.raises(ValueError):
            t.remove_rows((0, 2, 4))

    def test_col_get_slice(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', engine=engine)

        # get slice
        t2 = t[1:3]  # table slice
        assert_col_equal(t2['a'], [2, 3])
        assert np.all(t2.indices[0].sorted_data() == [0, 1])

        col_slice = t['a'][1:3]
        assert_col_equal(col_slice, [2, 3])
        # true column slices discard indices
        if isinstance(t['a'], BaseColumn):
            assert len(col_slice.info.indices) == 0

        # take slice of slice
        t2 = t[::2]
        assert_col_equal(t2['a'], np.array([1, 3, 5]))
        t3 = t2[::-1]
        assert_col_equal(t3['a'], np.array([5, 3, 1]))
        assert np.all(t3.indices[0].sorted_data() == [2, 1, 0])
        t3 = t2[:2]
        assert_col_equal(t3['a'], np.array([1, 3]))
        assert np.all(t3.indices[0].sorted_data() == [0, 1])
        # out-of-bound slices
        for t_empty in (t2[3:], t2[2:1], t3[2:]):
            assert len(t_empty['a']) == 0
            assert np.all(t_empty.indices[0].sorted_data() == [])

        if self.mutable:
            # get boolean mask
            mask = t['a'] % 2 == 1
            t2 = t[mask]
            assert_col_equal(t2['a'], [1, 3, 5])
            assert np.all(t2.indices[0].sorted_data() == [0, 1, 2])

    def test_col_set_slice(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        if not self.mutable:
            return
        t = self.t
        t.add_index('a', engine=engine)

        # set slice
        t2 = t.copy()
        t2['a'][1:3] = np.array([6, 7])
        assert_col_equal(t2['a'], np.array([1, 6, 7, 4, 5]))
        assert np.all(t2.indices[0].sorted_data() == [0, 3, 4, 1, 2])

        # change original table via slice reference
        t2 = t.copy()
        t3 = t2[1:3]
        assert_col_equal(t3['a'], np.array([2, 3]))
        assert np.all(t3.indices[0].sorted_data() == [0, 1])
        t3['a'][0] = 5
        assert_col_equal(t3['a'], np.array([5, 3]))
        assert_col_equal(t2['a'], np.array([1, 5, 3, 4, 5]))
        assert np.all(t3.indices[0].sorted_data() == [1, 0])
        assert np.all(t2.indices[0].sorted_data() == [0, 2, 3, 1, 4])

        # set boolean mask
        t2 = t.copy()
        mask = t['a'] % 2 == 1
        t2['a'][mask] = 0.
        assert_col_equal(t2['a'], [0, 2, 0, 4, 0])
        assert np.all(t2.indices[0].sorted_data() == [0, 2, 4, 1, 3])

    def test_multiple_slices(self, main_col, table_types, engine):
        self._setup(main_col, table_types)

        if not self.mutable:
            return

        t = self.t
        t.add_index('a', engine=engine)

        for i in range(6, 51):
            t.add_row((i, 1.0, 'A'))

        assert_col_equal(t['a'], [i for i in range(1, 51)])
        assert np.all(t.indices[0].sorted_data() == [i for i in range(50)])

        evens = t[::2]
        assert np.all(evens.indices[0].sorted_data() == [i for i in range(25)])
        reverse = evens[::-1]
        index = reverse.indices[0]
        assert (index.start, index.stop, index.step) == (48, -2, -2)
        assert np.all(index.sorted_data() == [i for i in range(24, -1, -1)])

        # modify slice of slice
        reverse[-10:] = 0
        expected = np.array([i for i in range(1, 51)])
        expected[:20][expected[:20] % 2 == 1] = 0
        assert_col_equal(t['a'], expected)
        assert_col_equal(evens['a'], expected[::2])
        assert_col_equal(reverse['a'], expected[::2][::-1])
        # first ten evens are now zero
        assert np.all(t.indices[0].sorted_data() ==
                      [0, 2, 4, 6, 8, 10, 12, 14, 16, 18,
                       1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
                      + [i for i in range(20, 50)])
        assert np.all(evens.indices[0].sorted_data() == [i for i in range(25)])
        assert np.all(reverse.indices[0].sorted_data() ==
                      [i for i in range(24, -1, -1)])

        # try different step sizes of slice
        t2 = t[1:20:2]
        assert_col_equal(t2['a'], [2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
        assert np.all(t2.indices[0].sorted_data() == [i for i in range(10)])
        t3 = t2[::3]
        assert_col_equal(t3['a'], [2, 8, 14, 20])
        assert np.all(t3.indices[0].sorted_data() == [0, 1, 2, 3])
        t4 = t3[2::-1]
        assert_col_equal(t4['a'], [14, 8, 2])
        assert np.all(t4.indices[0].sorted_data() == [2, 1, 0])

    def test_sort(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t[::-1]  # reverse table
        assert_col_equal(t['a'], [5, 4, 3, 2, 1])
        t.add_index('a', engine=engine)
        assert np.all(t.indices[0].sorted_data() == [4, 3, 2, 1, 0])

        if not self.mutable:
            return

        # sort table by column a
        t2 = t.copy()
        t2.sort('a')
        assert_col_equal(t2['a'], [1, 2, 3, 4, 5])
        assert np.all(t2.indices[0].sorted_data() == [0, 1, 2, 3, 4])

        # sort table by primary key
        t2 = t.copy()
        t2.sort()
        assert_col_equal(t2['a'], [1, 2, 3, 4, 5])
        assert np.all(t2.indices[0].sorted_data() == [0, 1, 2, 3, 4])

    def test_insert_row(self, main_col, table_types, engine):
        self._setup(main_col, table_types)

        if not self.mutable:
            return

        t = self.t
        t.add_index('a', engine=engine)
        t.insert_row(2, (6, 1.0, '12'))
        assert_col_equal(t['a'], [1, 2, 6, 3, 4, 5])
        assert np.all(t.indices[0].sorted_data() == [0, 1, 3, 4, 5, 2])
        t.insert_row(1, (0, 4.0, '13'))
        assert_col_equal(t['a'], [1, 0, 2, 6, 3, 4, 5])
        assert np.all(t.indices[0].sorted_data() == [1, 0, 2, 4, 5, 6, 3])

    def test_index_modes(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', engine=engine)

        # first, no special mode
        assert len(t[[1, 3]].indices) == 1
        assert len(t[::-1].indices) == 1
        assert len(self._table_type(t).indices) == 1
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])
        t2 = t.copy()

        # non-copy mode
        with t.index_mode('discard_on_copy'):
            assert len(t[[1, 3]].indices) == 0
            assert len(t[::-1].indices) == 0
            assert len(self._table_type(t).indices) == 0
            assert len(t2.copy().indices) == 1  # mode should only affect t

        # make sure non-copy mode is exited correctly
        assert len(t[[1, 3]].indices) == 1

        if not self.mutable:
            return

        # non-modify mode
        with t.index_mode('freeze'):
            assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])
            t['a'][0] = 6
            assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])
            t.add_row((2, 1.5, '12'))
            assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])
            t.remove_rows([1, 3])
            assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])
            assert_col_equal(t['a'], [6, 3, 5, 2])
            # mode should only affect t
            assert np.all(t2.indices[0].sorted_data() == [0, 1, 2, 3, 4])
            t2['a'][0] = 6
            assert np.all(t2.indices[0].sorted_data() == [1, 2, 3, 4, 0])

        # make sure non-modify mode is exited correctly
        assert np.all(t.indices[0].sorted_data() == [3, 1, 2, 0])

        if isinstance(t['a'], BaseColumn):
            assert len(t['a'][::-1].info.indices) == 0
            with t.index_mode('copy_on_getitem'):
                assert len(t['a'][[1, 2]].info.indices) == 1
                # mode should only affect t
                assert len(t2['a'][[1, 2]].info.indices) == 0

            assert len(t['a'][::-1].info.indices) == 0
            assert len(t2['a'][::-1].info.indices) == 0

    def test_index_retrieval(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', engine=engine)
        t.add_index(['a', 'c'], engine=engine)
        assert len(t.indices) == 2
        assert len(t.indices['a'].columns) == 1
        assert len(t.indices['a', 'c'].columns) == 2

        with pytest.raises(IndexError):
            t.indices['b']

    def test_col_rename(self, main_col, table_types, engine):
        '''
        Checks for a previous bug in which copying a Table
        with different column names raised an exception.
        '''
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', engine=engine)
        t2 = self._table_type(self.t, names=['d', 'e', 'f'])
        assert len(t2.indices) == 1

    def test_table_loc(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t

        t.add_index('a', engine=engine)
        t.add_index('b', engine=engine)

        t2 = t.loc[self.make_val(3)]  # single label, with primary key 'a'
        assert_col_equal(t2['a'], [3])
        assert isinstance(t2, Row)

        # list search
        t2 = t.loc[[self.make_val(1), self.make_val(4), self.make_val(2)]]
        assert_col_equal(t2['a'], [1, 4, 2])  # same order as input list
        if not isinstance(main_col, Time):
            # ndarray search
            t2 = t.loc[np.array([1, 4, 2])]
            assert_col_equal(t2['a'], [1, 4, 2])
        assert_col_equal(t2['a'], [1, 4, 2])
        t2 = t.loc[self.make_val(3): self.make_val(5)]  # range search
        assert_col_equal(t2['a'], [3, 4, 5])
        t2 = t.loc['b', 5.0:7.0]
        assert_col_equal(t2['b'], [5.1, 6.2, 7.0])
        # search by sorted index
        t2 = t.iloc[0:2]  # two smallest rows by column 'a'
        assert_col_equal(t2['a'], [1, 2])
        t2 = t.iloc['b', 2:]  # exclude two smallest rows in column 'b'
        assert_col_equal(t2['b'], [5.1, 6.2, 7.0])

        for t2 in (t.loc[:], t.iloc[:]):
            assert_col_equal(t2['a'], [1, 2, 3, 4, 5])

    def test_table_loc_indices(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t

        t.add_index('a', engine=engine)
        t.add_index('b', engine=engine)

        t2 = t.loc_indices[self.make_val(3)]  # single label, with primary key 'a'
        assert t2 == 2

        # list search
        t2 = t.loc_indices[[self.make_val(1), self.make_val(4), self.make_val(2)]]
        for i, p in zip(t2,[1,4,2]):  # same order as input list
            assert i == p-1

    def test_invalid_search(self, main_col, table_types, engine):
        # using .loc and .loc_indices with a value not present should raise an exception
        self._setup(main_col, table_types)
        t = self.t

        t.add_index('a')
        with pytest.raises(KeyError):
            t.loc[self.make_val(6)]
        with pytest.raises(KeyError):
            t.loc_indices[self.make_val(6)]


    def test_copy_index_references(self, main_col, table_types, engine):
        # check against a bug in which indices were given an incorrect
        # column reference when copied
        self._setup(main_col, table_types)
        t = self.t

        t.add_index('a')
        t.add_index('b')
        t2 = t.copy()
        assert t2.indices['a'].columns[0] is t2['a']
        assert t2.indices['b'].columns[0] is t2['b']

    def test_unique_index(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = self.t

        t.add_index('a', engine=engine, unique=True)
        assert np.all(t.indices['a'].sorted_data() == [0, 1, 2, 3, 4])

        if self.mutable:
            with pytest.raises(ValueError):
                t.add_row((5, 5.0, '9'))

    def test_copy_indexed_table(self, table_types):
        self._setup(_col, table_types)
        t = self.t
        t.add_index('a')
        t.add_index(['a', 'b'])
        for tp in (self._table_type(t), t.copy()):
            assert len(t.indices) == len(tp.indices)
            for index, indexp in zip(t.indices, tp.indices):
                assert np.all(index.data.data == indexp.data.data)
                assert index.data.data.colnames == indexp.data.data.colnames

    def test_updating_row_byindex(self, main_col, table_types, engine):
        self._setup(main_col, table_types)
        t = Table([['a', 'b', 'c', 'd'], [2, 3, 4, 5], [3, 4, 5, 6]], names=('a', 'b', 'c'), meta={'name': 'first table'})

        t.add_index('a', engine=engine)
        t.add_index('b', engine=engine)

        t.loc['c'] = ['g', 40, 50] # single label, with primary key 'a'
        t2 = t[2]
        assert list(t2) == ['g', 40, 50]

        # list search
        t.loc[['a', 'd', 'b']] = [['a', 20, 30], ['d', 50, 60], ['b', 30, 40]]
        t2 = [['a', 20, 30], ['d', 50, 60], ['b', 30, 40]]
        for i, p in zip(t2, [1, 4, 2]):  # same order as input list
            assert list(t[p-1]) == i

    def test_invalid_updates(self, main_col, table_types, engine):
        # using .loc and .loc_indices with a value not present should raise an exception
        self._setup(main_col, table_types)
        t = Table([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]], names=('a', 'b', 'c'), meta={'name': 'first table'})

        t.add_index('a')
        with pytest.raises(ValueError):
            t.loc[3] = [[1,2,3]]
        with pytest.raises(ValueError):
            t.loc[[1, 4, 2]] = [[1, 2, 3], [4, 5, 6]]
        with pytest.raises(ValueError):
            t.loc[[1, 4, 2]] = [[1, 2, 3], [4, 5, 6], [2, 3]]
        with pytest.raises(ValueError):
            t.loc[[1, 4, 2]] = [[1, 2, 3], [4, 5], [2, 3]]
