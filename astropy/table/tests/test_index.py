import pytest
import numpy as np
from .test_table import SetupData
from ..bst import BST, FastBST, FastRBT
from ..sorted_array import SortedArray
from ..table import Table, QTable, NdarrayMixin
from ...table import table_helpers
from ..index import index_mode
from ... import units as u
from ...coordinates import SkyCoord

@pytest.fixture(params=[BST, FastBST, FastRBT, SortedArray])
def impl(request):
    return request.param

@pytest.fixture(params=[
    [1, 2, 3, 4, 5],
    u.Quantity([1, 2, 3, 4, 5]),
    SkyCoord(ra=[1, 2, 3, 4, 5] * u.deg,
             dec=[1, 2, 3, 4, 5] * u.deg)])
def main_col(request):
    return request.param

def assert_col_equal(col, array):
    assert np.all(col == col.__class__(array))

@pytest.mark.usefixtures('table_types')
class TestIndex(SetupData):
    def _setup(self, main_col, table_types):
        super(TestIndex, self)._setup(table_types)
        self.main_col = main_col
        if isinstance(main_col, u.Quantity):
            self._table_type = QTable
        if not isinstance(main_col, list):
            self._column_type = lambda x: x

    def make_col(self, name, lst):
        return self._column_type(lst, name=name)

    @property
    def t(self):
        if not hasattr(self, '_t'):
            self._t = self._table_type()
            self._t['a'] = self._column_type(self.main_col)
            self._t['b'] = self._column_type([4.0, 5.1, 6.2, 7.0, 1.1])
            self._t['c'] = self._column_type(['7', '8', '9', '10', '11'])
        return self._t

    @pytest.mark.parametrize("composite", [False, True])
    def test_table_index(self, main_col, table_types, composite, impl):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index(('a', 'b') if composite else 'a', impl=impl)
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
            print(l)
            assert np.all(l == [((2,), [1]),
                                ((4,), [0, 5]),
                                ((5,), [3]),
                                ((6,), [4]),
                                ((10,), [2])])
        t.remove_indices('a')
        assert len(t.indices) == 0

    def test_table_slicing(self, main_col, table_types, impl):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', impl=impl)
        assert_col_equal(t['a'],  np.array([1, 2, 3, 4, 5]))
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

        t2 = t[[0, 2]]
        # t2 should retain an index on column 'a'
        assert len(t2.indices) == 1
        assert np.all(t2['a'].data == np.array([1, 3]))
        # the index in t2 should reorder row numbers after slicing
        assert np.all(t2.indices[0].sorted_data() == [0, 1])
        # however, this index should be a deep copy of t1's index
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

    def test_remove_rows(self, main_col, table_types, impl):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', impl=impl)
        
        # remove individual row
        t2 = t.copy()
        t2.remove_rows(2)
        assert np.all(t2['a'].data == np.array([1, 2, 4, 5]))
        assert np.all(t2.indices[0].sorted_data() == [0, 1, 2, 3])

        # remove by list, ndarray, or slice
        for cut in ([0, 2, 4], np.array([0, 2, 4]), slice(0, 5, 2)):
            t2 = t.copy()
            t2.remove_rows(cut)
            assert np.all(t2['a'].data == np.array([2, 4]))
            assert np.all(t2.indices[0].sorted_data() == [0, 1])

        with pytest.raises(ValueError):
            t.remove_rows((0, 2, 4))

    def test_col_slicing(self, main_col, table_types, impl):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', impl=impl)

        # get slice
        t2 = t[1:3] # table slice
        assert np.all(t2['a'].data == np.array([2, 3]))
        assert np.all(t2.indices[0].sorted_data() == [0, 1])

        col_slice = t['a'][1:3] # column slice should discard indices
        assert np.all(col_slice.data == np.array([2, 3]))
        assert np.all(col_slice.indices == [])

        # take slice of slice
        t2 = t[::2]
        assert np.all(t2['a'].data == np.array([1, 3, 5]))
        t3 = t2[::-1]
        assert np.all(t3['a'].data == np.array([5, 3, 1]))
        assert np.all(t3.indices[0].sorted_data() == [2, 1, 0])
        t3 = t2[:2]
        assert np.all(t3['a'].data == np.array([1, 3]))
        assert np.all(t3.indices[0].sorted_data() == [0, 1])
        # out-of-bound slices
        for t_empty in (t2[3:], t2[2:1], t3[2:]):
            assert len(t_empty['a'].data) == 0
            assert np.all(t_empty.indices[0].sorted_data() == [])

        # set slice
        t2 = t.copy()
        t2['a'][1:3] = np.array([6, 7])
        assert np.all(t2['a'].data == np.array([1, 6, 7, 4, 5]))
        assert np.all(t2.indices[0].sorted_data() == [0, 3, 4, 1, 2])

        # change original table via slice reference
        t2 = t.copy()
        t3 = t2[1:3]
        assert np.all(t3['a'].data == np.array([2, 3]))
        assert np.all(t3.indices[0].sorted_data() == [0, 1])
        t3['a'][0] = 5
        assert np.all(t3['a'].data == np.array([5, 3]))
        assert np.all(t2['a'].data == np.array([1, 5, 3, 4, 5]))
        assert np.all(t3.indices[0].sorted_data() == [1, 0])
        assert np.all(t2.indices[0].sorted_data() == [0, 2, 3, 1, 4])

        # get boolean mask
        mask = t['a'] % 2 == 1
        t2 = t[mask]
        assert np.all(t2['a'].data == [1, 3, 5])
        assert np.all(t2.indices[0].sorted_data() == [0, 1, 2])

        # set boolean mask
        t2 = t.copy()
        t2['a'][mask] = 0.
        assert np.all(t2['a'].data == [0, 2, 0, 4, 0])
        assert np.all(t2.indices[0].sorted_data() == [0, 2, 4, 1, 3])

    def test_sort(self, main_col, table_types, impl):
        self._setup(main_col, table_types)
        t = self.t[::-1] # reverse table
        assert_col_equal(t['a'], [5, 4, 3, 2, 1])
        t.add_index('a', impl=impl)
        assert np.all(t.indices[0].sorted_data() == [4, 3, 2, 1, 0])

        # sort table by column a
        t.sort('a')
        assert_col_equal(t['a'], [1, 2, 3, 4, 5])
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

    def test_insert_row(self, main_col, table_types, impl):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', impl=impl)
        t.insert_row(2, (6, 1.0, '12'))
        assert_col_equal(t['a'], [1, 2, 6, 3, 4, 5])
        assert np.all(t.indices[0].sorted_data() == [0, 1, 3, 4, 5, 2])
        t.insert_row(1, (0, 4.0, '13'))
        assert_col_equal(t['a'], [1, 0, 2, 6, 3, 4, 5])
        assert np.all(t.indices[0].sorted_data() == [1, 0, 2, 4, 5, 6, 3])

    def test_index_modes(self, main_col, table_types, impl):
        self._setup(main_col, table_types)
        t = self.t
        t.add_index('a', impl=impl)

        # first, no special mode
        assert len(t[[1, 3]].indices) == 1
        assert len(t[::-1].indices) == 1
        assert len(Table(t).indices) == 1
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])
        t['a'][0] = 6
        assert_col_equal(t['a'], [6, 2, 3, 4, 5])
        assert np.all(t.indices[0].sorted_data() == [1, 2, 3, 4, 0])
        t2 = t.copy()

        # non-copy mode
        with index_mode(t, 'discard_on_copy'):
            assert len(t[[1, 3]].indices) == 0
            assert len(t[::-1].indices) == 0
            assert len(Table(t).indices) == 0
            assert len(t2.copy().indices) == 1 # mode should only affect t

        # make sure non-copy mode is exited correctly
        assert len(t[[1, 3]].indices) == 1

        # non-modify mode
        with index_mode(t, 'freeze'):
            assert np.all(t.indices[0].sorted_data() == [1, 2, 3, 4, 0])
            t['a'][0] = 1
            assert np.all(t.indices[0].sorted_data() == [1, 2, 3, 4, 0])
            t.add_row((2, 1.5, '12'))
            assert np.all(t.indices[0].sorted_data() == [1, 2, 3, 4, 0])
            t.remove_rows([1, 3])
            assert np.all(t.indices[0].sorted_data() == [1, 2, 3, 4, 0])
            assert_col_equal(t['a'], [1, 3, 5, 2])
            # mode should only affect t
            assert np.all(t2.indices[0].sorted_data() == [1, 2, 3, 4, 0])
            t2['a'][0] = 1
            assert np.all(t2.indices[0].sorted_data() == [0, 1, 2, 3, 4])

        # make sure non-modify mode is exited correctly
        assert np.all(t.indices[0].sorted_data() == [0, 3, 1, 2])

        assert len(t['a'][::-1].indices) == 0
        with index_mode(t, 'copy_on_getitem'):
            assert len(t['a'][[1, 2]].indices) == 1
            # mode affects t2 as well as t
            assert len(t2['a'][[1, 2]].indices) == 1

        assert len(t['a'][::-1].indices) == 0
        assert len(t2['a'][::-1].indices) == 0
