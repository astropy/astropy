import pytest
import numpy as np
from .test_table import SetupData
from ..bst import BST, FastBST, FastRBT
from ..array import SortedArray

@pytest.fixture(params=[BST, FastBST, FastRBT, SortedArray])
def impl(request):
    return request.param

@pytest.mark.usefixtures('table_types')
class TestIndex(SetupData):
    def make_col(self, name, lst):
        if self._column_type is not None:
            return self._column_type(lst, name=name)

    @property
    def t(self):
        if self._table_type is not None:
            if not hasattr(self, '_t'):
                cols = [self.make_col('a', [1, 2, 3, 4, 5]),
                        self.make_col('b', [4.0, 5.1, 6.2, 7.0, 1.1]),
                        self.make_col('c', ['7', '8', '9', '10', '11'])]
                self._t = self._table_type(cols)
            return self._t

    @pytest.mark.parametrize("composite", [False, True])
    def test_table_index(self, table_types, composite, impl):
        self._setup(table_types)
        t = self.t
        t.add_index(('a', 'b') if composite else 'a', impl=impl)
        t['a'][0] = 4
        t.add_row((6, 6.0, '7'))
        t['a'][3] = 10
        t.remove_row(2)
        t.add_row((4, 5.0, '9'))
        assert np.all(t['a'].data == np.array([4, 2, 10, 5, 6, 4]))
        assert np.allclose(t['b'].data, np.array([4.0, 5.1, 7.0, 1.1, 6.0, 5.0]))
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

    def test_table_where(self, table_types, impl):
        ##TODO: test for expected errors
        self._setup(table_types)
        t = self.t
        t.add_index(['b', 'a'], impl=impl)
        t.add_row((2, 5.3, '7'))
        t.add_row((4, 4.0, '1'))
        '''
        a   b   c 
        --- --- ---
        1 4.0   7
        2 5.1   8
        3 6.2   9
        4 7.0  10
        5 1.1  11
        2 5.3   7
        4 4.0   1
        '''
        assert np.all(t['a'].data == np.array([1, 2, 3, 4, 5, 2, 4]))
        assert np.all(t['b'].data == np.array([4.0, 5.1, 6.2, 7.0, 1.1, 5.3, 4.0]))
        assert np.all(t['c'].data == np.array(['7', '8', '9', '10', '11', '7', '1']))
        
        # only leftmost column
        assert t.where('[b] = {0}', 4.0) == [0, 6]
        assert t.where('[b] = {0}', 5.1) == [1]
        # range query
        assert t.where('[b] in ({0}, {1})', 5.0, 5.5) == [1, 5]
        assert t.where('[b] in ({0}, {1})', 6.5, 7.0) == [] # exclusive range
        assert t.where('[b] in ({0}, {1}]', 6.5, 7.0) == [3] # inclusive range
        # both columns
        assert t.where('[a] = {0} and [b] = {1}', 4, 4.0) == [6]
        assert t.where('[b] = {1} and [a] = {0}', 1, 4.0) == [0]
        assert t.where('[a] = {0} and [b] = {1}', 4, 10.0) == []
        # range on both columns
        assert t.where('[b] = {0} and [a] in ({1}, {2})', 4.0, 3, 5) == [6]
        # query without index
        t.remove_indices('a')
        assert t.where('[b] = {0}', 4.0) == [0, 6]
        # range query without index
        assert t.where('[b] in ({0}, {1})', 5.0, 5.5) == [1, 5]

    def test_table_slicing(self, table_types, impl):
        self._setup(table_types)
        t = self.t
        t.add_index('a', impl=impl)
        assert np.all(t['a'].data == np.array([1, 2, 3, 4, 5]))
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

        t2 = t[[0, 2]]
        # t2 should retain an index on column 'a'
        assert len(t2.indices) == 1
        assert np.all(t2['a'].data == np.array([1, 3]))
        # the index in t2 should reorder row numbers after slicing
        assert np.all(t2.indices[0].sorted_data() == [0, 1])
        # however, this index should be a deep copy of t1's index
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

    def test_remove_rows(self, table_types, impl):
        self._setup(table_types)
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

    def test_col_slicing(self, table_types, impl):
        self._setup(table_types)
        t = self.t
        t.add_index('a', impl=impl)

        # get slice
        t2 = t[1:3] # table slice
        assert np.all(t2['a'].data == np.array([2, 3]))
        assert np.all(t2.indices[0].sorted_data() == [0, 1])
        col_slice = t['a'][1:3] # column slice
        assert np.all(col_slice == [2, 3])
        assert np.all(col_slice.indices[0].sorted_data() == [0, 1])

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

    def test_sort(self, table_types, impl):
        self._setup(table_types)
        t = self.t[::-1] # reverse table
        assert np.all(t['a'].data == [5, 4, 3, 2, 1])
        t.add_index('a', impl=impl)
        assert np.all(t.indices[0].sorted_data() == [4, 3, 2, 1, 0])

        # sort table by column a
        t.sort('a')
        assert np.all(t['a'].data == [1, 2, 3, 4, 5])
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2, 3, 4])

