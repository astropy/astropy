import pytest
import numpy as np
from .test_table import SetupData

@pytest.mark.usefixtures('table_types')
class TestIndex(SetupData):
    @property
    def b(self):
        if self._column_type is not None:
            if not hasattr(self, '_b'):
                self._b = self._column_type(name='b', data=[4.0, 5.1, 6.2])
            return self._b

    @property
    def c(self):
        if self._column_type is not None:
            if not hasattr(self, '_c'):
                self._c = self._column_type(name='c', data=['7', '8', '9'])
            return self._c

    @property
    def t(self):
        if self._table_type is not None:
            if not hasattr(self, '_t'):
                self._t = self._table_type([self.a, self.b, self.c])
            return self._t

    @pytest.mark.parametrize("composite", [False, True])
    def test_table_index(self, table_types, composite):
        self._setup(table_types)
        t = self.t
        t.add_index(('a', 'b') if composite else 'a')
        t['a'][0] = 4
        t.add_row((5, 6.0, '7'))
        t['a'][3] = 10
        t.remove_row(2)
        t.add_row((4, 5.0, '9'))
        assert np.all(t['a'].data == np.array([4, 2, 10, 4]))
        assert np.allclose(t['b'].data, np.array([4.0, 5.1, 6.0, 5.0]))
        assert np.all(t['c'].data == np.array(['7', '8', '7', '9']))
        index = t.indices[0]
        l = [(x.key, x.data) for x in index.data.traverse()]

        if composite:
            assert np.all(l == [((2, 5.1), [1]),
                                ((4, 4.0), [0]),
                                ((4, 5.0), [3]),
                                ((10, 6.0), [2])])
        else:
            assert np.all(l == [((2,), [1]),
                                ((4,), [0, 3]),
                                ((10,), [2])])
        t.remove_indices('a')
        assert len(t.indices) == 0

    def test_table_where(self, table_types):
        ##TODO: test for expected errors
        self._setup(table_types)
        t = self.t
        t.add_index(['b', 'a'])
        t.add_row((2, 5.3, '7'))
        t.add_row((4, 4.0, '1'))
        assert np.all(t['a'].data == np.array([1, 2, 3, 2, 4]))
        assert np.all(t['b'].data == np.array([4.0, 5.1, 6.2, 5.3, 4.0]))
        assert np.all(t['c'].data == np.array(['7', '8', '9', '7', '1']))
        
        # only leftmost column
        assert t.where('[b] = {0}', 4.0) == [0, 4]
        assert t.where('[b] = {0}', 5.1) == [1]
        # range query
        assert t.where('[b] in ({0}, {1})', 5.0, 5.5) == [1, 3]
        assert t.where('[b] in ({0}, {1})', 6.5, 7.0) == []
        # both columns
        assert t.where('[a] = {0} and [b] = {1}', 4, 4.0) == [4]
        assert t.where('[b] = {1} and [a] = {0}', 1, 4.0) == [0]
        # range on both columns
        assert t.where('[b] = {0} and [a] in ({1}, {2})', 4.0, 3, 5) == [4]

    def test_slicing(self, table_types):
        self._setup(table_types)
        t = self.t
        t.add_index('a')
        assert np.all(t['a'].data == np.array([1, 2, 3]))
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2])

        t2 = t[[0, 2]]
        # t2 should retain an index on column 'a'
        assert len(t2.indices) == 1
        assert np.all(t2['a'].data == np.array([1, 3]))
        # the index in t2 should reorder row numbers after slicing
        assert np.all(t2.indices[0].sorted_data() == [0, 1])
        # however, this index should be a deep copy of t1's index
        assert np.all(t.indices[0].sorted_data() == [0, 1, 2])

    def test_remove_rows(self, table_types):
        self._setup(table_types)
        t = self.t
        t.add_index('a')
        t.add_row((4, 9.0, 'A'))
        t.add_row((5, 1.4, 'B'))
        
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
