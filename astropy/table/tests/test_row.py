import pytest
import numpy as np
from .. import Column, Row, Table


class TestRow():

    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3], dtype=np.int64)
        self.b = Column('b', [4, 5, 6], dtype=np.int64)
        self.t = Table([self.a, self.b])

    def test_subclass(self):
        """Row is subclass of ndarray and Row"""
        c = Row(self.t, 2)
        assert isinstance(c, Row)

    def test_values(self):
        """Row accurately reflects table values and attributes"""
        table = self.t
        row = table[1]
        assert row['a'] == 2
        assert row['b'] == 5
        assert row[0] == 2
        assert row[1] == 5
        assert row.meta is table.meta
        assert row.colnames == table.colnames
        assert row.columns is table.columns
        with pytest.raises(IndexError):
            row[2]
        assert str(row.dtype) == "[('a', '<i8'), ('b', '<i8')]"

    def test_ref(self):
        """Row is a reference into original table data"""
        table = self.t
        row = table[1]
        row['a'] = 10
        assert table['a'][1] == 10

    def test_left_equal(self):
        """Compare a table row to the corresponding structured array row"""
        np_t = self.t._data.copy()
        for row, np_row in zip(self.t, np_t):
            assert np.all(row == np_row)

    def test_left_not_equal(self):
        """Compare a table row to the corresponding structured array row"""
        np_t = self.t._data.copy()
        np_t['a'] = [0, 0, 0]
        for row, np_row in zip(self.t, np_t):
            assert np.all(row != np_row)

    def test_right_equal(self):
        """Test right equal"""
        np_t = self.t._data.copy()
        for row, np_row in zip(self.t, np_t):
            assert np.all(np_row == row)

    @pytest.mark.xfail
    def test_set_slice(self):
        """Set row elements with a slice

        This currently fails because the underlying np.void object
        row.data = table._data[index] does not support slice assignment.
        """
        table = self.t
        row = table[0]
        row[:] = [-1, -1]
        row[:1] = np.array([-2])
        assert np.all(table._data == np.array([[-1, -1],
                                               [-2, 5],
                                               [3, 6]]))

    def test_convert_numpy_array(self):
        d = self.t[1]

        np_data = np.array(d)
        assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[('c', 'i8'), ('d', 'i8')])
