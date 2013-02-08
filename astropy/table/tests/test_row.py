# Licensed under a 3-clause BSD style license - see LICENSE.rst
from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table
from ...table import Row

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')

# Dummy init of Table, DATA for pyflakes and to be sure test fixture is working
Table = None
Column = None


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all tests for both an unmasked (ndarray) and masked (MaskedArray) column.
@pytest.fixture(params=[False] if numpy_lt_1p5 else [False, True])
def set_global_Table(request):
    global Table
    global Column
    Table = MaskedTable if request.param else table.Table
    Column = table.MaskedColumn if request.param else table.Column


@pytest.mark.usefixtures('set_global_Table')
class TestRow():

    @property
    def t(self):
        # py.test wants to run this method once before set_global_Table is run
        # to set Table and Column.  In this case just return None, which would
        # cause any downstream test to fail if this happened in any other context.
        if Column is None:
            return None
        if not hasattr(self, '_t'):
            a = Column(name='a', data=[1, 2, 3], dtype='i8')
            b = Column(name='b', data=[4, 5, 6], dtype='i8')
            self._t = Table([a, b])
        return self._t

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
        if Table is not MaskedTable:
            assert table['a'][1] == 10

    def test_left_equal(self):
        """Compare a table row to the corresponding structured array row"""
        np_t = self.t._data.copy()
        if Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(row == np_row)

    def test_left_not_equal(self):
        """Compare a table row to the corresponding structured array row"""
        np_t = self.t._data.copy()
        np_t['a'] = [0, 0, 0]
        if Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(row != np_row)

    def test_right_equal(self):
        """Test right equal"""
        np_t = self.t._data.copy()
        if Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(np_row == row)

    def test_convert_numpy_array(self):
        d = self.t[1]

        np_data = np.array(d)
        if Table is not MaskedTable:
            assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        if Table is not MaskedTable:
            assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[('c', 'i8'), ('d', 'i8')])

    def test_format_row(self):
        """Test formatting row"""
        table = self.t
        row = table[0]
        assert format(row, "").startswith("<Row 0 of table")
