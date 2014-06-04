# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import sys

from distutils import version
import numpy as np

from ...tests.helper import pytest, catch_warnings
from ... import table
from ...table import Row
from ...utils.exceptions import AstropyDeprecationWarning
from .conftest import MaskedTable

numpy_lt_1p8 = version.LooseVersion(np.__version__) < version.LooseVersion('1.8')


def test_masked_row_with_object_col():
    """
    Numpy < 1.8 has a bug in masked array that prevents access a row if there is
    a column with object type.
    """
    t = table.Table([[1]], dtype=['O'], masked=True)
    if numpy_lt_1p8:
        with pytest.raises(ValueError):
            t['col0'].mask = False
            t[0]
        with pytest.raises(ValueError):
            t['col0'].mask = True
            t[0]
    else:
        t['col0'].mask = False
        assert t[0]['col0'] == 1
        t['col0'].mask = True
        assert t[0]['col0'] is np.ma.masked


@pytest.mark.usefixtures('table_types')
class TestRow():
    def _setup(self, table_types):
        self._table_type = table_types.Table
        self._column_type = table_types.Column

    @property
    def t(self):
        # py.test wants to run this method once before table_types is run
        # to set Table and Column.  In this case just return None, which would
        # cause any downstream test to fail if this happened in any other context.
        if self._column_type is None:
            return None
        if not hasattr(self, '_t'):
            a = self._column_type(name='a', data=[1, 2, 3], dtype='i8')
            b = self._column_type(name='b', data=[4, 5, 6], dtype='i8')
            self._t = self._table_type([a, b])
        return self._t

    def test_subclass(self, table_types):
        """Row is subclass of ndarray and Row"""
        self._setup(table_types)
        c = Row(self.t, 2)
        assert isinstance(c, Row)

    def test_values(self, table_types):
        """Row accurately reflects table values and attributes"""
        self._setup(table_types)
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
        if sys.byteorder == 'little':
            assert str(row.dtype) == "[('a', '<i8'), ('b', '<i8')]"
        else:
            assert str(row.dtype) == "[('a', '>i8'), ('b', '>i8')]"

    def test_ref(self, table_types):
        """Row is a reference into original table data"""
        self._setup(table_types)
        table = self.t
        row = table[1]
        row['a'] = 10
        if table_types.Table is not MaskedTable:
            assert table['a'][1] == 10

    def test_left_equal(self, table_types):
        """Compare a table row to the corresponding structured array row"""
        self._setup(table_types)
        np_t = self.t._data.copy()
        if table_types.Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(row == np_row)

    def test_left_not_equal(self, table_types):
        """Compare a table row to the corresponding structured array row"""
        self._setup(table_types)
        np_t = self.t._data.copy()
        np_t['a'] = [0, 0, 0]
        if table_types.Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(row != np_row)

    def test_right_equal(self, table_types):
        """Test right equal"""
        self._setup(table_types)
        np_t = self.t._data.copy()
        if table_types.Table is MaskedTable:
            with pytest.raises(ValueError):
                self.t[0] == np_t[0]
        else:
            for row, np_row in zip(self.t, np_t):
                assert np.all(np_row == row)

    def test_convert_numpy_array(self, table_types):
        self._setup(table_types)
        d = self.t[1]

        np_data = np.array(d)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        np_data = np.array(d, copy=False)
        if table_types.Table is not MaskedTable:
            assert np.all(np_data == d._data)
        assert not np_data is d._data
        assert d.colnames == list(np_data.dtype.names)

        with pytest.raises(ValueError):
            np_data = np.array(d, dtype=[(str('c'), 'i8'), (str('d'), 'i8')])

    def test_format_row(self, table_types):
        """Test formatting row"""
        self._setup(table_types)
        table = self.t
        row = table[0]
        assert format(row, "").startswith("<{0} 0 of table".format(row.__class__.__name__))
