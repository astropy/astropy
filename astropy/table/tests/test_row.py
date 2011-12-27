import pytest
import numpy as np
from .. import Column, Row, Table


class TestRow():

    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])
        self.b = Column('b', [4, 5, 6])

    def test_subclass(self):
        """Row is subclass of ndarray and Row"""
        table = Table([self.a, self.b])
        c = Row(table, 2)
        assert isinstance(c, Row)

    def test_values(self):
        """Row accurately reflects table values and attributes"""
        table = Table([self.a, self.b], meta={'x': 1})
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
        table = Table([self.a, self.b])
        row = table[1]
        row['a'] = 10
        assert table['a'][1] == 10
