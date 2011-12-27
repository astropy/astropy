import pytest
import numpy as np
from .. import Column, Row, Table


class TestRow():

    def setup_method(self, method):
        self.a = Column('a', [1, 2, 3])
        self.b = Column('b', [4, 5, 6])

    def test_subclass(self):
        table = Table([self.a, self.b])
        c = Row(table, 2)
        assert isinstance(c, Row)
        assert isinstance(c, np.ndarray)

    def test_values(self):
        table = Table([self.a, self.b], meta={'x': 1})
        row = table[1]
        assert row.ndim == 0
        assert row['a'] == 2
        assert row['b'] == 5
        assert row.meta['x'] == 1
        assert row.colnames == ['a', 'b']
        assert row.columns['b'].name == 'b'

    def test_index_0d_array(self):
        table = Table([self.a, self.b], meta={'x': 1})
        row = table[1]
        with pytest.raises(IndexError):
            row[0]
