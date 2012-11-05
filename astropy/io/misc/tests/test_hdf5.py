import numpy as np

from ....table import Table, Column


def test_write_simple(tmpdir):
    t = Table()
    t.add_column(Column('a', [1,2,3]))
    t.add_column(Column('b', [4,5,6]))
    t.meta['a'] = 1
    t.meta['b'] = 'hello'
    t.meta['c'] = np.array([1,2,3])
    t.write(str(tmpdir.join('test.hdf5')), name='the_table')