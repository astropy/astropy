import re
import glob
import numpy as np
from ... import ascii as asciitable

from .common import (raises,
                     assert_equal, assert_almost_equal, assert_true,
                     setup_function, teardown_function)


def _test_values_equal(data, mem_data):
    for colname in data.dtype.names:
        matches = data[colname] == mem_data[colname]
        assert(matches.all())


def test_memory_from_table():
    table = asciitable.get_reader(Reader=asciitable.Daophot)
    data = table.read('t/daophot.dat')

    mem_table = asciitable.get_reader(Reader=asciitable.Memory)
    mem_data = mem_table.read(data)
    assert(data.dtype.names == mem_data.dtype.names)
    _test_values_equal(data, mem_data)

    mem_data = mem_table.read(mem_table)
    assert(data.dtype.names == mem_data.dtype.names)
    _test_values_equal(data, mem_data)


def test_memory_from_LOL():
    data = [[1, 2, 3], [4, 5.2, 6.1], [8, 9, 'hello']]
    mem_table = asciitable.get_reader(Reader=asciitable.Memory)
    mem_data = mem_table.read(data)
    print(mem_data.dtype.names)
    assert(mem_data.dtype.names == ('col1', 'col2', 'col3'))
    assert(mem_data[0][0] == 1)
    assert(mem_data[0][1] == 2)
    assert(mem_data[0][2] == '3')
    assert((mem_data['col2'] == np.array([2, 5.2, 9])).all())
    assert((mem_data['col3'] == np.array([3, 6.1, 'hello'])).all())


def test_memory_from_LOL2():
    data = [[1, 2, 3], [4, 5.2, 6.1], [8, 9, 'hello']]
    mem_table = asciitable.get_reader(Reader=asciitable.Memory,
                                      names=('c1', 'c2', 'c3'))
    mem_data = mem_table.read(data)
    print(mem_data.dtype.names)
    assert(mem_data.dtype.names == ('c1', 'c2', 'c3'))
    assert(mem_data[0][0] == 1)
    assert(mem_data[0][1] == 2)
    assert(mem_data[0][2] == '3')
    assert((mem_data['c2'] == np.array([2, 5.2, 9])).all())
    assert((mem_data['c3'] == np.array([3, 6.1, 'hello'])).all())


def test_memory_from_DOL():
    data = {'c1': [1, 2, 3],
            'c2': [4, 5.2, 6.1],
            'c3': [8, 9, 'hello']}
    mem_table = asciitable.get_reader(Reader=asciitable.Memory,
                                      names=sorted(data.keys()))
    mem_data = mem_table.read(data)
    assert(mem_data.dtype.names == ('c1', 'c2', 'c3'))
    assert(mem_data[0][0] == 1)
    assert(mem_data[0][1] == 4)
    assert(mem_data[0][2] == '8')
    assert((mem_data['c2'] == np.array([4, 5.2, 6.1])).all())
    assert((mem_data['c3'] == np.array([8, 9, 'hello'])).all())
