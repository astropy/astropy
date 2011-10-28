import re
import glob

try:
    from .. import ascii as asciitable
except ValueError:
    import asciitable

if asciitable.has_numpy:
    import numpy as np

from .common import *

def _test_values_equal(data, mem_data, numpy):
    for colname in data.dtype.names:
        matches = data[colname] == mem_data[colname]
        if numpy:
            assert(matches.all())
        else:
            assert(matches)

@has_numpy_and_not_has_numpy
def test_memory_from_table(numpy):
    table = asciitable.get_reader(numpy=numpy, Reader=asciitable.Daophot)
    data = table.read('t/daophot.dat')

    mem_table = asciitable.get_reader(Reader=asciitable.Memory, numpy=numpy)
    mem_data = mem_table.read(data)
    assert(data.dtype.names == mem_data.dtype.names)
    _test_values_equal(data, mem_data, numpy)

    mem_data = mem_table.read(mem_table)
    assert(data.dtype.names == mem_data.dtype.names)
    _test_values_equal(data, mem_data, numpy)

@has_numpy_and_not_has_numpy
def test_memory_from_LOL(numpy):
    data = [[1, 2, 3], [4, 5.2, 6.1], [8, 9, 'hello']]
    mem_table = asciitable.get_reader(Reader=asciitable.Memory, numpy=numpy)
    mem_data = mem_table.read(data)
    print(mem_data.dtype.names)
    assert(mem_data.dtype.names == ('col1', 'col2', 'col3'))
    if numpy:
        assert(mem_data[0][0] == 1)
        assert(mem_data[0][1] == 2)
        assert(mem_data[0][2] == '3')
        assert((mem_data['col2'] == np.array([2, 5.2, 9])).all())
        assert((mem_data['col3'] == np.array([3, 6.1, 'hello'])).all())
    else:
        assert(mem_data[0] == [1, 2, 3])
        assert(mem_data['col2'] == [2, 5.2, 9])
        assert(mem_data['col3'] == [3, 6.1, 'hello'])

@has_numpy_and_not_has_numpy
def test_memory_from_LOL2(numpy):
    data = [[1, 2, 3], [4, 5.2, 6.1], [8, 9, 'hello']]
    mem_table = asciitable.get_reader(Reader=asciitable.Memory, numpy=numpy, names=('c1','c2','c3'))
    mem_data = mem_table.read(data)
    print(mem_data.dtype.names)
    assert(mem_data.dtype.names == ('c1', 'c2', 'c3'))
    if numpy:
        assert(mem_data[0][0] == 1)
        assert(mem_data[0][1] == 2)
        assert(mem_data[0][2] == '3')
        assert((mem_data['c2'] == np.array([2, 5.2, 9])).all())
        assert((mem_data['c3'] == np.array([3, 6.1, 'hello'])).all())
    else:
        assert(mem_data[0] == [1, 2, 3])
        assert(mem_data['c2'] == [2, 5.2, 9])
        assert(mem_data['c3'] == [3, 6.1, 'hello'])

@has_numpy_and_not_has_numpy
def test_memory_from_DOL(numpy):
    data = {'c1': [1, 2, 3],
            'c2': [4, 5.2, 6.1],
            'c3': [8, 9, 'hello']}
    mem_table = asciitable.get_reader(Reader=asciitable.Memory, numpy=numpy,
                                      names=sorted(data.keys()))
    mem_data = mem_table.read(data)
    assert(mem_data.dtype.names == ('c1', 'c2', 'c3'))
    if numpy:
        assert(mem_data[0][0] == 1)
        assert(mem_data[0][1] == 4)
        assert(mem_data[0][2] == '8')
        assert((mem_data['c2'] == np.array([4, 5.2, 6.1])).all())
        assert((mem_data['c3'] == np.array([8, 9, 'hello'])).all())
    else:
        assert(mem_data[0] == [1, 4, 8])
        assert(mem_data['c2'] == [4, 5.2, 6.1])
        assert(mem_data['c3'] == [8, 9, 'hello'])
