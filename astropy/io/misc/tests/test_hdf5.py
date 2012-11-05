from __future__ import print_function

import pytest
import numpy as np

from ....table import Table, Column


ALL_DTYPES = [np.uint8, np.uint16, np.uint32, np.uint64, np.int8,
              np.int16, np.int32, np.int64, np.float32, np.float64,
              np.bool, '|S3']


def _default_values(dtype):
    if dtype == np.bool:
        return [0, 1, 1]
    elif dtype == '|S3':
        return ['abc', 'def', 'ghi']
    else:
        return [1, 2, 3]


def test_read_write_simple(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, name='the_table')
    t2 = Table.read(test_file, name='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


def test_read_fileobj(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, name='the_table')

    import h5py
    input_file = h5py.File(test_file, 'r')

    t2 = Table.read(input_file, name='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


def test_read_filobj_group(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, name='the_table', group='path/to/data')

    import h5py
    input_file = h5py.File(test_file, 'r')

    t2 = Table.read(input_file, name='the_table', group='path/to/data')
    assert np.all(t2['a'] == [1, 2, 3])


def test_write_fileobj(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    import h5py
    output_file = h5py.File(test_file, 'w')

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(output_file, name='the_table')
    output_file.close()

    t2 = Table.read(test_file, name='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


def test_write_filobj_group(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    import h5py
    output_file = h5py.File(test_file, 'w')

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(output_file, name='the_table', group='path/to/data')
    output_file.close()

    t2 = Table.read(test_file, name='the_table', group='path/to/data')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.parametrize(('dtype'), ALL_DTYPES)
def test_preserve_single_dtypes(tmpdir, dtype):

    test_file = str(tmpdir.join('test.hdf5'))

    values = _default_values(dtype)

    t1 = Table()
    t1.add_column(Column('a', np.array(values, dtype=dtype)))
    t1.write(test_file, name='the_table')

    t2 = Table.read(test_file, name='the_table')

    assert np.all(t2['a'] == values)
    assert t2['a'].dtype == dtype


def test_preserve_all_dtypes(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        t1.add_column(Column(str(dtype), np.array(values, dtype=dtype)))

    t1.write(test_file, name='the_table')

    t2 = Table.read(test_file, name='the_table')

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        assert np.all(t2[str(dtype)] == values)
        assert t2[str(dtype)].dtype == dtype


def test_preserve_meta(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))

    t1.meta['a'] = 1
    t1.meta['b'] = 'hello'
    t1.meta['c'] = 3.14159
    t1.meta['d'] = True
    t1.meta['e'] = np.array([1, 2, 3])

    t1.write(test_file, name='the_table')

    t2 = Table.read(test_file, name='the_table')

    for key in t1.meta:
        assert np.all(t1.meta[key] == t2.meta[key])
