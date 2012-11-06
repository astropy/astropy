# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import pytest
import numpy as np

from ....table import Table, Column


try:
    import h5py
except ImportError:
    HAS_H5PY = False
else:
    HAS_H5PY = True


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


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_simple(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, path='the_table')
    t2 = Table.read(test_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_memory(tmpdir):
    output_file = h5py.File('test', driver='core', backing_store=False)
    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(output_file, path='the_table')
    t2 = Table.read(output_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    f = h5py.File(test_file, 'w')
    f.close()
    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    with pytest.raises(IOError) as exc:
        t1.write(test_file, path='the_table')
    assert exc.value.args[0].startswith("File exists:")


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing_overwrite(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    f = h5py.File(test_file, 'w')
    f.close()
    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, path='the_table', overwrite=True)
    t2 = Table.read(test_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing_append(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    f = h5py.File(test_file, 'w')
    f.close()
    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, path='the_table_1', append=True)
    t1.write(test_file, path='the_table_2', append=True)
    t2 = Table.read(test_file, path='the_table_1')
    assert np.all(t2['a'] == [1, 2, 3])
    t3 = Table.read(test_file, path='the_table_2')
    assert np.all(t3['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing_append_groups(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    f = h5py.File(test_file, 'w')
    f.create_group('test_1')
    f.close()
    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, path='test_1/the_table_1', append=True)
    t1.write(test_file, path='test_2/the_table_2', append=True)
    t2 = Table.read(test_file, path='test_1/the_table_1')
    assert np.all(t2['a'] == [1, 2, 3])
    t3 = Table.read(test_file, path='test_2/the_table_2')
    assert np.all(t3['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_fileobj(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, path='the_table')

    import h5py
    input_file = h5py.File(test_file, 'r')

    t2 = Table.read(input_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_filobj_group(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(test_file, path='path/to/data/the_table')

    import h5py
    input_file = h5py.File(test_file, 'r')

    t2 = Table.read(input_file, path='path/to/data/the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_write_fileobj(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    import h5py
    output_file = h5py.File(test_file, 'w')

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(output_file, path='the_table')
    output_file.close()

    t2 = Table.read(test_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_write_filobj_group(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    import h5py
    output_file = h5py.File(test_file, 'w')

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))
    t1.write(output_file, path='path/to/data/the_table')
    output_file.close()

    t2 = Table.read(test_file, path='path/to/data/the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
@pytest.mark.parametrize(('dtype'), ALL_DTYPES)
def test_preserve_single_dtypes(tmpdir, dtype):

    test_file = str(tmpdir.join('test.hdf5'))

    values = _default_values(dtype)

    t1 = Table()
    t1.add_column(Column('a', np.array(values, dtype=dtype)))
    t1.write(test_file, path='the_table')

    t2 = Table.read(test_file, path='the_table')

    assert np.all(t2['a'] == values)
    assert t2['a'].dtype == dtype


@pytest.mark.skipif('not HAS_H5PY')
def test_preserve_all_dtypes(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        t1.add_column(Column(str(dtype), np.array(values, dtype=dtype)))

    t1.write(test_file, path='the_table')

    t2 = Table.read(test_file, path='the_table')

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        assert np.all(t2[str(dtype)] == values)
        assert t2[str(dtype)].dtype == dtype


@pytest.mark.skipif('not HAS_H5PY')
def test_preserve_meta(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column('a', [1, 2, 3]))

    t1.meta['a'] = 1
    t1.meta['b'] = 'hello'
    t1.meta['c'] = 3.14159
    t1.meta['d'] = True
    t1.meta['e'] = np.array([1, 2, 3])

    t1.write(test_file, path='the_table')

    t2 = Table.read(test_file, path='the_table')

    for key in t1.meta:
        assert np.all(t1.meta[key] == t2.meta[key])
