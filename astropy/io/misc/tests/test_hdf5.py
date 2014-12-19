# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os

import numpy as np

from ....tests.helper import pytest, catch_warnings
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
        return [b'abc', b'def', b'ghi']
    else:
        return [1, 2, 3]


@pytest.mark.skipif('not HAS_H5PY')
def test_write_nopath(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    with pytest.raises(ValueError) as exc:
        t1.write(test_file)
    assert exc.value.args[0] == "table path should be set via the path= argument"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_notable_nopath(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    h5py.File(test_file, 'w').close()  # create empty file
    with pytest.raises(ValueError) as exc:
        t1 = Table.read(test_file, path='/', format='hdf5')
    assert exc.value.args[0] == 'no table found in HDF5 group /'


@pytest.mark.skipif('not HAS_H5PY')
def test_read_nopath(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table')
    t2 = Table.read(test_file)
    assert np.all(t1['a'] == t2['a'])


@pytest.mark.skipif('not HAS_H5PY')
def test_write_invalid_path(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    with pytest.raises(ValueError) as exc:
        t1.write(test_file, path='test/')
    assert exc.value.args[0] == "table path should end with table name, not /"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_invalid_path(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table')
    with pytest.raises(IOError) as exc:
        Table.read(test_file, path='test/')
    assert exc.value.args[0] == "Path test/ does not exist"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_missing_group(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    h5py.File(test_file, 'w').close()  # create empty file
    with pytest.raises(IOError) as exc:
        Table.read(test_file, path='test/path/table')
    assert exc.value.args[0] == "Path test/path/table does not exist"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_missing_table(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    with h5py.File(test_file, 'w') as f:
        f.create_group('test').create_group('path')
    with pytest.raises(IOError) as exc:
        Table.read(test_file, path='test/path/table')
    assert exc.value.args[0] == "Path test/path/table does not exist"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_missing_group_fileobj(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    with h5py.File(test_file, 'w') as f:
        with pytest.raises(IOError) as exc:
            Table.read(f, path='test/path/table')
        assert exc.value.args[0] == "Path test/path/table does not exist"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_simple(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table')
    t2 = Table.read(test_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
@pytest.mark.skipif(os.environ.get('APPVEYOR'), reason="fails on AppVeyor")
def test_read_write_existing_table(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table')
    with pytest.raises(IOError) as exc:
        t1.write(test_file, path='the_table', append=True)
    assert exc.value.args[0] == "Table the_table already exists"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_memory(tmpdir):
    with h5py.File('test', 'w', driver='core', backing_store=False) as output_file:
        t1 = Table()
        t1.add_column(Column(name='a', data=[1, 2, 3]))
        t1.write(output_file, path='the_table')
        t2 = Table.read(output_file, path='the_table')
        assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    h5py.File(test_file, 'w').close()  # create empty file
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    with pytest.raises(IOError) as exc:
        t1.write(test_file, path='the_table')
    assert exc.value.args[0].startswith("File exists:")


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing_overwrite(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    h5py.File(test_file, 'w').close()  # create empty file
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table', overwrite=True)
    t2 = Table.read(test_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing_append(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    h5py.File(test_file, 'w').close()  # create empty file
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table_1', append=True)
    t1.write(test_file, path='the_table_2', append=True)
    t2 = Table.read(test_file, path='the_table_1')
    assert np.all(t2['a'] == [1, 2, 3])
    t3 = Table.read(test_file, path='the_table_2')
    assert np.all(t3['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_write_existing_append_groups(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    with h5py.File(test_file, 'w') as f:
        f.create_group('test_1')
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='test_1/the_table_1', append=True)
    t1.write(test_file, path='test_2/the_table_2', append=True)
    t2 = Table.read(test_file, path='test_1/the_table_1')
    assert np.all(t2['a'] == [1, 2, 3])
    t3 = Table.read(test_file, path='test_2/the_table_2')
    assert np.all(t3['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
@pytest.mark.skipif(os.environ.get('APPVEYOR'), reason="fails on AppVeyor")
def test_read_write_existing_append_overwrite(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='table1')
    t1.write(test_file, path='table2', append=True)
    t1v2 = Table()
    t1v2.add_column(Column(name='a', data=[4, 5, 6]))
    with pytest.raises(IOError) as exc:
        t1v2.write(test_file, path='table1', append=True)
    assert exc.value.args[0] == 'Table table1 already exists'
    t1v2.write(test_file, path='table1', append=True, overwrite=True)
    t2 = Table.read(test_file, path='table1')
    assert np.all(t2['a'] == [4, 5, 6])
    t3 = Table.read(test_file, path='table2')
    assert np.all(t3['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_fileobj(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table')

    import h5py
    with h5py.File(test_file, 'r') as input_file:
        t2 = Table.read(input_file, path='the_table')
        assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_filobj_path(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='path/to/data/the_table')

    import h5py
    with h5py.File(test_file, 'r') as input_file:
        t2 = Table.read(input_file, path='path/to/data/the_table')
        assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_filobj_group_path(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='path/to/data/the_table')

    import h5py
    with h5py.File(test_file, 'r') as input_file:
        t2 = Table.read(input_file['path/to'], path='data/the_table')
        assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_read_wrong_fileobj():

    class FakeFile(object):
        def read(self):
            pass

    f = FakeFile()

    with pytest.raises(TypeError) as exc:
        t1 = Table.read(f, format='hdf5')
    assert exc.value.args[0] == 'h5py can only open regular files'


@pytest.mark.skipif('not HAS_H5PY')
def test_write_fileobj(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    import h5py
    with h5py.File(test_file, 'w') as output_file:
        t1 = Table()
        t1.add_column(Column(name='a', data=[1, 2, 3]))
        t1.write(output_file, path='the_table')

    t2 = Table.read(test_file, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_write_filobj_group(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    import h5py
    with h5py.File(test_file, 'w') as output_file:
        t1 = Table()
        t1.add_column(Column(name='a', data=[1, 2, 3]))
        t1.write(output_file, path='path/to/data/the_table')

    t2 = Table.read(test_file, path='path/to/data/the_table')
    assert np.all(t2['a'] == [1, 2, 3])


@pytest.mark.skipif('not HAS_H5PY')
def test_write_wrong_type():

    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    with pytest.raises(TypeError) as exc:
        t1.write(1212, path='path/to/data/the_table', format='hdf5')
    assert exc.value.args[0] == ('output should be a string '
                                 'or an h5py File or Group object')


@pytest.mark.skipif('not HAS_H5PY')
@pytest.mark.parametrize(('dtype'), ALL_DTYPES)
def test_preserve_single_dtypes(tmpdir, dtype):

    test_file = str(tmpdir.join('test.hdf5'))

    values = _default_values(dtype)

    t1 = Table()
    t1.add_column(Column(name='a', data=np.array(values, dtype=dtype)))
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
        t1.add_column(Column(name=str(dtype), data=np.array(values, dtype=dtype)))

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
    t1.add_column(Column(name='a', data=[1, 2, 3]))

    t1.meta['a'] = 1
    t1.meta['b'] = 'hello'
    t1.meta['c'] = 3.14159
    t1.meta['d'] = True
    t1.meta['e'] = np.array([1, 2, 3])

    t1.write(test_file, path='the_table')

    t2 = Table.read(test_file, path='the_table')

    for key in t1.meta:
        assert np.all(t1.meta[key] == t2.meta[key])


@pytest.mark.skipif('not HAS_H5PY')
def test_skip_meta(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))

    t1.meta['a'] = 1
    t1.meta['b'] = 'hello'
    t1.meta['c'] = 3.14159
    t1.meta['d'] = True
    t1.meta['e'] = np.array([1, 2, 3])
    t1.meta['f'] = str

    with catch_warnings() as w:
        t1.write(test_file, path='the_table')
    assert len(w) == 1
    assert str(w[0].message).startswith(
        "Attribute `f` of type {0} cannot be written to HDF5 files - skipping".format(type(t1.meta['f'])))

@pytest.mark.skipif('not HAS_H5PY')
def test_read_h5py_objects(tmpdir):

    # Regression test - ensure that Datasets are recognized automatically

    test_file = str(tmpdir.join('test.hdf5'))

    import h5py
    with h5py.File(test_file, 'w') as output_file:
        t1 = Table()
        t1.add_column(Column(name='a', data=[1, 2, 3]))
        t1.write(output_file, path='the_table')

    f = h5py.File(test_file)

    t2 = Table.read(f, path='the_table')
    assert np.all(t2['a'] == [1, 2, 3])

    t3 = Table.read(f['/'], path='the_table')
    assert np.all(t3['a'] == [1, 2, 3])

    t4 = Table.read(f['the_table'])
    assert np.all(t4['a'] == [1, 2, 3])
