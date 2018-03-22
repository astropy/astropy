# Licensed under a 3-clause BSD style license - see LICENSE.rst


import pytest
import numpy as np

from ....tests.helper import catch_warnings
from ....table import Table, QTable, NdarrayMixin, Column

from .... import units as u

from ....coordinates import SkyCoord, Latitude, Longitude, Angle, EarthLocation
from ....time import Time, TimeDelta
from ....units.quantity import QuantityInfo

try:
    import h5py
except ImportError:
    HAS_H5PY = False
else:
    HAS_H5PY = True

try:
    import yaml
except ImportError:
    HAS_YAML = False
else:
    HAS_YAML = True

ALL_DTYPES = [np.uint8, np.uint16, np.uint32, np.uint64, np.int8,
              np.int16, np.int32, np.int64, np.float32, np.float64,
              np.bool_, '|S3']


def _default_values(dtype):
    if dtype == np.bool_:
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
    with pytest.raises(OSError) as exc:
        Table.read(test_file, path='test/')
    assert exc.value.args[0] == "Path test/ does not exist"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_missing_group(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    h5py.File(test_file, 'w').close()  # create empty file
    with pytest.raises(OSError) as exc:
        Table.read(test_file, path='test/path/table')
    assert exc.value.args[0] == "Path test/path/table does not exist"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_missing_table(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    with h5py.File(test_file, 'w') as f:
        f.create_group('test').create_group('path')
    with pytest.raises(OSError) as exc:
        Table.read(test_file, path='test/path/table')
    assert exc.value.args[0] == "Path test/path/table does not exist"


@pytest.mark.skipif('not HAS_H5PY')
def test_read_missing_group_fileobj(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    with h5py.File(test_file, 'w') as f:
        with pytest.raises(OSError) as exc:
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
def test_read_write_existing_table(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='the_table')
    with pytest.raises(OSError) as exc:
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
    with pytest.raises(OSError) as exc:
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
def test_read_write_existing_append_overwrite(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.write(test_file, path='table1')
    t1.write(test_file, path='table2', append=True)
    t1v2 = Table()
    t1v2.add_column(Column(name='a', data=[4, 5, 6]))
    with pytest.raises(OSError) as exc:
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

    class FakeFile:
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


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
def test_preserve_serialized(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1['a'] = Column(data=[1, 2, 3], unit="s")
    t1['a'].meta['a0'] = "A0"
    t1['a'].meta['a1'] = {"a1": [0, 1]}
    t1['a'].format = '7.3f'
    t1['a'].description = 'A column'
    t1.meta['b'] = 1
    t1.meta['c'] = {"c0": [0, 1]}

    t1.write(test_file, path='the_table', serialize_meta=True, overwrite=True)

    t2 = Table.read(test_file, path='the_table')

    assert t1['a'].unit == t2['a'].unit
    assert t1['a'].format == t2['a'].format
    assert t1['a'].description == t2['a'].description
    assert t1['a'].meta == t2['a'].meta
    assert t1.meta == t2.meta


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
def test_preserve_serialized_compatibility_mode(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1['a'] = Column(data=[1, 2, 3], unit="s")
    t1['a'].meta['a0'] = "A0"
    t1['a'].meta['a1'] = {"a1": [0, 1]}
    t1['a'].format = '7.3f'
    t1['a'].description = 'A column'
    t1.meta['b'] = 1
    t1.meta['c'] = {"c0": [0, 1]}

    with catch_warnings() as w:
        t1.write(test_file, path='the_table', serialize_meta=True,
                 overwrite=True, compatibility_mode=True)

    assert str(w[0].message).startswith(
        "compatibility mode for writing is deprecated")

    t2 = Table.read(test_file, path='the_table')

    assert t1['a'].unit == t2['a'].unit
    assert t1['a'].format == t2['a'].format
    assert t1['a'].description == t2['a'].description
    assert t1['a'].meta == t2['a'].meta
    assert t1.meta == t2.meta


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
def test_preserve_serialized_in_complicated_path(tmpdir):
    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1['a'] = Column(data=[1, 2, 3], unit="s")
    t1['a'].meta['a0'] = "A0"
    t1['a'].meta['a1'] = {"a1": [0, 1]}
    t1['a'].format = '7.3f'
    t1['a'].description = 'A column'
    t1.meta['b'] = 1
    t1.meta['c'] = {"c0": [0, 1]}

    t1.write(test_file, path='the_table/complicated/path', serialize_meta=True,
             overwrite=True)

    t2 = Table.read(test_file, path='the_table/complicated/path')

    assert t1['a'].format == t2['a'].format
    assert t1['a'].unit == t2['a'].unit
    assert t1['a'].description == t2['a'].description
    assert t1['a'].meta == t2['a'].meta
    assert t1.meta == t2.meta

@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
def test_metadata_very_large(tmpdir):
    """Test that very large datasets work, now!"""
    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1['a'] = Column(data=[1, 2, 3], unit="s")
    t1['a'].meta['a0'] = "A0"
    t1['a'].meta['a1'] = {"a1": [0, 1]}
    t1['a'].format = '7.3f'
    t1['a'].description = 'A column'
    t1.meta['b'] = 1
    t1.meta['c'] = {"c0": [0, 1]}
    t1.meta["meta_big"] = "0" * (2 ** 16 + 1)
    t1.meta["meta_biggerstill"] = "0" * (2 ** 18)

    t1.write(test_file, path='the_table', serialize_meta=True, overwrite=True)

    t2 = Table.read(test_file, path='the_table')

    assert t1['a'].unit == t2['a'].unit
    assert t1['a'].format == t2['a'].format
    assert t1['a'].description == t2['a'].description
    assert t1['a'].meta == t2['a'].meta
    assert t1.meta == t2.meta


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
def test_metadata_very_large_fails_compatibility_mode(tmpdir):
    """Test that very large metadata do not work in compatibility mode."""
    test_file = str(tmpdir.join('test.hdf5'))
    t1 = Table()
    t1['a'] = Column(data=[1, 2, 3])
    t1.meta["meta"] = "0" * (2 ** 16 + 1)
    with catch_warnings() as w:
        t1.write(test_file, path='the_table', serialize_meta=True,
                 overwrite=True, compatibility_mode=True)
    assert len(w) == 2

    # Error message slightly changed in h5py 2.7.1, thus the 2part assert
    assert str(w[1].message).startswith(
        "Attributes could not be written to the output HDF5 "
        "file: Unable to create attribute ")
    assert "bject header message is too large" in str(w[1].message)


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


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
def test_fail_meta_serialize(tmpdir):

    test_file = str(tmpdir.join('test.hdf5'))

    t1 = Table()
    t1.add_column(Column(name='a', data=[1, 2, 3]))
    t1.meta['f'] = str

    with pytest.raises(Exception) as err:
        t1.write(test_file, path='the_table', serialize_meta=True)
    assert "cannot represent an object: <class 'str'>" in str(err)


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

    f.close()         # don't raise an error in 'test --open-files'


def assert_objects_equal(obj1, obj2, attrs, compare_class=True):
    if compare_class:
        assert obj1.__class__ is obj2.__class__

    info_attrs = ['info.name', 'info.format', 'info.unit', 'info.description', 'info.meta']
    for attr in attrs + info_attrs:
        a1 = obj1
        a2 = obj2
        for subattr in attr.split('.'):
            try:
                a1 = getattr(a1, subattr)
                a2 = getattr(a2, subattr)
            except AttributeError:
                a1 = a1[subattr]
                a2 = a2[subattr]

        # Mixin info.meta can None instead of empty OrderedDict(), #6720 would
        # fix this.
        if attr == 'info.meta':
            if a1 is None:
                a1 = {}
            if a2 is None:
                a2 = {}

        assert np.all(a1 == a2)

# Testing HDF5 table read/write with mixins.  This is mostly
# copied from FITS mixin testing.

el = EarthLocation(x=1 * u.km, y=3 * u.km, z=5 * u.km)
el2 = EarthLocation(x=[1, 2] * u.km, y=[3, 4] * u.km, z=[5, 6] * u.km)
sc = SkyCoord([1, 2], [3, 4], unit='deg,deg', frame='fk4',
              obstime='J1990.5')
scc = sc.copy()
scc.representation = 'cartesian'
tm = Time([2450814.5, 2450815.5], format='jd', scale='tai', location=el)


mixin_cols = {
    'tm': tm,
    'dt': TimeDelta([1, 2] * u.day),
    'sc': sc,
    'scc': scc,
    'scd': SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,m', frame='fk4',
                    obstime=['J1990.5', 'J1991.5']),
    'q': [1, 2] * u.m,
    'lat': Latitude([1, 2] * u.deg),
    'lon': Longitude([1, 2] * u.deg, wrap_angle=180. * u.deg),
    'ang': Angle([1, 2] * u.deg),
    'el2': el2,
}

time_attrs = ['value', 'shape', 'format', 'scale', 'location']
compare_attrs = {
    'c1': ['data'],
    'c2': ['data'],
    'tm': time_attrs,
    'dt': ['shape', 'value', 'format', 'scale'],
    'sc': ['ra', 'dec', 'representation', 'frame.name'],
    'scc': ['x', 'y', 'z', 'representation', 'frame.name'],
    'scd': ['ra', 'dec', 'distance', 'representation', 'frame.name'],
    'q': ['value', 'unit'],
    'lon': ['value', 'unit', 'wrap_angle'],
    'lat': ['value', 'unit'],
    'ang': ['value', 'unit'],
    'el2': ['x', 'y', 'z', 'ellipsoid'],
    'nd': ['x', 'y', 'z'],
}


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
def test_hdf5_mixins_qtable_to_table(tmpdir):
    """Test writing as QTable and reading as Table.  Ensure correct classes
    come out.
    """
    filename = str(tmpdir.join('test_simple.hdf5'))

    names = sorted(mixin_cols)

    t = QTable([mixin_cols[name] for name in names], names=names)
    t.write(filename, format='hdf5', path='root', serialize_meta=True)
    t2 = Table.read(filename, format='hdf5', path='root')

    assert t.colnames == t2.colnames

    for name, col in t.columns.items():
        col2 = t2[name]

        # Special-case Time, which does not yet support round-tripping
        # the format.
        if isinstance(col2, Time):
            col2.format = col.format

        attrs = compare_attrs[name]
        compare_class = True

        if isinstance(col.info, QuantityInfo):
            # Downgrade Quantity to Column + unit
            assert type(col2) is Column
            # Class-specific attributes like `value` or `wrap_angle` are lost.
            attrs = ['unit']
            compare_class = False
            # Compare data values here (assert_objects_equal doesn't know how in this case)
            assert np.all(col.value == col2)

        assert_objects_equal(col, col2, attrs, compare_class)


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
@pytest.mark.parametrize('table_cls', (Table, QTable))
def test_hdf5_mixins_as_one(table_cls, tmpdir):
    """Test write/read all cols at once and validate intermediate column names"""
    filename = str(tmpdir.join('test_simple.hdf5'))
    names = sorted(mixin_cols)

    serialized_names = ['ang',
                        'dt.jd1', 'dt.jd2',
                        'el2.x', 'el2.y', 'el2.z',
                        'lat',
                        'lon',
                        'q',
                        'sc.ra', 'sc.dec',
                        'scc.x', 'scc.y', 'scc.z',
                        'scd.ra', 'scd.dec', 'scd.distance',
                        'scd.obstime.jd1', 'scd.obstime.jd2',
                        'tm.jd1', 'tm.jd2',
                        ]

    t = table_cls([mixin_cols[name] for name in names], names=names)
    t.meta['C'] = 'spam'
    t.meta['comments'] = ['this', 'is', 'a', 'comment']
    t.meta['history'] = ['first', 'second', 'third']

    t.write(filename, format="hdf5", path='root', serialize_meta=True)

    t2 = table_cls.read(filename, format='hdf5', path='root')
    assert t2.meta['C'] == 'spam'
    assert t2.meta['comments'] == ['this', 'is', 'a', 'comment']
    assert t2.meta['history'] == ['first', 'second', 'third']

    assert t.colnames == t2.colnames

    # Read directly via hdf5 and confirm column names
    h5 = h5py.File(filename, 'r')
    assert list(h5['root'].dtype.names) == serialized_names
    h5.close()


@pytest.mark.skipif('not HAS_H5PY or not HAS_YAML')
@pytest.mark.parametrize('name_col', list(mixin_cols.items()))
@pytest.mark.parametrize('table_cls', (Table, QTable))
def test_hdf5_mixins_per_column(table_cls, name_col, tmpdir):
    """Test write/read one col at a time and do detailed validation"""
    filename = str(tmpdir.join('test_simple.hdf5'))
    name, col = name_col

    c = [1.0, 2.0]
    t = table_cls([c, col, c], names=['c1', name, 'c2'])
    t[name].info.description = 'my description'
    t[name].info.meta = {'list': list(range(50)), 'dict': {'a': 'b' * 200}}

    if not t.has_mixin_columns:
        pytest.skip('column is not a mixin (e.g. Quantity subclass in Table)')

    if isinstance(t[name], NdarrayMixin):
        pytest.xfail('NdarrayMixin not supported')

    t.write(filename, format="hdf5", path='root', serialize_meta=True)
    t2 = table_cls.read(filename, format='hdf5', path='root')

    assert t.colnames == t2.colnames

    for colname in t.colnames:
        assert_objects_equal(t[colname], t2[colname], compare_attrs[colname])

    # Special case to make sure Column type doesn't leak into Time class data
    if name.startswith('tm'):
        assert t2[name]._time.jd1.__class__ is np.ndarray
        assert t2[name]._time.jd2.__class__ is np.ndarray


@pytest.mark.skipif('HAS_YAML or not HAS_H5PY')
def test_warn_for_dropped_info_attributes(tmpdir):
    filename = str(tmpdir.join('test.hdf5'))
    t = Table([[1, 2]])
    t['col0'].info.description = 'hello'
    with catch_warnings() as warns:
        t.write(filename, path='root', serialize_meta=False)
    assert len(warns) == 1
    assert str(warns[0].message).startswith(
        "table contains column(s) with defined 'unit'")


@pytest.mark.skipif('HAS_YAML or not HAS_H5PY')
def test_error_for_mixins_but_no_yaml(tmpdir):
    filename = str(tmpdir.join('test.hdf5'))
    t = Table([mixin_cols['sc']])
    with pytest.raises(TypeError) as err:
        t.write(filename, path='root', serialize_meta=True)
    assert "cannot write type SkyCoord column 'col0' to HDF5 without PyYAML" in str(err)
