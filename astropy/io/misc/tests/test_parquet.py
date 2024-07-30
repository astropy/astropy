# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import (
    Angle,
    CartesianRepresentation,
    EarthLocation,
    Latitude,
    Longitude,
    SkyCoord,
    SphericalCosLatDifferential,
    SphericalRepresentation,
)
from astropy.io.misc.parquet import get_pyarrow, parquet_identify
from astropy.table import Column, NdarrayMixin, QTable, Table
from astropy.table.table_helpers import simple_table
from astropy.time import Time, TimeDelta
from astropy.units import allclose as quantity_allclose
from astropy.units.quantity import QuantityInfo
from astropy.utils.compat.optional_deps import HAS_PANDAS
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

# Skip all tests in this file if we cannot import pyarrow
pyarrow = pytest.importorskip("pyarrow")

ALL_DTYPES = [
    np.uint8,
    np.uint16,
    np.uint32,
    np.uint64,
    np.int8,
    np.int16,
    np.int32,
    np.int64,
    np.float32,
    np.float64,
    np.bool_,
    "|S3",
    "U3",
]

BIGENDIAN_DTYPES = [
    np.dtype(">i4"),
    np.dtype(">i8"),
    np.dtype(">f4"),
    np.dtype(">f8"),
]


def _default_values(dtype):
    if dtype == np.bool_:
        return [0, 1, 1]
    elif dtype == "|S3":
        return [b"abc", b"def", b"ghi"]
    elif dtype == "U3":
        return ["abc", "def", "ghi"]
    else:
        return [1, 2, 3]


def _default_array_values(dtype):
    values = _default_values(dtype)
    return [values for i in range(3)]


def _default_var_length_array_values(dtype):
    values = _default_values(dtype)
    return [
        [
            values[0],
        ],
        [
            values[0],
            values[1],
        ],
        [
            values[0],
            values[1],
            values[2],
        ],
    ]


def test_read_write_simple(tmp_path):
    """Test writing/reading a simple parquet file."""
    test_file = tmp_path / "test.parquet"
    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))
    t1.write(test_file)
    t2 = Table.read(test_file)
    assert np.all(t2["a"] == [1, 2, 3])


def test_read_write_existing(tmp_path):
    """Test writing an existing file without overwriting."""
    test_file = tmp_path / "test.parquet"
    with open(test_file, "w") as f:  # create empty file
        pass
    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))

    with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
        t1.write(test_file)


def test_read_write_existing_overwrite(tmp_path):
    """Test overwriting an existing file."""

    test_file = tmp_path / "test.parquet"
    with open(test_file, "w") as f:  # create empty file
        pass
    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))
    t1.write(test_file, overwrite=True)
    t2 = Table.read(test_file)
    assert np.all(t2["a"] == [1, 2, 3])


def test_read_fileobj(tmp_path):
    """Test reading a file object."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))
    t1.write(test_file)

    import io

    with io.FileIO(test_file, mode="r") as input_file:
        t2 = Table.read(input_file)
        assert np.all(t2["a"] == [1, 2, 3])


def test_read_pathlikeobj(tmp_path):
    """Test reading a path-like object."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))
    t1.write(test_file)

    import pathlib

    p = pathlib.Path(test_file)
    t2 = Table.read(p)
    assert np.all(t2["a"] == [1, 2, 3])


def test_read_wrong_fileobj():
    """Test reading an incorrect fileobject type."""

    class FakeFile:
        def not_read(self):
            pass

    f = FakeFile()

    with pytest.raises(
        TypeError, match="pyarrow can only open path-like or file-like objects."
    ):
        Table.read(f, format="parquet")


def test_identify_wrong_fileobj():
    """Test identifying an incorrect fileobj."""

    class FakeFile:
        def not_read(self):
            pass

    f = FakeFile()

    assert not parquet_identify("test", "test", f)


def test_identify_file_wrong_extension():
    """Test identifying an incorrect extension."""

    assert not parquet_identify("test", "test.notparquet", None)


def test_identify_file_correct_extension():
    """Test identifying an incorrect extension."""

    assert parquet_identify("test", "test.parquet", None)
    assert parquet_identify("test", "test.parq", None)


def test_identify_file_noobject_nopath():
    """Test running identify with no object or path."""

    assert not parquet_identify("test", None, None)


def test_write_wrong_type():
    """Test writing to a filename of the wrong type."""

    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))
    with pytest.raises(TypeError, match="should be a string"):
        t1.write(1212, format="parquet")


def test_empty_roundtrip(tmp_path):
    """Test writing and reading an empty Table."""
    test_file = tmp_path / "test.parquet"
    t1 = Table()
    t1.write(test_file)
    t2 = Table.read(test_file)
    assert len(t2) == 0
    assert t1.colnames == t2.colnames


@pytest.mark.parametrize("dtype", ALL_DTYPES)
def test_preserve_single_dtypes(tmp_path, dtype):
    """Test that round-tripping a single column preserves datatypes."""

    test_file = tmp_path / "test.parquet"

    values = _default_values(dtype)

    t1 = Table()
    t1.add_column(Column(name="a", data=np.array(values, dtype=dtype)))
    t1.write(test_file)

    t2 = Table.read(test_file)

    assert np.all(t2["a"] == values)
    assert t2["a"].dtype == dtype


@pytest.mark.parametrize("dtype", BIGENDIAN_DTYPES)
def test_preserve_single_bigendian_dtypes(tmp_path, dtype):
    """Test that round-tripping a single big-endian column preserves data."""

    test_file = tmp_path / "test.parquet"

    values = _default_values(dtype)

    t1 = Table()
    t1.add_column(Column(name="a", data=np.array(values, dtype=dtype)))
    t1.write(test_file)

    t2 = Table.read(test_file)

    assert np.all(t2["a"] == values)
    # The parquet serialization will turn all arrays into little-endian.
    assert t2["a"].dtype == dtype.newbyteorder("<")


@pytest.mark.parametrize("dtype", ALL_DTYPES)
def test_preserve_single_array_dtypes(tmp_path, dtype):
    """Test that round-tripping a single array column preserves datatypes."""

    test_file = tmp_path / "test.parquet"

    values = _default_array_values(dtype)

    t1 = Table()
    t1.add_column(Column(name="a", data=np.array(values, dtype=dtype)))
    t1.write(test_file)

    t2 = Table.read(test_file)

    assert np.all(t2["a"] == t1["a"])
    assert np.all(t2["a"].shape == np.array(values).shape)
    assert t2["a"].dtype == dtype


@pytest.mark.parametrize("dtype", BIGENDIAN_DTYPES)
def test_preserve_single_bigendian_array_dtypes(tmp_path, dtype):
    """Test that round-tripping a single array column (big-endian) preserves data."""

    test_file = tmp_path / "test.parquet"

    values = _default_array_values(dtype)

    t1 = Table()
    t1.add_column(Column(name="a", data=np.array(values, dtype=dtype)))
    t1.write(test_file)

    t2 = Table.read(test_file)

    assert np.all(t2["a"] == t1["a"])
    assert np.all(t2["a"].shape == np.array(values).shape)
    assert t2["a"].dtype == dtype.newbyteorder("<")


@pytest.mark.parametrize("dtype", ALL_DTYPES)
def test_preserve_single_var_length_array_dtypes(tmp_path, dtype):
    """
    Test that round-tripping a single variable length array column preserves
    datatypes.
    """

    test_file = tmp_path / "test.parquet"

    values = _default_var_length_array_values(dtype)

    t1 = Table()
    data = np.array([np.array(val, dtype=dtype) for val in values], dtype=np.object_)
    t1.add_column(Column(name="a", data=data))
    t1.write(test_file)

    t2 = Table.read(test_file)

    for row1, row2 in zip(t1["a"], t2["a"]):
        assert np.all(row1 == row2)
        assert row1.dtype == row2.dtype


@pytest.mark.parametrize("dtype", BIGENDIAN_DTYPES)
def test_preserve_single_bigendian_var_length_array_dtypes(tmp_path, dtype):
    """
    Test that round-tripping a single big-endian variable length array column preserves
    datatypes.
    """

    test_file = tmp_path / "test.parquet"

    values = _default_var_length_array_values(dtype)

    t1 = Table()
    data = np.array([np.array(val, dtype=dtype) for val in values], dtype=np.object_)
    t1.add_column(Column(name="a", data=data))
    t1.write(test_file)

    t2 = Table.read(test_file)

    for row1, row2 in zip(t1["a"], t2["a"]):
        assert np.all(row1 == row2)
        assert row1.dtype.newbyteorder(">") == row2.dtype.newbyteorder(">")


def test_preserve_all_dtypes(tmp_path):
    """Test that round-tripping preserves a table with all the datatypes."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        t1.add_column(Column(name=str(dtype), data=np.array(values, dtype=dtype)))

        arr_values = _default_array_values(dtype)
        t1.add_column(
            Column(name=str(dtype) + "_arr", data=np.array(arr_values, dtype=dtype))
        )

    t1.write(test_file)

    t2 = Table.read(test_file)

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        assert np.all(t2[str(dtype)] == values)
        assert t2[str(dtype)].dtype == dtype

        arr_values = _default_array_values(dtype)
        assert np.all(t2[str(dtype) + "_arr"] == values)
        assert t2[str(dtype)].dtype == dtype
        assert np.all(t2[str(dtype) + "_arr"].shape == np.array(arr_values).shape)

    # Test just reading the schema
    schema2 = Table.read(test_file, schema_only=True)
    assert len(schema2) == 0
    assert schema2.dtype == t2.dtype


def test_preserve_all_var_length_dtypes(tmp_path):
    """Test that round-tripping preserves a table with all the var length datatypes."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()

    for dtype in ALL_DTYPES:
        varr_values = _default_var_length_array_values(dtype)
        data = np.array(
            [np.array(val, dtype=dtype) for val in varr_values], dtype=np.object_
        )
        t1.add_column(Column(name=str(dtype) + "_varr", data=data))

    t1.write(test_file)

    t2 = Table.read(test_file)

    for dtype in ALL_DTYPES:
        varr_values = _default_var_length_array_values(dtype)
        colname = str(dtype) + "_varr"
        for row1, row2 in zip(t1[colname], t2[colname]):
            assert np.all(row1 == row2)
            assert row1.dtype == row2.dtype


def test_write_empty_tables(tmp_path):
    """Test that we can save an empty table with var length datatypes."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        t1.add_column(Column(name=str(dtype), data=np.array(values, dtype=dtype)))

        arr_values = _default_array_values(dtype)
        t1.add_column(
            Column(name=str(dtype) + "_arr", data=np.array(arr_values, dtype=dtype))
        )

    # Write an empty table with values and arrays, and confirm it works.
    data = np.zeros(0, dtype=t1.dtype)
    t2 = Table(data=data)

    t2.write(test_file)

    t3 = Table.read(test_file)

    assert t3.dtype == t2.dtype

    test_file2 = tmp_path / "test2.parquet"

    t4 = Table()
    for dtype in ALL_DTYPES:
        varr_values = _default_var_length_array_values(dtype)
        data = np.array(
            [np.array(val, dtype=dtype) for val in varr_values], dtype=np.object_
        )
        t4.add_column(Column(name=str(dtype) + "_varr", data=data))

    # Write an empty table with variable-length arrays, and confirm this
    # raises an exception. (The datatype of an np.object_ type column
    # cannot be inferred from an empty table.)
    data = np.zeros(0, dtype=t4.dtype)
    t5 = Table(data=data)
    with pytest.raises(ValueError, match="Cannot serialize zero-length table") as err:
        t5.write(test_file2)


def test_heterogeneous_var_array_table(tmp_path):
    """Test exception when trying to serialize a mixed-type variable-length column."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()

    data = np.array(
        [
            np.array([0, 1, 2], dtype=np.int32),
            np.array([0, 1, 2, 3, 4], dtype=np.float64),
        ],
        dtype=np.object_,
    )
    t1.add_column(Column(name="a", data=data))

    with pytest.raises(ValueError, match="Cannot serialize mixed-type column") as err:
        t1.write(test_file)


def test_preserve_meta(tmp_path):
    """Test that writing/reading preserves metadata."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))

    t1.meta["a"] = 1
    t1.meta["b"] = "hello"
    t1.meta["c"] = 3.14159
    t1.meta["d"] = True
    t1.meta["e"] = np.array([1, 2, 3])

    t1.write(test_file)

    t2 = Table.read(test_file)

    for key in t1.meta:
        assert np.all(t1.meta[key] == t2.meta[key])


def test_preserve_serialized(tmp_path):
    """Test that writing/reading preserves unit/format/description."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()
    t1["a"] = Column(data=[1, 2, 3], unit="s")
    t1["a"].meta["a0"] = "A0"
    t1["a"].meta["a1"] = {"a1": [0, 1]}
    t1["a"].format = "7.3f"
    t1["a"].description = "A column"
    t1.meta["b"] = 1
    t1.meta["c"] = {"c0": [0, 1]}

    t1.write(test_file, overwrite=True)

    t2 = Table.read(test_file)

    assert t1["a"].unit == t2["a"].unit
    assert t1["a"].format == t2["a"].format
    assert t1["a"].description == t2["a"].description
    assert t1["a"].meta == t2["a"].meta
    assert t1.meta == t2.meta


def test_metadata_very_large(tmp_path):
    """Test that very large datasets work"""

    test_file = tmp_path / "test.parquet"

    t1 = Table()
    t1["a"] = Column(data=[1, 2, 3], unit="s")
    t1["a"].meta["a0"] = "A0"
    t1["a"].meta["a1"] = {"a1": [0, 1]}
    t1["a"].format = "7.3f"
    t1["a"].description = "A column"
    t1.meta["b"] = 1
    t1.meta["c"] = {"c0": [0, 1]}
    t1.meta["meta_big"] = "0" * (2**16 + 1)
    t1.meta["meta_biggerstill"] = "0" * (2**18)

    t1.write(test_file, overwrite=True)

    t2 = Table.read(test_file)

    assert t1["a"].unit == t2["a"].unit
    assert t1["a"].format == t2["a"].format
    assert t1["a"].description == t2["a"].description
    assert t1["a"].meta == t2["a"].meta
    assert t1.meta == t2.meta


def test_fail_meta_serialize(tmp_path):
    """Test that we cannot preserve objects in metadata."""

    test_file = tmp_path / "test.parquet"

    t1 = Table()
    t1.add_column(Column(name="a", data=[1, 2, 3]))
    t1.meta["f"] = str

    with pytest.raises(Exception) as err:
        t1.write(test_file)
    assert "cannot represent an object" in str(err.value)
    assert "<class 'str'>" in str(err.value)


def assert_objects_equal(obj1, obj2, attrs, compare_class=True):
    """Convenient routine to check objects and attributes match."""

    if compare_class:
        assert obj1.__class__ is obj2.__class__

    info_attrs = [
        "info.name",
        "info.format",
        "info.unit",
        "info.description",
        "info.meta",
    ]
    for attr in attrs + info_attrs:
        a1 = obj1
        a2 = obj2
        for subattr in attr.split("."):
            try:
                a1 = getattr(a1, subattr)
                a2 = getattr(a2, subattr)
            except AttributeError:
                a1 = a1[subattr]
                a2 = a2[subattr]

        # Mixin info.meta can None instead of empty OrderedDict(), #6720 would
        # fix this.
        if attr == "info.meta":
            if a1 is None:
                a1 = {}
            if a2 is None:
                a2 = {}

        if isinstance(a1, np.ndarray) and a1.dtype.kind == "f":
            assert quantity_allclose(a1, a2, rtol=1e-15)
        else:
            assert np.all(a1 == a2)


# Testing Parquet table read/write with mixins. This is mostly
# copied from HDF5/FITS mixin testing, and it might be good to unify it.
# Analogous tests also exist for ECSV.


el = EarthLocation(x=1 * u.km, y=3 * u.km, z=5 * u.km)
el2 = EarthLocation(x=[1, 2] * u.km, y=[3, 4] * u.km, z=[5, 6] * u.km)
sr = SphericalRepresentation([0, 1] * u.deg, [2, 3] * u.deg, 1 * u.kpc)
cr = CartesianRepresentation([0, 1] * u.pc, [4, 5] * u.pc, [8, 6] * u.pc)
sd = SphericalCosLatDifferential(
    [0, 1] * u.mas / u.yr, [0, 1] * u.mas / u.yr, 10 * u.km / u.s
)
srd = SphericalRepresentation(sr, differentials=sd)
sc = SkyCoord([1, 2], [3, 4], unit="deg,deg", frame="fk4", obstime="J1990.5")
scd = SkyCoord(
    [1, 2],
    [3, 4],
    [5, 6],
    unit="deg,deg,m",
    frame="fk4",
    obstime=["J1990.5", "J1991.5"],
)
scdc = scd.copy()
scdc.representation_type = "cartesian"
scpm = SkyCoord(
    [1, 2],
    [3, 4],
    [5, 6],
    unit="deg,deg,pc",
    pm_ra_cosdec=[7, 8] * u.mas / u.yr,
    pm_dec=[9, 10] * u.mas / u.yr,
)
scpmrv = SkyCoord(
    [1, 2],
    [3, 4],
    [5, 6],
    unit="deg,deg,pc",
    pm_ra_cosdec=[7, 8] * u.mas / u.yr,
    pm_dec=[9, 10] * u.mas / u.yr,
    radial_velocity=[11, 12] * u.km / u.s,
)
scrv = SkyCoord(
    [1, 2], [3, 4], [5, 6], unit="deg,deg,pc", radial_velocity=[11, 12] * u.km / u.s
)
tm = Time([2450814.5, 2450815.5], format="jd", scale="tai", location=el)

# NOTE: in the test below the name of the column "x" for the Quantity is
# important since it tests the fix for #10215 (namespace clash, where "x"
# clashes with "el2.x").
mixin_cols = {
    "tm": tm,
    "dt": TimeDelta([1, 2] * u.day),
    "sc": sc,
    "scd": scd,
    "scdc": scdc,
    "scpm": scpm,
    "scpmrv": scpmrv,
    "scrv": scrv,
    "x": [1, 2] * u.m,
    "qdb": [10, 20] * u.dB(u.mW),
    "qdex": [4.5, 5.5] * u.dex(u.cm / u.s**2),
    "qmag": [21, 22] * u.ABmag,
    "lat": Latitude([1, 2] * u.deg),
    "lon": Longitude([1, 2] * u.deg, wrap_angle=180.0 * u.deg),
    "ang": Angle([1, 2] * u.deg),
    "el2": el2,
    "sr": sr,
    "cr": cr,
    "sd": sd,
    "srd": srd,
}

time_attrs = ["value", "shape", "format", "scale", "location"]
compare_attrs = {
    "c1": ["data"],
    "c2": ["data"],
    "tm": time_attrs,
    "dt": ["shape", "value", "format", "scale"],
    "sc": ["ra", "dec", "representation_type", "frame.name"],
    "scd": ["ra", "dec", "distance", "representation_type", "frame.name"],
    "scdc": ["x", "y", "z", "representation_type", "frame.name"],
    "scpm": [
        "ra",
        "dec",
        "distance",
        "pm_ra_cosdec",
        "pm_dec",
        "representation_type",
        "frame.name",
    ],
    "scpmrv": [
        "ra",
        "dec",
        "distance",
        "pm_ra_cosdec",
        "pm_dec",
        "radial_velocity",
        "representation_type",
        "frame.name",
    ],
    "scrv": [
        "ra",
        "dec",
        "distance",
        "radial_velocity",
        "representation_type",
        "frame.name",
    ],
    "x": ["value", "unit"],
    "qdb": ["value", "unit"],
    "qdex": ["value", "unit"],
    "qmag": ["value", "unit"],
    "lon": ["value", "unit", "wrap_angle"],
    "lat": ["value", "unit"],
    "ang": ["value", "unit"],
    "el2": ["x", "y", "z", "ellipsoid"],
    "nd": ["x", "y", "z"],
    "sr": ["lon", "lat", "distance"],
    "cr": ["x", "y", "z"],
    "sd": ["d_lon_coslat", "d_lat", "d_distance"],
    "srd": [
        "lon",
        "lat",
        "distance",
        "differentials.s.d_lon_coslat",
        "differentials.s.d_lat",
        "differentials.s.d_distance",
    ],
}


def test_parquet_mixins_qtable_to_table(tmp_path):
    """Test writing as QTable and reading as Table. Ensure correct classes
    come out.
    """
    filename = tmp_path / "test_simple.parquet"

    names = sorted(mixin_cols)

    t = QTable([mixin_cols[name] for name in names], names=names)
    t.write(filename, format="parquet")
    t2 = Table.read(filename, format="parquet")

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
            attrs = ["unit"]
            compare_class = False
            # Compare data values here (assert_objects_equal doesn't know how in this case)
            assert np.all(col.value == col2)

        assert_objects_equal(col, col2, attrs, compare_class)


@pytest.mark.parametrize("table_cls", (Table, QTable))
def test_parquet_mixins_as_one(table_cls, tmp_path):
    """Test write/read all cols at once and validate intermediate column names"""
    filename = tmp_path / "test_simple.parquet"
    names = sorted(mixin_cols)

    t = table_cls([mixin_cols[name] for name in names], names=names)
    t.meta["C"] = "spam"
    t.meta["comments"] = ["this", "is", "a", "comment"]
    t.meta["history"] = ["first", "second", "third"]

    t.write(filename, format="parquet")

    t2 = table_cls.read(filename, format="parquet")
    assert t2.meta["C"] == "spam"
    assert t2.meta["comments"] == ["this", "is", "a", "comment"]
    assert t2.meta["history"] == ["first", "second", "third"]

    assert t.colnames == t2.colnames


@pytest.mark.parametrize("name_col", list(mixin_cols.items()))
@pytest.mark.parametrize("table_cls", (Table, QTable))
def test_parquet_mixins_per_column(table_cls, name_col, tmp_path):
    """Test write/read one col at a time and do detailed validation"""
    filename = tmp_path / "test_simple.parquet"
    name, col = name_col

    c = [1.0, 2.0]
    t = table_cls([c, col, c], names=["c1", name, "c2"])
    t[name].info.description = "my description"
    t[name].info.meta = {"list": list(range(50)), "dict": {"a": "b" * 200}}

    if not t.has_mixin_columns:
        pytest.skip("column is not a mixin (e.g. Quantity subclass in Table)")

    if isinstance(t[name], NdarrayMixin):
        pytest.xfail("NdarrayMixin not supported")

    t.write(filename, format="parquet")
    t2 = table_cls.read(filename, format="parquet")

    assert t.colnames == t2.colnames

    for colname in t.colnames:
        assert_objects_equal(t[colname], t2[colname], compare_attrs[colname])

    # Special case to make sure Column type doesn't leak into Time class data
    if name.startswith("tm"):
        assert t2[name]._time.jd1.__class__ is np.ndarray
        assert t2[name]._time.jd2.__class__ is np.ndarray


def test_round_trip_masked_table_default(tmp_path):
    """Test round-trip of MaskedColumn through Parquet using default serialization
    that writes a separate mask column.  Note:

    >>> simple_table(masked=True)
    <Table masked=True length=3>
      a      b     c
    int64 float64 str1
    ----- ------- ----
       --     1.0    c
        2     2.0   --
        3      --    e
    """
    filename = tmp_path / "test.parquet"

    t = simple_table(masked=True)  # int, float, and str cols with one masked element
    t["c"] = [b"c", b"d", b"e"]
    t["c"].mask[1] = True
    t.write(filename, format="parquet")

    t2 = Table.read(filename)
    assert t2.masked is False
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        assert np.all(t2[name].mask == t[name].mask)
        assert np.all(t2[name] == t[name])

        # Data under the mask round-trips also (unmask data to show this).
        t[name].mask = False
        t2[name].mask = False
        assert np.all(t2[name] == t[name])


@pytest.mark.parametrize("table_cls", (Table, QTable))
def test_parquet_mixins_read_one_name(table_cls, tmp_path):
    """Test write all cols at once, and read one at a time."""
    filename = tmp_path / "test_simple.parquet"
    names = sorted(mixin_cols)

    t = table_cls([mixin_cols[name] for name in names], names=names)
    t.meta["C"] = "spam"
    t.meta["comments"] = ["this", "is", "a", "comment"]
    t.meta["history"] = ["first", "second", "third"]

    t.write(filename, format="parquet")

    for name in names:
        t2 = table_cls.read(filename, format="parquet", include_names=[name])
        assert t2.meta["C"] == "spam"
        assert t2.meta["comments"] == ["this", "is", "a", "comment"]
        assert t2.meta["history"] == ["first", "second", "third"]

        assert t2.colnames == [name]


@pytest.mark.parametrize("table_cls", (Table, QTable))
def test_parquet_mixins_read_exclude_names(table_cls, tmp_path):
    """Test write all cols at once, and read all but one at a time."""
    filename = tmp_path / "test_simple.parquet"
    names = sorted(mixin_cols)

    t = table_cls([mixin_cols[name] for name in names], names=names)
    t.meta["C"] = "spam"
    t.meta["comments"] = ["this", "is", "a", "comment"]
    t.meta["history"] = ["first", "second", "third"]

    t.write(filename, format="parquet")

    t2 = table_cls.read(filename, format="parquet", exclude_names=names[0:5])
    assert t.colnames[5:] == t2.colnames


@pytest.mark.parametrize("table_cls", (Table, QTable))
def test_parquet_mixins_read_no_columns(table_cls, tmp_path):
    """Test write all cols at once, and try to read no valid columns."""
    filename = tmp_path / "test_simple.parquet"
    names = sorted(mixin_cols)

    t = table_cls([mixin_cols[name] for name in names], names=names)
    t.meta["C"] = "spam"
    t.meta["comments"] = ["this", "is", "a", "comment"]
    t.meta["history"] = ["first", "second", "third"]

    t.write(filename, format="parquet")

    with pytest.raises(ValueError, match="No include_names specified"):
        t2 = table_cls.read(
            filename,
            format="parquet",
            include_names=["not_a_column", "also_not_a_column"],
        )


@pytest.mark.parametrize("table_cls", (Table, QTable))
def test_parquet_mixins_read_schema(table_cls, tmp_path):
    """Test write all cols at once, and read the schema."""
    filename = tmp_path / "test_simple.parquet"
    names = sorted(mixin_cols)

    t = table_cls([mixin_cols[name] for name in names], names=names)
    t.meta["C"] = "spam"
    t.meta["comments"] = ["this", "is", "a", "comment"]
    t.meta["history"] = ["first", "second", "third"]

    t.write(filename, format="parquet")

    t2 = table_cls.read(filename, format="parquet", schema_only=True)

    assert t2.meta["C"] == "spam"
    assert t2.meta["comments"] == ["this", "is", "a", "comment"]
    assert t2.meta["history"] == ["first", "second", "third"]

    assert t.colnames == t2.colnames

    assert len(t2) == 0


def test_parquet_filter(tmp_path):
    """Test reading a parquet file with a filter."""
    filename = tmp_path / "test_simple.parquet"

    t1 = Table()
    t1["a"] = Column(data=np.arange(100), dtype=np.int32)
    t1["b"] = Column(data=np.arange(100, 0, -1), dtype=np.float64)

    t1.write(filename, overwrite=True)

    t2 = Table.read(filename, filters=[("a", "<", 50)])

    assert t2["a"].max() < 50

    t2 = Table.read(filename, filters=[("b", "<", 50)])

    assert t2["b"].max() < 50


def test_parquet_read_generic(tmp_path):
    """Test reading a generic parquet file."""
    filename = tmp_path / "test_generic.parq"

    t1 = Table()

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        t1.add_column(Column(name=str(dtype), data=np.array(values, dtype=dtype)))

    # Write the table generically via pyarrow.parquet
    names = t1.dtype.names
    type_list = [
        (name, pyarrow.from_numpy_dtype(t1[name].dtype.type)) for name in names
    ]
    schema = pyarrow.schema(type_list)

    _, parquet = get_pyarrow()

    with parquet.ParquetWriter(filename, schema, version="2.4") as writer:
        arrays = [pyarrow.array(t1[name].data) for name in names]
        writer.write_table(pyarrow.Table.from_arrays(arrays, schema=schema))

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        t2 = Table.read(filename)

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        assert np.all(t2[str(dtype)] == values)
        assert t2[str(dtype)].dtype == dtype


@pytest.mark.skipif(not HAS_PANDAS, reason="requires pandas")
def test_parquet_read_pandas(tmp_path):
    """Test reading a pandas parquet file."""
    filename = tmp_path / "test_pandas.parq"

    t1 = Table()

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        t1.add_column(Column(name=str(dtype), data=np.array(values, dtype=dtype)))

    df = t1.to_pandas()
    df.to_parquet(filename, version="2.4")

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        t2 = Table.read(filename)

    for dtype in ALL_DTYPES:
        values = _default_values(dtype)
        assert np.all(t2[str(dtype)] == values)
        assert t2[str(dtype)].dtype == dtype
