# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to YAML serialization.
"""

from io import StringIO

import numpy as np
import pytest
from yaml import SafeDumper

import astropy.coordinates as coords
import astropy.units as u
from astropy.coordinates import (
    Angle,
    CartesianDifferential,
    CartesianRepresentation,
    EarthLocation,
    Latitude,
    Longitude,
    SkyCoord,
    SphericalCosLatDifferential,
    SphericalDifferential,
    SphericalRepresentation,
    UnitSphericalRepresentation,
    frame_transform_graph,
)
from astropy.coordinates.tests.test_representation import representation_equal
from astropy.io.ascii import InconsistentTableError
from astropy.io.misc.yaml import AstropyDumper, dump, load, load_all
from astropy.table import QTable, SerializedColumn, Table
from astropy.time import Time


@pytest.mark.parametrize(
    "c",
    [
        True,
        np.uint8(8),
        np.int16(4),
        np.int32(1),
        np.int64(3),
        np.int64(2**63 - 1),
        2.0,
        np.float64(),
        3 + 4j,
        np.complex64(3 + 4j),
        np.complex128(1.0 - 2**-52 + 1j * (1.0 - 2**-52)),
    ],
)
def test_numpy_types(c):
    cy = load(dump(c))
    assert c == cy


@pytest.mark.parametrize("c", [float("inf"), float("-inf"), np.inf, -np.inf])
def test_astropydumper_represent_float_override(c):
    # AstropyDumper overrides SafeDumper.represent_float, which has multiple
    # branches, not all of which are intended to deviate.
    # Check that the subclass behave as its parent in these cases.
    d1 = SafeDumper(stream=StringIO())
    d2 = AstropyDumper(stream=StringIO())
    assert d2.represent_float(c).value == d1.represent_float(c).value


@pytest.mark.parametrize(
    "c", [u.m, u.m / u.s, u.hPa, u.dimensionless_unscaled, u.Unit("m, (cm, um)")]
)
def test_unit(c):
    cy = load(dump(c))
    if isinstance(c, (u.CompositeUnit, u.StructuredUnit)):
        assert c == cy
    else:
        assert c is cy


@pytest.mark.parametrize("c", [u.Unit("bakers_dozen", 13 * u.one), u.def_unit("magic")])
def test_custom_unit(c):
    s = dump(c)
    with pytest.warns(u.UnitsWarning, match=f"'{c!s}' did not parse") as w:
        cy = load(s)
    assert len(w) == 1
    assert isinstance(cy, u.UnrecognizedUnit)
    assert str(cy) == str(c)

    with u.add_enabled_units(c):
        cy2 = load(s)
        assert cy2 is c


@pytest.mark.parametrize(
    "c",
    [
        Angle("1 2 3", unit="deg"),
        Longitude("1 2 3", unit="deg"),
        Latitude("1 2 3", unit="deg"),
        [[1], [3]] * u.m,
        np.array([[1, 2, 3], [4, 5, 6]], order="F"),
        np.arange(24).reshape(2, 3, 4, order="F"),
        np.array([[1, 2], [3, 4]], order="C"),
        np.array([1, 2, 3, 4])[::2],
        np.array([(1.0, 2), (3.0, 4)], dtype="f8,i4"),  # array with structured dtype.
        np.array((1.0, 2), dtype="f8,i4"),  # array scalar with structured dtype.
        np.array((1.0, 2), dtype="f8,i4")[()],  # numpy void.
        np.array((1.0, 2.0), dtype="f8,f8") * u.s,  # Quantity structured scalar.
        [
            ((1.0, 2.0, 3.0), (4.0, 5.0, 6.0)),  # Quantity with structured unit.
            ((11.0, 12.0, 13.0), (14.0, 15.0, 16.0)),
        ]
        * u.Unit("m, m/s"),
        np.array(
            [
                ((1.0, 2.0, 3.0), (4.0, 5.0, 6.0)),
                ((11.0, 12.0, 13.0), (14.0, 15.0, 16.0)),
            ],
            dtype=[("p", "3f8"), ("v", "3f8")],
        )
        * u.Unit("m, m/s"),
    ],
)
def test_ndarray_subclasses(c):
    cy = load(dump(c))

    assert np.all(c == cy)
    assert c.shape == cy.shape
    assert c.dtype == cy.dtype
    assert type(c) is type(cy)

    cc = "C_CONTIGUOUS"
    fc = "F_CONTIGUOUS"
    if c.flags[cc] or c.flags[fc]:
        assert c.flags[cc] == cy.flags[cc]
        assert c.flags[fc] == cy.flags[fc]
    else:
        # Original was not contiguous but round-trip version
        # should be c-contig.
        assert cy.flags[cc]

    if hasattr(c, "unit"):
        assert c.unit == cy.unit


def test_masked_constant():
    values = {"masked": np.ma.masked, "nested": [np.ma.masked], "none": None}
    result = load(dump(values))
    assert result["masked"] is np.ma.masked
    assert result["nested"][0] is np.ma.masked
    assert result["none"] is None


@pytest.mark.parametrize(
    "obj",
    [
        np.ma.array([1, 2, 3]),
        np.ma.array([1, 2, 3], mask=False),
        np.ma.array(
            [[1, 2], [3, 4]], mask=[[False, True], [True, False]], fill_value=-99
        ),
        np.ma.array(
            np.arange(6).reshape(2, 3, order="F"),
            mask=[[False, True, False], [True, False, True]],
        ),
        np.ma.array(3.5, mask=True),
        np.ma.array([], dtype=np.int16),
        np.ma.array([1.0, 2.0], dtype=np.float16),
        np.ma.array([1, 2], dtype=np.int32, mask=[False, True], fill_value=999999),
        np.ma.array(["a", "bc"], mask=[True, False], fill_value="--"),
        np.ma.array(
            ["2000-01-01", "2001-01-01"],
            dtype="datetime64[D]",
            mask=[False, True],
            fill_value=np.datetime64("1999-01-01"),
        ),
        np.ma.array(
            [(1, 2.5), (3, 4.5)],
            dtype=[("a", "i4"), ("b", "f8")],
            mask=[(False, True), (True, False)],
            fill_value=(7, 8.5),
        ),
        np.ma.array(np.arange(12), mask=np.arange(12) % 3 == 0)[::2],
    ],
)
def test_masked_array(obj):
    result = load(dump(obj))
    assert type(result) is np.ma.MaskedArray
    assert result.shape == obj.shape
    assert result.dtype == obj.dtype
    np.testing.assert_array_equal(result.data, obj.data)
    np.testing.assert_array_equal(result.mask, obj.mask)
    assert (result.mask is np.ma.nomask) == (obj.mask is np.ma.nomask)
    np.testing.assert_array_equal(result.fill_value, obj.fill_value)
    assert np.asarray(result.fill_value).dtype == np.asarray(obj.fill_value).dtype
    assert result.hardmask == obj.hardmask


def test_masked_array_default_fill_survives_astype():
    obj = np.ma.array([1], dtype=np.int16)
    result = load(dump(obj))
    assert result.astype(float).fill_value == obj.astype(float).fill_value


def test_masked_array_initialized_fill_value():
    obj = np.ma.array([1, 2], dtype=np.int16, mask=[False, True])
    obj.get_fill_value()
    obj.fill_value = -9

    result = load(dump(obj))

    assert result.fill_value == obj.fill_value
    assert np.asarray(result.fill_value).dtype == np.asarray(obj.fill_value).dtype


@pytest.mark.parametrize("hard_mask", [False, True])
def test_masked_array_assignment(hard_mask):
    obj = np.ma.array([1, 2], mask=[False, True], hard_mask=hard_mask)
    result = load(dump(obj))
    result[0] = 3
    result[1] = 4
    assert result[0] == 3
    assert result.hardmask == hard_mask
    assert result.mask[1] == hard_mask
    assert result.data[1] == (2 if hard_mask else 4)


@pytest.mark.parametrize("table_cls", [Table, QTable])
def test_ecsv_masked_metadata(table_cls):
    table = table_cls({"a": [1, 2]})
    values = np.ma.array([1.0, 2.0], mask=[False, True], fill_value=-99.0)
    table.meta = {"missing": np.ma.masked, "values": values}
    table["a"].info.meta = {"missing": np.ma.masked, "values": values}
    out = StringIO()
    table.write(out, format="ascii.ecsv")
    result = table_cls.read(out.getvalue(), format="ascii.ecsv")

    for metadata in (result.meta, result["a"].info.meta):
        assert metadata["missing"] is np.ma.masked
        np.testing.assert_array_equal(metadata["values"].data, values.data)
        np.testing.assert_array_equal(metadata["values"].mask, values.mask)
        assert metadata["values"].fill_value == values.fill_value


@pytest.mark.parametrize(
    "dtype",
    [
        np.dtype([("a", "i1"), ("b", "i8")], align=True),
        np.dtype(
            {
                "names": ["a", "b"],
                "formats": ["<i2", ">f8"],
                "offsets": [2, 16],
                "titles": [None, "sample"],
                "itemsize": 32,
            }
        ),
        np.dtype(
            [
                ("record", np.dtype([("a", "i1"), ("b", "i8")])),
                ("samples", "f8", (2,)),
            ],
            align=True,
        ),
        np.dtype(
            [("records", np.dtype([("a", "i1"), ("b", "i8")]), (2,))],
            align=True,
        ),
        np.dtype(
            {
                "names": ["a", "b"],
                "formats": ["i4", "i2"],
                "offsets": [0, 0],
                "itemsize": 4,
            }
        ),
        np.dtype(
            {
                "names": ["b", "a"],
                "formats": ["i4", "i2"],
                "offsets": [4, 0],
                "itemsize": 8,
            }
        ),
    ],
)
@pytest.mark.parametrize("kind", ["array", "scalar", "masked"])
def test_structured_dtype_layout(dtype, kind):
    data = np.ones(2, dtype=dtype)
    data[1] = 2
    if kind == "scalar":
        obj = data[0]
    elif kind == "masked":
        mask = np.zeros(data.shape, dtype=np.ma.make_mask_descr(dtype))
        mask[dtype.names[-1]][0] = True
        fill_value = np.zeros((), dtype=dtype)
        fill_value[()] = -9
        obj = np.ma.array(data, mask=mask, fill_value=fill_value)
    else:
        obj = data

    result = load(dump(obj))

    assert result.dtype == obj.dtype
    assert result.dtype.isalignedstruct == obj.dtype.isalignedstruct
    assert result.shape == obj.shape
    if kind == "masked":
        np.testing.assert_array_equal(result.data, obj.data)
        np.testing.assert_array_equal(result.mask, obj.mask)
        assert result.mask.dtype == obj.mask.dtype
        np.testing.assert_array_equal(result.fill_value, obj.fill_value)
        assert result.fill_value.dtype == obj.fill_value.dtype
    else:
        np.testing.assert_array_equal(result, obj)


@pytest.mark.parametrize("tag", ["ndarray", "void"])
def test_legacy_structured_dtype_yaml(tag):
    # The descriptor-list representation written before field layouts were
    # stored explicitly must remain readable.
    serialized = f"""!numpy.{tag}
buffer: !!binary |
  QVFBQ0FBPT0=
dtype:
- !!python/tuple [a, <i2]
- !!python/tuple [b, <i2]
"""
    if tag == "ndarray":
        serialized += "shape: !!python/tuple [1]\norder: C\n"

    result = load(serialized)

    expected = np.array([(1, 2)], dtype=[("a", "<i2"), ("b", "<i2")])
    if tag == "void":
        expected = expected[0]
    assert result.dtype == expected.dtype
    np.testing.assert_array_equal(result, expected)


def compare_coord(c, cy, has_data=True):
    assert c.shape == cy.shape
    if hasattr(c, "frame"):
        assert c.frame.name == cy.frame.name
    else:
        assert c.name == cy.name

    assert list(c.frame_attributes) == list(cy.frame_attributes)
    for attr in c.frame_attributes:
        assert getattr(c, attr) == getattr(cy, attr)

    assert list(c.representation_component_names) == list(
        cy.representation_component_names
    )
    if has_data:
        for attr in c.representation_component_names:
            assert np.all(getattr(c, attr) == getattr(cy, attr))


def test_skycoord_galactocentric():
    """Integration test for #18526 from https://stackoverflow.com/q/79735659/323631"""
    sc = SkyCoord(
        x=1 * u.kpc,
        y=2 * u.kpc,
        z=3 * u.kpc,
        v_x=4 * u.km / u.s,
        v_y=5 * u.km / u.s,
        v_z=6 * u.km / u.s,
        frame="galactocentric",
    )
    scr = load(dump(sc))
    compare_coord(sc, scr)


@pytest.mark.parametrize("frame", ["fk4", "altaz"])
def test_skycoord(frame):
    c = SkyCoord(
        [[1, 2], [3, 4]],
        [[5, 6], [7, 8]],
        unit="deg",
        frame=frame,
        obstime=Time("2016-01-02"),
        location=EarthLocation(1000, 2000, 3000, unit=u.km),
    )
    cy = load(dump(c))
    compare_coord(c, cy)


@pytest.mark.parametrize(
    "rep",
    [
        CartesianRepresentation(1 * u.m, 2.0 * u.m, 3.0 * u.m),
        SphericalRepresentation(
            [[1, 2], [3, 4]] * u.deg, [[5, 6], [7, 8]] * u.deg, 10 * u.pc
        ),
        UnitSphericalRepresentation(0 * u.deg, 10 * u.deg),
        SphericalCosLatDifferential(
            [[1.0], [2.0]] * u.mas / u.yr,
            [4.0, 5.0] * u.mas / u.yr,
            [[[10]], [[20]]] * u.km / u.s,
        ),
        CartesianDifferential([10, 20, 30] * u.km / u.s),
        CartesianRepresentation(
            [1, 2, 3] * u.m,
            differentials=CartesianDifferential([10, 20, 30] * u.km / u.s),
        ),
        SphericalRepresentation(
            [[1, 2], [3, 4]] * u.deg,
            [[5, 6], [7, 8]] * u.deg,
            10 * u.pc,
            differentials={
                "s": SphericalDifferential(
                    [[0.0, 1.0], [2.0, 3.0]] * u.mas / u.yr,
                    [[4.0, 5.0], [6.0, 7.0]] * u.mas / u.yr,
                    10 * u.km / u.s,
                )
            },
        ),
    ],
)
def test_representations(rep):
    rrep = load(dump(rep))
    assert np.all(representation_equal(rrep, rep))


# Note: consistent test order (e.g. sorting) needed for parallel testing.
frame_clss = sorted(
    [cls for cls in frame_transform_graph.frame_set if cls.__name__ in dir(coords)],
    key=lambda x: x.__name__,
)


@pytest.mark.parametrize("frame_cls", frame_clss)
def test_frames(frame_cls):
    """Test that bare instances of built-in frames round-trip through YAML

    This excludes dynamically-created SkyOffset frames like ``abc.SkyOffsetAltAz``.
    This one in particular fails tests since it tries to construct an ``AltAz``
    instance. Possibly related to https://github.com/astropy/astropy/issues/10157.
    """
    frame = frame_cls()
    frame_dump = dump(frame)
    frame_rt = load(frame_dump)
    compare_coord(frame, frame_rt, has_data=False)


def _get_time():
    t = Time(
        [[1], [2]], format="cxcsec", location=EarthLocation(1000, 2000, 3000, unit=u.km)
    )
    t.format = "iso"
    t.precision = 5
    t.delta_ut1_utc = np.array([[3.0], [4.0]])
    t.delta_tdb_tt = np.array([[5.0], [6.0]])
    t.out_subfmt = "date_hm"

    return t


def compare_time(t, ty):
    assert type(t) is type(ty)
    assert np.all(t == ty)
    for attr in (
        "shape",
        "jd1",
        "jd2",
        "format",
        "scale",
        "precision",
        "in_subfmt",
        "out_subfmt",
        "location",
        "delta_ut1_utc",
        "delta_tdb_tt",
    ):
        assert np.all(getattr(t, attr) == getattr(ty, attr))


def test_time():
    t = _get_time()
    ty = load(dump(t))
    compare_time(t, ty)


def test_timedelta():
    t = _get_time()
    dt = t - t + 0.1234556 * u.s
    dty = load(dump(dt))

    assert type(dt) is type(dty)
    for attr in ("shape", "jd1", "jd2", "format", "scale"):
        assert np.all(getattr(dt, attr) == getattr(dty, attr))


def test_serialized_column():
    sc = SerializedColumn({"name": "hello", "other": 1, "other2": 2.0})
    scy = load(dump(sc))

    assert sc == scy


def test_load_all():
    t = _get_time()
    unit = u.m / u.s
    c = SkyCoord(
        [[1, 2], [3, 4]],
        [[5, 6], [7, 8]],
        unit="deg",
        frame="fk4",
        obstime=Time("2016-01-02"),
        location=EarthLocation(1000, 2000, 3000, unit=u.km),
    )

    # Make a multi-document stream
    out = "---\n" + dump(t) + "---\n" + dump(unit) + "---\n" + dump(c)

    ty, unity, cy = list(load_all(out))

    compare_time(t, ty)
    compare_coord(c, cy)
    assert unity == unit


def test_ecsv_astropy_objects_in_meta():
    """
    Test that astropy core objects in ``meta`` are serialized.
    """
    t = QTable([[1, 2] * u.m, [4, 5]], names=["a", "b"])
    tm = _get_time()
    c = SkyCoord(
        [[1, 2], [3, 4]],
        [[5, 6], [7, 8]],
        unit="deg",
        frame="fk4",
        obstime=Time("2016-01-02"),
        location=EarthLocation(1000, 2000, 3000, unit=u.km),
    )
    unit = u.m / u.s

    t.meta = {"tm": tm, "c": c, "unit": unit}
    out = StringIO()
    t.write(out, format="ascii.ecsv")
    t2 = QTable.read(out.getvalue(), format="ascii.ecsv")

    compare_time(tm, t2.meta["tm"])
    compare_coord(c, t2.meta["c"])
    assert t2.meta["unit"] == unit


@pytest.mark.parametrize(
    "obj",
    [
        np.array([1, 2, 3], dtype=object),
        np.array([(1, 2)], dtype=[("a", "i4"), ("b", "O")]),
        np.array((1, 2), dtype=[("a", "i4"), ("b", "O")])[()],  # numpy void.
    ],
)
def test_yaml_dump_of_object_arrays_fail(obj):
    """Test that dumping arrays containing objects fails."""
    with pytest.raises(TypeError, match="cannot serialize"):
        dump(obj)


def test_yaml_load_of_object_arrays_fail():
    """Test that dumping and loading object arrays fails.

    The string to load was obtained by suppressing the exception and dumping
    ``np.array([1, 2, 3], dtype=object)`` to a yaml file.
    """
    with pytest.raises(TypeError, match="cannot load numpy array"):
        load(
            """!numpy.ndarray
            buffer: !!binary |
              WndBQUFISUFBQUJwQUFBQQ==
            dtype: object
            order: C
            shape: !!python/tuple [3]"""
        )


# A plain and a nested structured object dtype for each of the ndarray and void
# constructors. Each buffer holds an invalid pointer, so if a guard is missing, these
# can crash the interpreter instead of failing the test.
BAD_DTYPE_YAML = [
    "!numpy.void {dtype: O, buffer: QUFBQUFBQUE=}",
    "!numpy.void {dtype: [!!python/tuple [a, [!!python/tuple [b, O]]]], buffer: QUFBQUFBQUE=}",
    "!numpy.ndarray {dtype: O, shape: [1], buffer: QUFBQUFBQUE=, order: C}",
    (
        "!numpy.ndarray {dtype: [!!python/tuple [a, [!!python/tuple [b, O]]]], shape: [1], "
        "buffer: QUFBQUFBQUE=, order: C}"
    ),
]


@pytest.mark.parametrize("text", BAD_DTYPE_YAML)
def test_yaml_load_refuses_object_dtype(text):
    """Object arrays must not be built from a buffer supplied by the YAML."""
    with pytest.raises(TypeError, match="cannot load numpy array with object dtype"):
        load(text)


def test_ecsv_read_refuses_object_dtype_in_meta():
    """Table metadata is the route by which an untrusted file reaches the YAML loader."""
    ecsv = """\
# %ECSV 1.0
# ---
# datatype:
# - {name: a, datatype: int64}
# meta: {bad: !numpy.void {dtype: O, buffer: QUFBQUFBQUE=}}
# schema: astropy-2.0
a
1
"""
    with pytest.raises(InconsistentTableError, match="unable to parse yaml") as exc:
        Table.read(ecsv, format="ascii.ecsv")
    assert isinstance(exc.value.__cause__.__cause__, TypeError)


@pytest.mark.parametrize(
    "obj",
    [
        np.array([1.5, 2.5]).reshape(2, 1),
        np.array([1.5, 2.5, 3.5, 4.5]).reshape(2, 2, order="F"),
        np.array(["a", "bc"]),
        np.array([(1, 2.5)], dtype=[("a", "i4"), ("b", "f8")]),
        np.array(
            [(1, (2.5, 3))], dtype=[("a", "i4"), ("b", [("c", "f8"), ("d", "i2")])]
        ),
        np.array([(1, [2, 3])], dtype=[("a", "i4"), ("b", "f8", (2,))]),
        np.array((1, 2.5), dtype=[("a", "i4"), ("b", "f8")])[()],  # numpy void.
    ],
)
def test_yaml_round_trip_without_objects(obj):
    """Dtypes with no objects, in particular structured ones, must still load."""
    result = load(dump(obj))
    assert result.dtype == obj.dtype
    assert result.shape == obj.shape
    assert result.tobytes() == obj.tobytes()
