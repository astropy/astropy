# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to YAML serialization.
"""

from io import StringIO

import numpy as np
import pytest
from yaml import SafeDumper

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
)
from astropy.coordinates.tests.test_representation import representation_equal
from astropy.io.misc.yaml import AstropyDumper, dump, load, load_all
from astropy.table import QTable, SerializedColumn
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
        np.array([[1, 2], [3, 4]], order="F"),
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


def compare_coord(c, cy):
    assert c.shape == cy.shape
    assert c.frame.name == cy.frame.name

    assert list(c.frame_attributes) == list(cy.frame_attributes)
    for attr in c.frame_attributes:
        assert getattr(c, attr) == getattr(cy, attr)

    assert list(c.representation_component_names) == list(
        cy.representation_component_names
    )
    for name in c.representation_component_names:
        assert np.all(getattr(c, attr) == getattr(cy, attr))


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


def test_yaml_dump_of_object_arrays_fail():
    """Test that dumping and loading object arrays fails."""
    with pytest.raises(TypeError, match="cannot serialize"):
        dump(np.array([1, 2, 3], dtype=object))


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
