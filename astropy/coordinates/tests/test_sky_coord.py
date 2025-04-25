# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for the SkyCoord class.  Note that there are also SkyCoord tests in
test_api_ape5.py
"""

import copy
from copy import deepcopy

import numpy as np
import numpy.testing as npt
import pytest
from erfa import ErfaWarning

from astropy import units as u
from astropy.coordinates import (
    FK4,
    FK5,
    GCRS,
    ICRS,
    AltAz,
    Angle,
    Attribute,
    BaseCoordinateFrame,
    CartesianRepresentation,
    EarthLocation,
    Galactic,
    Latitude,
    RepresentationMapping,
    SkyCoord,
    SphericalRepresentation,
    UnitSphericalRepresentation,
    frame_transform_graph,
)
from astropy.coordinates.representation import (
    DUPLICATE_REPRESENTATIONS,
    REPRESENTATION_CLASSES,
)
from astropy.coordinates.tests.helper import skycoord_equal
from astropy.coordinates.transformations import FunctionTransform
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time
from astropy.units import allclose as quantity_allclose
from astropy.utils.compat import NUMPY_LT_2_0
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.wcs import WCS

RA = 1.0 * u.deg
DEC = 2.0 * u.deg
C_ICRS = ICRS(RA, DEC)
C_FK5 = C_ICRS.transform_to(FK5())
J2001 = Time("J2001")


def allclose(a, b, rtol=0.0, atol=None):
    if atol is None:
        atol = 1.0e-8 * getattr(a, "unit", 1.0)
    return quantity_allclose(a, b, rtol, atol)


def setup_function(func):
    func.REPRESENTATION_CLASSES_ORIG = deepcopy(REPRESENTATION_CLASSES)
    func.DUPLICATE_REPRESENTATIONS_ORIG = deepcopy(DUPLICATE_REPRESENTATIONS)


def teardown_function(func):
    REPRESENTATION_CLASSES.clear()
    REPRESENTATION_CLASSES.update(func.REPRESENTATION_CLASSES_ORIG)
    DUPLICATE_REPRESENTATIONS.clear()
    DUPLICATE_REPRESENTATIONS.update(func.DUPLICATE_REPRESENTATIONS_ORIG)


def test_is_transformable_to_str_input():
    """Test method ``is_transformable_to`` with string input.

    The only difference from the frame method of the same name is that
    strings are allowed. As the frame tests cover ``is_transform_to``, here
    we only test the added string option.

    """
    # make example SkyCoord
    c = SkyCoord(90 * u.deg, -11 * u.deg)

    # iterate through some frames, checking consistency
    names = frame_transform_graph.get_names()
    for name in names:
        frame = frame_transform_graph.lookup_name(name)()
        assert c.is_transformable_to(name) == c.is_transformable_to(frame)


def test_transform_to():
    for frame in (
        FK5(),
        FK5(equinox=Time("J1975.0")),
        FK4(),
        FK4(equinox=Time("J1975.0")),
        SkyCoord(RA, DEC, frame="fk4", equinox="J1980"),
    ):
        c_frame = C_ICRS.transform_to(frame)
        s_icrs = SkyCoord(RA, DEC, frame="icrs")
        s_frame = s_icrs.transform_to(frame)
        assert allclose(c_frame.ra, s_frame.ra)
        assert allclose(c_frame.dec, s_frame.dec)
        assert allclose(c_frame.distance, s_frame.distance)


# set up for parametrized test
rt_sets = []
rt_frames = [ICRS, FK4, FK5, Galactic]
for rt_frame0 in rt_frames:
    for rt_frame1 in rt_frames:
        for equinox0 in (None, "J1975.0"):
            for obstime0 in (None, "J1980.0"):
                for equinox1 in (None, "J1975.0"):
                    for obstime1 in (None, "J1980.0"):
                        rt_sets.append(
                            (
                                rt_frame0,
                                rt_frame1,
                                equinox0,
                                equinox1,
                                obstime0,
                                obstime1,
                            )
                        )
rt_args = ("frame0", "frame1", "equinox0", "equinox1", "obstime0", "obstime1")


@pytest.mark.parametrize(rt_args, rt_sets)
def test_round_tripping(frame0, frame1, equinox0, equinox1, obstime0, obstime1):
    """
    Test round tripping out and back using transform_to in every combination.
    """
    attrs0 = {"equinox": equinox0, "obstime": obstime0}
    attrs1 = {"equinox": equinox1, "obstime": obstime1}

    # Remove None values
    attrs0 = {k: v for k, v in attrs0.items() if v is not None}
    attrs1 = {k: v for k, v in attrs1.items() if v is not None}

    # Go out and back
    sc = SkyCoord(RA, DEC, frame=frame0, **attrs0)

    # Keep only frame attributes for frame1
    attrs1 = {
        attr: val for attr, val in attrs1.items() if attr in frame1.frame_attributes
    }
    sc2 = sc.transform_to(frame1(**attrs1))

    # When coming back only keep frame0 attributes for transform_to
    attrs0 = {
        attr: val for attr, val in attrs0.items() if attr in frame0.frame_attributes
    }
    # also, if any are None, fill in with defaults
    for attrnm in frame0.frame_attributes:
        if attrs0.get(attrnm, None) is None:
            if attrnm == "obstime" and frame0.get_frame_attr_defaults()[attrnm] is None:
                if "equinox" in attrs0:
                    attrs0[attrnm] = attrs0["equinox"]
            else:
                attrs0[attrnm] = frame0.get_frame_attr_defaults()[attrnm]
    sc_rt = sc2.transform_to(frame0(**attrs0))

    if frame0 is Galactic:
        assert allclose(sc.l, sc_rt.l)
        assert allclose(sc.b, sc_rt.b)
    else:
        assert allclose(sc.ra, sc_rt.ra)
        assert allclose(sc.dec, sc_rt.dec)
    if equinox0:
        assert type(sc.equinox) is Time and sc.equinox == sc_rt.equinox
    if obstime0:
        assert type(sc.obstime) is Time and sc.obstime == sc_rt.obstime


def test_coord_init_string():
    """
    Spherical or Cartesian representation input coordinates.
    """
    sc = SkyCoord("1d 2d")
    assert allclose(sc.ra, 1 * u.deg)
    assert allclose(sc.dec, 2 * u.deg)

    sc = SkyCoord("1d", "2d")
    assert allclose(sc.ra, 1 * u.deg)
    assert allclose(sc.dec, 2 * u.deg)

    sc = SkyCoord("1°2′3″", "2°3′4″")
    assert allclose(sc.ra, Angle("1°2′3″"))
    assert allclose(sc.dec, Angle("2°3′4″"))

    sc = SkyCoord("1°2′3″ 2°3′4″")
    assert allclose(sc.ra, Angle("1°2′3″"))
    assert allclose(sc.dec, Angle("2°3′4″"))

    with pytest.raises(ValueError) as err:
        SkyCoord("1d 2d 3d")
    assert "Cannot parse first argument data" in str(err.value)

    sc1 = SkyCoord("8 00 00 +5 00 00.0", unit=(u.hour, u.deg), frame="icrs")
    assert isinstance(sc1, SkyCoord)
    assert allclose(sc1.ra, Angle(120 * u.deg))
    assert allclose(sc1.dec, Angle(5 * u.deg))

    sc11 = SkyCoord("8h00m00s+5d00m00.0s", unit=(u.hour, u.deg), frame="icrs")
    assert isinstance(sc11, SkyCoord)
    assert allclose(sc1.ra, Angle(120 * u.deg))
    assert allclose(sc1.dec, Angle(5 * u.deg))

    sc2 = SkyCoord("8 00 -5 00 00.0", unit=(u.hour, u.deg), frame="icrs")
    assert isinstance(sc2, SkyCoord)
    assert allclose(sc2.ra, Angle(120 * u.deg))
    assert allclose(sc2.dec, Angle(-5 * u.deg))

    sc3 = SkyCoord("8 00 -5 00.6", unit=(u.hour, u.deg), frame="icrs")
    assert isinstance(sc3, SkyCoord)
    assert allclose(sc3.ra, Angle(120 * u.deg))
    assert allclose(sc3.dec, Angle(-5.01 * u.deg))

    sc4 = SkyCoord("J080000.00-050036.00", unit=(u.hour, u.deg), frame="icrs")
    assert isinstance(sc4, SkyCoord)
    assert allclose(sc4.ra, Angle(120 * u.deg))
    assert allclose(sc4.dec, Angle(-5.01 * u.deg))

    sc41 = SkyCoord("J080000+050036", unit=(u.hour, u.deg), frame="icrs")
    assert isinstance(sc41, SkyCoord)
    assert allclose(sc41.ra, Angle(120 * u.deg))
    assert allclose(sc41.dec, Angle(+5.01 * u.deg))

    sc5 = SkyCoord("8h00.6m -5d00.6m", unit=(u.hour, u.deg), frame="icrs")
    assert isinstance(sc5, SkyCoord)
    assert allclose(sc5.ra, Angle(120.15 * u.deg))
    assert allclose(sc5.dec, Angle(-5.01 * u.deg))

    sc6 = SkyCoord("8h00.6m -5d00.6m", unit=(u.hour, u.deg), frame="fk4")
    assert isinstance(sc6, SkyCoord)
    assert allclose(sc6.ra, Angle(120.15 * u.deg))
    assert allclose(sc6.dec, Angle(-5.01 * u.deg))

    sc61 = SkyCoord("8h00.6m-5d00.6m", unit=(u.hour, u.deg), frame="fk4")
    assert isinstance(sc61, SkyCoord)
    assert allclose(sc6.ra, Angle(120.15 * u.deg))
    assert allclose(sc6.dec, Angle(-5.01 * u.deg))

    sc61 = SkyCoord("8h00.6-5d00.6", unit=(u.hour, u.deg), frame="fk4")
    assert isinstance(sc61, SkyCoord)
    assert allclose(sc6.ra, Angle(120.15 * u.deg))
    assert allclose(sc6.dec, Angle(-5.01 * u.deg))

    sc7 = SkyCoord("J1874221.60+122421.6", unit=u.deg)
    assert isinstance(sc7, SkyCoord)
    assert allclose(sc7.ra, Angle(187.706 * u.deg))
    assert allclose(sc7.dec, Angle(12.406 * u.deg))

    with pytest.raises(ValueError):
        SkyCoord("8 00 -5 00.6", unit=(u.deg, u.deg), frame="galactic")


def test_coord_init_unit():
    """
    Test variations of the unit keyword.
    """
    for unit in (
        "deg",
        "deg,deg",
        " deg , deg ",
        u.deg,
        (u.deg, u.deg),
        np.array(["deg", "deg"]),
    ):
        sc = SkyCoord(1, 2, unit=unit)
        assert allclose(sc.ra, Angle(1 * u.deg))
        assert allclose(sc.dec, Angle(2 * u.deg))

    for unit in (
        "hourangle",
        "hourangle,hourangle",
        " hourangle , hourangle ",
        u.hourangle,
        [u.hourangle, u.hourangle],
    ):
        sc = SkyCoord(1, 2, unit=unit)
        assert allclose(sc.ra, Angle(15 * u.deg))
        assert allclose(sc.dec, Angle(30 * u.deg))

    for unit in ("hourangle,deg", (u.hourangle, u.deg)):
        sc = SkyCoord(1, 2, unit=unit)
        assert allclose(sc.ra, Angle(15 * u.deg))
        assert allclose(sc.dec, Angle(2 * u.deg))

    for unit in ("deg,deg,deg,deg", [u.deg, u.deg, u.deg, u.deg], None):
        with pytest.raises(ValueError) as err:
            SkyCoord(1, 2, unit=unit)
        assert "Unit keyword must have one to three unit values" in str(err.value)

    for unit in ("m", (u.m, u.deg), ""):
        with pytest.raises(u.UnitsError) as err:
            SkyCoord(1, 2, unit=unit)


def test_coord_init_list():
    """
    Spherical or Cartesian representation input coordinates.
    """
    sc = SkyCoord(
        [("1d", "2d"), (1 * u.deg, 2 * u.deg), "1d 2d", ("1°", "2°"), "1° 2°"],
        unit="deg",
    )
    assert allclose(sc.ra, Angle("1d"))
    assert allclose(sc.dec, Angle("2d"))

    with pytest.raises(ValueError) as err:
        SkyCoord(["1d 2d 3d"])
    assert "Cannot parse first argument data" in str(err.value)

    with pytest.raises(ValueError) as err:
        SkyCoord([("1d", "2d", "3d")])
    assert "Cannot parse first argument data" in str(err.value)

    sc = SkyCoord([1 * u.deg, 1 * u.deg], [2 * u.deg, 2 * u.deg])
    assert allclose(sc.ra, Angle("1d"))
    assert allclose(sc.dec, Angle("2d"))

    with pytest.raises(
        ValueError,
        match="One or more elements of input sequence does not have a length",
    ):
        SkyCoord([1 * u.deg, 2 * u.deg])  # this list is taken as RA w/ missing dec


def test_coord_init_array():
    """
    Input in the form of a list array or numpy array
    """
    for a in (["1 2", "3 4"], [["1", "2"], ["3", "4"]], [[1, 2], [3, 4]]):
        sc = SkyCoord(a, unit="deg")
        assert allclose(sc.ra - [1, 3] * u.deg, 0 * u.deg)
        assert allclose(sc.dec - [2, 4] * u.deg, 0 * u.deg)

        sc = SkyCoord(np.array(a), unit="deg")
        assert allclose(sc.ra - [1, 3] * u.deg, 0 * u.deg)
        assert allclose(sc.dec - [2, 4] * u.deg, 0 * u.deg)


def test_coord_init_representation():
    """
    Spherical or Cartesian representation input coordinates.
    """
    coord = SphericalRepresentation(lon=8 * u.deg, lat=5 * u.deg, distance=1 * u.kpc)
    sc = SkyCoord(coord, frame="icrs")
    assert allclose(sc.ra, coord.lon)
    assert allclose(sc.dec, coord.lat)
    assert allclose(sc.distance, coord.distance)

    with pytest.raises(ValueError) as err:
        SkyCoord(coord, frame="icrs", ra="1d")
    assert "conflicts with keyword argument 'ra'" in str(err.value)

    coord = CartesianRepresentation(1 * u.one, 2 * u.one, 3 * u.one)
    sc = SkyCoord(coord, frame="icrs")
    sc_cart = sc.represent_as(CartesianRepresentation)
    assert allclose(sc_cart.x, 1.0)
    assert allclose(sc_cart.y, 2.0)
    assert allclose(sc_cart.z, 3.0)


def test_frame_init():
    """
    Different ways of providing the frame.
    """

    sc = SkyCoord(RA, DEC, frame="icrs")
    assert sc.frame.name == "icrs"

    sc = SkyCoord(RA, DEC, frame=ICRS)
    assert sc.frame.name == "icrs"

    sc = SkyCoord(sc)
    assert sc.frame.name == "icrs"

    sc = SkyCoord(C_ICRS)
    assert sc.frame.name == "icrs"

    SkyCoord(C_ICRS, frame="icrs")
    assert sc.frame.name == "icrs"

    with pytest.raises(ValueError) as err:
        SkyCoord(C_ICRS, frame="galactic")
    assert "Cannot override frame=" in str(err.value)


def test_equal():
    obstime = "B1955"
    sc1 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg, obstime=obstime)
    sc2 = SkyCoord([1, 20] * u.deg, [3, 4] * u.deg, obstime=obstime)

    # Compare arrays and scalars
    eq = sc1 == sc2
    ne = sc1 != sc2
    assert np.all(eq == [True, False])
    assert np.all(ne == [False, True])
    v = sc1[0] == sc2[0]
    assert isinstance(v, (bool, np.bool_))
    assert v
    v = sc1[0] != sc2[0]
    assert isinstance(v, (bool, np.bool_))
    assert not v

    # Broadcasting
    eq = sc1[0] == sc2
    ne = sc1[0] != sc2
    assert np.all(eq == [True, False])
    assert np.all(ne == [False, True])

    # With diff only in velocity
    sc1 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg, radial_velocity=[1, 2] * u.km / u.s)
    sc2 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg, radial_velocity=[1, 20] * u.km / u.s)

    eq = sc1 == sc2
    ne = sc1 != sc2
    assert np.all(eq == [True, False])
    assert np.all(ne == [False, True])
    v = sc1[0] == sc2[0]
    assert isinstance(v, (bool, np.bool_))
    assert v
    v = sc1[0] != sc2[0]
    assert isinstance(v, (bool, np.bool_))
    assert not v


def test_equal_different_type():
    sc1 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg, obstime="B1955")
    # Test equals and not equals operators against different types
    assert sc1 != "a string"
    assert not (sc1 == "a string")


def test_equal_exceptions():
    sc1 = SkyCoord(1 * u.deg, 2 * u.deg, obstime="B1955")
    sc2 = SkyCoord(1 * u.deg, 2 * u.deg)
    with pytest.raises(
        ValueError,
        match=(
            "cannot compare: extra frame attribute 'obstime' is not equivalent"
            r" \(perhaps compare the frames directly to avoid this exception\)"
        ),
    ):
        sc1 == sc2  # noqa: B015
    # Note that this exception is the only one raised directly in SkyCoord.
    # All others come from lower-level classes and are tested in test_frames.py.


def test_attr_inheritance():
    """
    When initializing from an existing coord the representation attrs like
    equinox should be inherited to the SkyCoord.  If there is a conflict
    then raise an exception.
    """
    sc = SkyCoord(1, 2, frame="icrs", unit="deg", equinox="J1999", obstime="J2001")
    sc2 = SkyCoord(sc)
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc2 = SkyCoord(sc.frame)  # Doesn't have equinox there so we get FK4 defaults
    assert sc2.equinox != sc.equinox
    assert sc2.obstime != sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc = SkyCoord(1, 2, frame="fk4", unit="deg", equinox="J1999", obstime="J2001")
    sc2 = SkyCoord(sc)
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc2 = SkyCoord(sc.frame)  # sc.frame has equinox, obstime
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)


@pytest.mark.parametrize("frame", ["fk4", "fk5", "icrs"])
def test_setitem_no_velocity(frame):
    """Test different flavors of item setting for a SkyCoord without a velocity
    for different frames.  Include a frame attribute that is sometimes an
    actual frame attribute and sometimes an extra frame attribute.
    """
    sc0 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg, obstime="B1955", frame=frame)
    sc2 = SkyCoord([10, 20] * u.deg, [30, 40] * u.deg, obstime="B1955", frame=frame)

    sc1 = sc0.copy()
    sc1[1] = sc2[0]
    assert np.allclose(sc1.ra.to_value(u.deg), [1, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [3, 30])
    assert sc1.obstime == Time("B1955")
    assert sc1.frame.name == frame

    sc1 = sc0.copy()
    sc1[:] = sc2[0]
    assert np.allclose(sc1.ra.to_value(u.deg), [10, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [30, 30])

    sc1 = sc0.copy()
    sc1[:] = sc2[:]
    assert np.allclose(sc1.ra.to_value(u.deg), [10, 20])
    assert np.allclose(sc1.dec.to_value(u.deg), [30, 40])

    sc1 = sc0.copy()
    sc1[[1, 0]] = sc2[:]
    assert np.allclose(sc1.ra.to_value(u.deg), [20, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [40, 30])


def test_setitem_initially_broadcast():
    sc = SkyCoord(np.ones((2, 1)) * u.deg, np.ones((1, 3)) * u.deg)
    sc[1, 1] = SkyCoord(0 * u.deg, 0 * u.deg)
    expected = np.ones((2, 3)) * u.deg
    expected[1, 1] = 0.0
    assert np.all(sc.ra == expected)
    assert np.all(sc.dec == expected)


def test_setitem_velocities():
    """Test different flavors of item setting for a SkyCoord with a velocity."""
    sc0 = SkyCoord(
        [1, 2] * u.deg,
        [3, 4] * u.deg,
        radial_velocity=[1, 2] * u.km / u.s,
        obstime="B1950",
        frame="fk4",
    )
    sc2 = SkyCoord(
        [10, 20] * u.deg,
        [30, 40] * u.deg,
        radial_velocity=[10, 20] * u.km / u.s,
        obstime="B1950",
        frame="fk4",
    )

    sc1 = sc0.copy()
    sc1[1] = sc2[0]
    assert np.allclose(sc1.ra.to_value(u.deg), [1, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [3, 30])
    assert np.allclose(sc1.radial_velocity.to_value(u.km / u.s), [1, 10])
    assert sc1.obstime == Time("B1950")
    assert sc1.frame.name == "fk4"

    sc1 = sc0.copy()
    sc1[:] = sc2[0]
    assert np.allclose(sc1.ra.to_value(u.deg), [10, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [30, 30])
    assert np.allclose(sc1.radial_velocity.to_value(u.km / u.s), [10, 10])

    sc1 = sc0.copy()
    sc1[:] = sc2[:]
    assert np.allclose(sc1.ra.to_value(u.deg), [10, 20])
    assert np.allclose(sc1.dec.to_value(u.deg), [30, 40])
    assert np.allclose(sc1.radial_velocity.to_value(u.km / u.s), [10, 20])

    sc1 = sc0.copy()
    sc1[[1, 0]] = sc2[:]
    assert np.allclose(sc1.ra.to_value(u.deg), [20, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [40, 30])
    assert np.allclose(sc1.radial_velocity.to_value(u.km / u.s), [20, 10])


def test_setitem_exceptions():
    class SkyCoordSub(SkyCoord):
        pass

    obstime = "B1955"
    sc0 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg, frame="fk4")
    sc2 = SkyCoord([10, 20] * u.deg, [30, 40] * u.deg, frame="fk4", obstime=obstime)

    sc1 = SkyCoordSub(sc0)
    with pytest.raises(
        TypeError,
        match="an only set from object of same class: SkyCoordSub vs. SkyCoord",
    ):
        sc1[0] = sc2[0]

    sc1 = SkyCoord(sc0.ra, sc0.dec, frame="fk4", obstime="B2001")
    with pytest.raises(
        ValueError, match="can only set frame item from an equivalent frame"
    ):
        sc1.frame[0] = sc2.frame[0]

    sc1 = SkyCoord(sc0.ra[0], sc0.dec[0], frame="fk4", obstime=obstime)
    with pytest.raises(
        TypeError, match="scalar 'FK4' frame object does not support item assignment"
    ):
        sc1[0] = sc2[0]

    # Different differentials
    sc1 = SkyCoord(
        [1, 2] * u.deg,
        [3, 4] * u.deg,
        pm_ra_cosdec=[1, 2] * u.mas / u.yr,
        pm_dec=[3, 4] * u.mas / u.yr,
    )
    sc2 = SkyCoord(
        [10, 20] * u.deg, [30, 40] * u.deg, radial_velocity=[10, 20] * u.km / u.s
    )
    with pytest.raises(
        TypeError,
        match=(
            "can only set from object of same class: "
            "UnitSphericalCosLatDifferential vs. RadialDifferential"
        ),
    ):
        sc1[0] = sc2[0]


def test_insert():
    sc0 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
    sc1 = SkyCoord(5 * u.deg, 6 * u.deg)
    sc3 = SkyCoord([10, 20] * u.deg, [30, 40] * u.deg)
    sc4 = SkyCoord([[1, 2], [3, 4]] * u.deg, [[5, 6], [7, 8]] * u.deg)
    sc5 = SkyCoord([[10, 2], [30, 4]] * u.deg, [[50, 6], [70, 8]] * u.deg)

    # Insert a scalar
    sc = sc0.insert(1, sc1)
    assert skycoord_equal(sc, SkyCoord([1, 5, 2] * u.deg, [3, 6, 4] * u.deg))

    # Insert length=2 array at start of array
    sc = sc0.insert(0, sc3)
    assert skycoord_equal(sc, SkyCoord([10, 20, 1, 2] * u.deg, [30, 40, 3, 4] * u.deg))

    # Insert length=2 array at end of array
    sc = sc0.insert(2, sc3)
    assert skycoord_equal(sc, SkyCoord([1, 2, 10, 20] * u.deg, [3, 4, 30, 40] * u.deg))

    # Multidimensional
    sc = sc4.insert(1, sc5)
    assert skycoord_equal(
        sc,
        SkyCoord(
            [[1, 2], [10, 2], [30, 4], [3, 4]] * u.deg,
            [[5, 6], [50, 6], [70, 8], [7, 8]] * u.deg,
        ),
    )


def test_insert_exceptions():
    sc0 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
    sc1 = SkyCoord(5 * u.deg, 6 * u.deg)
    # sc3 = SkyCoord([10, 20]*u.deg, [30, 40]*u.deg)
    sc4 = SkyCoord([[1, 2], [3, 4]] * u.deg, [[5, 6], [7, 8]] * u.deg)

    with pytest.raises(TypeError, match="cannot insert into scalar"):
        sc1.insert(0, sc0)

    with pytest.raises(ValueError, match="axis must be 0"):
        sc0.insert(0, sc1, axis=1)

    with pytest.raises(TypeError, match="obj arg must be an integer"):
        sc0.insert(slice(None), sc0)

    with pytest.raises(
        IndexError, match="index -100 is out of bounds for axis 0 with size 2"
    ):
        sc0.insert(-100, sc0)

    # Bad shape
    with pytest.raises(
        ValueError,
        match=r"could not broadcast input array from shape \(2,2\) into shape \(2,?\)",
    ):
        sc0.insert(0, sc4)


def test_attr_conflicts():
    """
    Check conflicts resolution between coordinate attributes and init kwargs.
    """
    sc = SkyCoord(1, 2, frame="icrs", unit="deg", equinox="J1999", obstime="J2001")

    # OK if attrs both specified but with identical values
    SkyCoord(sc, equinox="J1999", obstime="J2001")

    # OK because sc.frame doesn't have obstime
    SkyCoord(sc.frame, equinox="J1999", obstime="J2100")

    # Not OK if attrs don't match
    with pytest.raises(ValueError) as err:
        SkyCoord(sc, equinox="J1999", obstime="J2002")
    assert "Coordinate attribute 'obstime'=" in str(err.value)

    # Same game but with fk4 which has equinox and obstime frame attrs
    sc = SkyCoord(1, 2, frame="fk4", unit="deg", equinox="J1999", obstime="J2001")

    # OK if attrs both specified but with identical values
    SkyCoord(sc, equinox="J1999", obstime="J2001")

    # Not OK if SkyCoord attrs don't match
    with pytest.raises(ValueError) as err:
        SkyCoord(sc, equinox="J1999", obstime="J2002")
    assert "Frame attribute 'obstime' has conflicting" in str(err.value)

    # Not OK because sc.frame has different attrs
    with pytest.raises(ValueError) as err:
        SkyCoord(sc.frame, equinox="J1999", obstime="J2002")
    assert "Frame attribute 'obstime' has conflicting" in str(err.value)


def test_frame_attr_getattr():
    """
    When accessing frame attributes like equinox, the value should come
    from self.frame when that object has the relevant attribute, otherwise
    from self.
    """
    sc = SkyCoord(1, 2, frame="icrs", unit="deg", equinox="J1999", obstime="J2001")
    assert sc.equinox == "J1999"  # Just the raw value (not validated)
    assert sc.obstime == "J2001"

    sc = SkyCoord(1, 2, frame="fk4", unit="deg", equinox="J1999", obstime="J2001")
    assert sc.equinox == Time("J1999")  # Coming from the self.frame object
    assert sc.obstime == Time("J2001")

    sc = SkyCoord(1, 2, frame="fk4", unit="deg", equinox="J1999")
    assert sc.equinox == Time("J1999")
    assert sc.obstime == Time("J1999")


def test_to_string():
    """
    Basic testing of converting SkyCoord to strings.  This just tests
    for a single input coordinate and and 1-element list.  It does not
    test the underlying `Angle.to_string` method itself.
    """
    coord = "1h2m3s 1d2m3s"
    for wrap in (lambda x: x, lambda x: [x]):
        sc = SkyCoord(wrap(coord))
        assert sc.to_string() == wrap("15.5125 1.03417")
        assert sc.to_string("dms") == wrap("15d30m45s 1d02m03s")
        assert sc.to_string("hmsdms") == wrap("01h02m03s +01d02m03s")
        with_kwargs = sc.to_string("hmsdms", precision=3, pad=True, alwayssign=True)
        assert with_kwargs == wrap("+01h02m03.000s +01d02m03.000s")


def test_repr():
    sc1 = SkyCoord(0 * u.deg, 1 * u.deg, frame="icrs")
    sc2 = SkyCoord(1 * u.deg, 1 * u.deg, frame="icrs", distance=1 * u.kpc)

    assert repr(sc1) == "<SkyCoord (ICRS): (ra, dec) in deg\n    (0., 1.)>"
    assert (
        repr(sc2)
        == "<SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, kpc)\n    (1., 1., 1.)>"
    )

    sc3 = SkyCoord(0.25 * u.deg, [1, 2.5] * u.deg, frame="icrs")
    assert repr(sc3).startswith("<SkyCoord (ICRS): (ra, dec) in deg\n")

    sc_default = SkyCoord(0 * u.deg, 1 * u.deg)
    assert repr(sc_default) == "<SkyCoord (ICRS): (ra, dec) in deg\n    (0., 1.)>"


def test_repr_altaz():
    sc2 = SkyCoord(1 * u.deg, 1 * u.deg, frame="icrs", distance=1 * u.kpc)

    if NUMPY_LT_2_0:
        expected_el_repr = "(-2309223., -3695529., -4641767.)"
    else:
        expected_el_repr = "(-2309223.0, -3695529.0, -4641767.0)"

    loc = EarthLocation(-2309223 * u.m, -3695529 * u.m, -4641767 * u.m)
    time = Time("2005-03-21 00:00:00")
    sc4 = sc2.transform_to(AltAz(location=loc, obstime=time))
    assert repr(sc4).startswith(
        "<SkyCoord (AltAz: obstime=2005-03-21 00:00:00.000, "
        f"location={expected_el_repr} m, pressure=0.0 hPa, "
        "temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron):"
        " (az, alt, distance) in (deg, deg, kpc)\n"
    )


def test_ops():
    """
    Tests miscellaneous operations like `len`
    """
    sc = SkyCoord(0 * u.deg, 1 * u.deg, frame="icrs")
    sc_arr = SkyCoord(0 * u.deg, [1, 2] * u.deg, frame="icrs")
    sc_empty = SkyCoord([] * u.deg, [] * u.deg, frame="icrs")

    assert sc.isscalar
    assert not sc_arr.isscalar
    assert not sc_empty.isscalar

    with pytest.raises(TypeError):
        len(sc)
    assert len(sc_arr) == 2
    assert len(sc_empty) == 0

    assert bool(sc)
    assert bool(sc_arr)
    assert not bool(sc_empty)

    assert sc_arr[0].isscalar
    assert len(sc_arr[:1]) == 1
    # A scalar shouldn't be indexable
    with pytest.raises(TypeError):
        sc[0:]
    # but it should be possible to just get an item
    sc_item = sc[()]
    assert sc_item.shape == ()
    # and to turn it into an array
    sc_1d = sc[np.newaxis]
    assert sc_1d.shape == (1,)

    with pytest.raises(TypeError):
        iter(sc)
    assert not np.iterable(sc)
    assert np.iterable(sc_arr)
    assert np.iterable(sc_empty)
    it = iter(sc_arr)
    assert next(it).dec == sc_arr[0].dec
    assert next(it).dec == sc_arr[1].dec
    with pytest.raises(StopIteration):
        next(it)


def test_none_transform():
    """
    Ensure that transforming from a SkyCoord with no frame provided works like
    ICRS
    """
    sc = SkyCoord(0 * u.deg, 1 * u.deg)
    sc_arr = SkyCoord(0 * u.deg, [1, 2] * u.deg)

    sc2 = sc.transform_to(ICRS)
    assert sc.ra == sc2.ra and sc.dec == sc2.dec

    sc5 = sc.transform_to("fk5")
    assert sc5.ra == sc2.transform_to("fk5").ra

    sc_arr2 = sc_arr.transform_to(ICRS)
    sc_arr5 = sc_arr.transform_to("fk5")
    npt.assert_array_equal(sc_arr5.ra, sc_arr2.transform_to("fk5").ra)


def test_position_angle_directly():
    """Regression check for #3800: position_angle should accept floats."""
    from astropy.coordinates import position_angle

    result = position_angle(10.0, 20.0, 10.0, 20.0)
    assert result.unit is u.radian
    assert result.value == 0.0


def test_directional_offset_by():
    # Round-trip tests: where is sc2 from sc1?
    # Use those offsets from sc1 and verify you get to sc2.
    npoints = 7  # How many points when doing vectors of SkyCoords
    for sc1 in [
        SkyCoord(0 * u.deg, -90 * u.deg),  # South pole
        SkyCoord(0 * u.deg, 90 * u.deg),  # North pole
        SkyCoord(1 * u.deg, 2 * u.deg),
        SkyCoord(
            np.linspace(0, 359, npoints),
            np.linspace(-90, 90, npoints),
            unit=u.deg,
            frame="fk4",
        ),
        SkyCoord(
            np.linspace(359, 0, npoints),
            np.linspace(-90, 90, npoints),
            unit=u.deg,
            frame="icrs",
        ),
        SkyCoord(
            np.linspace(-3, 3, npoints),
            np.linspace(-90, 90, npoints),
            unit=(u.rad, u.deg),
            frame="barycentricmeanecliptic",
        ),
    ]:
        for sc2 in [
            SkyCoord(5 * u.deg, 10 * u.deg),
            SkyCoord(
                np.linspace(0, 359, npoints),
                np.linspace(-90, 90, npoints),
                unit=u.deg,
                frame="galactic",
            ),
        ]:
            # Find the displacement from sc1 to sc2,
            posang = sc1.position_angle(sc2)
            sep = sc1.separation(sc2, origin_mismatch="ignore")

            # then do the offset from sc1 and verify that you are at sc2
            sc2a = sc1.directional_offset_by(position_angle=posang, separation=sep)
            assert (
                np.max(np.abs(sc2.separation(sc2a, origin_mismatch="ignore").arcsec))
                < 1e-3
            )

    # Specific test cases
    # Go over the North pole a little way, and
    # over the South pole a long way, to get to same spot
    sc1 = SkyCoord(0 * u.deg, 89 * u.deg)
    for posang, sep in [(0 * u.deg, 2 * u.deg), (180 * u.deg, 358 * u.deg)]:
        sc2 = sc1.directional_offset_by(posang, sep)
        assert allclose([sc2.ra.degree, sc2.dec.degree], [180, 89])
        # Go twice as far to ensure that dec is actually changing
        # and that >360deg is supported
        sc2 = sc1.directional_offset_by(posang, 2 * sep)
        assert allclose([sc2.ra.degree, sc2.dec.degree], [180, 87])

    # Verify that a separation of 180 deg in any direction gets to the antipode
    # and 360 deg returns to start
    sc1 = SkyCoord(10 * u.deg, 47 * u.deg)
    for posang in np.linspace(0, 377, npoints):
        sc2 = sc1.directional_offset_by(posang, 180 * u.deg)
        assert allclose([sc2.ra.degree, sc2.dec.degree], [190, -47])
        sc2 = sc1.directional_offset_by(posang, 360 * u.deg)
        assert allclose([sc2.ra.degree, sc2.dec.degree], [10, 47])

    # Verify that a 90 degree posang, which means East
    # corresponds to an increase in RA, by ~separation/cos(dec) and
    # a slight convergence to equator
    sc1 = SkyCoord(10 * u.deg, 60 * u.deg)
    sc2 = sc1.directional_offset_by(90 * u.deg, 1.0 * u.deg)
    assert 11.9 < sc2.ra.degree < 12.0
    assert 59.9 < sc2.dec.degree < 60.0


def test_table_to_coord():
    """
    Checks "end-to-end" use of `Table` with `SkyCoord` - the `Quantity`
    initializer is the intermediary that translate the table columns into
    something coordinates understands.

    (Regression test for #1762 )
    """
    from astropy.table import Column, Table

    t = Table()
    t.add_column(Column(data=[1, 2, 3], name="ra", unit=u.deg))
    t.add_column(Column(data=[4, 5, 6], name="dec", unit=u.deg))

    c = SkyCoord(t["ra"], t["dec"])

    assert allclose(c.ra.to(u.deg), [1, 2, 3] * u.deg)
    assert allclose(c.dec.to(u.deg), [4, 5, 6] * u.deg)


def assert_quantities_allclose(coord, q1s, attrs):
    """
    Compare two tuples of quantities.  This assumes that the values in q1 are of
    order(1) and uses atol=1e-13, rtol=0.  It also asserts that the units of the
    two quantities are the *same*, in order to check that the representation
    output has the expected units.
    """
    q2s = [getattr(coord, attr) for attr in attrs]
    assert len(q1s) == len(q2s)
    for q1, q2 in zip(q1s, q2s):
        assert q1.shape == q2.shape
        assert allclose(q1, q2, rtol=0, atol=1e-13 * q1.unit)


# Sets of inputs corresponding to Galactic frame
base_unit_attr_sets = [
    ("spherical", u.karcsec, u.karcsec, u.kpc, Latitude, "l", "b", "distance"),
    ("unitspherical", u.karcsec, u.karcsec, None, Latitude, "l", "b", None),
    ("physicsspherical", u.karcsec, u.karcsec, u.kpc, Angle, "phi", "theta", "r"),
    ("cartesian", u.km, u.km, u.km, u.Quantity, "u", "v", "w"),
    ("cylindrical", u.km, u.karcsec, u.km, Angle, "rho", "phi", "z"),
]

units_attr_sets = []
for base_unit_attr_set in base_unit_attr_sets:
    repr_name = base_unit_attr_set[0]
    for representation in (repr_name, REPRESENTATION_CLASSES[repr_name]):
        for c1, c2, c3 in ((1, 2, 3), ([1], [2], [3])):
            for arrayify in True, False:
                if arrayify:
                    c1 = np.array(c1)
                    c2 = np.array(c2)
                    c3 = np.array(c3)
                units_attr_sets.append(
                    base_unit_attr_set + (representation, c1, c2, c3)
                )
units_attr_args = (
    "repr_name",
    "unit1",
    "unit2",
    "unit3",
    "cls2",
    "attr1",
    "attr2",
    "attr3",
    "representation",
    "c1",
    "c2",
    "c3",
)


@pytest.mark.parametrize(
    units_attr_args, [x for x in units_attr_sets if x[0] != "unitspherical"]
)
def test_skycoord_three_components(
    repr_name,
    unit1,
    unit2,
    unit3,
    cls2,
    attr1,
    attr2,
    attr3,
    representation,
    c1,
    c2,
    c3,
):
    """
    Tests positional inputs using components (COMP1, COMP2, COMP3)
    and various representations.  Use weird units and Galactic frame.
    """
    sc = SkyCoord(
        c1,
        c2,
        c3,
        unit=(unit1, unit2, unit3),
        representation_type=representation,
        frame=Galactic,
    )
    assert_quantities_allclose(
        sc, (c1 * unit1, c2 * unit2, c3 * unit3), (attr1, attr2, attr3)
    )

    sc = SkyCoord(
        1000 * c1 * u.Unit(unit1 / 1000),
        cls2(c2, unit=unit2),
        1000 * c3 * u.Unit(unit3 / 1000),
        frame=Galactic,
        unit=(unit1, unit2, unit3),
        representation_type=representation,
    )
    assert_quantities_allclose(
        sc, (c1 * unit1, c2 * unit2, c3 * unit3), (attr1, attr2, attr3)
    )

    kwargs = {attr3: c3}
    sc = SkyCoord(
        c1,
        c2,
        unit=(unit1, unit2, unit3),
        frame=Galactic,
        representation_type=representation,
        **kwargs,
    )
    assert_quantities_allclose(
        sc, (c1 * unit1, c2 * unit2, c3 * unit3), (attr1, attr2, attr3)
    )

    kwargs = {attr1: c1, attr2: c2, attr3: c3}
    sc = SkyCoord(
        frame=Galactic,
        unit=(unit1, unit2, unit3),
        representation_type=representation,
        **kwargs,
    )
    assert_quantities_allclose(
        sc, (c1 * unit1, c2 * unit2, c3 * unit3), (attr1, attr2, attr3)
    )


@pytest.mark.parametrize(
    units_attr_args,
    [x for x in units_attr_sets if x[0] in ("spherical", "unitspherical")],
)
def test_skycoord_spherical_two_components(
    repr_name,
    unit1,
    unit2,
    unit3,
    cls2,
    attr1,
    attr2,
    attr3,
    representation,
    c1,
    c2,
    c3,
):
    """
    Tests positional inputs using components (COMP1, COMP2) for spherical
    representations.  Use weird units and Galactic frame.
    """
    sc = SkyCoord(
        c1, c2, unit=(unit1, unit2), frame=Galactic, representation_type=representation
    )
    assert_quantities_allclose(sc, (c1 * unit1, c2 * unit2), (attr1, attr2))

    sc = SkyCoord(
        1000 * c1 * u.Unit(unit1 / 1000),
        cls2(c2, unit=unit2),
        frame=Galactic,
        unit=(unit1, unit2, unit3),
        representation_type=representation,
    )
    assert_quantities_allclose(sc, (c1 * unit1, c2 * unit2), (attr1, attr2))

    kwargs = {attr1: c1, attr2: c2}
    sc = SkyCoord(
        frame=Galactic,
        unit=(unit1, unit2),
        representation_type=representation,
        **kwargs,
    )
    assert_quantities_allclose(sc, (c1 * unit1, c2 * unit2), (attr1, attr2))


@pytest.mark.parametrize(
    units_attr_args, [x for x in units_attr_sets if x[0] != "unitspherical"]
)
def test_galactic_three_components(
    repr_name,
    unit1,
    unit2,
    unit3,
    cls2,
    attr1,
    attr2,
    attr3,
    representation,
    c1,
    c2,
    c3,
):
    """
    Tests positional inputs using components (COMP1, COMP2, COMP3)
    and various representations.  Use weird units and Galactic frame.
    """
    sc = Galactic(
        1000 * c1 * u.Unit(unit1 / 1000),
        cls2(c2, unit=unit2),
        1000 * c3 * u.Unit(unit3 / 1000),
        representation_type=representation,
    )
    assert_quantities_allclose(
        sc, (c1 * unit1, c2 * unit2, c3 * unit3), (attr1, attr2, attr3)
    )

    kwargs = {attr3: c3 * unit3}
    sc = Galactic(c1 * unit1, c2 * unit2, representation_type=representation, **kwargs)
    assert_quantities_allclose(
        sc, (c1 * unit1, c2 * unit2, c3 * unit3), (attr1, attr2, attr3)
    )

    kwargs = {attr1: c1 * unit1, attr2: c2 * unit2, attr3: c3 * unit3}
    sc = Galactic(representation_type=representation, **kwargs)
    assert_quantities_allclose(
        sc, (c1 * unit1, c2 * unit2, c3 * unit3), (attr1, attr2, attr3)
    )


@pytest.mark.parametrize(
    units_attr_args,
    [x for x in units_attr_sets if x[0] in ("spherical", "unitspherical")],
)
def test_galactic_spherical_two_components(
    repr_name,
    unit1,
    unit2,
    unit3,
    cls2,
    attr1,
    attr2,
    attr3,
    representation,
    c1,
    c2,
    c3,
):
    """
    Tests positional inputs using components (COMP1, COMP2) for spherical
    representations.  Use weird units and Galactic frame.
    """

    sc = Galactic(
        1000 * c1 * u.Unit(unit1 / 1000),
        cls2(c2, unit=unit2),
        representation_type=representation,
    )
    assert_quantities_allclose(sc, (c1 * unit1, c2 * unit2), (attr1, attr2))

    sc = Galactic(c1 * unit1, c2 * unit2, representation_type=representation)
    assert_quantities_allclose(sc, (c1 * unit1, c2 * unit2), (attr1, attr2))

    kwargs = {attr1: c1 * unit1, attr2: c2 * unit2}
    sc = Galactic(representation_type=representation, **kwargs)
    assert_quantities_allclose(sc, (c1 * unit1, c2 * unit2), (attr1, attr2))


@pytest.mark.parametrize(
    ("repr_name", "unit1", "unit2", "unit3", "cls2", "attr1", "attr2", "attr3"),
    [x for x in base_unit_attr_sets if x[0] != "unitspherical"],
)
def test_skycoord_coordinate_input(
    repr_name, unit1, unit2, unit3, cls2, attr1, attr2, attr3
):
    c1, c2, c3 = 1, 2, 3
    sc = SkyCoord(
        [(c1, c2, c3)],
        unit=(unit1, unit2, unit3),
        representation_type=repr_name,
        frame="galactic",
    )
    assert_quantities_allclose(
        sc, ([c1] * unit1, [c2] * unit2, [c3] * unit3), (attr1, attr2, attr3)
    )

    c1, c2, c3 = 1 * unit1, 2 * unit2, 3 * unit3
    sc = SkyCoord([(c1, c2, c3)], representation_type=repr_name, frame="galactic")
    assert_quantities_allclose(
        sc, ([1] * unit1, [2] * unit2, [3] * unit3), (attr1, attr2, attr3)
    )


def test_skycoord_string_coordinate_input():
    sc = SkyCoord("01 02 03 +02 03 04", unit="deg", representation_type="unitspherical")
    assert_quantities_allclose(
        sc,
        (Angle("01:02:03", unit="deg"), Angle("02:03:04", unit="deg")),
        ("ra", "dec"),
    )
    sc = SkyCoord(
        ["01 02 03 +02 03 04"], unit="deg", representation_type="unitspherical"
    )
    assert_quantities_allclose(
        sc,
        (Angle(["01:02:03"], unit="deg"), Angle(["02:03:04"], unit="deg")),
        ("ra", "dec"),
    )


def test_units():
    sc = SkyCoord(1, 2, 3, unit="m", representation_type="cartesian")  # All get meters
    assert sc.x.unit is u.m
    assert sc.y.unit is u.m
    assert sc.z.unit is u.m

    # All get u.m
    sc = SkyCoord(1, 2 * u.km, 3, unit="m", representation_type="cartesian")
    assert sc.x.unit is u.m
    assert sc.y.unit is u.m
    assert sc.z.unit is u.m

    sc = SkyCoord(1, 2, 3, unit=u.m, representation_type="cartesian")  # All get u.m
    assert sc.x.unit is u.m
    assert sc.y.unit is u.m
    assert sc.z.unit is u.m

    sc = SkyCoord(1, 2, 3, unit="m, km, pc", representation_type="cartesian")
    assert_quantities_allclose(sc, (1 * u.m, 2 * u.km, 3 * u.pc), ("x", "y", "z"))

    with pytest.raises(u.UnitsError) as err:
        SkyCoord(1, 2, 3, unit=(u.m, u.m), representation_type="cartesian")
    assert "should have matching physical types" in str(err.value)

    SkyCoord(1, 2, 3, unit=(u.m, u.km, u.pc), representation_type="cartesian")
    assert_quantities_allclose(sc, (1 * u.m, 2 * u.km, 3 * u.pc), ("x", "y", "z"))


@pytest.mark.xfail
def test_units_known_fail():
    # should fail but doesn't => corner case oddity
    with pytest.raises(u.UnitsError):
        SkyCoord(1, 2, 3, unit=u.deg, representation_type="spherical")


def test_nodata_failure():
    with pytest.raises(ValueError):
        SkyCoord()


@pytest.mark.parametrize(("mode", "origin"), [("wcs", 0), ("all", 0), ("all", 1)])
def test_wcs_methods(mode, origin):
    from astropy.utils.data import get_pkg_data_contents
    from astropy.wcs import WCS
    from astropy.wcs.utils import pixel_to_skycoord

    header = get_pkg_data_contents(
        "../../wcs/tests/data/maps/1904-66_TAN.hdr", encoding="binary"
    )
    wcs = WCS(header)

    ref = SkyCoord(0.1 * u.deg, -89.0 * u.deg, frame="icrs")

    xp, yp = ref.to_pixel(wcs, mode=mode, origin=origin)

    # WCS is in FK5 so we need to transform back to ICRS
    new = pixel_to_skycoord(xp, yp, wcs, mode=mode, origin=origin).transform_to("icrs")

    assert_allclose(new.ra.degree, ref.ra.degree)
    assert_allclose(new.dec.degree, ref.dec.degree)

    # also try to round-trip with `from_pixel`
    scnew = SkyCoord.from_pixel(xp, yp, wcs, mode=mode, origin=origin).transform_to(
        "icrs"
    )
    assert_allclose(scnew.ra.degree, ref.ra.degree)
    assert_allclose(scnew.dec.degree, ref.dec.degree)

    # Also make sure the right type comes out
    class SkyCoord2(SkyCoord):
        pass

    scnew2 = SkyCoord2.from_pixel(xp, yp, wcs, mode=mode, origin=origin)
    assert scnew.__class__ is SkyCoord
    assert scnew2.__class__ is SkyCoord2


def test_frame_attr_transform_inherit():
    """
    Test that frame attributes get inherited as expected during transform.
    Driven by #3106.
    """
    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK5)
    c2 = c.transform_to(FK4)
    assert c2.equinox.value == "B1950.000"
    assert c2.obstime.value == "B1950.000"

    c2 = c.transform_to(FK4(equinox="J1975", obstime="J1980"))
    assert c2.equinox.value == "J1975.000"
    assert c2.obstime.value == "J1980.000"

    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK4)
    c2 = c.transform_to(FK5)
    assert c2.equinox.value == "J2000.000"
    assert c2.obstime is None

    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK4, obstime="J1980")
    c2 = c.transform_to(FK5)
    assert c2.equinox.value == "J2000.000"
    assert c2.obstime.value == "J1980.000"

    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK4, equinox="J1975", obstime="J1980")
    c2 = c.transform_to(FK5)
    assert c2.equinox.value == "J1975.000"
    assert c2.obstime.value == "J1980.000"

    c2 = c.transform_to(FK5(equinox="J1990"))
    assert c2.equinox.value == "J1990.000"
    assert c2.obstime.value == "J1980.000"

    # The work-around for #5722
    c = SkyCoord(1 * u.deg, 2 * u.deg, frame="fk5")
    c1 = SkyCoord(1 * u.deg, 2 * u.deg, frame="fk5", equinox="B1950.000")
    c2 = c1.transform_to(c)
    assert not c2.is_equivalent_frame(c)  # counterintuitive, but documented
    assert c2.equinox.value == "B1950.000"
    c3 = c1.transform_to(c, merge_attributes=False)
    assert c3.equinox.value == "J2000.000"
    assert c3.is_equivalent_frame(c)


def test_deepcopy():
    c1 = SkyCoord(1 * u.deg, 2 * u.deg)
    c2 = copy.copy(c1)
    c3 = copy.deepcopy(c1)

    c4 = SkyCoord(
        [1, 2] * u.m,
        [2, 3] * u.m,
        [3, 4] * u.m,
        representation_type="cartesian",
        frame="fk5",
        obstime="J1999.9",
        equinox="J1988.8",
    )
    c5 = copy.deepcopy(c4)
    assert np.all(c5.x == c4.x)  # and y and z
    assert c5.frame.name == c4.frame.name
    assert c5.obstime == c4.obstime
    assert c5.equinox == c4.equinox
    assert c5.representation_type == c4.representation_type


def test_no_copy():
    c1 = SkyCoord(np.arange(10.0) * u.hourangle, np.arange(20.0, 30.0) * u.deg)
    c2 = SkyCoord(c1, copy=False)
    # Note: c1.ra and c2.ra will *not* share memory, as these are recalculated
    # to be in "preferred" units.  See discussion in #4883.
    assert np.may_share_memory(c1.data.lon, c2.data.lon)
    c3 = SkyCoord(c1, copy=True)
    assert not np.may_share_memory(c1.data.lon, c3.data.lon)


def test_immutable():
    c1 = SkyCoord(1 * u.deg, 2 * u.deg)
    with pytest.raises(AttributeError):
        c1.ra = 3.0

    c1.foo = 42
    assert c1.foo == 42


@pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
def test_search_around():
    """
    Test the search_around_* methods

    Here we don't actually test the values are right, just that the methods of
    SkyCoord work.  The accuracy tests are in ``test_matching.py``
    """
    from astropy.utils import NumpyRNGContext

    with NumpyRNGContext(987654321):
        sc1 = SkyCoord(
            np.random.rand(20) * 360.0 * u.degree,
            (np.random.rand(20) * 180.0 - 90.0) * u.degree,
        )
        sc2 = SkyCoord(
            np.random.rand(100) * 360.0 * u.degree,
            (np.random.rand(100) * 180.0 - 90.0) * u.degree,
        )

        sc1ds = SkyCoord(ra=sc1.ra, dec=sc1.dec, distance=np.random.rand(20) * u.kpc)
        sc2ds = SkyCoord(ra=sc2.ra, dec=sc2.dec, distance=np.random.rand(100) * u.kpc)

    idx1_sky, idx2_sky, d2d_sky, d3d_sky = sc1.search_around_sky(sc2, 10 * u.deg)
    idx1_3d, idx2_3d, d2d_3d, d3d_3d = sc1ds.search_around_3d(sc2ds, 250 * u.pc)


def test_init_with_frame_instance_keyword():
    # Frame instance
    c1 = SkyCoord(3 * u.deg, 4 * u.deg, frame=FK5(equinox="J2010"))
    assert c1.equinox == Time("J2010")

    # Frame instance with data (data gets ignored)
    c2 = SkyCoord(
        3 * u.deg, 4 * u.deg, frame=FK5(1.0 * u.deg, 2 * u.deg, equinox="J2010")
    )
    assert c2.equinox == Time("J2010")
    assert allclose(c2.ra.degree, 3)
    assert allclose(c2.dec.degree, 4)

    # SkyCoord instance
    c3 = SkyCoord(3 * u.deg, 4 * u.deg, frame=c1)
    assert c3.equinox == Time("J2010")

    # Check duplicate arguments
    with pytest.raises(ValueError) as err:
        c = SkyCoord(3 * u.deg, 4 * u.deg, frame=FK5(equinox="J2010"), equinox="J2001")
    assert "Cannot specify frame attribute 'equinox'" in str(err.value)


def test_guess_from_table():
    from astropy.table import Column, Table
    from astropy.utils import NumpyRNGContext

    tab = Table()
    with NumpyRNGContext(987654321):
        tab.add_column(Column(data=np.random.rand(10), unit="deg", name="RA[J2000]"))
        tab.add_column(Column(data=np.random.rand(10), unit="deg", name="DEC[J2000]"))

    sc = SkyCoord.guess_from_table(tab)
    npt.assert_array_equal(sc.ra.deg, tab["RA[J2000]"])
    npt.assert_array_equal(sc.dec.deg, tab["DEC[J2000]"])

    # try without units in the table
    tab["RA[J2000]"].unit = None
    tab["DEC[J2000]"].unit = None
    # should fail if not given explicitly
    with pytest.raises(u.UnitsError):
        sc2 = SkyCoord.guess_from_table(tab)

    # but should work if provided
    sc2 = SkyCoord.guess_from_table(tab, unit=u.deg)
    npt.assert_array_equal(sc2.ra.deg, tab["RA[J2000]"])
    npt.assert_array_equal(sc2.dec.deg, tab["DEC[J2000]"])

    # should fail if two options are available - ambiguity bad!
    tab.add_column(Column(data=np.random.rand(10), name="RA_J1900"))
    with pytest.raises(ValueError) as excinfo:
        SkyCoord.guess_from_table(tab, unit=u.deg)
    assert "J1900" in excinfo.value.args[0] and "J2000" in excinfo.value.args[0]

    tab.remove_column("RA_J1900")
    tab["RA[J2000]"].unit = u.deg
    tab["DEC[J2000]"].unit = u.deg

    # but should succeed if the ambiguity can be broken b/c one of the matches
    # is the name of a different component
    tab.add_column(Column(data=np.random.rand(10) * u.mas / u.yr, name="pm_ra_cosdec"))
    tab.add_column(Column(data=np.random.rand(10) * u.mas / u.yr, name="pm_dec"))
    sc3 = SkyCoord.guess_from_table(tab)
    assert u.allclose(sc3.ra, tab["RA[J2000]"])
    assert u.allclose(sc3.dec, tab["DEC[J2000]"])
    assert u.allclose(sc3.pm_ra_cosdec, tab["pm_ra_cosdec"])
    assert u.allclose(sc3.pm_dec, tab["pm_dec"])

    # should fail if stuff doesn't have proper units
    tab["RA[J2000]"].unit = None
    tab["DEC[J2000]"].unit = None
    with pytest.raises(u.UnitTypeError, match="no unit was given."):
        SkyCoord.guess_from_table(tab)

    tab.remove_column("pm_ra_cosdec")
    tab.remove_column("pm_dec")

    # should also fail if user specifies something already in the table, but
    # should succeed even if the user has to give one of the components
    with pytest.raises(ValueError):
        SkyCoord.guess_from_table(tab, ra=tab["RA[J2000]"], unit=u.deg)

    oldra = tab["RA[J2000]"]
    tab.remove_column("RA[J2000]")
    sc3 = SkyCoord.guess_from_table(tab, ra=oldra, unit=u.deg)
    npt.assert_array_equal(sc3.ra.deg, oldra)
    npt.assert_array_equal(sc3.dec.deg, tab["DEC[J2000]"])

    # check a few non-ICRS/spherical systems
    x, y, z = np.arange(3).reshape(3, 1) * u.pc
    l, b = np.arange(2).reshape(2, 1) * u.deg

    tabcart = Table([x, y, z], names=("x", "y", "z"))
    tabgal = Table([b, l], names=("b", "l"))

    sc_cart = SkyCoord.guess_from_table(tabcart, representation_type="cartesian")
    npt.assert_array_equal(sc_cart.x, x)
    npt.assert_array_equal(sc_cart.y, y)
    npt.assert_array_equal(sc_cart.z, z)

    sc_gal = SkyCoord.guess_from_table(tabgal, frame="galactic")
    npt.assert_array_equal(sc_gal.l, l)
    npt.assert_array_equal(sc_gal.b, b)

    # also try some column names that *end* with the attribute name
    tabgal["b"].name = "gal_b"
    tabgal["l"].name = "gal_l"
    SkyCoord.guess_from_table(tabgal, frame="galactic")

    tabgal["gal_b"].name = "blob"
    tabgal["gal_l"].name = "central"
    with pytest.raises(ValueError):
        SkyCoord.guess_from_table(tabgal, frame="galactic")


def test_skycoord_list_creation():
    """
    Test that SkyCoord can be created in a reasonable way with lists of SkyCoords
    (regression for #2702)
    """
    sc = SkyCoord(ra=[1, 2, 3] * u.deg, dec=[4, 5, 6] * u.deg)
    sc0 = sc[0]
    sc2 = sc[2]
    scnew = SkyCoord([sc0, sc2])
    assert np.all(scnew.ra == [1, 3] * u.deg)
    assert np.all(scnew.dec == [4, 6] * u.deg)

    # also check ranges
    sc01 = sc[:2]
    scnew2 = SkyCoord([sc01, sc2])
    assert np.all(scnew2.ra == sc.ra)
    assert np.all(scnew2.dec == sc.dec)

    # now try with a mix of skycoord, frame, and repr objects
    frobj = ICRS(2 * u.deg, 5 * u.deg)
    reprobj = UnitSphericalRepresentation(3 * u.deg, 6 * u.deg)
    scnew3 = SkyCoord([sc0, frobj, reprobj])
    assert np.all(scnew3.ra == sc.ra)
    assert np.all(scnew3.dec == sc.dec)

    # should *fail* if different frame attributes or types are passed in
    scfk5_j2000 = SkyCoord(1 * u.deg, 4 * u.deg, frame="fk5")
    with pytest.raises(ValueError):
        SkyCoord([sc0, scfk5_j2000])
    scfk5_j2010 = SkyCoord(1 * u.deg, 4 * u.deg, frame="fk5", equinox="J2010")
    with pytest.raises(ValueError):
        SkyCoord([scfk5_j2000, scfk5_j2010])

    # but they should inherit if they're all consistent
    scfk5_2_j2010 = SkyCoord(2 * u.deg, 5 * u.deg, frame="fk5", equinox="J2010")
    scfk5_3_j2010 = SkyCoord(3 * u.deg, 6 * u.deg, frame="fk5", equinox="J2010")

    scnew4 = SkyCoord([scfk5_j2010, scfk5_2_j2010, scfk5_3_j2010])
    assert np.all(scnew4.ra == sc.ra)
    assert np.all(scnew4.dec == sc.dec)
    assert scnew4.equinox == Time("J2010")


def test_nd_skycoord_to_string():
    c = SkyCoord(np.ones((2, 2)), 1, unit=("deg", "deg"))
    ts = c.to_string()
    assert np.all(ts.shape == c.shape)
    assert np.all(ts == "1 1")


def test_equiv_skycoord():
    sci1 = SkyCoord(1 * u.deg, 2 * u.deg, frame="icrs")
    sci2 = SkyCoord(1 * u.deg, 3 * u.deg, frame="icrs")
    assert sci1.is_equivalent_frame(sci1)
    assert sci1.is_equivalent_frame(sci2)

    assert sci1.is_equivalent_frame(ICRS())
    assert not sci1.is_equivalent_frame(FK5())
    with pytest.raises(TypeError):
        sci1.is_equivalent_frame(10)

    scf1 = SkyCoord(1 * u.deg, 2 * u.deg, frame="fk5")
    scf2 = SkyCoord(1 * u.deg, 2 * u.deg, frame="fk5", equinox="J2005")
    # obstime is *not* an FK5 attribute, but we still want scf1 and scf3 to come
    # to come out different because they're part of SkyCoord
    scf3 = SkyCoord(1 * u.deg, 2 * u.deg, frame="fk5", obstime="J2005")

    assert scf1.is_equivalent_frame(scf1)
    assert not scf1.is_equivalent_frame(sci1)
    assert scf1.is_equivalent_frame(FK5())

    assert not scf1.is_equivalent_frame(scf2)
    assert scf2.is_equivalent_frame(FK5(equinox="J2005"))
    assert not scf3.is_equivalent_frame(scf1)
    assert not scf3.is_equivalent_frame(FK5(equinox="J2005"))


def test_equiv_skycoord_with_extra_attrs():
    """Regression test for #10658."""
    # GCRS has a CartesianRepresentationAttribute called obsgeoloc
    gcrs = GCRS(
        1 * u.deg, 2 * u.deg, obsgeoloc=CartesianRepresentation([1, 2, 3], unit=u.m)
    )
    # Create a SkyCoord where obsgeoloc tags along as an extra attribute
    sc1 = SkyCoord(gcrs).transform_to(ICRS)
    # Now create a SkyCoord with an equivalent frame but without the extra attribute
    sc2 = SkyCoord(sc1.frame)
    # The SkyCoords are therefore not equivalent, but check both directions
    assert not sc1.is_equivalent_frame(sc2)
    # This way around raised a TypeError which is fixed by #10658
    assert not sc2.is_equivalent_frame(sc1)


def test_constellations():
    # the actual test for accuracy is in test_funcs - this is just meant to make
    # sure we get sensible answers
    sc = SkyCoord(135 * u.deg, 65 * u.deg)
    assert sc.get_constellation() == "Ursa Major"
    assert sc.get_constellation(short_name=True) == "UMa"

    scs = SkyCoord([135] * 2 * u.deg, [65] * 2 * u.deg)
    npt.assert_equal(scs.get_constellation(), ["Ursa Major"] * 2)
    npt.assert_equal(scs.get_constellation(short_name=True), ["UMa"] * 2)


@pytest.mark.remote_data
def test_constellations_with_nameresolve():
    assert SkyCoord.from_name("And I").get_constellation(short_name=True) == "And"

    # you'd think "And ..." should be in Andromeda.  But you'd be wrong.
    assert SkyCoord.from_name("And VI").get_constellation() == "Pegasus"

    # maybe it's because And VI isn't really a galaxy?
    assert SkyCoord.from_name("And XXII").get_constellation() == "Pisces"
    assert SkyCoord.from_name("And XXX").get_constellation() == "Cassiopeia"
    # ok maybe not

    # ok, but at least some of the others do make sense...
    assert (
        SkyCoord.from_name("Coma Cluster").get_constellation(short_name=True) == "Com"
    )
    assert SkyCoord.from_name("Orion Nebula").get_constellation() == "Orion"
    assert SkyCoord.from_name("Triangulum Galaxy").get_constellation() == "Triangulum"


def test_getitem_representation():
    """
    Make sure current representation survives __getitem__ even if different
    from data representation.
    """
    sc = SkyCoord([1, 1] * u.deg, [2, 2] * u.deg)
    sc.representation_type = "cartesian"
    assert sc[0].representation_type is CartesianRepresentation


def test_spherical_offsets_to_api():
    i00 = SkyCoord(0 * u.arcmin, 0 * u.arcmin, frame="icrs")

    fk5 = SkyCoord(0 * u.arcmin, 0 * u.arcmin, frame="fk5")
    with pytest.raises(ValueError):
        # different frames should fail
        i00.spherical_offsets_to(fk5)

    i1deg = ICRS(1 * u.deg, 1 * u.deg)
    dra, ddec = i00.spherical_offsets_to(i1deg)
    assert_allclose(dra, 1 * u.deg)
    assert_allclose(ddec, 1 * u.deg)

    # make sure an abbreviated array-based version of the above also works
    i00s = SkyCoord([0] * 4 * u.arcmin, [0] * 4 * u.arcmin, frame="icrs")
    i01s = SkyCoord([0] * 4 * u.arcmin, np.arange(4) * u.arcmin, frame="icrs")
    dra, ddec = i00s.spherical_offsets_to(i01s)
    assert_allclose(dra, 0 * u.arcmin)
    assert_allclose(ddec, np.arange(4) * u.arcmin)


@pytest.mark.parametrize("frame", ["icrs", "galactic"])
@pytest.mark.parametrize(
    "comparison_data",
    [
        (0 * u.arcmin, 1 * u.arcmin),
        (1 * u.arcmin, 0 * u.arcmin),
        (1 * u.arcmin, 1 * u.arcmin),
    ],
)
def test_spherical_offsets_roundtrip(frame, comparison_data):
    i00 = SkyCoord(0 * u.arcmin, 0 * u.arcmin, frame=frame)
    comparison = SkyCoord(*comparison_data, frame=frame)

    dlon, dlat = i00.spherical_offsets_to(comparison)
    assert_allclose(dlon, comparison.data.lon)
    assert_allclose(dlat, comparison.data.lat)

    i00_back = comparison.spherical_offsets_by(-dlon, -dlat)

    # This reaches machine precision when only one component is changed, but for
    # the third parametrized case (both lon and lat change), the transformation
    # will have finite accuracy:
    assert_allclose(i00_back.data.lon, i00.data.lon, atol=1e-10 * u.rad)
    assert_allclose(i00_back.data.lat, i00.data.lat, atol=1e-10 * u.rad)

    # Test roundtripping the other direction:
    init_c = SkyCoord(40.0 * u.deg, 40.0 * u.deg, frame=frame)
    new_c = init_c.spherical_offsets_by(3.534 * u.deg, 2.2134 * u.deg)
    dlon, dlat = new_c.spherical_offsets_to(init_c)
    back_c = new_c.spherical_offsets_by(dlon, dlat)
    assert init_c.separation(back_c) < 1e-10 * u.deg


def test_frame_attr_changes():
    """
    This tests the case where a frame is added with a new frame attribute after
    a SkyCoord has been created.  This is necessary because SkyCoords get the
    attributes set at creation time, but the set of attributes can change as
    frames are added or removed from the transform graph.  This makes sure that
    everything continues to work consistently.
    """
    sc_before = SkyCoord(1 * u.deg, 2 * u.deg, frame="icrs")

    assert "fakeattr" not in dir(sc_before)

    class FakeFrame(BaseCoordinateFrame):
        fakeattr = Attribute()

    # doesn't matter what this does as long as it just puts the frame in the
    # transform graph
    transset = (ICRS, FakeFrame, lambda c, f: c)
    frame_transform_graph.add_transform(*transset)
    try:
        assert "fakeattr" in dir(sc_before)
        assert sc_before.fakeattr is None

        sc_after1 = SkyCoord(1 * u.deg, 2 * u.deg, frame="icrs")
        assert "fakeattr" in dir(sc_after1)
        assert sc_after1.fakeattr is None

        sc_after2 = SkyCoord(1 * u.deg, 2 * u.deg, frame="icrs", fakeattr=1)
        assert sc_after2.fakeattr == 1
    finally:
        frame_transform_graph.remove_transform(*transset)

    assert "fakeattr" not in dir(sc_before)
    assert "fakeattr" not in dir(sc_after1)
    assert "fakeattr" not in dir(sc_after2)


def test_cache_clear_sc():
    from astropy.coordinates import SkyCoord

    i = SkyCoord(1 * u.deg, 2 * u.deg)

    # Add an in frame units version of the rep to the cache.
    repr(i)

    assert len(i.cache["representation"]) == 2

    i.cache.clear()

    assert len(i.cache["representation"]) == 0


def test_set_attribute_exceptions():
    """Ensure no attribute for any frame can be set directly.

    Though it is fine if the current frame does not have it."""
    sc = SkyCoord(1.0 * u.deg, 2.0 * u.deg, frame="fk5")
    assert hasattr(sc.frame, "equinox")
    with pytest.raises(AttributeError):
        sc.equinox = "B1950"

    assert sc.relative_humidity is None
    sc.relative_humidity = 0.5
    assert sc.relative_humidity == 0.5
    assert not hasattr(sc.frame, "relative_humidity")


def test_extra_attributes():
    """Ensure any extra attributes are dealt with correctly.

    Regression test against #5743.
    """
    obstime_string = ["2017-01-01T00:00", "2017-01-01T00:10"]
    obstime = Time(obstime_string)
    sc = SkyCoord([5, 10], [20, 30], unit=u.deg, obstime=obstime_string)
    assert not hasattr(sc.frame, "obstime")
    assert type(sc.obstime) is Time
    assert sc.obstime.shape == (2,)
    assert np.all(sc.obstime == obstime)
    # ensure equivalency still works for more than one obstime.
    assert sc.is_equivalent_frame(sc)
    sc_1 = sc[1]
    assert sc_1.obstime == obstime[1]
    # Transforming to FK4 should use sc.obstime.
    sc_fk4 = sc.transform_to("fk4")
    assert np.all(sc_fk4.frame.obstime == obstime)
    # And transforming back should not loose it.
    sc2 = sc_fk4.transform_to("icrs")
    assert not hasattr(sc2.frame, "obstime")
    assert np.all(sc2.obstime == obstime)
    # Ensure obstime get taken from the SkyCoord if passed in directly.
    # (regression test for #5749).
    sc3 = SkyCoord([0.0, 1.0], [2.0, 3.0], unit="deg", frame=sc)
    assert np.all(sc3.obstime == obstime)
    # Finally, check that we can delete such attributes.
    del sc3.obstime
    assert sc3.obstime is None


def test_apply_space_motion():
    # use this 12 year period because it's a multiple of 4 to avoid the quirks
    # of leap years while having 2 leap seconds in it
    t1 = Time("2000-01-01T00:00")
    t2 = Time("2012-01-01T00:00")

    # Check a very simple case first:
    frame = ICRS(
        ra=10.0 * u.deg,
        dec=0 * u.deg,
        distance=10.0 * u.pc,
        pm_ra_cosdec=0.1 * u.deg / u.yr,
        pm_dec=0 * u.mas / u.yr,
        radial_velocity=0 * u.km / u.s,
    )

    # Cases that should work (just testing input for now):
    c1 = SkyCoord(frame, obstime=t1, pressure=101 * u.kPa)
    with pytest.warns(ErfaWarning, match='ERFA function "pmsafe" yielded .*'):
        # warning raised due to high PM chosen above
        applied1 = c1.apply_space_motion(new_obstime=t2)
        applied2 = c1.apply_space_motion(dt=12 * u.year)

    assert isinstance(applied1.frame, c1.frame.__class__)
    assert isinstance(applied2.frame, c1.frame.__class__)
    assert_allclose(applied1.ra, applied2.ra)
    assert_allclose(applied1.pm_ra_cosdec, applied2.pm_ra_cosdec)
    assert_allclose(applied1.dec, applied2.dec)
    assert_allclose(applied1.distance, applied2.distance)

    # ensure any frame attributes that were there before get passed through
    assert applied1.pressure == c1.pressure

    # there were 2 leap seconds between 2000 and 2010, so the difference in
    # the two forms of time evolution should be ~2 sec
    adt = np.abs(applied2.obstime - applied1.obstime)
    assert 1.9 * u.second < adt.to(u.second) < 2.1 * u.second

    c2 = SkyCoord(frame)
    with pytest.warns(ErfaWarning, match='ERFA function "pmsafe" yielded .*'):
        # warning raised due to high PM chosen above
        applied3 = c2.apply_space_motion(dt=6 * u.year)
    assert isinstance(applied3.frame, c1.frame.__class__)
    assert applied3.obstime is None

    # this should *not* be .6 deg due to space-motion on a sphere, but it
    # should be fairly close
    assert 0.5 * u.deg < applied3.ra - c1.ra < 0.7 * u.deg

    # the two cases should only match somewhat due to it being space motion, but
    # they should be at least this close
    assert quantity_allclose(
        applied1.ra - c1.ra, (applied3.ra - c1.ra) * 2, atol=1e-3 * u.deg
    )
    # but *not* this close
    assert not quantity_allclose(
        applied1.ra - c1.ra, (applied3.ra - c1.ra) * 2, atol=1e-4 * u.deg
    )

    with pytest.raises(ValueError):
        c2.apply_space_motion(new_obstime=t2)


def test_custom_frame_skycoord():
    # also regression check for the case from #7069

    class BlahBleeBlopFrame(BaseCoordinateFrame):
        default_representation = SphericalRepresentation
        # without a differential, SkyCoord creation fails
        # default_differential = SphericalDifferential

        _frame_specific_representation_info = {
            "spherical": [
                RepresentationMapping("lon", "lon", "recommended"),
                RepresentationMapping("lat", "lat", "recommended"),
                RepresentationMapping("distance", "radius", "recommended"),
            ]
        }

    SkyCoord(lat=1 * u.deg, lon=2 * u.deg, frame=BlahBleeBlopFrame)


def test_user_friendly_pm_error():
    """
    This checks that a more user-friendly error message is raised for the user
    if they pass, e.g., pm_ra instead of pm_ra_cosdec
    """

    with pytest.raises(ValueError) as e:
        SkyCoord(
            ra=150 * u.deg,
            dec=-11 * u.deg,
            pm_ra=100 * u.mas / u.yr,
            pm_dec=10 * u.mas / u.yr,
        )
    assert "pm_ra_cosdec" in str(e.value)

    with pytest.raises(ValueError) as e:
        SkyCoord(
            l=150 * u.deg,
            b=-11 * u.deg,
            pm_l=100 * u.mas / u.yr,
            pm_b=10 * u.mas / u.yr,
            frame="galactic",
        )
    assert "pm_l_cosb" in str(e.value)

    # The special error should not turn on here:
    with pytest.raises(ValueError) as e:
        SkyCoord(
            x=1 * u.pc,
            y=2 * u.pc,
            z=3 * u.pc,
            pm_ra=100 * u.mas / u.yr,
            pm_dec=10 * u.mas / u.yr,
            representation_type="cartesian",
        )
    assert "pm_ra_cosdec" not in str(e.value)


def test_contained_by():
    """
    Test Skycoord.contained(wcs,image)
    """

    header = """
WCSAXES =                    2 / Number of coordinate axes
CRPIX1  =               1045.0 / Pixel coordinate of reference point
CRPIX2  =               1001.0 / Pixel coordinate of reference point
PC1_1   =    -0.00556448550786 / Coordinate transformation matrix element
PC1_2   =   -0.001042120133257 / Coordinate transformation matrix element
PC2_1   =    0.001181477028705 / Coordinate transformation matrix element
PC2_2   =   -0.005590809742987 / Coordinate transformation matrix element
CDELT1  =                  1.0 / [deg] Coordinate increment at reference point
CDELT2  =                  1.0 / [deg] Coordinate increment at reference point
CUNIT1  = 'deg'                / Units of coordinate increment and value
CUNIT2  = 'deg'                / Units of coordinate increment and value
CTYPE1  = 'RA---TAN'           / TAN (gnomonic) projection + SIP distortions
CTYPE2  = 'DEC--TAN'           / TAN (gnomonic) projection + SIP distortions
CRVAL1  =      250.34971683647 / [deg] Coordinate value at reference point
CRVAL2  =      2.2808772582495 / [deg] Coordinate value at reference point
LONPOLE =                180.0 / [deg] Native longitude of celestial pole
LATPOLE =      2.2808772582495 / [deg] Native latitude of celestial pole
RADESYS = 'ICRS'               / Equatorial coordinate system
MJD-OBS =      58612.339199259 / [d] MJD of observation matching DATE-OBS
DATE-OBS= '2019-05-09T08:08:26.816Z' / ISO-8601 observation date matching MJD-OB
NAXIS   =                    2 / NAXIS
NAXIS1  =                 2136 / length of first array dimension
NAXIS2  =                 2078 / length of second array dimension
    """

    test_wcs = WCS(fits.Header.fromstring(header.strip(), "\n"))
    assert SkyCoord(254, 2, unit="deg").contained_by(test_wcs)
    assert not SkyCoord(240, 2, unit="deg").contained_by(test_wcs)

    img = np.zeros((2136, 2078))
    assert SkyCoord(250, 2, unit="deg").contained_by(test_wcs, img)
    assert not SkyCoord(240, 2, unit="deg").contained_by(test_wcs, img)

    ra = np.array([254.2, 254.1])
    dec = np.array([2, 12.1])
    coords = SkyCoord(ra, dec, unit="deg")
    assert np.all(test_wcs.footprint_contains(coords) == np.array([True, False]))


def test_none_differential_type():
    """
    This is a regression test for #8021
    """
    from astropy.coordinates import BaseCoordinateFrame

    class MockHeliographicStonyhurst(BaseCoordinateFrame):
        default_representation = SphericalRepresentation

        frame_specific_representation_info = {
            SphericalRepresentation: [
                RepresentationMapping(
                    reprname="lon", framename="lon", defaultunit=u.deg
                ),
                RepresentationMapping(
                    reprname="lat", framename="lat", defaultunit=u.deg
                ),
                RepresentationMapping(
                    reprname="distance", framename="radius", defaultunit=None
                ),
            ]
        }

    fr = MockHeliographicStonyhurst(lon=1 * u.deg, lat=2 * u.deg, radius=10 * u.au)
    SkyCoord(0 * u.deg, fr.lat, fr.radius, frame=fr)  # this was the failure


def test_multiple_aliases():
    # Define a frame with multiple aliases
    class MultipleAliasesFrame(BaseCoordinateFrame):
        name = ["alias_1", "alias_2"]
        default_representation = SphericalRepresentation

    # Register a transform, which adds the aliases to the transform graph
    tfun = lambda c, f: f.__class__(lon=c.lon, lat=c.lat)
    ftrans = FunctionTransform(
        tfun,
        MultipleAliasesFrame,
        MultipleAliasesFrame,
        register_graph=frame_transform_graph,
    )

    coord = SkyCoord(lon=1 * u.deg, lat=2 * u.deg, frame=MultipleAliasesFrame)

    # Test attribute-style access returns self (not a copy)
    assert coord.alias_1 is coord
    assert coord.alias_2 is coord

    # Test for aliases in __dir__()
    assert "alias_1" in coord.__dir__()
    assert "alias_2" in coord.__dir__()

    # Test transform_to() calls
    assert isinstance(coord.transform_to("alias_1").frame, MultipleAliasesFrame)
    assert isinstance(coord.transform_to("alias_2").frame, MultipleAliasesFrame)

    ftrans.unregister(frame_transform_graph)


@pytest.mark.parametrize(
    "kwargs, error_message",
    [
        (
            {"ra": 1, "dec": 1, "distance": 1 * u.pc, "unit": "deg"},
            r"Unit 'deg' \(angle\) could not be applied to 'distance'. ",
        ),
        (
            {
                "rho": 1 * u.m,
                "phi": 1,
                "z": 1 * u.m,
                "unit": "deg",
                "representation_type": "cylindrical",
            },
            r"Unit 'deg' \(angle\) could not be applied to 'rho'. ",
        ),
    ],
)
def test_passing_inconsistent_coordinates_and_units_raises_helpful_error(
    kwargs, error_message
):
    # https://github.com/astropy/astropy/issues/10725
    with pytest.raises(ValueError, match=error_message):
        SkyCoord(**kwargs)


@pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy.")
def test_match_to_catalog_3d_and_sky():
    # Test for issue #5857. See PR #11449
    cfk5_default = SkyCoord(
        [1, 2, 3, 4] * u.degree,
        [0, 0, 0, 0] * u.degree,
        distance=[1, 1, 1.5, 1] * u.kpc,
        frame="fk5",
    )
    cfk5_J1950 = cfk5_default.transform_to(FK5(equinox="J1950"))

    idx, angle, quantity = cfk5_J1950.match_to_catalog_3d(cfk5_default)
    npt.assert_array_equal(idx, [0, 1, 2, 3])
    assert_allclose(angle, 0 * u.deg, atol=1e-14 * u.deg, rtol=0)
    assert_allclose(quantity, 0 * u.kpc, atol=1e-14 * u.kpc, rtol=0)

    idx, angle, distance = cfk5_J1950.match_to_catalog_sky(cfk5_default)
    npt.assert_array_equal(idx, [0, 1, 2, 3])
    assert_allclose(angle, 0 * u.deg, atol=1e-14 * u.deg, rtol=0)
    assert_allclose(distance, 0 * u.kpc, atol=1e-14 * u.kpc, rtol=0)


def test_subclass_property_exception_error():
    """Regression test for gh-8340.

    Non-existing attribute access inside a property should give attribute
    error for the attribute, not for the property.
    """

    class custom_coord(SkyCoord):
        @property
        def prop(self):
            return self.random_attr

    c = custom_coord("00h42m30s", "+41d12m00s", frame="icrs")
    with pytest.raises(AttributeError, match="random_attr"):
        # Before this matched "prop" rather than "random_attr"
        c.prop
