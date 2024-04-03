# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
from copy import deepcopy

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import (
    EarthLocation,
    SkyCoord,
    galactocentric_frame_defaults,
)
from astropy.coordinates import representation as r
from astropy.coordinates.attributes import (
    Attribute,
    CoordinateAttribute,
    DifferentialAttribute,
    EarthLocationAttribute,
    QuantityAttribute,
    TimeAttribute,
)
from astropy.coordinates.baseframe import BaseCoordinateFrame, RepresentationMapping
from astropy.coordinates.builtin_frames import (
    FK4,
    FK5,
    GCRS,
    HCRS,
    ICRS,
    ITRS,
    AltAz,
    Galactic,
    Galactocentric,
    HADec,
)
from astropy.coordinates.representation import (
    REPRESENTATION_CLASSES,
    CartesianDifferential,
)
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time
from astropy.units import allclose
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning

from .test_representation import unitphysics  # this fixture is used below  # noqa: F401


def setup_function(func):
    """Copy original 'REPRESENTATIONCLASSES' as attribute in function."""
    func.REPRESENTATION_CLASSES_ORIG = deepcopy(REPRESENTATION_CLASSES)


def teardown_function(func):
    """Reset REPRESENTATION_CLASSES to original value."""
    REPRESENTATION_CLASSES.clear()
    REPRESENTATION_CLASSES.update(func.REPRESENTATION_CLASSES_ORIG)


def test_frame_attribute_descriptor():
    """Unit tests of the Attribute descriptor."""

    class TestAttributes:
        attr_none = Attribute()
        attr_2 = Attribute(default=2)
        attr_3_attr2 = Attribute(default=3, secondary_attribute="attr_2")
        attr_none_attr2 = Attribute(default=None, secondary_attribute="attr_2")
        attr_none_nonexist = Attribute(default=None, secondary_attribute="nonexist")

    t = TestAttributes()

    # Defaults
    assert t.attr_none is None
    assert t.attr_2 == 2
    assert t.attr_3_attr2 == 3
    assert t.attr_none_attr2 == t.attr_2
    assert t.attr_none_nonexist is None  # No default and non-existent secondary attr

    # Setting values via '_'-prefixed internal vars
    # (as would normally done in __init__)
    t._attr_none = 10
    assert t.attr_none == 10

    t._attr_2 = 20
    assert t.attr_2 == 20
    assert t.attr_3_attr2 == 3
    assert t.attr_none_attr2 == t.attr_2

    t._attr_none_attr2 = 40
    assert t.attr_none_attr2 == 40

    # Make sure setting values via public attribute fails
    with pytest.raises(AttributeError) as err:
        t.attr_none = 5
    assert "Cannot set frame attribute" in str(err.value)


def test_frame_subclass_attribute_descriptor():
    """Unit test of the attribute descriptors in subclasses."""
    _EQUINOX_B1980 = Time("B1980", scale="tai")

    class MyFK4(FK4):
        # equinox inherited from FK4, obstime overridden, and newattr is new
        obstime = TimeAttribute(default=_EQUINOX_B1980)
        newattr = Attribute(default="newattr")

    mfk4 = MyFK4()
    assert mfk4.equinox.value == "B1950.000"
    assert mfk4.obstime.value == "B1980.000"
    assert mfk4.newattr == "newattr"

    with pytest.warns(AstropyDeprecationWarning):
        assert set(mfk4.get_frame_attr_names()) == {"equinox", "obstime", "newattr"}

    mfk4 = MyFK4(equinox="J1980.0", obstime="J1990.0", newattr="world")
    assert mfk4.equinox.value == "J1980.000"
    assert mfk4.obstime.value == "J1990.000"
    assert mfk4.newattr == "world"


def test_frame_multiple_inheritance_attribute_descriptor():
    """
    Ensure that all attributes are accumulated in case of inheritance from
    multiple BaseCoordinateFrames.  See
    https://github.com/astropy/astropy/pull/11099#issuecomment-735829157
    """

    class Frame1(BaseCoordinateFrame):
        attr1 = Attribute()

    class Frame2(BaseCoordinateFrame):
        attr2 = Attribute()

    class Frame3(Frame1, Frame2):
        pass

    assert len(Frame3.frame_attributes) == 2
    assert "attr1" in Frame3.frame_attributes
    assert "attr2" in Frame3.frame_attributes

    # In case the same attribute exists in both frames, the one from the
    # left-most class in the MRO should take precedence
    class Frame4(BaseCoordinateFrame):
        attr1 = Attribute()
        attr2 = Attribute()

    class Frame5(Frame1, Frame4):
        pass

    assert Frame5.frame_attributes["attr1"] is Frame1.frame_attributes["attr1"]
    assert Frame5.frame_attributes["attr2"] is Frame4.frame_attributes["attr2"]


def test_differentialattribute():
    # Test logic of passing input through to allowed class
    vel = [1, 2, 3] * u.km / u.s
    dif = r.CartesianDifferential(vel)

    class TestFrame(BaseCoordinateFrame):
        attrtest = DifferentialAttribute(
            default=dif, allowed_classes=[r.CartesianDifferential]
        )

    frame1 = TestFrame()
    frame2 = TestFrame(attrtest=dif)
    frame3 = TestFrame(attrtest=vel)

    assert np.all(frame1.attrtest.d_xyz == frame2.attrtest.d_xyz)
    assert np.all(frame1.attrtest.d_xyz == frame3.attrtest.d_xyz)

    # This shouldn't work if there is more than one allowed class:
    class TestFrame2(BaseCoordinateFrame):
        attrtest = DifferentialAttribute(
            default=dif,
            allowed_classes=[r.CartesianDifferential, r.CylindricalDifferential],
        )

    frame1 = TestFrame2()
    frame2 = TestFrame2(attrtest=dif)
    with pytest.raises(TypeError):
        TestFrame2(attrtest=vel)


def test_create_data_frames():
    # from repr
    i1 = ICRS(r.SphericalRepresentation(1 * u.deg, 2 * u.deg, 3 * u.kpc))
    i2 = ICRS(r.UnitSphericalRepresentation(lon=1 * u.deg, lat=2 * u.deg))

    # from preferred name
    i3 = ICRS(ra=1 * u.deg, dec=2 * u.deg, distance=3 * u.kpc)
    i4 = ICRS(ra=1 * u.deg, dec=2 * u.deg)

    assert i1.data.lat == i3.data.lat
    assert i1.data.lon == i3.data.lon
    assert i1.data.distance == i3.data.distance

    assert i2.data.lat == i4.data.lat
    assert i2.data.lon == i4.data.lon

    # now make sure the preferred names work as properties
    assert_allclose(i1.ra, i3.ra)
    assert_allclose(i2.ra, i4.ra)
    assert_allclose(i1.distance, i3.distance)

    with pytest.raises(AttributeError):
        i1.ra = [11.0] * u.deg


def test_create_orderered_data():
    TOL = 1e-10 * u.deg

    i = ICRS(1 * u.deg, 2 * u.deg)
    assert (i.ra - 1 * u.deg) < TOL
    assert (i.dec - 2 * u.deg) < TOL

    g = Galactic(1 * u.deg, 2 * u.deg)
    assert (g.l - 1 * u.deg) < TOL
    assert (g.b - 2 * u.deg) < TOL

    a = AltAz(1 * u.deg, 2 * u.deg)
    assert (a.az - 1 * u.deg) < TOL
    assert (a.alt - 2 * u.deg) < TOL

    with pytest.raises(TypeError):
        ICRS(1 * u.deg, 2 * u.deg, 1 * u.deg, 2 * u.deg)

    with pytest.raises(TypeError):
        sph = r.SphericalRepresentation(1 * u.deg, 2 * u.deg, 3 * u.kpc)
        ICRS(sph, 1 * u.deg, 2 * u.deg)


def test_create_nodata_frames():
    i = ICRS()
    assert len(i.frame_attributes) == 0

    f5 = FK5()
    assert f5.equinox == FK5.get_frame_attr_defaults()["equinox"]

    f4 = FK4()
    assert f4.equinox == FK4.get_frame_attr_defaults()["equinox"]

    # obstime is special because it's a property that uses equinox if obstime is not set
    assert f4.obstime in (
        FK4.get_frame_attr_defaults()["obstime"],
        FK4.get_frame_attr_defaults()["equinox"],
    )


def test_no_data_nonscalar_frames():
    a1 = AltAz(
        obstime=Time("2012-01-01") + np.arange(10.0) * u.day,
        temperature=np.ones((3, 1)) * u.deg_C,
    )
    assert a1.obstime.shape == (3, 10)
    assert a1.temperature.shape == (3, 10)
    assert a1.shape == (3, 10)
    with pytest.raises(ValueError) as exc:
        AltAz(
            obstime=Time("2012-01-01") + np.arange(10.0) * u.day,
            temperature=np.ones((3,)) * u.deg_C,
        )
    assert "inconsistent shapes" in str(exc.value)


def test_frame_repr():
    i = ICRS()
    assert repr(i) == "<ICRS Frame>"

    f5 = FK5()
    assert repr(f5).startswith("<FK5 Frame (equinox=")

    i2 = ICRS(ra=1 * u.deg, dec=2 * u.deg)
    i3 = ICRS(ra=1 * u.deg, dec=2 * u.deg, distance=3 * u.kpc)

    assert repr(i2) == "<ICRS Coordinate: (ra, dec) in deg\n    (1., 2.)>"
    assert (
        repr(i3)
        == "<ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)\n    (1., 2., 3.)>"
    )

    # try with arrays
    i2 = ICRS(ra=[1.1, 2.1] * u.deg, dec=[2.1, 3.1] * u.deg)
    i3 = ICRS(
        ra=[1.1, 2.1] * u.deg, dec=[-15.6, 17.1] * u.deg, distance=[11.0, 21.0] * u.kpc
    )

    assert (
        repr(i2) == "<ICRS Coordinate: (ra, dec) in deg\n    [(1.1, 2.1), (2.1, 3.1)]>"
    )

    assert (
        repr(i3) == "<ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)\n"
        "    [(1.1, -15.6, 11.), (2.1,  17.1, 21.)]>"
    )


def test_frame_repr_vels():
    i = ICRS(
        ra=1 * u.deg,
        dec=2 * u.deg,
        pm_ra_cosdec=1 * u.marcsec / u.yr,
        pm_dec=2 * u.marcsec / u.yr,
    )

    # unit comes out as mas/yr because of the preferred units defined in the
    # frame RepresentationMapping
    assert (
        repr(i) == "<ICRS Coordinate: (ra, dec) in deg\n"
        "    (1., 2.)\n"
        " (pm_ra_cosdec, pm_dec) in mas / yr\n"
        "    (1., 2.)>"
    )


def test_converting_units():
    # this is a regular expression that with split (see below) removes what's
    # the decimal point  to fix rounding problems
    rexrepr = re.compile(r"(.*?=\d\.).*?( .*?=\d\.).*?( .*)")

    # Use values that aren't subject to rounding down to X.9999...
    i2 = ICRS(ra=2.0 * u.deg, dec=2.0 * u.deg)
    i2_many = ICRS(ra=[2.0, 4.0] * u.deg, dec=[2.0, -8.1] * u.deg)

    # converting from FK5 to ICRS and back changes the *internal* representation,
    # but it should still come out in the preferred form

    i4 = i2.transform_to(FK5()).transform_to(ICRS())
    i4_many = i2_many.transform_to(FK5()).transform_to(ICRS())

    ri2 = "".join(rexrepr.split(repr(i2)))
    ri4 = "".join(rexrepr.split(repr(i4)))
    assert ri2 == ri4
    assert i2.data.lon.unit != i4.data.lon.unit  # Internal repr changed

    ri2_many = "".join(rexrepr.split(repr(i2_many)))
    ri4_many = "".join(rexrepr.split(repr(i4_many)))

    assert ri2_many == ri4_many
    assert i2_many.data.lon.unit != i4_many.data.lon.unit  # Internal repr changed

    # but that *shouldn't* hold if we turn off units for the representation
    class FakeICRS(ICRS):
        frame_specific_representation_info = {
            "spherical": [
                RepresentationMapping("lon", "ra", u.hourangle),
                RepresentationMapping("lat", "dec", None),
                RepresentationMapping("distance", "distance"),
            ]  # should fall back to default of None unit
        }

    fi = FakeICRS(i4.data)
    ri2 = "".join(rexrepr.split(repr(i2)))
    rfi = "".join(rexrepr.split(repr(fi)))
    rfi = re.sub("FakeICRS", "ICRS", rfi)  # Force frame name to match
    assert ri2 != rfi

    # the attributes should also get the right units
    assert i2.dec.unit == i4.dec.unit
    # unless no/explicitly given units
    assert i2.dec.unit != fi.dec.unit
    assert i2.ra.unit != fi.ra.unit
    assert fi.ra.unit == u.hourangle


def test_representation_info():
    class NewICRS1(ICRS):
        frame_specific_representation_info = {
            r.SphericalRepresentation: [
                RepresentationMapping("lon", "rara", u.hourangle),
                RepresentationMapping("lat", "decdec", u.degree),
                RepresentationMapping("distance", "distance", u.kpc),
            ]
        }

    i1 = NewICRS1(
        rara=10 * u.degree,
        decdec=-12 * u.deg,
        distance=1000 * u.pc,
        pm_rara_cosdecdec=100 * u.mas / u.yr,
        pm_decdec=17 * u.mas / u.yr,
        radial_velocity=10 * u.km / u.s,
    )
    assert allclose(i1.rara, 10 * u.deg)
    assert i1.rara.unit == u.hourangle
    assert allclose(i1.decdec, -12 * u.deg)
    assert allclose(i1.distance, 1000 * u.pc)
    assert i1.distance.unit == u.kpc
    assert allclose(i1.pm_rara_cosdecdec, 100 * u.mas / u.yr)
    assert allclose(i1.pm_decdec, 17 * u.mas / u.yr)

    # this should auto-set the names of UnitSpherical:
    i1.set_representation_cls(
        r.UnitSphericalRepresentation, s=r.UnitSphericalCosLatDifferential
    )
    assert allclose(i1.rara, 10 * u.deg)
    assert allclose(i1.decdec, -12 * u.deg)
    assert allclose(i1.pm_rara_cosdecdec, 100 * u.mas / u.yr)
    assert allclose(i1.pm_decdec, 17 * u.mas / u.yr)

    # For backwards compatibility, we also support the string name in the
    # representation info dictionary:
    class NewICRS2(ICRS):
        frame_specific_representation_info = {
            "spherical": [
                RepresentationMapping("lon", "ang1", u.hourangle),
                RepresentationMapping("lat", "ang2", u.degree),
                RepresentationMapping("distance", "howfar", u.kpc),
            ]
        }

    i2 = NewICRS2(ang1=10 * u.degree, ang2=-12 * u.deg, howfar=1000 * u.pc)
    assert allclose(i2.ang1, 10 * u.deg)
    assert i2.ang1.unit == u.hourangle
    assert allclose(i2.ang2, -12 * u.deg)
    assert allclose(i2.howfar, 1000 * u.pc)
    assert i2.howfar.unit == u.kpc

    # Test that the differential kwargs get overridden
    class NewICRS3(ICRS):
        frame_specific_representation_info = {
            r.SphericalCosLatDifferential: [
                RepresentationMapping("d_lon_coslat", "pm_ang1", u.hourangle / u.year),
                RepresentationMapping("d_lat", "pm_ang2"),
                RepresentationMapping("d_distance", "vlos", u.kpc / u.Myr),
            ]
        }

    i3 = NewICRS3(
        lon=10 * u.degree,
        lat=-12 * u.deg,
        distance=1000 * u.pc,
        pm_ang1=1 * u.mas / u.yr,
        pm_ang2=2 * u.mas / u.yr,
        vlos=100 * u.km / u.s,
    )
    assert allclose(i3.pm_ang1, 1 * u.mas / u.yr)
    assert i3.pm_ang1.unit == u.hourangle / u.year
    assert allclose(i3.pm_ang2, 2 * u.mas / u.yr)
    assert allclose(i3.vlos, 100 * u.km / u.s)
    assert i3.vlos.unit == u.kpc / u.Myr


def test_realizing():
    rep = r.SphericalRepresentation(1 * u.deg, 2 * u.deg, 3 * u.kpc)

    i = ICRS()
    i2 = i.realize_frame(rep)

    assert not i.has_data
    assert i2.has_data

    f = FK5(equinox=Time("J2001"))
    f2 = f.realize_frame(rep)

    assert not f.has_data
    assert f2.has_data

    assert f2.equinox == f.equinox
    assert f2.equinox != FK5.get_frame_attr_defaults()["equinox"]

    # Check that a nicer error message is returned:
    with pytest.raises(
        TypeError, match="Class passed as data instead of a representation"
    ):
        f.realize_frame(f.representation_type)


def test_replicating():
    i = ICRS(ra=[1] * u.deg, dec=[2] * u.deg)

    icopy = i.replicate(copy=True)
    irepl = i.replicate(copy=False)
    i.data._lat[:] = 0 * u.deg
    assert np.all(i.data.lat == irepl.data.lat)
    assert np.all(i.data.lat != icopy.data.lat)

    iclone = i.replicate_without_data()
    assert i.has_data
    assert not iclone.has_data

    aa = AltAz(alt=1 * u.deg, az=2 * u.deg, obstime=Time("J2000"))
    aaclone = aa.replicate_without_data(obstime=Time("J2001"))
    assert not aaclone.has_data
    assert aa.obstime != aaclone.obstime
    assert aa.pressure == aaclone.pressure
    assert aa.obswl == aaclone.obswl


def test_getitem():
    rep = r.SphericalRepresentation(
        [1, 2, 3] * u.deg, [4, 5, 6] * u.deg, [7, 8, 9] * u.kpc
    )

    i = ICRS(rep)
    assert len(i.ra) == 3

    iidx = i[1:]
    assert len(iidx.ra) == 2

    iidx2 = i[0]
    assert iidx2.ra.isscalar


def test_transform():
    """
    This test just makes sure the transform architecture works, but does *not*
    actually test all the builtin transforms themselves are accurate.

    """
    i = ICRS(ra=[1, 2] * u.deg, dec=[3, 4] * u.deg)
    f = i.transform_to(FK5())
    i2 = f.transform_to(ICRS())

    assert i2.data.__class__ == r.UnitSphericalRepresentation

    assert_allclose(i.ra, i2.ra)
    assert_allclose(i.dec, i2.dec)

    i = ICRS(ra=[1, 2] * u.deg, dec=[3, 4] * u.deg, distance=[5, 6] * u.kpc)
    f = i.transform_to(FK5())
    i2 = f.transform_to(ICRS())

    assert i2.data.__class__ != r.UnitSphericalRepresentation

    f = FK5(ra=1 * u.deg, dec=2 * u.deg, equinox=Time("J2001"))
    f4 = f.transform_to(FK4())
    f4_2 = f.transform_to(FK4(equinox=f.equinox))

    # make sure attributes are copied over correctly
    assert f4.equinox == FK4().equinox
    assert f4_2.equinox == f.equinox

    # make sure self-transforms also work
    i = ICRS(ra=[1, 2] * u.deg, dec=[3, 4] * u.deg)
    i2 = i.transform_to(ICRS())

    assert_allclose(i.ra, i2.ra)
    assert_allclose(i.dec, i2.dec)

    f = FK5(ra=1 * u.deg, dec=2 * u.deg, equinox=Time("J2001"))
    f2 = f.transform_to(FK5())  # default equinox, so should be *different*
    assert f2.equinox == FK5().equinox
    with pytest.raises(AssertionError):
        assert_allclose(f.ra, f2.ra)
    with pytest.raises(AssertionError):
        assert_allclose(f.dec, f2.dec)

    # finally, check Galactic round-tripping
    i1 = ICRS(ra=[1, 2] * u.deg, dec=[3, 4] * u.deg)
    i2 = i1.transform_to(Galactic()).transform_to(ICRS())

    assert_allclose(i1.ra, i2.ra)
    assert_allclose(i1.dec, i2.dec)


def test_transform_to_nonscalar_nodata_frame():
    # https://github.com/astropy/astropy/pull/5254#issuecomment-241592353
    times = Time("2016-08-23") + np.linspace(0, 10, 12) * u.day
    coo1 = ICRS(
        ra=[[0.0], [10.0], [20.0]] * u.deg, dec=[[-30.0], [30.0], [60.0]] * u.deg
    )
    coo2 = coo1.transform_to(FK5(equinox=times))
    assert coo2.shape == (3, 12)


def test_setitem_no_velocity():
    """Test different flavors of item setting for a Frame without a velocity."""
    obstime = "B1955"
    sc0 = FK4([1, 2] * u.deg, [3, 4] * u.deg, obstime=obstime)
    sc2 = FK4([10, 20] * u.deg, [30, 40] * u.deg, obstime=obstime)

    sc1 = sc0.copy()
    sc1_repr = repr(sc1)
    assert "representation" in sc1.cache
    sc1[1] = sc2[0]
    assert sc1.cache == {}
    assert repr(sc2) != sc1_repr

    assert np.allclose(sc1.ra.to_value(u.deg), [1, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [3, 30])
    assert sc1.obstime == sc2.obstime
    assert sc1.name == "fk4"

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

    # Works for array-valued obstime so long as they are considered equivalent
    sc1 = FK4(sc0.ra, sc0.dec, obstime=[obstime, obstime])
    sc1[0] = sc2[0]

    # Multidimensional coordinates
    sc1 = FK4([[1, 2], [3, 4]] * u.deg, [[5, 6], [7, 8]] * u.deg)
    sc2 = FK4([[10, 20], [30, 40]] * u.deg, [[50, 60], [70, 80]] * u.deg)
    sc1[0] = sc2[0]
    assert np.allclose(sc1.ra.to_value(u.deg), [[10, 20], [3, 4]])
    assert np.allclose(sc1.dec.to_value(u.deg), [[50, 60], [7, 8]])


def test_setitem_velocities():
    """Test different flavors of item setting for a Frame with a velocity."""
    sc0 = FK4(
        [1, 2] * u.deg,
        [3, 4] * u.deg,
        radial_velocity=[1, 2] * u.km / u.s,
        obstime="B1950",
    )
    sc2 = FK4(
        [10, 20] * u.deg,
        [30, 40] * u.deg,
        radial_velocity=[10, 20] * u.km / u.s,
        obstime="B1950",
    )

    sc1 = sc0.copy()
    sc1[1] = sc2[0]
    assert np.allclose(sc1.ra.to_value(u.deg), [1, 10])
    assert np.allclose(sc1.dec.to_value(u.deg), [3, 30])
    assert np.allclose(sc1.radial_velocity.to_value(u.km / u.s), [1, 10])
    assert sc1.obstime == sc2.obstime
    assert sc1.name == "fk4"

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
    obstime = "B1950"
    sc0 = FK4([1, 2] * u.deg, [3, 4] * u.deg)
    sc2 = FK4([10, 20] * u.deg, [30, 40] * u.deg, obstime=obstime)

    sc1 = Galactic(sc0.ra, sc0.dec)
    with pytest.raises(
        TypeError, match="can only set from object of same class: Galactic vs. FK4"
    ):
        sc1[0] = sc2[0]

    sc1 = FK4(sc0.ra, sc0.dec, obstime="B2001")
    with pytest.raises(
        ValueError, match="can only set frame item from an equivalent frame"
    ):
        sc1[0] = sc2[0]

    sc1 = FK4(sc0.ra[0], sc0.dec[0], obstime=obstime)
    with pytest.raises(
        TypeError, match="scalar 'FK4' frame object does not support item assignment"
    ):
        sc1[0] = sc2[0]

    sc1 = FK4(obstime=obstime)
    with pytest.raises(ValueError, match="cannot set frame which has no data"):
        sc1[0] = sc2[0]

    sc1 = FK4(sc0.ra, sc0.dec, obstime=[obstime, "B1980"])
    with pytest.raises(
        ValueError, match="can only set frame item from an equivalent frame"
    ):
        sc1[0] = sc2[0]

    # Wrong shape
    sc1 = FK4([sc0.ra], [sc0.dec], obstime=[obstime, "B1980"])
    with pytest.raises(
        ValueError, match="can only set frame item from an equivalent frame"
    ):
        sc1[0] = sc2[0]


def test_time_inputs():
    """
    Test validation and conversion of inputs for equinox and obstime attributes.

    """
    c = FK4(1 * u.deg, 2 * u.deg, equinox="J2001.5", obstime="2000-01-01 12:00:00")
    assert c.equinox == Time("J2001.5")
    assert c.obstime == Time("2000-01-01 12:00:00")

    with pytest.raises(ValueError) as err:
        c = FK4(1 * u.deg, 2 * u.deg, equinox=1.5)
    assert "Invalid time input" in str(err.value)

    with pytest.raises(ValueError) as err:
        c = FK4(1 * u.deg, 2 * u.deg, obstime="hello")
    assert "Invalid time input" in str(err.value)

    # A vector time should work if the shapes match, and we automatically
    # broadcast the basic data.
    c = FK4([1, 2] * u.deg, [2, 3] * u.deg, obstime=["J2000", "J2001"])
    assert c.shape == (2,)
    c = FK4(1 * u.deg, 2 * u.deg, obstime=["J2000", "J2001"])
    assert c.shape == (2,)

    # If the shapes are not broadcastable, then we should raise an exception.
    with pytest.raises(ValueError, match="inconsistent shapes"):
        FK4([1, 2, 3] * u.deg, [4, 5, 6] * u.deg, obstime=["J2000", "J2001"])


def test_is_frame_attr_default():
    """
    Check that the `is_frame_attr_default` machinery works as expected

    """
    c1 = FK5(ra=1 * u.deg, dec=1 * u.deg)
    c2 = FK5(
        ra=1 * u.deg, dec=1 * u.deg, equinox=FK5.get_frame_attr_defaults()["equinox"]
    )
    c3 = FK5(ra=1 * u.deg, dec=1 * u.deg, equinox=Time("J2001.5"))

    assert c1.equinox == c2.equinox
    assert c1.equinox != c3.equinox

    assert c1.is_frame_attr_default("equinox")
    assert not c2.is_frame_attr_default("equinox")
    assert not c3.is_frame_attr_default("equinox")

    c4 = c1.realize_frame(r.UnitSphericalRepresentation(3 * u.deg, 4 * u.deg))
    c5 = c2.realize_frame(r.UnitSphericalRepresentation(3 * u.deg, 4 * u.deg))

    assert c4.is_frame_attr_default("equinox")
    assert not c5.is_frame_attr_default("equinox")


def test_altaz_attributes():
    aa = AltAz(1 * u.deg, 2 * u.deg)
    assert aa.obstime is None
    assert aa.location is None

    aa2 = AltAz(1 * u.deg, 2 * u.deg, obstime="J2000")
    assert aa2.obstime == Time("J2000")

    aa3 = AltAz(
        1 * u.deg, 2 * u.deg, location=EarthLocation(0 * u.deg, 0 * u.deg, 0 * u.m)
    )
    assert isinstance(aa3.location, EarthLocation)


def test_hadec_attributes():
    hd = HADec(1 * u.hourangle, 2 * u.deg)
    assert hd.ha == 1.0 * u.hourangle
    assert hd.dec == 2 * u.deg
    assert hd.obstime is None
    assert hd.location is None

    hd2 = HADec(
        23 * u.hourangle,
        -2 * u.deg,
        obstime="J2000",
        location=EarthLocation(0 * u.deg, 0 * u.deg, 0 * u.m),
    )
    assert_allclose(hd2.ha, -1 * u.hourangle)
    assert hd2.dec == -2 * u.deg
    assert hd2.obstime == Time("J2000")
    assert isinstance(hd2.location, EarthLocation)

    sr = hd2.represent_as(r.SphericalRepresentation)
    assert_allclose(sr.lon, -1 * u.hourangle)


def test_itrs_earth_location():
    loc = EarthLocation(lat=0 * u.deg, lon=0 * u.deg, height=0 * u.m)
    sat = EarthLocation(
        lat=-24.6609379 * u.deg, lon=160.34199789 * u.deg, height=420.17927591 * u.km
    )

    itrs_geo = sat.get_itrs()
    eloc = itrs_geo.earth_location
    assert_allclose(sat.lon, eloc.lon)
    assert_allclose(sat.lat, eloc.lat)
    assert_allclose(sat.height, eloc.height)

    topo_itrs_repr = itrs_geo.cartesian - loc.get_itrs().cartesian
    itrs_topo = ITRS(topo_itrs_repr, location=loc)
    eloc = itrs_topo.earth_location
    assert_allclose(sat.lon, eloc.lon)
    assert_allclose(sat.lat, eloc.lat)
    assert_allclose(sat.height, eloc.height)

    obstime = Time("J2010")  # Anything different from default
    topo_itrs_repr2 = sat.get_itrs(obstime).cartesian - loc.get_itrs(obstime).cartesian
    itrs_topo2 = ITRS(topo_itrs_repr2, location=loc, obstime=obstime)
    eloc2 = itrs_topo2.earth_location
    assert_allclose(sat.lon, eloc2.lon)
    assert_allclose(sat.lat, eloc2.lat)
    assert_allclose(sat.height, eloc2.height)

    wgs84 = ITRS(325 * u.deg, 2 * u.deg, representation_type="wgs84geodetic")
    assert wgs84.lon == 325 * u.deg
    assert wgs84.lat == 2 * u.deg
    assert wgs84.height == 0.0 * u.m


def test_representation():
    """
    Test the getter and setter properties for `representation`

    """
    # Create the frame object.
    icrs = ICRS(ra=1 * u.deg, dec=1 * u.deg)
    data = icrs.data

    # Create some representation objects.
    icrs_cart = icrs.cartesian
    icrs_spher = icrs.spherical
    icrs_cyl = icrs.cylindrical

    # Testing when `_representation` set to `CartesianRepresentation`.
    icrs.representation_type = r.CartesianRepresentation

    assert icrs.representation_type == r.CartesianRepresentation
    assert icrs_cart.x == icrs.x
    assert icrs_cart.y == icrs.y
    assert icrs_cart.z == icrs.z
    assert icrs.data == data

    # Testing that an ICRS object in CartesianRepresentation must not have spherical attributes.
    for attr in ("ra", "dec", "distance"):
        with pytest.raises(AttributeError) as err:
            getattr(icrs, attr)
        assert "object has no attribute" in str(err.value)

    # Testing when `_representation` set to `CylindricalRepresentation`.
    icrs.representation_type = r.CylindricalRepresentation

    assert icrs.representation_type == r.CylindricalRepresentation
    assert icrs.data == data

    # Testing setter input using text argument for spherical.
    icrs.representation_type = "spherical"

    assert icrs.representation_type is r.SphericalRepresentation
    assert icrs_spher.lat == icrs.dec
    assert icrs_spher.lon == icrs.ra
    assert icrs_spher.distance == icrs.distance
    assert icrs.data == data

    # Testing that an ICRS object in SphericalRepresentation must not have cartesian attributes.
    for attr in ("x", "y", "z"):
        with pytest.raises(AttributeError) as err:
            getattr(icrs, attr)
        assert "object has no attribute" in str(err.value)

    # Testing setter input using text argument for cylindrical.
    icrs.representation_type = "cylindrical"

    assert icrs.representation_type is r.CylindricalRepresentation
    assert icrs_cyl.rho == icrs.rho
    assert icrs_cyl.phi == icrs.phi
    assert icrs_cyl.z == icrs.z
    assert icrs.data == data

    # Testing that an ICRS object in CylindricalRepresentation must not have spherical attributes.
    for attr in ("ra", "dec", "distance"):
        with pytest.raises(AttributeError) as err:
            getattr(icrs, attr)
        assert "object has no attribute" in str(err.value)

    with pytest.raises(ValueError) as err:
        icrs.representation_type = "WRONG"
    assert "but must be a BaseRepresentation class" in str(err.value)

    with pytest.raises(ValueError) as err:
        icrs.representation_type = ICRS
    assert "but must be a BaseRepresentation class" in str(err.value)


def test_represent_as():
    icrs = ICRS(ra=1 * u.deg, dec=1 * u.deg)

    cart1 = icrs.represent_as("cartesian")
    cart2 = icrs.represent_as(r.CartesianRepresentation)

    assert cart1.x == cart2.x
    assert cart1.y == cart2.y
    assert cart1.z == cart2.z

    # now try with velocities
    icrs = ICRS(
        ra=0 * u.deg,
        dec=0 * u.deg,
        distance=10 * u.kpc,
        pm_ra_cosdec=0 * u.mas / u.yr,
        pm_dec=0 * u.mas / u.yr,
        radial_velocity=1 * u.km / u.s,
    )

    # single string
    rep2 = icrs.represent_as("cylindrical")
    assert isinstance(rep2, r.CylindricalRepresentation)
    assert isinstance(rep2.differentials["s"], r.CylindricalDifferential)

    # single class with positional in_frame_units, verify that warning raised
    with pytest.warns(AstropyWarning, match="argument position") as w:
        icrs.represent_as(r.CylindricalRepresentation, False)
    assert len(w) == 1

    # TODO: this should probably fail in the future once we figure out a better
    # workaround for dealing with UnitSphericalRepresentation's with
    # RadialDifferential's
    # two classes
    # rep2 = icrs.represent_as(r.CartesianRepresentation,
    #                          r.SphericalCosLatDifferential)
    # assert isinstance(rep2, r.CartesianRepresentation)
    # assert isinstance(rep2.differentials['s'], r.SphericalCosLatDifferential)

    with pytest.raises(ValueError):
        icrs.represent_as("odaigahara")


def test_shorthand_representations():
    rep = r.CartesianRepresentation([1, 2, 3] * u.pc)
    dif = r.CartesianDifferential([1, 2, 3] * u.km / u.s)
    rep = rep.with_differentials(dif)

    icrs = ICRS(rep)

    cyl = icrs.cylindrical
    assert isinstance(cyl, r.CylindricalRepresentation)
    assert isinstance(cyl.differentials["s"], r.CylindricalDifferential)

    sph = icrs.spherical
    assert isinstance(sph, r.SphericalRepresentation)
    assert isinstance(sph.differentials["s"], r.SphericalDifferential)

    sph = icrs.sphericalcoslat
    assert isinstance(sph, r.SphericalRepresentation)
    assert isinstance(sph.differentials["s"], r.SphericalCosLatDifferential)


def test_equal():
    obstime = "B1955"
    sc1 = FK4([1, 2] * u.deg, [3, 4] * u.deg, obstime=obstime)
    sc2 = FK4([1, 20] * u.deg, [3, 4] * u.deg, obstime=obstime)

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
    sc1 = FK4([1, 2] * u.deg, [3, 4] * u.deg, radial_velocity=[1, 2] * u.km / u.s)
    sc2 = FK4([1, 2] * u.deg, [3, 4] * u.deg, radial_velocity=[1, 20] * u.km / u.s)

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

    assert (FK4() == ICRS()) is False
    assert (FK4() == FK4(obstime="J1999")) is False


def test_equal_exceptions():
    # Shape mismatch
    sc1 = FK4([1, 2, 3] * u.deg, [3, 4, 5] * u.deg)
    with pytest.raises(ValueError, match="cannot compare: shape mismatch"):
        sc1 == sc1[:2]  # noqa: B015

    # Different representation_type
    sc1 = FK4(1, 2, 3, representation_type="cartesian")
    sc2 = FK4(1 * u.deg, 2 * u.deg, 2, representation_type="spherical")
    with pytest.raises(
        TypeError,
        match=(
            "cannot compare: objects must have same "
            "class: CartesianRepresentation vs. SphericalRepresentation"
        ),
    ):
        sc1 == sc2  # noqa: B015

    # Different differential type
    sc1 = FK4(1 * u.deg, 2 * u.deg, radial_velocity=1 * u.km / u.s)
    sc2 = FK4(
        1 * u.deg, 2 * u.deg, pm_ra_cosdec=1 * u.mas / u.yr, pm_dec=1 * u.mas / u.yr
    )
    with pytest.raises(
        TypeError,
        match=(
            "cannot compare: objects must have same "
            "class: RadialDifferential vs. UnitSphericalCosLatDifferential"
        ),
    ):
        sc1 == sc2  # noqa: B015

    # Different frame attribute
    sc1 = FK5(1 * u.deg, 2 * u.deg)
    sc2 = FK5(1 * u.deg, 2 * u.deg, equinox="J1999")
    with pytest.raises(
        TypeError,
        match=r"cannot compare: objects must have equivalent "
        r"frames: <FK5 Frame \(equinox=J2000.000\)> "
        r"vs. <FK5 Frame \(equinox=J1999.000\)>",
    ):
        sc1 == sc2  # noqa: B015

    # Different frame
    sc1 = FK4(1 * u.deg, 2 * u.deg)
    sc2 = FK5(1 * u.deg, 2 * u.deg, equinox="J2000")
    with pytest.raises(
        TypeError,
        match="cannot compare: objects must have equivalent "
        r"frames: <FK4 Frame \(equinox=B1950.000, obstime=B1950.000\)> "
        r"vs. <FK5 Frame \(equinox=J2000.000\)>",
    ):
        sc1 == sc2  # noqa: B015

    sc1 = FK4(1 * u.deg, 2 * u.deg)
    sc2 = FK4()
    with pytest.raises(
        ValueError, match="cannot compare: one frame has data and the other does not"
    ):
        sc1 == sc2  # noqa: B015
    with pytest.raises(
        ValueError, match="cannot compare: one frame has data and the other does not"
    ):
        sc2 == sc1  # noqa: B015


def test_dynamic_attrs():
    c = ICRS(1 * u.deg, 2 * u.deg)
    assert "ra" in dir(c)
    assert "dec" in dir(c)

    with pytest.raises(AttributeError) as err:
        c.blahblah
    assert "object has no attribute 'blahblah'" in str(err.value)

    with pytest.raises(AttributeError) as err:
        c.ra = 1
    assert "Cannot set any frame attribute" in str(err.value)

    c.blahblah = 1
    assert c.blahblah == 1


def test_nodata_error():
    i = ICRS()
    with pytest.raises(ValueError) as excinfo:
        i.data

    assert "does not have associated data" in str(excinfo.value)


def test_len0_data():
    i = ICRS([] * u.deg, [] * u.deg)
    assert i.has_data
    repr(i)


def test_quantity_attributes():
    # make sure we can create a GCRS frame with valid inputs
    GCRS(obstime="J2002", obsgeoloc=[1, 2, 3] * u.km, obsgeovel=[4, 5, 6] * u.km / u.s)

    # make sure it fails for invalid lovs or vels
    with pytest.raises(TypeError):
        GCRS(obsgeoloc=[1, 2, 3])  # no unit
    with pytest.raises(u.UnitsError):
        GCRS(obsgeoloc=[1, 2, 3] * u.km / u.s)  # incorrect unit
    with pytest.raises(ValueError):
        GCRS(obsgeoloc=[1, 3] * u.km)  # incorrect shape


def test_quantity_attribute_default():
    # The default default (yes) is None:
    class MyCoord(BaseCoordinateFrame):
        someval = QuantityAttribute(unit=u.deg)

    frame = MyCoord()
    assert frame.someval is None

    frame = MyCoord(someval=15 * u.deg)
    assert u.isclose(frame.someval, 15 * u.deg)

    # This should work if we don't explicitly pass in a unit, but we pass in a
    # default value with a unit
    class MyCoord2(BaseCoordinateFrame):
        someval = QuantityAttribute(15 * u.deg)

    frame = MyCoord2()
    assert u.isclose(frame.someval, 15 * u.deg)

    # Since here no shape was given, we can set to any shape we like.
    frame = MyCoord2(someval=np.ones(3) * u.deg)
    assert frame.someval.shape == (3,)
    assert np.all(frame.someval == 1 * u.deg)

    # We should also be able to insist on a given shape.
    class MyCoord3(BaseCoordinateFrame):
        someval = QuantityAttribute(unit=u.arcsec, shape=(3,))

    frame = MyCoord3(someval=np.ones(3) * u.deg)
    assert frame.someval.shape == (3,)
    assert frame.someval.unit == u.arcsec
    assert u.allclose(frame.someval.value, 3600.0)

    # The wrong shape raises.
    with pytest.raises(ValueError, match="shape"):
        MyCoord3(someval=1.0 * u.deg)

    # As does the wrong unit.
    with pytest.raises(u.UnitsError):
        MyCoord3(someval=np.ones(3) * u.m)

    # We are allowed a short-cut for zero.
    frame0 = MyCoord3(someval=0)
    assert frame0.someval.shape == (3,)
    assert frame0.someval.unit == u.arcsec
    assert np.all(frame0.someval.value == 0.0)

    # But not if it has the wrong shape.
    with pytest.raises(ValueError, match="shape"):
        MyCoord3(someval=np.zeros(2))

    # This should fail, if we don't pass in a default or a unit
    with pytest.raises(ValueError):

        class MyCoord(BaseCoordinateFrame):
            someval = QuantityAttribute()


def test_eloc_attributes():
    el = EarthLocation(lon=12.3 * u.deg, lat=45.6 * u.deg, height=1 * u.km)
    it = ITRS(
        r.SphericalRepresentation(lon=12.3 * u.deg, lat=45.6 * u.deg, distance=1 * u.km)
    )
    gc = GCRS(ra=12.3 * u.deg, dec=45.6 * u.deg, distance=6375 * u.km)

    el1 = AltAz(location=el).location
    assert isinstance(el1, EarthLocation)
    # these should match *exactly* because the EarthLocation
    assert el1.lat == el.lat
    assert el1.lon == el.lon
    assert el1.height == el.height

    el2 = AltAz(location=it).location
    assert isinstance(el2, EarthLocation)
    # these should *not* match because giving something in Spherical ITRS is
    # *not* the same as giving it as an EarthLocation: EarthLocation is on an
    # elliptical geoid. So the longitude should match (because flattening is
    # only along the z-axis), but latitude should not. Also, height is relative
    # to the *surface* in EarthLocation, but the ITRS distance is relative to
    # the center of the Earth
    assert not allclose(el2.lat, it.spherical.lat)
    assert allclose(el2.lon, it.spherical.lon)
    assert el2.height < -6000 * u.km

    el3 = AltAz(location=gc).location
    # GCRS inputs implicitly get transformed to ITRS and then onto
    # EarthLocation's elliptical geoid. So both lat and lon shouldn't match
    assert isinstance(el3, EarthLocation)
    assert not allclose(el3.lat, gc.dec)
    assert not allclose(el3.lon, gc.ra)
    assert np.abs(el3.height) < 500 * u.km


def test_equivalent_frames():
    i = ICRS()
    i2 = ICRS(1 * u.deg, 2 * u.deg)
    assert i.is_equivalent_frame(i)
    assert i.is_equivalent_frame(i2)
    with pytest.raises(TypeError):
        assert i.is_equivalent_frame(10)
    with pytest.raises(TypeError):
        assert i2.is_equivalent_frame(SkyCoord(i2))

    f0 = FK5()  # this J2000 is TT
    f1 = FK5(equinox="J2000")
    f2 = FK5(1 * u.deg, 2 * u.deg, equinox="J2000")
    f3 = FK5(equinox="J2010")
    f4 = FK4(equinox="J2010")

    assert f1.is_equivalent_frame(f1)
    assert not i.is_equivalent_frame(f1)
    assert f0.is_equivalent_frame(f1)
    assert f1.is_equivalent_frame(f2)
    assert not f1.is_equivalent_frame(f3)
    assert not f3.is_equivalent_frame(f4)

    aa1 = AltAz()
    aa2 = AltAz(obstime="J2010")

    assert aa2.is_equivalent_frame(aa2)
    assert not aa1.is_equivalent_frame(i)
    assert not aa1.is_equivalent_frame(aa2)


def test_equivalent_frame_coordinateattribute():
    class FrameWithCoordinateAttribute(BaseCoordinateFrame):
        coord_attr = CoordinateAttribute(HCRS)

    # These frames should not be considered equivalent
    f0 = FrameWithCoordinateAttribute()
    f1 = FrameWithCoordinateAttribute(
        coord_attr=HCRS(1 * u.deg, 2 * u.deg, obstime="J2000")
    )
    f2 = FrameWithCoordinateAttribute(
        coord_attr=HCRS(3 * u.deg, 4 * u.deg, obstime="J2000")
    )
    f3 = FrameWithCoordinateAttribute(
        coord_attr=HCRS(1 * u.deg, 2 * u.deg, obstime="J2001")
    )

    assert not f0.is_equivalent_frame(f1)
    assert not f1.is_equivalent_frame(f0)
    assert not f1.is_equivalent_frame(f2)
    assert not f1.is_equivalent_frame(f3)
    assert not f2.is_equivalent_frame(f3)

    # They each should still be equivalent to a deep copy of themselves
    assert f0.is_equivalent_frame(deepcopy(f0))
    assert f1.is_equivalent_frame(deepcopy(f1))
    assert f2.is_equivalent_frame(deepcopy(f2))
    assert f3.is_equivalent_frame(deepcopy(f3))


def test_equivalent_frame_locationattribute():
    class FrameWithLocationAttribute(BaseCoordinateFrame):
        loc_attr = EarthLocationAttribute()

    # These frames should not be considered equivalent
    f0 = FrameWithLocationAttribute()
    location = EarthLocation(lat=-34, lon=19, height=300)
    f1 = FrameWithLocationAttribute(loc_attr=location)

    assert not f0.is_equivalent_frame(f1)
    assert not f1.is_equivalent_frame(f0)

    # They each should still be equivalent to a deep copy of themselves
    assert f0.is_equivalent_frame(deepcopy(f0))
    assert f1.is_equivalent_frame(deepcopy(f1))


def test_representation_subclass():
    # Regression test for #3354

    # Normally when instantiating a frame without a distance the frame will try
    # and use UnitSphericalRepresentation internally instead of
    # SphericalRepresentation.
    frame = FK5(
        representation_type=r.SphericalRepresentation, ra=32 * u.deg, dec=20 * u.deg
    )
    assert type(frame._data) == r.UnitSphericalRepresentation
    assert frame.representation_type == r.SphericalRepresentation

    # If using a SphericalRepresentation class this used to not work, so we
    # test here that this is now fixed.
    class NewSphericalRepresentation(r.SphericalRepresentation):
        attr_classes = r.SphericalRepresentation.attr_classes

    frame = FK5(
        representation_type=NewSphericalRepresentation, lon=32 * u.deg, lat=20 * u.deg
    )
    assert type(frame._data) == r.UnitSphericalRepresentation
    assert frame.representation_type == NewSphericalRepresentation

    # A similar issue then happened in __repr__ with subclasses of
    # SphericalRepresentation.
    assert (
        repr(frame)
        == "<FK5 Coordinate (equinox=J2000.000): (lon, lat) in deg\n    (32., 20.)>"
    )

    # A more subtle issue is when specifying a custom
    # UnitSphericalRepresentation subclass for the data and
    # SphericalRepresentation or a subclass for the representation.

    class NewUnitSphericalRepresentation(r.UnitSphericalRepresentation):
        attr_classes = r.UnitSphericalRepresentation.attr_classes

        def __repr__(self):
            return "<NewUnitSphericalRepresentation: spam spam spam>"

    frame = FK5(
        NewUnitSphericalRepresentation(lon=32 * u.deg, lat=20 * u.deg),
        representation_type=NewSphericalRepresentation,
    )

    assert repr(frame) == "<FK5 Coordinate (equinox=J2000.000):  spam spam spam>"


def test_getitem_representation():
    """
    Make sure current representation survives __getitem__ even if different
    from data representation.
    """
    c = ICRS([1, 1] * u.deg, [2, 2] * u.deg)
    c.representation_type = "cartesian"
    assert c[0].representation_type is r.CartesianRepresentation


def test_component_error_useful():
    """
    Check that a data-less frame gives useful error messages about not having
    data when the attributes asked for are possible coordinate components
    """

    i = ICRS()

    with pytest.raises(ValueError) as excinfo:
        i.ra
    assert "does not have associated data" in str(excinfo.value)

    with pytest.raises(AttributeError) as excinfo1:
        i.foobar
    with pytest.raises(AttributeError) as excinfo2:
        i.lon  # lon is *not* the component name despite being the underlying representation's name
    assert "object has no attribute 'foobar'" in str(excinfo1.value)
    assert "object has no attribute 'lon'" in str(excinfo2.value)


def test_cache_clear():
    i = ICRS(1 * u.deg, 2 * u.deg)

    # Add an in frame units version of the rep to the cache.
    repr(i)

    assert len(i.cache["representation"]) == 2

    i.cache.clear()

    assert len(i.cache["representation"]) == 0


def test_inplace_array():
    i = ICRS([[1, 2], [3, 4]] * u.deg, [[10, 20], [30, 40]] * u.deg)

    # Add an in frame units version of the rep to the cache.
    repr(i)

    # Check that repr() has added a rep to the cache
    assert len(i.cache["representation"]) == 2

    # Modify the data
    i.data.lon[:, 0] = [100, 200] * u.deg

    # Clear the cache
    i.cache.clear()

    # This will use a second (potentially cached rep)
    assert_allclose(i.ra, [[100, 2], [200, 4]] * u.deg)
    assert_allclose(i.dec, [[10, 20], [30, 40]] * u.deg)


def test_inplace_change():
    i = ICRS(1 * u.deg, 2 * u.deg)

    # Add an in frame units version of the rep to the cache.
    repr(i)

    # Check that repr() has added a rep to the cache
    assert len(i.cache["representation"]) == 2

    # Modify the data
    i.data.lon[()] = 10 * u.deg

    # Clear the cache
    i.cache.clear()

    # This will use a second (potentially cached rep)
    assert i.ra == 10 * u.deg
    assert i.dec == 2 * u.deg


def test_representation_with_multiple_differentials():
    dif1 = r.CartesianDifferential([1, 2, 3] * u.km / u.s)
    dif2 = r.CartesianDifferential([1, 2, 3] * u.km / u.s**2)
    rep = r.CartesianRepresentation(
        [1, 2, 3] * u.pc, differentials={"s": dif1, "s2": dif2}
    )

    # check warning is raised for a scalar
    with pytest.raises(ValueError):
        ICRS(rep)


def test_missing_component_error_names():
    """
    This test checks that the component names are frame component names, not
    representation or differential names, when referenced in an exception raised
    when not passing in enough data. For example:

    ICRS(ra=10*u.deg)

    should state:

    TypeError: __init__() missing 1 required positional argument: 'dec'

    """
    with pytest.raises(TypeError) as e:
        ICRS(ra=150 * u.deg)
    assert "missing 1 required positional argument: 'dec'" in str(e.value)

    with pytest.raises(TypeError) as e:
        ICRS(
            ra=150 * u.deg,
            dec=-11 * u.deg,
            pm_ra=100 * u.mas / u.yr,
            pm_dec=10 * u.mas / u.yr,
        )
    assert "pm_ra_cosdec" in str(e.value)


def test_non_spherical_representation_unit_creation(unitphysics):  # noqa: F811
    class PhysicsICRS(ICRS):
        default_representation = r.PhysicsSphericalRepresentation

    pic = PhysicsICRS(phi=1 * u.deg, theta=25 * u.deg, r=1 * u.kpc)
    assert isinstance(pic.data, r.PhysicsSphericalRepresentation)

    picu = PhysicsICRS(phi=1 * u.deg, theta=25 * u.deg)
    assert isinstance(picu.data, unitphysics)


def test_attribute_repr():
    class Spam:
        def _astropy_repr_in_frame(self):
            return "TEST REPR"

    class TestFrame(BaseCoordinateFrame):
        attrtest = Attribute(default=Spam())

    assert "TEST REPR" in repr(TestFrame())


def test_component_names_repr():
    # Frame class with new component names that includes a name swap
    class NameChangeFrame(BaseCoordinateFrame):
        default_representation = r.PhysicsSphericalRepresentation

        frame_specific_representation_info = {
            r.PhysicsSphericalRepresentation: [
                RepresentationMapping("phi", "theta", u.deg),
                RepresentationMapping("theta", "phi", u.arcsec),
                RepresentationMapping("r", "JUSTONCE", u.AU),
            ]
        }

    frame = NameChangeFrame(0 * u.deg, 0 * u.arcsec, 0 * u.AU)

    # Check for the new names in the Frame repr
    assert "(theta, phi, JUSTONCE)" in repr(frame)

    # Check that the letter "r" has not been replaced more than once in the Frame repr
    assert repr(frame).count("JUSTONCE") == 1


def test_galactocentric_defaults():
    with galactocentric_frame_defaults.set("pre-v4.0"):
        galcen_pre40 = Galactocentric()

    with galactocentric_frame_defaults.set("v4.0"):
        galcen_40 = Galactocentric()

    with galactocentric_frame_defaults.set("latest"):
        galcen_latest = Galactocentric()

    # parameters that changed
    assert not u.allclose(galcen_pre40.galcen_distance, galcen_40.galcen_distance)
    assert not u.allclose(galcen_pre40.z_sun, galcen_40.z_sun)

    for k in galcen_40.frame_attributes:
        if isinstance(getattr(galcen_40, k), BaseCoordinateFrame):
            continue  # skip coordinate comparison...

        elif isinstance(getattr(galcen_40, k), CartesianDifferential):
            assert u.allclose(
                getattr(galcen_40, k).d_xyz, getattr(galcen_latest, k).d_xyz
            )
        else:
            assert getattr(galcen_40, k) == getattr(galcen_latest, k)

    # test validate Galactocentric
    with galactocentric_frame_defaults.set("latest"):
        params = galactocentric_frame_defaults.validate(galcen_latest)
        references = galcen_latest.frame_attribute_references
        state = {"parameters": params, "references": references}

        assert galactocentric_frame_defaults.parameters == params
        assert galactocentric_frame_defaults.references == references
        assert galactocentric_frame_defaults._state == state

    # Test not one of accepted parameter types
    with pytest.raises(ValueError):
        galactocentric_frame_defaults.validate(ValueError)

    # test parameters property
    assert (
        galactocentric_frame_defaults.parameters
        == galactocentric_frame_defaults.parameters
    )


def test_galactocentric_references():
    # references in the "scientific paper"-sense

    with galactocentric_frame_defaults.set("pre-v4.0"):
        galcen_pre40 = Galactocentric()

        for k in galcen_pre40.frame_attributes:
            if k == "roll":  # no reference for this parameter
                continue

            assert k in galcen_pre40.frame_attribute_references

    with galactocentric_frame_defaults.set("v4.0"):
        galcen_40 = Galactocentric()

        for k in galcen_40.frame_attributes:
            if k == "roll":  # no reference for this parameter
                continue

            assert k in galcen_40.frame_attribute_references

    with galactocentric_frame_defaults.set("v4.0"):
        galcen_custom = Galactocentric(z_sun=15 * u.pc)

        for k in galcen_custom.frame_attributes:
            if k == "roll":  # no reference for this parameter
                continue

            if k == "z_sun":
                assert k not in galcen_custom.frame_attribute_references
            else:
                assert k in galcen_custom.frame_attribute_references


def test_coordinateattribute_transformation():
    class FrameWithCoordinateAttribute(BaseCoordinateFrame):
        coord_attr = CoordinateAttribute(HCRS)

    hcrs = HCRS(1 * u.deg, 2 * u.deg, 3 * u.AU, obstime="2001-02-03")
    f1_frame = FrameWithCoordinateAttribute(coord_attr=hcrs)
    f1_skycoord = FrameWithCoordinateAttribute(coord_attr=SkyCoord(hcrs))

    # The input is already HCRS, so the frame attribute should not change it
    assert f1_frame.coord_attr == hcrs
    # The output should not be different if a SkyCoord is provided
    assert f1_skycoord.coord_attr == f1_frame.coord_attr

    gcrs = GCRS(4 * u.deg, 5 * u.deg, 6 * u.AU, obstime="2004-05-06")
    f2_frame = FrameWithCoordinateAttribute(coord_attr=gcrs)
    f2_skycoord = FrameWithCoordinateAttribute(coord_attr=SkyCoord(gcrs))

    # The input needs to be converted from GCRS to HCRS
    assert isinstance(f2_frame.coord_attr, HCRS)
    # The `obstime` frame attribute should have been "merged" in a SkyCoord-style transformation
    assert f2_frame.coord_attr.obstime == gcrs.obstime
    # The output should not be different if a SkyCoord is provided
    assert f2_skycoord.coord_attr == f2_frame.coord_attr


def test_realize_frame_accepts_kwargs():
    c1 = ICRS(
        x=1 * u.pc,
        y=2 * u.pc,
        z=3 * u.pc,
        representation_type=r.CartesianRepresentation,
    )
    new_data = r.CartesianRepresentation(x=11 * u.pc, y=12 * u.pc, z=13 * u.pc)

    c2 = c1.realize_frame(new_data, representation_type="cartesian")
    c3 = c1.realize_frame(new_data, representation_type="cylindrical")

    assert c2.representation_type == r.CartesianRepresentation
    assert c3.representation_type == r.CylindricalRepresentation


def test_nameless_frame_subclass():
    """Note: this is a regression test for #11096"""

    class Test:
        pass

    # Subclass from a frame class and a non-frame class.
    # This subclassing is the test!
    class NewFrame(ICRS, Test):
        pass


def test_frame_coord_comparison():
    """Test that frame can be compared to a SkyCoord"""
    frame = ICRS(0 * u.deg, 0 * u.deg)
    coord = SkyCoord(frame)
    other = SkyCoord(ICRS(0 * u.deg, 1 * u.deg))

    assert frame == coord
    assert frame != other
    assert not (frame == other)
    error_msg = "objects must have equivalent frames"
    with pytest.raises(TypeError, match=error_msg):
        frame == SkyCoord(AltAz("0d", "1d"))  # noqa: B015

    coord = SkyCoord(ra=12 * u.hourangle, dec=5 * u.deg, frame=FK5(equinox="J1950"))
    frame = FK5(ra=12 * u.hourangle, dec=5 * u.deg, equinox="J2000")
    with pytest.raises(TypeError, match=error_msg):
        coord == frame  # noqa: B015

    frame = ICRS()
    coord = SkyCoord(0 * u.deg, 0 * u.deg, frame=frame)
    error_msg = "Can only compare SkyCoord to Frame with data"
    with pytest.raises(ValueError, match=error_msg):
        frame == coord  # noqa: B015


@pytest.mark.parametrize(
    ["s1", "s2"],
    (
        ((1,), (1,)),
        ((2,), (1,)),
        ((1,), (2,)),
        ((2,), (2,)),
        ((2, 1), (1,)),
        ((1,), (2, 1)),
        ((2, 1), (1, 3)),
    ),
)
def test_altaz_broadcast(s1, s2):
    """Note: Regression test for #5982"""
    where = EarthLocation.from_geodetic(lat=45 * u.deg, lon=30 * u.deg, height=0 * u.m)
    time = Time(np.full(s1, 58000.0), format="mjd")
    angle = np.full(s2, 45.0) * u.deg
    result = AltAz(alt=angle, az=angle, obstime=time, location=where)
    assert result.shape == np.broadcast_shapes(s1, s2)


def test_transform_altaz_array_obstime():
    """Note: Regression test for #12965"""
    obstime = Time("2010-01-01T00:00:00")
    location = EarthLocation(0 * u.deg, 0 * u.deg, 0 * u.m)

    frame1 = AltAz(location=location, obstime=obstime)
    coord1 = SkyCoord(alt=80 * u.deg, az=0 * u.deg, frame=frame1)

    obstimes = obstime + np.linspace(0, 15, 50) * u.min
    frame2 = AltAz(location=location, obstime=obstimes)
    coord2 = SkyCoord(alt=coord1.alt, az=coord1.az, frame=frame2)
    assert np.all(coord2.alt == 80 * u.deg)
    assert np.all(coord2.az == 0 * u.deg)
    assert coord2.shape == (50,)

    # test transformation to ICRS works
    assert len(coord2.icrs) == 50


def test_spherical_offsets_by_broadcast():
    """Note: Regression test for #14383"""
    assert SkyCoord(
        ra=np.array([123, 134, 145]), dec=np.array([45, 56, 67]), unit=u.deg
    ).spherical_offsets_by(2 * u.deg, 2 * u.deg).shape == (3,)


@pytest.mark.parametrize("shape", [(1,), (2,)])
def test_spherical_offsets_with_wrap(shape):
    # see https://github.com/astropy/astropy/issues/16219
    sc = SkyCoord(ra=np.broadcast_to(123.0, shape), dec=90.0, unit=u.deg)
    scop = sc.spherical_offsets_by(+2 * u.deg, 0 * u.deg)
    assert scop.shape == shape

    scom = sc.spherical_offsets_by(-2 * u.deg, 0 * u.deg)
    assert scom.shape == shape
