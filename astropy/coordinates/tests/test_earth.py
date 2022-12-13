# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test initialization of angles not already covered by the API tests"""

import pickle

import numpy as np
import pytest

from astropy import constants
from astropy import units as u
from astropy.coordinates.angles import Latitude, Longitude
from astropy.coordinates.earth import ELLIPSOIDS, EarthLocation
from astropy.coordinates.name_resolve import NameResolveError
from astropy.time import Time
from astropy.units import allclose as quantity_allclose


def allclose_m14(a, b, rtol=1.0e-14, atol=None):
    if atol is None:
        atol = 1.0e-14 * getattr(a, "unit", 1)
    return quantity_allclose(a, b, rtol, atol)


def allclose_m8(a, b, rtol=1.0e-8, atol=None):
    if atol is None:
        atol = 1.0e-8 * getattr(a, "unit", 1)
    return quantity_allclose(a, b, rtol, atol)


def isclose_m14(val, ref):
    return np.array([allclose_m14(v, r) for (v, r) in zip(val, ref)])


def isclose_m8(val, ref):
    return np.array([allclose_m8(v, r) for (v, r) in zip(val, ref)])


def vvd(val, valok, dval, func, test, status):
    """Mimic routine of erfa/src/t_erfa_c.c (to help copy & paste)"""
    assert quantity_allclose(val, valok * val.unit, atol=dval * val.unit)


def test_gc2gd():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""
    x, y, z = (2e6, 3e6, 5.244e6)

    status = 0  # help for copy & paste of vvd

    location = EarthLocation.from_geocentric(x, y, z, u.m)
    e, p, h = location.to_geodetic("WGS84")
    e, p, h = e.to(u.radian), p.to(u.radian), h.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
    vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
    vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

    e, p, h = location.to_geodetic("GRS80")
    e, p, h = e.to(u.radian), p.to(u.radian), h.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
    vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
    vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

    e, p, h = location.to_geodetic("WGS72")
    e, p, h = e.to(u.radian), p.to(u.radian), h.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e3", status)
    vvd(p, 0.97160181811015119, 1e-14, "eraGc2gd", "p3", status)
    vvd(h, 333.27707261303181, 1e-8, "eraGc2gd", "h3", status)


def test_gd2gc():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gd2gc"""
    e = 3.1 * u.rad
    p = -0.5 * u.rad
    h = 2500.0 * u.m

    status = 0  # help for copy & paste of vvd

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid="WGS84")
    xyz = tuple(v.to(u.m) for v in location.to_geocentric())
    vvd(xyz[0], -5599000.5577049947, 1e-7, "eraGd2gc", "0/1", status)
    vvd(xyz[1], 233011.67223479203, 1e-7, "eraGd2gc", "1/1", status)
    vvd(xyz[2], -3040909.4706983363, 1e-7, "eraGd2gc", "2/1", status)

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid="GRS80")
    xyz = tuple(v.to(u.m) for v in location.to_geocentric())
    vvd(xyz[0], -5599000.5577260984, 1e-7, "eraGd2gc", "0/2", status)
    vvd(xyz[1], 233011.6722356703, 1e-7, "eraGd2gc", "1/2", status)
    vvd(xyz[2], -3040909.4706095476, 1e-7, "eraGd2gc", "2/2", status)

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid="WGS72")
    xyz = tuple(v.to(u.m) for v in location.to_geocentric())
    vvd(xyz[0], -5598998.7626301490, 1e-7, "eraGd2gc", "0/3", status)
    vvd(xyz[1], 233011.5975297822, 1e-7, "eraGd2gc", "1/3", status)
    vvd(xyz[2], -3040908.6861467111, 1e-7, "eraGd2gc", "2/3", status)


class TestInput:
    def setup_method(self):
        self.lon = Longitude(
            [0.0, 45.0, 90.0, 135.0, 180.0, -180, -90, -45],
            u.deg,
            wrap_angle=180 * u.deg,
        )
        self.lat = Latitude([+0.0, 30.0, 60.0, +90.0, -90.0, -60.0, -30.0, 0.0], u.deg)
        self.h = u.Quantity([0.1, 0.5, 1.0, -0.5, -1.0, +4.2, -11.0, -0.1], u.m)
        self.location = EarthLocation.from_geodetic(self.lon, self.lat, self.h)
        self.x, self.y, self.z = self.location.to_geocentric()

    def test_default_ellipsoid(self):
        assert self.location.ellipsoid == EarthLocation._ellipsoid

    def test_geo_attributes(self):
        assert all(
            np.all(_1 == _2)
            for _1, _2 in zip(self.location.geodetic, self.location.to_geodetic())
        )
        assert all(
            np.all(_1 == _2)
            for _1, _2 in zip(self.location.geocentric, self.location.to_geocentric())
        )

    def test_attribute_classes(self):
        """Test that attribute classes are correct (and not EarthLocation)"""
        assert type(self.location.x) is u.Quantity
        assert type(self.location.y) is u.Quantity
        assert type(self.location.z) is u.Quantity
        assert type(self.location.lon) is Longitude
        assert type(self.location.lat) is Latitude
        assert type(self.location.height) is u.Quantity

    def test_input(self):
        """Check input is parsed correctly"""

        # units of length should be assumed geocentric
        geocentric = EarthLocation(self.x, self.y, self.z)
        assert np.all(geocentric == self.location)
        geocentric2 = EarthLocation(
            self.x.value, self.y.value, self.z.value, self.x.unit
        )
        assert np.all(geocentric2 == self.location)
        geodetic = EarthLocation(self.lon, self.lat, self.h)
        assert np.all(geodetic == self.location)
        geodetic2 = EarthLocation(
            self.lon.to_value(u.degree),
            self.lat.to_value(u.degree),
            self.h.to_value(u.m),
        )
        assert np.all(geodetic2 == self.location)
        geodetic3 = EarthLocation(self.lon, self.lat)
        assert allclose_m14(geodetic3.lon.value, self.location.lon.value)
        assert allclose_m14(geodetic3.lat.value, self.location.lat.value)
        assert not np.any(
            isclose_m14(geodetic3.height.value, self.location.height.value)
        )
        geodetic4 = EarthLocation(self.lon, self.lat, self.h[-1])
        assert allclose_m14(geodetic4.lon.value, self.location.lon.value)
        assert allclose_m14(geodetic4.lat.value, self.location.lat.value)
        assert allclose_m14(geodetic4.height[-1].value, self.location.height[-1].value)
        assert not np.any(
            isclose_m14(geodetic4.height[:-1].value, self.location.height[:-1].value)
        )
        # check length unit preservation
        geocentric5 = EarthLocation(self.x, self.y, self.z, u.pc)
        assert geocentric5.unit is u.pc
        assert geocentric5.x.unit is u.pc
        assert geocentric5.height.unit is u.pc
        assert allclose_m14(geocentric5.x.to_value(self.x.unit), self.x.value)
        geodetic5 = EarthLocation(self.lon, self.lat, self.h.to(u.pc))
        assert geodetic5.unit is u.pc
        assert geodetic5.x.unit is u.pc
        assert geodetic5.height.unit is u.pc
        assert allclose_m14(geodetic5.x.to_value(self.x.unit), self.x.value)

    def test_invalid_input(self):
        """Check invalid input raises exception"""
        # incomprehensible by either raises TypeError
        with pytest.raises(TypeError):
            EarthLocation(self.lon, self.y, self.z)

        # wrong units
        with pytest.raises(u.UnitsError):
            EarthLocation.from_geocentric(self.lon, self.lat, self.lat)
        # inconsistent units
        with pytest.raises(u.UnitsError):
            EarthLocation.from_geocentric(self.h, self.lon, self.lat)
        # floats without a unit
        with pytest.raises(TypeError):
            EarthLocation.from_geocentric(self.x.value, self.y.value, self.z.value)
        # inconsistent shape
        with pytest.raises(ValueError):
            EarthLocation.from_geocentric(self.x, self.y, self.z[:5])

        # inconsistent units
        with pytest.raises(u.UnitsError):
            EarthLocation.from_geodetic(self.x, self.y, self.z)
        # inconsistent shape
        with pytest.raises(ValueError):
            EarthLocation.from_geodetic(self.lon, self.lat, self.h[:5])

    def test_slicing(self):
        # test on WGS72 location, so we can check the ellipsoid is passed on
        locwgs72 = EarthLocation.from_geodetic(
            self.lon, self.lat, self.h, ellipsoid="WGS72"
        )
        loc_slice1 = locwgs72[4]
        assert isinstance(loc_slice1, EarthLocation)
        assert loc_slice1.unit is locwgs72.unit
        assert loc_slice1.ellipsoid == locwgs72.ellipsoid == "WGS72"
        assert not loc_slice1.shape
        with pytest.raises(TypeError):
            loc_slice1[0]
        with pytest.raises(IndexError):
            len(loc_slice1)

        loc_slice2 = locwgs72[4:6]
        assert isinstance(loc_slice2, EarthLocation)
        assert len(loc_slice2) == 2
        assert loc_slice2.unit is locwgs72.unit
        assert loc_slice2.ellipsoid == locwgs72.ellipsoid
        assert loc_slice2.shape == (2,)
        loc_x = locwgs72["x"]
        assert type(loc_x) is u.Quantity
        assert loc_x.shape == locwgs72.shape
        assert loc_x.unit is locwgs72.unit

    def test_invalid_ellipsoid(self):
        # unknown ellipsoid
        with pytest.raises(ValueError):
            EarthLocation.from_geodetic(self.lon, self.lat, self.h, ellipsoid="foo")
        with pytest.raises(TypeError):
            EarthLocation(self.lon, self.lat, self.h, ellipsoid="foo")

        with pytest.raises(ValueError):
            self.location.ellipsoid = "foo"

        with pytest.raises(ValueError):
            self.location.to_geodetic("foo")

    @pytest.mark.parametrize("ellipsoid", ELLIPSOIDS)
    def test_ellipsoid(self, ellipsoid):
        """Test that different ellipsoids are understood, and differ"""
        # check that heights differ for different ellipsoids
        # need different tolerance, since heights are relative to ~6000 km
        lon, lat, h = self.location.to_geodetic(ellipsoid)
        if ellipsoid == self.location.ellipsoid:
            assert allclose_m8(h.value, self.h.value)
        else:
            # Some heights are very similar for some; some lon, lat identical.
            assert not np.all(isclose_m8(h.value, self.h.value))

        # given lon, lat, height, check that x,y,z differ
        location = EarthLocation.from_geodetic(
            self.lon, self.lat, self.h, ellipsoid=ellipsoid
        )
        if ellipsoid == self.location.ellipsoid:
            assert allclose_m14(location.z.value, self.z.value)
        else:
            assert not np.all(isclose_m14(location.z.value, self.z.value))

        def test_to_value(self):
            loc = self.location
            loc_ndarray = loc.view(np.ndarray)
            assert np.all(loc.value == loc_ndarray)
            loc2 = self.location.to(u.km)
            loc2_ndarray = np.empty_like(loc_ndarray)
            for coo in "x", "y", "z":
                loc2_ndarray[coo] = loc_ndarray[coo] / 1000.0
            assert np.all(loc2.value == loc2_ndarray)
            loc2_value = self.location.to_value(u.km)
            assert np.all(loc2_value == loc2_ndarray)


def test_pickling():
    """Regression test against #4304."""
    el = EarthLocation(0.0 * u.m, 6000 * u.km, 6000 * u.km)
    s = pickle.dumps(el)
    el2 = pickle.loads(s)
    assert el == el2


def test_repr_latex():
    """
    Regression test for issue #4542
    """
    somelocation = EarthLocation(lon="149:3:57.9", lat="-31:16:37.3")
    somelocation._repr_latex_()
    somelocation2 = EarthLocation(lon=[1.0, 2.0] * u.deg, lat=[-1.0, 9.0] * u.deg)
    somelocation2._repr_latex_()


@pytest.mark.remote_data
# TODO: this parametrize should include a second option with a valid Google API
# key. For example, we should make an API key for Astropy, and add it to GitHub Actions
# as an environment variable (for security).
@pytest.mark.parametrize("google_api_key", [None])
def test_of_address(google_api_key):
    NYC_lon = -74.0 * u.deg
    NYC_lat = 40.7 * u.deg
    # ~10 km tolerance to address difference between OpenStreetMap and Google
    # for "New York, NY". This doesn't matter in practice because this test is
    # only used to verify that the query succeeded, not that the returned
    # position is precise.
    NYC_tol = 0.1 * u.deg

    # just a location
    try:
        loc = EarthLocation.of_address("New York, NY")
    except NameResolveError as e:
        # API limit might surface even here in CI.
        if "unknown failure with" not in str(e):
            pytest.xfail(str(e))
    else:
        assert quantity_allclose(loc.lat, NYC_lat, atol=NYC_tol)
        assert quantity_allclose(loc.lon, NYC_lon, atol=NYC_tol)
        assert np.allclose(loc.height.value, 0.0)

    # Put this one here as buffer to get around Google map API limit per sec.
    # no match: This always raises NameResolveError
    with pytest.raises(NameResolveError):
        EarthLocation.of_address("lkjasdflkja")

    if google_api_key is not None:
        # a location and height
        try:
            loc = EarthLocation.of_address("New York, NY", get_height=True)
        except NameResolveError as e:
            # Buffer above sometimes insufficient to get around API limit but
            # we also do not want to drag things out with time.sleep(0.195),
            # where 0.195 was empirically determined on some physical machine.
            pytest.xfail(str(e.value))
        else:
            assert quantity_allclose(loc.lat, NYC_lat, atol=NYC_tol)
            assert quantity_allclose(loc.lon, NYC_lon, atol=NYC_tol)
            assert quantity_allclose(loc.height, 10.438 * u.meter, atol=1.0 * u.cm)


def test_geodetic_tuple():
    lat = 2 * u.deg
    lon = 10 * u.deg
    height = 100 * u.m

    el = EarthLocation.from_geodetic(lat=lat, lon=lon, height=height)

    res1 = el.to_geodetic()
    res2 = el.geodetic

    assert res1.lat == res2.lat and quantity_allclose(res1.lat, lat)
    assert res1.lon == res2.lon and quantity_allclose(res1.lon, lon)
    assert res1.height == res2.height and quantity_allclose(res1.height, height)


def test_gravitational_redshift():
    someloc = EarthLocation(lon=-87.7 * u.deg, lat=37 * u.deg)
    sometime = Time("2017-8-21 18:26:40")
    zg0 = someloc.gravitational_redshift(sometime)

    # should be of order ~few mm/s change per week
    zg_week = someloc.gravitational_redshift(sometime + 7 * u.day)
    assert 1.0 * u.mm / u.s < abs(zg_week - zg0) < 1 * u.cm / u.s

    # ~cm/s over a half-year
    zg_halfyear = someloc.gravitational_redshift(sometime + 0.5 * u.yr)
    assert 1 * u.cm / u.s < abs(zg_halfyear - zg0) < 1 * u.dm / u.s

    # but when back to the same time in a year, should be tenths of mm
    # even over decades
    zg_year = someloc.gravitational_redshift(sometime - 20 * u.year)
    assert 0.1 * u.mm / u.s < abs(zg_year - zg0) < 1 * u.mm / u.s

    # Check mass adjustments.
    # If Jupiter and the moon are ignored, effect should be off by ~ .5 mm/s
    masses = {
        "sun": constants.G * constants.M_sun,
        "jupiter": 0 * constants.G * u.kg,
        "moon": 0 * constants.G * u.kg,
    }
    zg_moonjup = someloc.gravitational_redshift(sometime, masses=masses)
    assert 0.1 * u.mm / u.s < abs(zg_moonjup - zg0) < 1 * u.mm / u.s
    # Check that simply not including the bodies gives the same result.
    assert zg_moonjup == someloc.gravitational_redshift(sometime, bodies=("sun",))
    # And that earth can be given, even not as last argument
    assert zg_moonjup == someloc.gravitational_redshift(
        sometime, bodies=("earth", "sun")
    )

    # If the earth is also ignored, effect should be off by ~ 20 cm/s
    # This also tests the conversion of kg to gravitational units.
    masses["earth"] = 0 * u.kg
    zg_moonjupearth = someloc.gravitational_redshift(sometime, masses=masses)
    assert 1 * u.dm / u.s < abs(zg_moonjupearth - zg0) < 1 * u.m / u.s

    # If all masses are zero, redshift should be 0 as well.
    masses["sun"] = 0 * u.kg
    assert someloc.gravitational_redshift(sometime, masses=masses) == 0

    with pytest.raises(KeyError):
        someloc.gravitational_redshift(sometime, bodies=("saturn",))

    with pytest.raises(u.UnitsError):
        masses = {
            "sun": constants.G * constants.M_sun,
            "jupiter": constants.G * constants.M_jup,
            "moon": 1 * u.km,  # wrong units!
            "earth": constants.G * constants.M_earth,
        }
        someloc.gravitational_redshift(sometime, masses=masses)


def test_read_only_input():
    lon = np.array([80.0, 440.0]) * u.deg
    lat = np.array([45.0]) * u.deg
    lon.flags.writeable = lat.flags.writeable = False
    loc = EarthLocation.from_geodetic(lon=lon, lat=lat)
    assert quantity_allclose(loc[1].x, loc[0].x)


def test_info():
    EarthLocation._get_site_registry(force_builtin=True)
    greenwich = EarthLocation.of_site("greenwich")
    assert str(greenwich.info).startswith("name = Royal Observatory Greenwich")
