# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test that Distribution works with classes other than ndarray and Quantity."""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import astropy.units as u
from astropy.coordinates import (
    Angle,
    CartesianDifferential,
    CartesianRepresentation,
    EarthLocation,
    Latitude,
    Longitude,
    SkyCoord,
    SphericalDifferential,
    SphericalRepresentation,
)
from astropy.coordinates.representation import (
    DIFFERENTIAL_CLASSES,
    REPRESENTATION_CLASSES,
)
from astropy.coordinates.tests.test_representation import (
    components_allclose,
    representation_equal,
)
from astropy.uncertainty import Distribution


def assert_representation_equal(rep1, rep2):
    assert np.all(representation_equal(rep1, rep2))


def assert_representation_allclose(rep1, rep2):
    result = True
    if type(rep1) is not type(rep2):
        return False
    if getattr(rep1, "_differentials", False):
        if rep1._differentials.keys() != rep2._differentials.keys():
            return False
        for key, diff1 in rep1._differentials.items():
            result &= components_allclose(diff1, rep2._differentials[key])
    elif getattr(rep2, "_differentials", False):
        return False

    return np.all(result & components_allclose(rep1, rep2))


class TestAngles:
    @classmethod
    def setup_class(cls):
        cls.a = np.arange(27.0).reshape(3, 9)
        cls.d = Distribution(cls.a)
        cls.q = cls.a << u.deg
        cls.dq = Distribution(cls.q)

    @pytest.mark.parametrize("angle_cls", [Angle, Longitude, Latitude])
    def test_as_input_for_angle(self, angle_cls):
        da = angle_cls(self.dq)
        assert isinstance(da, angle_cls)
        assert isinstance(da, Distribution)
        assert_array_equal(da.distribution, angle_cls(self.q))

    @pytest.mark.parametrize("angle_cls", [Angle, Longitude, Latitude])
    def test_using_angle_as_input(self, angle_cls):
        a = angle_cls(self.q)
        da = Distribution(a)
        assert isinstance(da, angle_cls)
        assert isinstance(da, Distribution)

    # Parametrize the unit to check the various branches in Latitude._validate_angles
    @pytest.mark.parametrize("dtype", ["f8", "f4"])
    @pytest.mark.parametrize(
        "value", [90 * u.deg, np.pi / 2 * u.radian, 90 * 60 * u.arcmin]
    )
    def test_at_limit_for_latitude(self, value, dtype):
        q = u.Quantity(value, dtype=dtype).reshape(1)
        qd = Distribution(q)
        ld = Latitude(qd)
        assert_array_equal(ld.distribution, Latitude(q))

    # Parametrize the unit in case Longitude._wrap_at becomes unit-dependent.
    @pytest.mark.parametrize("dtype", ["f8", "f4"])
    @pytest.mark.parametrize(
        "value", [360 * u.deg, 2 * np.pi * u.radian, 360 * 60 * u.arcmin]
    )
    def test_at_wrap_angle_for_longitude(self, value, dtype):
        q = u.Quantity(value, dtype=dtype).reshape(1)
        qd = Distribution(q)
        ld = Longitude(qd)
        assert_array_equal(ld.distribution, Longitude(q))
        assert np.all(ld.distribution == 0)

    @pytest.mark.parametrize("angle_cls", [Longitude, Latitude])
    def test_operation_gives_correct_subclass(self, angle_cls):
        # Lon and Lat always fall back to Angle
        da = angle_cls(self.dq)
        da2 = da + da
        assert isinstance(da, Angle)
        assert isinstance(da, Distribution)

    @pytest.mark.parametrize("angle_cls", [Longitude, Latitude])
    def test_pdfstd_gives_correct_subclass(self, angle_cls):
        # Lon and Lat always fall back to Angle
        da = angle_cls(self.dq)
        std = da.pdf_std()
        assert isinstance(std, Angle)
        assert_array_equal(std, Angle(self.q.std(-1)))

    def test_earthlocation_geocentric_distribution(self):
        x = y = z = self.a << u.km

        eloc = EarthLocation.from_geocentric(x=x, y=y, z=z)

        xd = Distribution(x)
        yd = Distribution(y)
        zd = Distribution(z)
        deloc = EarthLocation.from_geocentric(x=xd, y=yd, z=zd)

        assert isinstance(deloc.x, Distribution)
        assert_array_equal(np.median(eloc.x, axis=1), deloc.x.pdf_median())

    def test_earthlocation_geodetic_distribution(self):
        h = self.a << u.km
        eloc = EarthLocation.from_geodetic(lon=self.q, lat=self.q, height=h)

        hq = Distribution(h)
        deloc = EarthLocation.from_geodetic(lon=self.dq, lat=self.dq, height=hq)

        assert isinstance(deloc.x, Distribution)
        assert_array_equal(np.median(eloc.x, axis=1), deloc.x.pdf_median())


class TestRepresentation:
    @classmethod
    def setup_class(cls):
        cls.lon = Distribution(np.linspace(0.0, 360 * u.deg, 10, endpoint=False))
        cls.lat = Angle([-45.0, 0.0, 45.0], u.deg)
        cls.r = 6000 * u.km  # Sort of OK for Geodetic representations.
        cls.sph = SphericalRepresentation(
            cls.lon.distribution, cls.lat[:, np.newaxis], cls.r
        )
        cls.dsph = SphericalRepresentation(cls.lon, cls.lat, cls.r)

    def get_distribution(self, rep):
        return rep._apply(lambda x: getattr(x, "distribution", x[..., np.newaxis]))

    def test_cartesian(self):
        dcart = self.dsph.to_cartesian()
        cart = self.sph.to_cartesian()
        assert isinstance(dcart.x, Distribution)
        assert_array_equal(dcart.x.distribution, cart.x)
        assert_array_equal(dcart.y.distribution, cart.y)
        assert_array_equal(dcart.z.distribution, cart.z)

    def test_cartesian_roundtrip(self):
        roundtrip = SphericalRepresentation.from_cartesian(self.dsph.to_cartesian())
        assert_representation_allclose(self.get_distribution(roundtrip), self.sph)

    @pytest.mark.parametrize("rep_cls", REPRESENTATION_CLASSES.values())
    def test_other_reps(self, rep_cls):
        drep = self.dsph.represent_as(rep_cls)
        rep = self.sph.represent_as(rep_cls)
        assert_representation_equal(self.get_distribution(drep), rep)


class TestRepresentationWithDifferential(TestRepresentation):
    @classmethod
    def setup_class(cls):
        super().setup_class()
        cls.d_lon = Distribution(np.ones(10) << u.deg / u.hour)
        cls.d_lat = np.ones(cls.lat.shape) << u.deg / u.hour
        cls.d_r = 1 * u.m / u.s
        cls.d_sph = SphericalDifferential(
            cls.d_lon.distribution, cls.d_lat[:, np.newaxis], cls.d_r
        )
        cls.dd_sph = SphericalDifferential(cls.d_lon, cls.d_lat, cls.d_r)
        cls.sph = cls.sph.with_differentials({"s": cls.d_sph})
        cls.dsph = cls.dsph.with_differentials({"s": cls.dd_sph})

    def test_cartesian(self):
        dcart = self.dsph.represent_as(CartesianRepresentation, CartesianDifferential)
        cart = self.sph.represent_as(CartesianRepresentation, CartesianDifferential)
        assert isinstance(dcart.x, Distribution)
        assert_array_equal(dcart.x.distribution, cart.x)
        assert_array_equal(dcart.y.distribution, cart.y)
        assert_array_equal(dcart.z.distribution, cart.z)
        assert "s" in dcart._differentials
        dd_cart = dcart._differentials["s"]
        d_cart = cart._differentials["s"]
        assert_array_equal(dd_cart.d_x.distribution, d_cart.d_x)
        assert_array_equal(dd_cart.d_y.distribution, d_cart.d_y)
        assert_array_equal(dd_cart.d_z.distribution, d_cart.d_z)

    def test_cartesian_roundtrip(self):
        roundtrip = self.dsph.represent_as(
            CartesianRepresentation, CartesianDifferential
        ).represent_as(SphericalRepresentation, SphericalDifferential)
        assert_representation_allclose(self.get_distribution(roundtrip), self.sph)

    @pytest.mark.parametrize("diff_cls", DIFFERENTIAL_CLASSES.values())
    def test_other_reps(self, diff_cls):
        repr_cls = diff_cls.base_representation
        drep = self.dsph.represent_as(repr_cls, diff_cls)
        rep = self.sph.represent_as(repr_cls, diff_cls)
        assert_representation_equal(self.get_distribution(drep), rep)


class TestSkyCoord:
    @classmethod
    def setup_class(cls):
        cls.ra = Distribution(np.linspace(0.0, 360 * u.deg, 10, endpoint=False))
        cls.dec = Angle([-45.0, 0.0, 45.0], u.deg)
        cls.dsc = SkyCoord(cls.ra, cls.dec)
        cls.sc = SkyCoord(cls.ra.distribution, cls.dec[:, np.newaxis])

    def test_init(self):
        assert isinstance(self.dsc.ra, Distribution)
        assert not isinstance(self.dsc.dec, Distribution)

    @pytest.mark.parametrize("frame", ["fk5", "galactic"])
    def test_convert(self, frame):
        dconv = self.dsc.transform_to(frame)
        conv = self.sc.transform_to(frame)
        # Coordinate names are different for FK5 and galactic, so
        # just get them from the representation.
        assert_array_equal(dconv.data.lon.distribution, conv.data.lon)
        assert_array_equal(dconv.data.lat.distribution, conv.data.lat)
