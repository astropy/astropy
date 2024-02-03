# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test geodetic representations"""
from copy import deepcopy

import pytest

from astropy import units as u
from astropy.coordinates.representation import (
    DUPLICATE_REPRESENTATIONS,
    REPRESENTATION_CLASSES,
    BaseBodycentricRepresentation,
    BaseGeodeticRepresentation,
    CartesianRepresentation,
    GRS80GeodeticRepresentation,
    WGS72GeodeticRepresentation,
    WGS84GeodeticRepresentation,
)
from astropy.coordinates.representation.geodetic import ELLIPSOIDS
from astropy.tests.helper import assert_quantity_allclose
from astropy.units.tests.test_quantity_erfa_ufuncs import vvd


class TestCustomGeodeticRepresentations:
    @classmethod
    def setup_class(self):
        # Preserve the original REPRESENTATION_CLASSES dict so that importing
        # the test file doesn't add a persistent test subclass (CustomGeodetic, etc.)
        self.REPRESENTATION_CLASSES_ORIG = deepcopy(REPRESENTATION_CLASSES)
        self.DUPLICATE_REPRESENTATIONS_ORIG = deepcopy(DUPLICATE_REPRESENTATIONS)

        class CustomGeodetic(BaseGeodeticRepresentation):
            _flattening = 0.01832
            _equatorial_radius = 4000000.0 * u.m

        class CustomSphericGeodetic(BaseGeodeticRepresentation):
            _flattening = 0.0
            _equatorial_radius = 4000000.0 * u.m

        class CustomSphericBodycentric(BaseBodycentricRepresentation):
            _flattening = 0.0
            _equatorial_radius = 4000000.0 * u.m

        class IAUMARS2000GeodeticRepresentation(BaseGeodeticRepresentation):
            _equatorial_radius = 3396190.0 * u.m
            _flattening = 0.5886007555512007 * u.percent

        class IAUMARS2000BodycentricRepresentation(BaseBodycentricRepresentation):
            _equatorial_radius = 3396190.0 * u.m
            _flattening = 0.5886007555512007 * u.percent

        self.CustomGeodetic = CustomGeodetic
        self.CustomSphericGeodetic = CustomSphericGeodetic
        self.CustomSphericBodycentric = CustomSphericBodycentric
        self.IAUMARS2000GeodeticRepresentation = IAUMARS2000GeodeticRepresentation
        self.IAUMARS2000BodycentricRepresentation = IAUMARS2000BodycentricRepresentation

    @classmethod
    def teardown_class(self):
        REPRESENTATION_CLASSES.clear()
        REPRESENTATION_CLASSES.update(self.REPRESENTATION_CLASSES_ORIG)
        DUPLICATE_REPRESENTATIONS.clear()
        DUPLICATE_REPRESENTATIONS.update(self.DUPLICATE_REPRESENTATIONS_ORIG)

    def get_representation(self, representation):
        if isinstance(representation, str):
            return getattr(self, representation)
        else:
            return representation

    def test_geodetic_bodycentric_equivalence_spherical_bodies(self):
        initial_cartesian = CartesianRepresentation(
            x=[1, 3000.0] * u.km, y=[7000.0, 4.0] * u.km, z=[5.0, 6000.0] * u.km
        )

        gd_transformed = self.CustomSphericGeodetic.from_representation(
            initial_cartesian
        )
        bc_transformed = self.CustomSphericBodycentric.from_representation(
            initial_cartesian
        )
        assert_quantity_allclose(gd_transformed.lon, bc_transformed.lon)
        assert_quantity_allclose(gd_transformed.lat, bc_transformed.lat)
        assert_quantity_allclose(gd_transformed.height, bc_transformed.height)

    @pytest.mark.parametrize(
        "geodeticrepresentation",
        [
            "CustomGeodetic",
            WGS84GeodeticRepresentation,
            "IAUMARS2000GeodeticRepresentation",
            "IAUMARS2000BodycentricRepresentation",
        ],
    )
    def test_cartesian_geodetic_roundtrip(self, geodeticrepresentation):
        geodeticrepresentation = self.get_representation(geodeticrepresentation)
        # Test array-valued input in the process.
        initial_cartesian = CartesianRepresentation(
            x=[1, 3000.0] * u.km, y=[7000.0, 4.0] * u.km, z=[5.0, 6000.0] * u.km
        )

        transformed = geodeticrepresentation.from_representation(initial_cartesian)

        roundtripped = CartesianRepresentation.from_representation(transformed)

        assert_quantity_allclose(initial_cartesian.x, roundtripped.x)
        assert_quantity_allclose(initial_cartesian.y, roundtripped.y)
        assert_quantity_allclose(initial_cartesian.z, roundtripped.z)

    @pytest.mark.parametrize(
        "geodeticrepresentation",
        [
            "CustomGeodetic",
            WGS84GeodeticRepresentation,
            "IAUMARS2000GeodeticRepresentation",
            "IAUMARS2000BodycentricRepresentation",
        ],
    )
    def test_geodetic_cartesian_roundtrip(self, geodeticrepresentation):
        geodeticrepresentation = self.get_representation(geodeticrepresentation)
        initial_geodetic = geodeticrepresentation(
            lon=[0.8, 1.3] * u.radian,
            lat=[0.3, 0.98] * u.radian,
            height=[100.0, 367.0] * u.m,
        )

        transformed = CartesianRepresentation.from_representation(initial_geodetic)

        roundtripped = geodeticrepresentation.from_representation(transformed)

        assert_quantity_allclose(initial_geodetic.lon, roundtripped.lon)
        assert_quantity_allclose(initial_geodetic.lat, roundtripped.lat)
        assert_quantity_allclose(initial_geodetic.height, roundtripped.height)

    def test_geocentric_to_geodetic(self):
        """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""
        # Here, test the chain.  Direct conversion from Cartesian to
        # various Geodetic representations is done indirectly in test_earth.
        x, y, z = (2e6, 3e6, 5.244e6)

        status = 0  # help for copy & paste of vvd

        gc = CartesianRepresentation(x, y, z, u.m)
        gd = WGS84GeodeticRepresentation.from_cartesian(gc)
        e, p, h = gd.lon.to(u.radian), gd.lat.to(u.radian), gd.height.to(u.m)
        vvd(e, 0.9827937232473290680, 1e-14, "eraGc2gd", "e1", status)
        vvd(p, 0.97160184819075459, 1e-14, "eraGc2gd", "p1", status)
        vvd(h, 331.4172461426059892, 1e-8, "eraGc2gd", "h1", status)

        gd = gd.represent_as(GRS80GeodeticRepresentation)
        e, p, h = gd.lon.to(u.radian), gd.lat.to(u.radian), gd.height.to(u.m)
        vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
        vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
        vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

        gd = gd.represent_as(WGS72GeodeticRepresentation)
        e, p, h = gd.lon.to(u.radian), gd.lat.to(u.radian), gd.height.to(u.m)
        vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e3", status)
        vvd(p, 0.97160181811015119, 1e-14, "eraGc2gd", "p3", status)
        vvd(h, 333.27707261303181, 1e-8, "eraGc2gd", "h3", status)

    def test_geodetic_to_geocentric(self):
        """Test that we reproduce erfa/src/t_erfa_c.c t_gd2gc"""
        # These tests are also done implicitly in test_earth.py.
        e = 3.1 * u.rad
        p = -0.5 * u.rad
        h = 2500.0 * u.m

        status = 0  # help for copy & paste of vvd

        gd = WGS84GeodeticRepresentation(e, p, h)
        xyz = gd.to_cartesian().get_xyz()
        vvd(xyz[0], -5599000.5577049947, 1e-7, "eraGd2gc", "0/1", status)
        vvd(xyz[1], 233011.67223479203, 1e-7, "eraGd2gc", "1/1", status)
        vvd(xyz[2], -3040909.4706983363, 1e-7, "eraGd2gc", "2/1", status)

        gd = GRS80GeodeticRepresentation(e, p, h)
        xyz = gd.to_cartesian().get_xyz()
        vvd(xyz[0], -5599000.5577260984, 1e-7, "eraGd2gc", "0/2", status)
        vvd(xyz[1], 233011.6722356703, 1e-7, "eraGd2gc", "1/2", status)
        vvd(xyz[2], -3040909.4706095476, 1e-7, "eraGd2gc", "2/2", status)

        gd = WGS72GeodeticRepresentation(e, p, h)
        xyz = gd.to_cartesian().get_xyz()
        vvd(xyz[0], -5598998.7626301490, 1e-7, "eraGd2gc", "0/3", status)
        vvd(xyz[1], 233011.5975297822, 1e-7, "eraGd2gc", "1/3", status)
        vvd(xyz[2], -3040908.6861467111, 1e-7, "eraGd2gc", "2/3", status)

    @pytest.mark.parametrize(
        "representation",
        [WGS84GeodeticRepresentation, "IAUMARS2000BodycentricRepresentation"],
    )
    def test_default_height_is_zero(self, representation):
        representation = self.get_representation(representation)
        gd = representation(10 * u.deg, 20 * u.deg)
        assert gd.lon == 10 * u.deg
        assert gd.lat == 20 * u.deg
        assert gd.height == 0 * u.m

    @pytest.mark.parametrize(
        "representation",
        [WGS84GeodeticRepresentation, "IAUMARS2000BodycentricRepresentation"],
    )
    def test_non_angle_error(self, representation):
        representation = self.get_representation(representation)
        with pytest.raises(u.UnitTypeError, match="require units equivalent to 'rad'"):
            representation(20 * u.m, 20 * u.deg, 20 * u.m)

    @pytest.mark.parametrize(
        "representation",
        [WGS84GeodeticRepresentation, "IAUMARS2000BodycentricRepresentation"],
    )
    def test_non_length_error(self, representation):
        representation = self.get_representation(representation)
        with pytest.raises(u.UnitTypeError, match="units of length"):
            representation(10 * u.deg, 20 * u.deg, 30)

    def test_subclass_bad_ellipsoid(self):
        # Test incomplete initialization.

        msg = "module 'erfa' has no attribute 'foo'"
        with pytest.raises(AttributeError, match=msg):

            class InvalidCustomEllipsoid(BaseGeodeticRepresentation):
                _ellipsoid = "foo"

        assert "foo" not in ELLIPSOIDS
        assert "invalidcustomellipsoid" not in REPRESENTATION_CLASSES

    @pytest.mark.parametrize(
        "baserepresentation",
        [BaseGeodeticRepresentation, BaseBodycentricRepresentation],
    )
    def test_geodetic_subclass_missing_equatorial_radius(self, baserepresentation):
        msg = "'_equatorial_radius' and '_flattening'."
        with pytest.raises(AttributeError, match=msg):

            class MissingCustomAttribute(baserepresentation):
                _flattening = 0.075 * u.dimensionless_unscaled

        assert "missingcustomattribute" not in REPRESENTATION_CLASSES
