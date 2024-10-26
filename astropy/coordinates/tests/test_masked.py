# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import astropy.coordinates.representation as r
import astropy.units as u
from astropy.coordinates import FK5, SkyCoord
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.coordinates.tests.helper import skycoord_equal
from astropy.coordinates.tests.test_representation import representation_equal
from astropy.utils.masked import Masked


class MaskedSphericalSetup:
    @classmethod
    def setup_class(cls):
        cls.lon = np.array([0.0, 3.0, 6.0, 12.0, 15.0, 18.0]) << u.hourangle
        cls.lat = np.array([-15.0, 30.0, 60.0, -60.0, 89.0, -80.0]) << u.deg
        cls.dis = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0]) << u.pc
        cls.mask_lon = np.array([False, True, False, False, True, True])
        cls.mask_lat = np.array([False, False, True, False, True, True])
        cls.mask_dis = np.array([False, False, False, True, False, True])
        cls.mlon = Masked(cls.lon, cls.mask_lon)
        cls.mlat = Masked(cls.lat, cls.mask_lat)
        cls.mdis = Masked(cls.dis, cls.mask_dis)
        cls.msph = r.SphericalRepresentation(cls.mlon, cls.mlat, cls.mdis)
        cls.mask_ang = cls.mask_lon | cls.mask_lat
        cls.mask = cls.mask_ang | cls.mask_dis


class TestSphericalRepresentationSeparateMasks(MaskedSphericalSetup):
    """Tests for mask propagation for Spherical with separate masks."""

    def test_initialization(self):
        assert_array_equal(self.msph.lon.mask, self.mask_lon)
        assert_array_equal(self.msph.lat.mask, self.mask_lat)
        assert_array_equal(self.msph.distance.mask, self.mask_dis)
        assert_array_equal(self.msph.unmasked.lon, self.lon)
        assert_array_equal(self.msph.unmasked.lat, self.lat)
        assert_array_equal(self.msph.unmasked.distance, self.dis)

        assert_array_equal(self.msph.mask, self.mask)
        assert_array_equal(self.msph.get_mask(), self.mask)
        assert_array_equal(self.msph.get_mask("lon", "lat"), self.mask_ang)
        assert repr(self.msph) == (
            "<SphericalRepresentation (lon, lat, distance) in (hourangle, deg, pc)\n"
            "    [( 0., -15., 10.), (———,  30., 20.), ( 6.,  ———, 30.),\n"
            "     (12., -60., ———), (———,  ———, 50.), (———,  ———, ———)]>"
        )
        assert str(self.msph) == (
            "[( 0., -15., 10.), (———,  30., 20.), ( 6.,  ———, 30.), (12., -60., ———),\n"
            " (———,  ———, 50.), (———,  ———, ———)] (hourangle, deg, pc)"
        )

    def test_convert_to_cartesian(self):
        mcart = self.msph.represent_as(r.CartesianRepresentation)
        assert_array_equal(mcart.mask, self.mask)

    def test_convert_to_unit_spherical(self):
        musph = self.msph.represent_as(r.UnitSphericalRepresentation)
        assert_array_equal(musph.lon.mask, self.mask_lon)
        assert_array_equal(musph.lat.mask, self.mask_lat)
        assert_array_equal(musph.mask, self.mask_ang)

    def test_convert_to_radial(self):
        mrad = r.RadialRepresentation.from_representation(self.msph)
        assert_array_equal(mrad.mask, self.mask_dis)

    def test_convert_to_physics_spherical(self):
        mpsph = self.msph.represent_as(r.PhysicsSphericalRepresentation)
        assert_array_equal(mpsph.phi.mask, self.mask_lon)
        assert_array_equal(mpsph.theta.mask, self.mask_lat)
        assert_array_equal(mpsph.r.mask, self.mask_dis)
        assert_array_equal(mpsph.mask, self.mask)
        assert_array_equal(mpsph.get_mask("phi", "theta"), self.mask_ang)

    def test_set_mask(self):
        msph = self.msph.copy()
        msph[0] = np.ma.masked
        assert_array_equal(msph.mask, np.concatenate(([True], self.mask[1:])))
        msph[0] = np.ma.nomask
        assert_array_equal(msph.mask, self.mask)

    def test_setitem(self):
        msph = self.msph.copy()
        msph[0] = self.msph[1]
        assert_array_equal(msph.mask, np.concatenate(([True], self.mask[1:])))
        assert_array_equal(
            msph.unmasked.lon, np.concatenate((msph.lon[1:2], self.lon[1:]))
        )
        msph[0] = self.msph[0]
        assert_array_equal(msph.mask, self.mask)
        assert_array_equal(msph.unmasked.lon, self.lon)

    def test_set_masked_item_on_unmasked_instance(self):
        # Currently, the mask on items is *ignored*, just as it is for ndarray,
        # Quantity, and Time. In principle, this could be changed for containers
        # like representations and Time. See
        # https://github.com/astropy/astropy/pull/17016#issuecomment-2439607869
        sph = self.msph.unmasked.copy()
        sph[0] = self.msph[1]
        assert not sph.masked
        assert_array_equal(sph.mask, np.zeros_like(self.mask))
        assert_array_equal(sph.lon, np.concatenate((sph.lon[1:2], self.lon[1:])))
        sph[0] = self.msph[0].unmasked
        assert not sph.masked
        assert_array_equal(sph.mask, np.zeros_like(self.mask))
        assert_array_equal(sph.unmasked.lon, self.lon)

    def test_set_np_ma_masked_on_unmasked_instance(self):
        sph = self.msph.copy().unmasked
        sph[0] = np.ma.masked
        assert sph.masked
        assert_array_equal(
            sph.mask, np.concatenate(([True], np.zeros_like(self.mask[1:])))
        )
        assert_array_equal(sph.unmasked.lon, self.lon)
        sph[0] = np.ma.nomask
        assert sph.masked  # Does not get reset
        assert_array_equal(sph.mask, np.zeros_like(self.mask))
        assert_array_equal(sph.unmasked.lon, self.lon)

    def test_set_np_ma_nomasked_on_unmasked_instance(self):
        sph = self.msph.copy().unmasked
        sph[0] = np.ma.nomask
        assert not sph.masked
        assert_array_equal(sph.unmasked.lon, self.lon)

    def test_filled(self):
        unmasked = self.msph.unmasked
        sph = self.msph.filled(unmasked[1])
        expected = unmasked.copy()
        expected[self.mask] = unmasked[1]
        assert np.all(representation_equal(sph, expected))

    def test_filled_with_masked_value(self):
        # Filled ignores the mask (this will be true as long as __setitem__
        # ignores it; it may be a logical choice to actually use it).
        sph = self.msph.filled(self.msph[1])
        assert not sph.masked
        expected = self.msph.unmasked.copy()
        expected[self.mask] = self.msph.unmasked[1]
        assert np.all(representation_equal(sph, expected))

    def test_transform_keeps_distance_angular_masks(self):
        m = rotation_matrix(30.0, "x") * 2.0  # rotation and scale
        sph = self.msph.transform(m)
        assert sph.masked
        # Distance now also masked if angular coordinates were masked.
        assert_array_equal(sph.distance.mask, self.mask)
        # But angular coordinates just depend on the angular mask.
        assert_array_equal(sph.lon.mask, self.mask_ang)
        assert_array_equal(sph.lat.mask, self.mask_ang)
        assert_array_equal(self.msph.get_mask(), self.mask)
        assert_array_equal(self.msph.get_mask("lon", "lat"), self.mask_ang)

    def test_unmasked_representation_masked_differential(self):
        rv = np.arange(6.0) << u.km / u.s
        mask_rv = [True, False] * 3
        mrv = Masked(rv, mask_rv)
        mdiff = r.RadialDifferential(mrv)
        msph = r.SphericalRepresentation(
            self.lon,
            self.lat,
            self.dis,
            differentials={"s": mdiff},
        )
        assert msph.masked
        assert_array_equal(msph.lon.mask, False)
        assert_array_equal(msph.lat.mask, False)
        assert_array_equal(msph.distance.mask, False)
        assert_array_equal(msph.differentials["s"].d_distance.mask, mask_rv)
        sph = msph.unmasked
        assert not sph.masked
        assert not sph.differentials["s"].masked
        # Sanity checks on "with[out]_differentials"
        assert msph.without_differentials().masked
        sph2 = r.SphericalRepresentation(self.lon, self.lat, self.dis)
        assert not sph2.masked
        sph3 = sph2.with_differentials({"s": mdiff})
        assert sph3.masked

    def test_masked_representation_unmasked_differential(self):
        diff = r.RadialDifferential(np.arange(6.0) << u.km / u.s)
        msph = r.SphericalRepresentation(
            self.mlon,
            self.mlat,
            self.mdis,
            differentials={"s": diff},
        )
        assert msph.masked
        assert msph.differentials["s"].masked
        assert_array_equal(msph.differentials["s"].d_distance.mask, False)
        # Sanity check on using with_differentials.
        msph2 = self.msph.with_differentials({"s": diff})
        assert msph2.masked
        assert msph2.differentials["s"].masked


class TestFrame(MaskedSphericalSetup):
    """Tests that mask is calculated properly for frames, using FK5."""

    @classmethod
    def setup_class(cls):
        super().setup_class()
        cls.fk5 = FK5(cls.msph)

    def test_initialization_directly(self):
        d = Masked([50, 1.0] * u.kpc, mask=[False, True])
        fk5 = FK5([0, 30] * u.deg, [-10, 10] * u.deg, distance=d)
        assert fk5.masked
        assert_array_equal(fk5.ra.mask, [False, False])
        assert_array_equal(fk5.dec.mask, [False, False])
        assert_array_equal(fk5.distance.mask, [False, True])
        assert_array_equal(fk5.mask, [False, True])
        assert_array_equal(fk5.get_mask("ra", "dec"), False)
        assert "—" in repr(fk5)

    def test_class_initialization(self):
        assert_array_equal(self.fk5.ra.mask, self.mask_lon)
        assert_array_equal(self.fk5.dec.mask, self.mask_lat)
        assert_array_equal(self.fk5.distance.mask, self.mask_dis)

        assert_array_equal(self.fk5.mask, self.mask)
        assert_array_equal(self.fk5.get_mask("ra", "dec"), self.mask_ang)

        unmasked = self.fk5.unmasked
        assert_array_equal(unmasked.ra, self.lon)
        assert_array_equal(unmasked.dec, self.lat)
        assert_array_equal(unmasked.distance, self.dis)

    def test_cache_clearing(self):
        mfk5 = self.fk5.copy()
        assert_array_equal(mfk5.ra.mask, self.mask_lon)
        assert "—" in repr(mfk5)
        mfk5[...] = np.ma.nomask
        assert_array_equal(mfk5.data.mask, np.zeros(mfk5.shape, bool))
        assert_array_equal(mfk5.ra.mask, np.zeros(mfk5.shape, bool))
        assert "—" not in repr(mfk5)
        mfk5[...] = self.fk5
        assert_array_equal(mfk5.ra.mask, self.mask_lon)
        assert "—" in repr(mfk5)

    def test_cartesian(self):
        mcart = FK5(self.msph.represent_as(r.CartesianRepresentation))
        assert_array_equal(mcart.mask, self.mask)

    def test_unit_spherical(self):
        musph = FK5(self.msph.represent_as(r.UnitSphericalRepresentation))
        assert_array_equal(musph.mask, self.mask_ang)
        assert_array_equal(musph.get_mask(), self.mask_ang)

    def test_physics_spherical(self):
        mpsph = FK5(self.msph.represent_as(r.PhysicsSphericalRepresentation))
        assert_array_equal(mpsph.mask, self.mask)
        assert_array_equal(mpsph.get_mask("data.phi", "data.theta"), self.mask_ang)

    def test_get_mask(self):
        assert_array_equal(self.fk5.get_mask("ra"), self.mask_lon)
        assert_array_equal(self.fk5.get_mask("ra", "dec"), self.mask_ang)
        assert_array_equal(self.fk5.get_mask("data.lat"), self.mask_lat)
        assert_array_equal(self.fk5.get_mask("data.lat", "data.lon"), self.mask_ang)
        assert_array_equal(self.fk5.get_mask("cartesian"), self.mask)
        assert_array_equal(self.fk5.get_mask("equinox"), np.zeros(self.fk5.shape, bool))

    def test_unmasked_frame(self):
        fk5 = self.fk5.unmasked
        assert not fk5.masked
        assert_array_equal(fk5.mask, np.zeros(fk5.shape, bool))
        assert_array_equal(fk5.get_mask(), np.zeros(fk5.shape, bool))


def test_frame_without_data():
    fk5_no_data = FK5(equinox=["J2000", "J2001"])
    with pytest.raises(ValueError, match="does not have associated data"):
        fk5_no_data.masked
    with pytest.raises(ValueError, match="does not have associated data"):
        fk5_no_data.mask
    with pytest.raises(ValueError, match="does not have associated data"):
        fk5_no_data.get_mask()
    assert_array_equal(fk5_no_data.get_mask("equinox"), np.zeros((2,), bool))


class TestSkyCoord(TestFrame):
    """Tests that mask is calculated properly for SkyCoord.

    Note that this does all the tests from TestFrame, as well as a few
    specific to SkyCoord, i.e., that use attributes the frame does not have.
    """

    @classmethod
    def setup_class(cls):
        super().setup_class()
        cls.p = np.linspace(900, 1000, cls.msph.size) << u.hPa
        cls.mask_p = np.array([True, False, False, False, False, False])
        cls.mp = Masked(cls.p, cls.mask_p)
        # Ensure we have an attribute not associated with the frame.
        cls.fk5 = SkyCoord(cls.msph, frame="fk5", pressure=cls.mp)

    def test_non_frame_attribute(self):
        assert_array_equal(self.fk5.get_mask("pressure"), self.mask_p)
        assert_array_equal(
            self.fk5.get_mask("pressure", "data"), self.mask | self.mask_p
        )


class TestSkyCoordWithDifferentials:
    @classmethod
    def setup_class(cls):
        cls.ra = [0.0, 3.0, 6.0, 12.0, 15.0, 18.0] << u.hourangle
        cls.dec = [-15.0, 30.0, 60.0, -60.0, 89.0, -80.0] << u.deg
        cls.dis = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0] << u.pc
        cls.mask_dis = np.array([False, True, False, True, False, True])
        cls.pm_ra_cosdec = [1.0, 2.0, 3.0, -4.0, -5.0, -6.0] << (u.mas / u.yr)
        cls.mask_pm_ra_cosdec = np.array([False, False, True, False, False, True])
        cls.pm_dec = [-9.0, -7.0, 5.0, 3.0, 1.0, 0.0] << (u.mas / u.yr)
        cls.mask_pm_dec = np.array([False, False, True, True, False, True])
        cls.rv = [40.0, 50.0, 0.0, 0.0, -30.0, -10.0] << (u.km / u.s)
        cls.mask_rv = np.array([False, False, False, False, True, True])
        cls.mdis = Masked(cls.dis, cls.mask_dis)
        cls.mpm_ra_cosdec = Masked(cls.pm_ra_cosdec, cls.mask_pm_ra_cosdec)
        cls.mpm_dec = Masked(cls.pm_dec, cls.mask_pm_dec)
        cls.mrv = Masked(cls.rv, cls.mask_rv)
        cls.sc = SkyCoord(
            ra=cls.ra,
            dec=cls.dec,
            distance=cls.mdis,
            pm_ra_cosdec=cls.mpm_ra_cosdec,
            pm_dec=cls.mpm_dec,
            radial_velocity=cls.mrv,
        )
        cls.mask = cls.mask_dis | cls.mask_pm_ra_cosdec | cls.mask_pm_dec | cls.mask_rv

    def test_setup(self):
        assert self.sc.masked
        assert_array_equal(self.sc.ra.mask, False)
        assert_array_equal(self.sc.dec.mask, False)
        assert_array_equal(self.sc.distance.mask, self.mask_dis)
        assert_array_equal(self.sc.pm_ra_cosdec.mask, self.mask_pm_ra_cosdec)
        assert_array_equal(self.sc.pm_dec.mask, self.mask_pm_dec)
        assert_array_equal(self.sc.radial_velocity.mask, self.mask_rv)
        assert_array_equal(self.sc.mask, self.mask)

    def test_get_mask(self):
        assert_array_equal(self.sc.get_mask(), self.mask)
        assert_array_equal(self.sc.get_mask("ra", "dec"), False)

    def test_filled(self):
        unmasked = self.sc.unmasked
        filled = self.sc.filled(unmasked[1])
        expected = unmasked.copy()
        expected[self.mask] = unmasked[1]
        assert skycoord_equal(filled, expected)

    def test_filled_with_masked_value(self):
        # Filled ignores the mask (this will be true as long as __setitem__
        # ignores it; it may be a logical choice to actually use it).
        filled = self.sc.filled(self.sc[1])
        expected = self.sc.unmasked.copy()
        expected[self.mask] = self.sc.unmasked[1]
        assert skycoord_equal(filled, expected)

    @pytest.mark.parametrize(
        "dt",
        [
            1 * u.yr,
            Masked([1, 2, 3] * u.yr, mask=[False, True, False])[:, np.newaxis],
        ],
    )
    def test_apply_space_motion(self, dt):
        sc = self.sc.apply_space_motion(dt=dt)
        # All parts of the coordinate influence the final positions.
        expected_mask = self.sc.get_mask() | getattr(dt, "mask", False)
        assert_array_equal(sc.mask, expected_mask)
        expected_unmasked = self.sc.unmasked.apply_space_motion(dt=dt)
        assert skycoord_equal(sc.unmasked, expected_unmasked)


class TestSkyCoordWithOnlyDifferentialsMasked(TestSkyCoordWithDifferentials):
    @classmethod
    def setup_class(cls):
        super().setup_class()
        # Overwrite SkyCoord using unmasked distance.
        cls.mask_dis = False
        cls.sc = SkyCoord(
            ra=cls.ra,
            dec=cls.dec,
            distance=cls.dis,
            pm_ra_cosdec=cls.mpm_ra_cosdec,
            pm_dec=cls.mpm_dec,
            radial_velocity=cls.mrv,
        )
        cls.mask = cls.mask_dis | cls.mask_pm_ra_cosdec | cls.mask_pm_dec | cls.mask_rv
