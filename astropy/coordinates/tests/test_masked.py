# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_array_equal

import astropy.coordinates.representation as r
import astropy.units as u
from astropy.coordinates.tests.test_representation import representation_equal
from astropy.utils.masked import Masked


class TestSphericalRepresentationSeparateMasks:
    """Tests for mask propagation for Spherical with separate masks."""

    def setup_class(self):
        self.lon = np.array([0.0, 3.0, 6.0, 12.0, 15.0, 18.0]) << u.hourangle
        self.lat = np.array([-15.0, 30.0, 60.0, -60.0, 89.0, -80.0]) << u.deg
        self.dis = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0]) << u.pc
        self.mask_lon = np.array([False, True, False, False, True, True])
        self.mask_lat = np.array([False, False, True, False, True, True])
        self.mask_dis = np.array([False, False, False, True, False, True])
        self.mlon = Masked(self.lon, self.mask_lon)
        self.mlat = Masked(self.lat, self.mask_lat)
        self.mdis = Masked(self.dis, self.mask_dis)
        self.msph = r.SphericalRepresentation(self.mlon, self.mlat, self.mdis)
        self.mask_ang = self.mask_lon | self.mask_lat
        self.mask = self.mask_ang | self.mask_dis

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
        sph = self.msph.copy().unmasked
        sph[0] = self.msph[1]
        assert sph.masked
        assert_array_equal(
            sph.mask, np.concatenate(([True], np.zeros_like(self.mask[1:])))
        )
        assert_array_equal(
            sph.unmasked.lon, np.concatenate((sph.lon[1:2], self.lon[1:]))
        )
        sph[0] = self.msph[0].unmasked
        assert sph.masked  # Does not get reset
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
