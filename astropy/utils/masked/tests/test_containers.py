# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord, representation as r

from .. import Masked


class TestRepresentations:
    def setup_class(self):
        self.x = np.array([3., 5., 0.]) << u.m
        self.y = np.array([4., 12., 1.]) << u.m
        self.z = np.array([0., 0., 1.]) << u.m
        self.c = r.CartesianRepresentation(self.x, self.y, self.z)
        self.mask = np.array([False, False, True])
        self.mx = Masked(self.x, self.mask)
        self.my = Masked(self.y, self.mask)
        self.mz = Masked(self.z, self.mask)
        self.mc = r.CartesianRepresentation(self.mx, self.my, self.mz)

    def test_initialization(self):
        check = self.mc.z == self.mz
        assert_array_equal(check.unmasked, np.ones(3, bool))
        assert_array_equal(check.mask, self.mask)
        assert_array_equal(self.mc.x, self.mx)
        assert_array_equal(self.mc.y, self.my)
        assert_array_equal(self.mc.z, self.mz)

    def test_norm(self):
        # Need stacking and erfa override.
        norm = self.mc.norm()
        assert_array_equal(norm.unmasked, self.c.norm())
        assert_array_equal(norm.mask, self.mask)

    def test_transformation(self):
        msr = self.mc.represent_as(r.SphericalRepresentation)
        sr = self.c.represent_as(r.SphericalRepresentation)
        for comp in msr.components:
            mc = getattr(msr, comp)
            c = getattr(sr, comp)
            assert_array_equal(mc.unmasked, c)
            assert_array_equal(mc.mask, self.mask)

        # Transformation back.  This also tests erfa.ufunc.s2p, which
        # is special in having a core dimension only in the output.
        cr2 = sr.represent_as(r.CartesianRepresentation)
        mcr2 = msr.represent_as(r.CartesianRepresentation)
        for comp in mcr2.components:
            mc = getattr(mcr2, comp)
            c = getattr(cr2, comp)
            assert_array_equal(mc.unmasked, c)
            assert_array_equal(mc.mask, self.mask)


class TestSkyCoord:
    def setup_class(self):
        self.ra = np.array([3., 5., 0.]) << u.hourangle
        self.dec = np.array([4., 12., 1.]) << u.deg
        self.sc = SkyCoord(self.ra, self.dec)
        self.mask = np.array([False, False, True])
        self.mra = Masked(self.ra, self.mask)
        self.mdec = Masked(self.dec, self.mask)
        self.msc = SkyCoord(self.mra, self.mdec)

    def test_initialization(self):
        check = self.msc.dec == self.mdec
        assert_array_equal(check.unmasked, np.ones(3, bool))
        assert_array_equal(check.mask, self.mask)
        assert_array_equal(self.msc.data.lon, self.mra)
        assert_array_equal(self.msc.data.lat, self.mdec)

    def test_transformation(self):
        gcrs = self.sc.gcrs
        mgcrs = self.msc.gcrs
        assert_array_equal(mgcrs.data.lon.mask, self.msc.data.lon.mask)
        assert_array_equal(mgcrs.data.lon.unmasked, gcrs.data.lon)
        assert_array_equal(mgcrs.data.lat.unmasked, gcrs.data.lat)
