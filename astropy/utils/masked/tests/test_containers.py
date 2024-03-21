# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import representation as r
from astropy.time import Time
from astropy.utils.masked import Masked


class TestRepresentations:
    def setup_class(self):
        self.x = np.array([3.0, 5.0, 0.0]) << u.m
        self.y = np.array([4.0, 12.0, 1.0]) << u.m
        self.z = np.array([0.0, 0.0, 1.0]) << u.m
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
        self.ra = np.array([3.0, 5.0, 0.0]) << u.hourangle
        self.dec = np.array([4.0, 12.0, 1.0]) << u.deg
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

    @pytest.mark.filterwarnings("ignore:.*ERFA.*distance overridden.*")
    def test_apply_space_motion(self):
        # Regression test for gh-13041. It is important that the
        # distance is missing here, so that a warning is raised.
        # Before gh-16224, the ERFA warning machinery let to an error.
        kwargs = dict(
            ra=[1, 2] * u.deg,
            dec=[3, 4] * u.deg,
            pm_ra_cosdec=[5, 6] * u.mas / u.yr,
            pm_dec=[1, 2] * u.mas / u.yr,
            obstime=Time(["2000-10-22T12:23:45", "2000-10-22T12:23:45"]),
        )
        sc = SkyCoord(**kwargs)
        mask = np.array([False, True])
        kwargs["pm_ra_cosdec"] = Masked(kwargs["pm_ra_cosdec"], mask=mask)
        kwargs["pm_dec"] = Masked(kwargs["pm_dec"], mask=mask)
        msc = SkyCoord(**kwargs)
        new_time = Time("2020-10-22T12:23:45")
        sc2 = sc.apply_space_motion(new_obstime=new_time)
        msc2 = msc.apply_space_motion(new_obstime=new_time)
        assert_array_equal(msc2.data.lon.mask, [False, True])
        assert_array_equal(msc2.data.lon.unmasked, sc2.data.lon)
        assert_array_equal(msc2.data.lat.unmasked, sc2.data.lat)


class TestTime:
    def setup_class(self):
        self.s = np.array(
            [
                "2010-11-12T13:14:15.160",
                "2010-11-12T13:14:15.161",
                "2011-12-13T14:15:16.170",
            ]
        )
        self.t = Time(self.s)
        # Time formats will currently strip any ndarray subtypes, so we cannot
        # initialize a Time with a Masked version of self.s yet. Instead, we
        # work around it, for now only testing that masked are preserved by
        # transformations.
        self.mask = np.array([False, False, True])
        self.mt = self.t._apply(Masked, self.mask)

    def test_initialization(self):
        assert_array_equal(self.mt.jd1.mask, self.mask)
        assert_array_equal(self.mt.jd2.mask, self.mask)
        assert_array_equal(self.mt.jd1.unmasked, self.t.jd1)
        assert_array_equal(self.mt.jd2.unmasked, self.t.jd2)

    @pytest.mark.parametrize("format_", ["jd", "cxcsec", "jyear"])
    def test_different_formats(self, format_):
        # Formats do not yet work with everything; e.g., isot is not supported
        # since the Masked class does not yet support structured arrays.
        tfmt = getattr(self.t, format_)
        mtfmt = getattr(self.mt, format_)
        check = mtfmt == tfmt
        assert_array_equal(check.unmasked, np.ones(3, bool))
        assert_array_equal(check.mask, self.mask)

    @pytest.mark.parametrize("scale", ["tai", "tcb", "ut1"])
    def test_transformation(self, scale):
        tscl = getattr(self.t, scale)
        mtscl = getattr(self.mt, scale)
        assert_array_equal(mtscl.jd1.mask, self.mask)
        assert_array_equal(mtscl.jd2.mask, self.mask)
        assert_array_equal(mtscl.jd1.unmasked, tscl.jd1)
        assert_array_equal(mtscl.jd2.unmasked, tscl.jd2)
