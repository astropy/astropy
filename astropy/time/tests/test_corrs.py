# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, solar_system_ephemeris
from astropy.time import Time, TimeDelta
from astropy.utils import iers
from astropy.utils.compat.optional_deps import HAS_JPLEPHEM


class TestHelioBaryCentric:
    """
    Verify time offsets to the solar system barycentre and the heliocentre.
    Uses the WHT observing site.

    Tests are against values returned at time of initial creation of these
    routines.  They agree to an independent SLALIB based implementation
    to 20 microseconds.
    """

    @classmethod
    def setup_class(cls):
        cls.orig_auto_download = iers.conf.auto_download
        iers.conf.auto_download = False

    @classmethod
    def teardown_class(cls):
        iers.conf.auto_download = cls.orig_auto_download

    def setup_method(self):
        wht = EarthLocation(342.12 * u.deg, 28.758333333333333 * u.deg, 2327 * u.m)
        self.obstime = Time("2013-02-02T23:00", location=wht)
        self.obstime2 = Time("2013-08-02T23:00", location=wht)
        self.obstimeArr = Time(["2013-02-02T23:00", "2013-08-02T23:00"], location=wht)
        self.star = SkyCoord(
            "08:08:08 +32:00:00", unit=(u.hour, u.degree), frame="icrs"
        )

    def test_heliocentric(self):
        hval = self.obstime.light_travel_time(self.star, "heliocentric")
        assert isinstance(hval, TimeDelta)
        assert hval.scale == "tdb"
        assert abs(hval - 461.43037870502235 * u.s) < 1.0 * u.us

    def test_barycentric(self):
        bval = self.obstime.light_travel_time(self.star, "barycentric")
        assert isinstance(bval, TimeDelta)
        assert bval.scale == "tdb"
        assert abs(bval - 460.58538779827836 * u.s) < 1.0 * u.us

    def test_arrays(self):
        bval1 = self.obstime.light_travel_time(self.star, "barycentric")
        bval2 = self.obstime2.light_travel_time(self.star, "barycentric")
        bval_arr = self.obstimeArr.light_travel_time(self.star, "barycentric")
        hval1 = self.obstime.light_travel_time(self.star, "heliocentric")
        hval2 = self.obstime2.light_travel_time(self.star, "heliocentric")
        hval_arr = self.obstimeArr.light_travel_time(self.star, "heliocentric")
        assert hval_arr[0] - hval1 < 1.0 * u.us
        assert hval_arr[1] - hval2 < 1.0 * u.us
        assert bval_arr[0] - bval1 < 1.0 * u.us
        assert bval_arr[1] - bval2 < 1.0 * u.us

    @pytest.mark.remote_data
    @pytest.mark.skipif(not HAS_JPLEPHEM, reason="requires jplephem")
    def test_ephemerides(self):
        bval1 = self.obstime.light_travel_time(self.star, "barycentric")
        with solar_system_ephemeris.set("jpl"):
            bval2 = self.obstime.light_travel_time(
                self.star, "barycentric", ephemeris="jpl"
            )
        # should differ by less than 0.1 ms, but not be the same
        assert abs(bval1 - bval2) < 1.0 * u.ms
        assert abs(bval1 - bval2) > 1.0 * u.us
