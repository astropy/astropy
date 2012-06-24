import pytest
import numpy as np

from astropy.time import Time

class TestBasic():
    """Basic tests stemming from initial example and API reference"""

    def test_simple(self):
        times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
        t = Time(times, format='iso', system='utc')
        assert (repr(t) == "<Time object: system='utc' format='iso' "
                "vals=['1999-01-01 00:00:00.123' '2010-01-01 00:00:00.000']>")
        assert np.allclose(t.jd1, np.array([2451179.5, 2455197.5]))
        assert np.allclose(t.jd2, np.array([1.42889802e-06, 0.00000000e+00]))

        # Set system to TAI
        t.set_system('tai')
        assert (repr(t) == "<Time object: system='tai' format='iso' "
                "vals=['1999-01-01 00:00:32.123' '2010-01-01 00:00:34.000']>")
        assert np.allclose(t.jd1, np.array([2451179.5, 2455197.5]))
        assert np.allclose(t.jd2, np.array([0.0003718, 0.00039352]))

        # Get a new ``Time`` object which is referenced to the TT system
        # (internal JD1 and JD1 are now with respect to TT system)"""

        assert (repr(t.tt) == "<Time object: system='tt' format='iso' "
                "vals=['1999-01-01 00:01:04.307' '2010-01-01 00:01:06.184']>")

        # Get the representation of the ``Time`` object in a particular format
        # (in this case seconds since 1998.0).  This returns either a scalar or
        # array, depending on whether the input was a scalar or array"""

        print t.cxcsec
        assert np.allclose(t.cxcsec, np.array([3.15360643e+07,
                                               3.78691266e+08]))

    def test_properties(self):
        """Use properties to convert systems and formats.  Note that the UT1 to
        UTC transformation requires a supplementary value (``delta_ut1_utc``)
        that can be obtained by interpolating from a table supplied by IERS.
        This will be included in the package later."""

        t = Time('2010-01-01 00:00:00', format='iso', system='utc')
        t.set_delta_ut1_utc(0.3341)  # Explicitly set one part of the xform
        assert np.allclose(t.jd, 2455197.5)
        assert t.iso == '2010-01-01 00:00:00.000'
        assert t.tt.iso == '2010-01-01 00:01:06.184'
        assert t.tai.iso == '2010-01-01 00:00:34.000'
        assert np.allclose(t.utc.jd, 2455197.5)
        assert np.allclose(t.ut1.jd, 2455197.500003867)
        assert t.tcg.isot == '2010-01-01T00:00:00.000'
        assert np.allclose(t.unix, 1262304000.0)
        assert np.allclose(t.cxcsec, 378691266.184)

    def test_precision(self):
        """Set the output precision which is used for some formats"""
        t = Time('2010-01-01 00:00:00', format='iso', system='utc')
        t.precision = 9
        assert t.iso == '2010-01-01 00:00:00.000000000'

    def test_transforms(self):
        """Transform from UTC to all supported time systems (TAI, TCB, TCG,
        TDB, TT, UT1, UTC).  This requires auxilliary information (latitude and
        longitude)."""

        lat = 19.48125
        lon = -155.933222
        t = Time('2006-01-15 21:24:37.5', format='iso', system='utc',
                 precision=6, lat=lat, lon=lon)
        t.set_delta_ut1_utc(0.3341)  # Explicitly set one part of the xform
        assert t.utc.iso == '2006-01-15 21:24:37.500000'
        assert t.ut1.iso == '2006-01-15 21:24:37.834100'
        assert t.tai.iso == '2006-01-15 21:25:10.500000'
        assert t.tt.iso == '2006-01-15 21:25:42.684000'
        assert t.tcg.iso == '2006-01-15 21:25:43.322690'
        assert t.tdb.iso == '2006-01-15 21:25:42.683799'
        assert t.tcb.iso == '2006-01-15 21:25:56.893378'
