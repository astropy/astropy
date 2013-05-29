# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np

from ...tests.helper import pytest
from .. import Time, ScaleValueError, sofa_time

allclose_jd = functools.partial(np.allclose, rtol=1e-15, atol=0)
allclose_jd2 = functools.partial(np.allclose, rtol=1e-15, atol=1e-11)  # 1 microsec atol
allclose_sec = functools.partial(np.allclose, rtol=1e-15, atol=1e-9)  # 1 nanosec atol
allclose_year = functools.partial(np.allclose, rtol=1e-15, atol=0)  # 60 microsec at current epoch


class TestBasic():
    """Basic tests stemming from initial example and API reference"""

    def test_simple(self):
        times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
        t = Time(times, format='iso', scale='utc')
        assert (repr(t) == "<Time object: scale='utc' format='iso' "
                "vals=['1999-01-01 00:00:00.123' '2010-01-01 00:00:00.000']>")
        assert allclose_jd(t.jd1, np.array([2451179.5, 2455197.5]))
        assert allclose_jd2(t.jd2, np.array([1.42889802e-06, 0.00000000e+00]))

        # Set scale to TAI
        t = t.tai
        assert (repr(t) == "<Time object: scale='tai' format='iso' "
                "vals=['1999-01-01 00:00:32.123' '2010-01-01 00:00:34.000']>")
        assert allclose_jd(t.jd1, np.array([2451179.5, 2455197.5]))
        assert allclose_jd2(t.jd2, np.array([0.00037179926839122024, 0.00039351851851851852]))

        # Get a new ``Time`` object which is referenced to the TT scale
        # (internal JD1 and JD1 are now with respect to TT scale)"""

        assert (repr(t.tt) == "<Time object: scale='tt' format='iso' "
                "vals=['1999-01-01 00:01:04.307' '2010-01-01 00:01:06.184']>")

        # Get the representation of the ``Time`` object in a particular format
        # (in this case seconds since 1998.0).  This returns either a scalar or
        # array, depending on whether the input was a scalar or array"""

        print t.cxcsec
        assert allclose_sec(t.cxcsec, np.array([31536064.307456788, 378691266.18400002]))

    def test_copy_time(self):
        """Test copying the values of a Time object by passing it into the
        Time initializer.
        """
        t = Time(2455197.5, format='jd', scale='utc')

        t2 = Time(t, copy=False)
        assert t.jd - t2.jd == 0
        assert (t - t2).jd == 0
        assert t._time.jd1 is t2._time.jd1
        assert t._time.jd2 is t2._time.jd2

        t2 = Time(t, copy=True)
        assert t.jd - t2.jd == 0
        assert (t - t2).jd == 0
        assert t._time.jd1 is not t2._time.jd1
        assert t._time.jd2 is not t2._time.jd2

        # Include initializers
        t2 = Time(t, format='iso', scale='tai', precision=1)
        assert t2.val == '2010-01-01 00:00:34.0'
        t2 = Time(t, format='iso', scale='tai', out_subfmt='date')
        assert t2.val == '2010-01-01'

    def test_properties(self):
        """Use properties to convert scales and formats.  Note that the UT1 to
        UTC transformation requires a supplementary value (``delta_ut1_utc``)
        that can be obtained by interpolating from a table supplied by IERS.
        This will be included in the package later."""

        t = Time('2010-01-01 00:00:00', format='iso', scale='utc')
        t.delta_ut1_utc = 0.3341  # Explicitly set one part of the xform
        assert allclose_jd(t.jd, 2455197.5)
        assert t.iso == '2010-01-01 00:00:00.000'
        assert t.tt.iso == '2010-01-01 00:01:06.184'
        assert t.tai.iso == '2010-01-01 00:00:34.000'
        assert allclose_jd(t.utc.jd, 2455197.5)
        assert allclose_jd(t.ut1.jd, 2455197.500003867)
        assert t.tcg.isot == '2010-01-01T00:01:06.910'
        assert allclose_sec(t.unix, 1262304000.0)
        assert allclose_sec(t.cxcsec, 378691266.184)

    def test_precision(self):
        """Set the output precision which is used for some formats.  This is
        also a test of the code that provides a dict for global and instance
        options."""

        t = Time('2010-01-01 00:00:00', format='iso', scale='utc')
        # Uses initial class-defined precision=3
        assert t.iso == '2010-01-01 00:00:00.000'

        # Set global precision = 5  XXX this uses private var, FIX THIS
        Time._precision = 5
        assert t.iso == '2010-01-01 00:00:00.00000'

        # Set instance precision to 9
        t.precision = 9
        assert t.iso == '2010-01-01 00:00:00.000000000'
        assert t.tai.utc.iso == '2010-01-01 00:00:00.000000000'

        # Restore global to original default of 3, instance is still at 9
        Time._precision = 3
        assert t.iso == '2010-01-01 00:00:00.000000000'

        # Make a new time instance and confirm precision = 3
        t = Time('2010-01-01 00:00:00', format='iso', scale='utc')
        assert t.iso == '2010-01-01 00:00:00.000'

    def test_transforms(self):
        """Transform from UTC to all supported time scales (TAI, TCB, TCG,
        TDB, TT, UT1, UTC).  This requires auxilliary information (latitude and
        longitude)."""

        lat = 19.48125
        lon = -155.933222
        t = Time('2006-01-15 21:24:37.5', format='iso', scale='utc',
                 precision=6, lat=lat, lon=lon)
        t.delta_ut1_utc = 0.3341  # Explicitly set one part of the xform
        assert t.utc.iso == '2006-01-15 21:24:37.500000'
        assert t.ut1.iso == '2006-01-15 21:24:37.834100'
        assert t.tai.iso == '2006-01-15 21:25:10.500000'
        assert t.tt.iso == '2006-01-15 21:25:42.684000'
        assert t.tcg.iso == '2006-01-15 21:25:43.322690'
        assert t.tdb.iso == '2006-01-15 21:25:42.683799'
        assert t.tcb.iso == '2006-01-15 21:25:56.893378'

    def test_creating_all_formats(self):
        """Create a time object using each defined format"""
        Time(100.0, format='cxcsec')
        Time(100.0, format='unix')
        Time(1950.0, format='byear', scale='tai')
        Time(2000.0, format='jyear', scale='tai')
        Time('B1950.0', format='byear_str', scale='tai')
        Time('J2000.0', format='jyear_str', scale='tai')
        Time('2000-01-01 12:23:34.0', format='iso', scale='tai')
        Time('2000-01-01T12:23:34.0', format='isot', scale='tai')
        Time(2400000.5, 51544.0333981, format='jd', scale='tai')
        Time(0.0, 51544.0333981, format='mjd', scale='tai')
        Time('2000:001:12:23:34.0', format='yday', scale='tai')

    def test_epoch_transform(self):
        """Besselian and julian epoch transforms"""
        jd = 2457073.05631
        t = Time(jd, format='jd', scale='tai', precision=6)
        assert allclose_year(t.byear, 2015.1365941020817)
        assert allclose_year(t.jyear, 2015.1349933196439)
        assert t.byear_str == 'B2015.136594'
        assert t.jyear_str == 'J2015.134993'
        t2 = Time(t.byear, format='byear', scale='tai')
        assert allclose_jd(t2.jd, jd)
        t2 = Time(t.jyear, format='jyear', scale='tai')
        assert allclose_jd(t2.jd, jd)

        t = Time('J2015.134993', scale='tai', precision=6)
        assert np.allclose(t.jd, jd, rtol=1e-10, atol=0)  # J2015.134993 has 10 digit precision
        assert t.byear_str == 'B2015.136594'

    def test_input_validation(self):
        """Wrong input type raises error"""
        times = [10, 20]
        with pytest.raises(ValueError):
            Time(times, format='iso', scale='utc')
        with pytest.raises(ValueError):
            Time('2000:001', format='jd', scale='utc')
        with pytest.raises(ValueError):
            Time([50000.0], ['bad'], format='mjd', scale='tai')
        with pytest.raises(ValueError):
            Time(50000.0, 'bad', format='mjd', scale='tai')

    def test_utc_leap_sec(self):
        """Time behaves properly near or in UTC leap second.  This
        uses the 2012-06-30 leap second for testing."""

        # Start with a day without a leap second and note rollover
        t1 = Time('2012-06-01 23:59:60.0', scale='utc')
        assert t1.iso == '2012-06-02 00:00:00.000'

        # Leap second is different
        t1 = Time('2012-06-30 23:59:59.900', scale='utc')
        assert t1.iso == '2012-06-30 23:59:59.900'

        t1 = Time('2012-06-30 23:59:60.000', scale='utc')
        assert t1.iso == '2012-06-30 23:59:60.000'

        t1 = Time('2012-06-30 23:59:60.999', scale='utc')
        assert t1.iso == '2012-06-30 23:59:60.999'

        t1 = Time('2012-06-30 23:59:61.0', scale='utc')
        assert t1.iso == '2012-07-01 00:00:00.000'

        # Delta time gives 2 seconds here as expected
        t0 = Time('2012-06-30 23:59:59', scale='utc')
        t1 = Time('2012-07-01 00:00:00', scale='utc')
        assert allclose_sec((t1 - t0).sec, 2.0)


class TestVal2():
    """Tests related to val2"""

    def test_val2_ignored(self):
        """Test that val2 is ignored for string input"""
        t = Time('2001:001', 'ignored', scale='utc')
        assert t.yday == '2001:001:00:00:00.000'

    def test_val2(self):
        """Various tests of the val2 input"""
        t = Time([0.0, 50000.0], [50000.0, 0.0], format='mjd', scale='tai')
        assert t.mjd[0] == t.mjd[1]
        assert t.jd[0] == t.jd[1]

    def test_val_matches_val2(self):
        with pytest.raises(ValueError):
            Time([0.0, 50000.0], [0.0], format='mjd', scale='tai')
        with pytest.raises(ValueError):
            Time([0.0], 0.0, format='mjd', scale='tai')


class TestSubFormat():
    """Test input and output subformat functionality"""

    def test_input_subformat(self):
        """Input subformat selection"""
        # Heterogeneous input formats with in_subfmt='*' (default)
        times = ['2000-01-01', '2000-01-01 01:01',
                 '2000-01-01 01:01:01', '2000-01-01 01:01:01.123']
        t = Time(times, format='iso', scale='tai')
        assert np.all(t.iso == np.array(['2000-01-01 00:00:00.000',
                                         '2000-01-01 01:01:00.000',
                                         '2000-01-01 01:01:01.000',
                                         '2000-01-01 01:01:01.123']))

        # Heterogeneous input formats with in_subfmt='date_*'
        times = ['2000-01-01 01:01',
                 '2000-01-01 01:01:01', '2000-01-01 01:01:01.123']
        t = Time(times, format='iso', scale='tai',
                     in_subfmt='date_*')
        assert np.all(t.iso == np.array(['2000-01-01 01:01:00.000',
                                         '2000-01-01 01:01:01.000',
                                         '2000-01-01 01:01:01.123']))

    def test_input_subformat_fail(self):
        """Failed format matching"""
        with pytest.raises(ValueError):
            Time('2000-01-01 01:01', format='iso', scale='tai',
                     in_subfmt='date')

    def test_bad_input_subformat(self):
        """Non-existent input subformat"""
        with pytest.raises(ValueError):
            Time('2000-01-01 01:01', format='iso', scale='tai',
                     in_subfmt='doesnt exist')

    def test_output_subformat(self):
        """Input subformat selection"""
        # Heterogeneous input formats with in_subfmt='*' (default)
        times = ['2000-01-01', '2000-01-01 01:01',
                 '2000-01-01 01:01:01', '2000-01-01 01:01:01.123']
        t = Time(times, format='iso', scale='tai',
                     out_subfmt='date_hm')
        assert np.all(t.iso == np.array(['2000-01-01 00:00',
                                         '2000-01-01 01:01',
                                         '2000-01-01 01:01',
                                         '2000-01-01 01:01']))

    def test_yday_format(self):
        """Year:Day_of_year format"""
        # Heterogeneous input formats with in_subfmt='*' (default)
        times = ['2000-12-01', '2001-12-01 01:01:01.123']
        t = Time(times, format='iso', scale='tai')
        t.out_subfmt = 'date_hm'
        assert np.all(t.yday == np.array(['2000:336:00:00',
                                          '2001:335:01:01']))
        t.out_subfmt = '*'
        assert np.all(t.yday == np.array(['2000:336:00:00:00.000',
                                          '2001:335:01:01:01.123']))

    def test_scale_input(self):
        """Test for issues related to scale input"""
        # Check case where required scale is defined by the TimeFormat.
        # All three should work.
        t = Time(100.0, format='cxcsec', scale='utc')
        assert t.scale == 'utc'
        t = Time(100.0, format='unix', scale='tai')
        assert t.scale == 'tai'

        # Check that bad scale is caught when format is specified
        with pytest.raises(ScaleValueError):
            Time(1950.0, format='byear')
        with pytest.raises(ScaleValueError):
            Time(1950.0, format='byear', scale='bad scale')

        # Check that bad scale is caught when format is auto-determined
        with pytest.raises(ScaleValueError):
            Time('2000:001:00:00:00')
        with pytest.raises(ScaleValueError):
            Time('2000:001:00:00:00', scale='bad scale')

    def test_epoch_times(self):
        """Test time formats derived from EpochFromTime"""
        t = Time(0.0, format='cxcsec', scale='tai')
        assert t.tt.iso == '1998-01-01 00:00:00.000'

        # Create new time object from this one and change scale, format
        t2 = Time(t, scale='tt', format='iso')
        assert t2.val == '1998-01-01 00:00:00.000'

        # Value take from Chandra.Time.DateTime('2010:001:00:00:00').secs
        t = Time(378691266.184, format='cxcsec', scale='utc')
        assert t.yday == '2010:001:00:00:00.000'
        t = Time('2010:001:00:00:00.000', scale='utc')
        assert allclose_sec(t.cxcsec, 378691266.184)
        assert allclose_sec(t.tt.cxcsec, 378691266.184)

        # Round trip through epoch time
        for scale in ('utc', 'tt'):
            t = Time('2000:001', scale=scale)
            t2 = Time(t.unix, scale=scale, format='unix')
            assert getattr(t2, scale).iso == '2000-01-01 00:00:00.000'

        # Test unix time.  Values taken from http://en.wikipedia.org/wiki/Unix_time
        t = Time('2013-05-20 21:18:46', scale='utc')
        assert allclose_sec(t.unix, 1369084726.0)
        assert allclose_sec(t.tt.unix, 1369084726.0)

        # Values from issue #1118
        t = Time('2004-09-16T23:59:59', scale='utc')
        assert allclose_sec(t.unix, 1095379199.0)


class TestSofaErrors():
    """Test that sofa_time.pyx handles sofa status return values correctly"""

    def test_bad_time(self):
        iy = np.array([2000], dtype=np.intc)
        im = np.array([2000], dtype=np.intc) # bad month
        id = np.array([2000], dtype=np.intc) # bad day
        djm0 = np.array([0], dtype=np.double)
        djm = np.array([0], dtype=np.double)
        with pytest.raises(ValueError):  # bad month, fatal error
            sofa_time.cal2jd(iy, im, id, djm0, djm)

        # Set month to a good value so now the bad day just gives a warning
        im[0] = 2
        sofa_time.cal2jd(iy, im, id, djm0, djm)
        assert allclose_jd(djm0, [2400000.5])
        assert allclose_jd(djm, [53574.])

        # How do you test for warnings in pytest?  Test that dubious year for
        # UTC works.


class TestCopyReplicate():
    """Test issues related to copying and replicating data"""

    def test_mutable_input(self):
        """By default if JDs are provided with copy=False then internals are
        mutable."""
        jds = np.array([2450000.5], dtype=np.double)
        t = Time(jds, format='jd', scale='tai')
        assert allclose_jd(t.jd, jds)
        jds[0] = 2459009.5
        assert allclose_jd(t.jd, jds)

        t = Time(jds, format='jd', scale='tai', copy=True)
        assert allclose_jd(t.jd, jds)
        jds[0] = 2458654
        assert not allclose_jd(t.jd, jds)

        # MJD does not suffer from this mutability
        mjds = np.array([50000.0], dtype=np.double)
        t = Time(mjds, format='mjd', scale='tai')
        assert allclose_jd(t.jd, [2450000.5])
        mjds[0] = 0.0
        assert allclose_jd(t.jd, [2450000.5])

    def test_replicate(self):
        """Test replicate method"""
        t = Time('2000:001', format='yday', scale='tai')
        t_yday = t.yday
        t2 = t.replicate()
        assert t.yday == t2.yday
        assert t.format == t2.format
        assert t.scale == t2.scale
        # This is not allowed publicly, but here we hack the internal time values
        # to show that t and t2 are sharing references.
        t2._time.jd1 += 100.0
        assert t.yday == t2.yday
        assert t.yday != t_yday  # prove that it changed

    def test_copy(self):
        """Test copy method"""
        t = Time('2000:001', format='yday', scale='tai')
        t_yday = t.yday
        t2 = t.copy()
        assert t.yday == t2.yday
        # This is not allowed publicly, but here we hack the internal time values
        # to show that t and t2 are not sharing references.
        t2._time.jd1 += 100.0
        assert t.yday != t2.yday
        assert t.yday == t_yday  # prove that it did not change

def test_python_builtin_copy():
    import copy

    t = Time('2000:001', format='yday', scale='tai')
    t2 = copy.copy(t)
    t3 = copy.deepcopy(t)

    assert t.jd == t2.jd
    assert t.jd == t3.jd
