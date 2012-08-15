import pytest
import numpy as np

import astropy.time as atm


class TestBasic():
    """Basic tests stemming from initial example and API reference"""

    def test_simple(self):
        times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
        t = atm.Time(times, format='iso', scale='utc')
        assert (repr(t) == "<Time object: scale='utc' format='iso' "
                "vals=['1999-01-01 00:00:00.123' '2010-01-01 00:00:00.000']>")
        assert np.allclose(t.jd1, np.array([2451179.5, 2455197.5]))
        assert np.allclose(t.jd2, np.array([1.42889802e-06, 0.00000000e+00]))

        # Set scale to TAI
        t = t.tai
        assert (repr(t) == "<Time object: scale='tai' format='iso' "
                "vals=['1999-01-01 00:00:32.123' '2010-01-01 00:00:34.000']>")
        assert np.allclose(t.jd1, np.array([2451179.5, 2455197.5]))
        assert np.allclose(t.jd2, np.array([0.0003718, 0.00039352]))

        # Get a new ``Time`` object which is referenced to the TT scale
        # (internal JD1 and JD1 are now with respect to TT scale)"""

        assert (repr(t.tt) == "<Time object: scale='tt' format='iso' "
                "vals=['1999-01-01 00:01:04.307' '2010-01-01 00:01:06.184']>")

        # Get the representation of the ``Time`` object in a particular format
        # (in this case seconds since 1998.0).  This returns either a scalar or
        # array, depending on whether the input was a scalar or array"""

        print t.cxcsec
        assert np.allclose(t.cxcsec, np.array([3.15360643e+07,
                                               3.78691266e+08]))

    def test_properties(self):
        """Use properties to convert scales and formats.  Note that the UT1 to
        UTC transformation requires a supplementary value (``delta_ut1_utc``)
        that can be obtained by interpolating from a table supplied by IERS.
        This will be included in the package later."""

        t = atm.Time('2010-01-01 00:00:00', format='iso', scale='utc')
        t.delta_ut1_utc = 0.3341  # Explicitly set one part of the xform
        assert np.allclose(t.jd, 2455197.5)
        assert t.iso == '2010-01-01 00:00:00.000'
        assert t.tt.iso == '2010-01-01 00:01:06.184'
        assert t.tai.iso == '2010-01-01 00:00:34.000'
        assert np.allclose(t.utc.jd, 2455197.5)
        assert np.allclose(t.ut1.jd, 2455197.500003867)
        assert t.tcg.isot == '2010-01-01T00:01:06.910'
        assert np.allclose(t.unix, 1262304000.0)
        assert np.allclose(t.cxcsec, 378691266.184)

    def test_precision(self):
        """Set the output precision which is used for some formats.  This is
        also a test of the code that provides a dict for global and instance
        options."""

        t = atm.Time('2010-01-01 00:00:00', format='iso', scale='utc')
        # Uses initial class-defined precision=3
        assert t.iso == '2010-01-01 00:00:00.000'

        # Set global precision = 5  XXX this uses private var, FIX THIS
        atm.Time._precision = 5
        assert t.iso == '2010-01-01 00:00:00.00000'

        # Set instance precision to 9
        t.precision = 9
        assert t.iso == '2010-01-01 00:00:00.000000000'
        assert t.tai.utc.iso == '2010-01-01 00:00:00.000000000'

        # Restore global to original default of 3, instance is still at 9
        atm.Time._precision = 3
        assert t.iso == '2010-01-01 00:00:00.000000000'

        # Make a new time instance and confirm precision = 3
        t = atm.Time('2010-01-01 00:00:00', format='iso', scale='utc')
        assert t.iso == '2010-01-01 00:00:00.000'

    def test_transforms(self):
        """Transform from UTC to all supported time scales (TAI, TCB, TCG,
        TDB, TT, UT1, UTC).  This requires auxilliary information (latitude and
        longitude)."""

        lat = 19.48125
        lon = -155.933222
        t = atm.Time('2006-01-15 21:24:37.5', format='iso', scale='utc',
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
        atm.Time(100.0, format='cxcsec')
        atm.Time(100.0, format='unix')
        atm.Time(1950.0, format='byear', scale='tai')
        atm.Time(2000.0, format='jyear', scale='tai')
        atm.Time('2000-01-01 12:23:34.0', format='iso', scale='tai')
        atm.Time('2000-01-01T12:23:34.0', format='isot', scale='tai')
        atm.Time(2400000.5, 51544.0333981, format='jd', scale='tai')
        atm.Time(0.0, 51544.0333981, format='mjd', scale='tai')
        atm.Time('2000:001:12:23:34.0', format='yday', scale='tai')

    def test_epoch_transform(self):
        """Besselian and julian epoch transforms"""
        jd = 2457073.05631
        t = atm.Time(jd, format='jd', scale='tai')
        assert np.allclose(t.byear, 2015.1365941021)
        assert np.allclose(t.jyear, 2015.1349933196)
        t2 = atm.Time(t.byear, format='byear', scale='tai')
        assert np.allclose(t2.jd, jd)
        t2 = atm.Time(t.jyear, format='jyear', scale='tai')
        assert np.allclose(t2.jd, jd)

    def test_input_validation(self):
        """Wrong input type raises error"""
        times = [10, 20]
        with pytest.raises(ValueError):
            atm.Time(times, format='iso', scale='utc')
        with pytest.raises(ValueError):
            atm.Time('2000:001', format='jd', scale='utc')
        with pytest.raises(ValueError):
            atm.Time([50000.0], ['bad'], format='mjd', scale='tai')
        with pytest.raises(ValueError):
            atm.Time(50000.0, 'bad', format='mjd', scale='tai')

class TestVal2():
    """Tests related to val2"""

    def test_val2_ignored(self):
        """Test that val2 is ignored for string input"""
        t = atm.Time('2001:001', 'ignored', scale='utc')
        assert t.yday == '2001:001:00:00:00.000'

    def test_val2(self):
        """Various tests of the val2 input"""
        t = atm.Time([0.0, 50000.0], [50000.0, 0.0], format='mjd', scale='tai')
        assert t.mjd[0] == t.mjd[1]
        assert t.jd[0] == t.jd[1]

    def test_val_matches_val2(self):
        with pytest.raises(ValueError):
            atm.Time([0.0, 50000.0], [0.0], format='mjd', scale='tai')
        with pytest.raises(ValueError):
            atm.Time([0.0], 0.0, format='mjd', scale='tai')

class TestSubFormat():
    """Test input and output subformat functionality"""

    def test_input_subformat(self):
        """Input subformat selection"""
        # Heterogeneous input formats with in_subfmt='*' (default)
        times = ['2000-01-01', '2000-01-01 01:01',
                 '2000-01-01 01:01:01', '2000-01-01 01:01:01.123']
        t = atm.Time(times, format='iso', scale='tai')
        assert np.all(t.iso == np.array(['2000-01-01 00:00:00.000',
                                         '2000-01-01 01:01:00.000',
                                         '2000-01-01 01:01:01.000',
                                         '2000-01-01 01:01:01.123']))

        # Heterogeneous input formats with in_subfmt='date_*'
        times = ['2000-01-01 01:01',
                 '2000-01-01 01:01:01', '2000-01-01 01:01:01.123']
        t = atm.Time(times, format='iso', scale='tai',
                     in_subfmt='date_*')
        assert np.all(t.iso == np.array(['2000-01-01 01:01:00.000',
                                         '2000-01-01 01:01:01.000',
                                         '2000-01-01 01:01:01.123']))

    def test_input_subformat_fail(self):
        """Failed format matching"""
        with pytest.raises(ValueError):
            atm.Time('2000-01-01 01:01', format='iso', scale='tai',
                     in_subfmt='date')

    def test_bad_input_subformat(self):
        """Non-existent input subformat"""
        with pytest.raises(ValueError):
            atm.Time('2000-01-01 01:01', format='iso', scale='tai',
                     in_subfmt='doesnt exist')

    def test_output_subformat(self):
        """Input subformat selection"""
        # Heterogeneous input formats with in_subfmt='*' (default)
        times = ['2000-01-01', '2000-01-01 01:01',
                 '2000-01-01 01:01:01', '2000-01-01 01:01:01.123']
        t = atm.Time(times, format='iso', scale='tai',
                     out_subfmt='date_hm')
        assert np.all(t.iso == np.array(['2000-01-01 00:00',
                                         '2000-01-01 01:01',
                                         '2000-01-01 01:01',
                                         '2000-01-01 01:01']))

    def test_yday_format(self):
        """Year:Day_of_year format"""
        # Heterogeneous input formats with in_subfmt='*' (default)
        times = ['2000-12-01', '2001-12-01 01:01:01.123']
        t = atm.Time(times, format='iso', scale='tai')
        t.out_subfmt = 'date_hm'
        assert np.all(t.yday == np.array(['2000:336:00:00',
                                          '2001:335:01:01']))
        t.out_subfmt = '*'
        assert np.all(t.yday == np.array(['2000:336:00:00:00.000',
                                          '2001:335:01:01:01.123']))

    def test_scale_input(self):
        """Test for issues related to scale input"""
        # Check case where required scale is defined by the TimeFormat.
        # The first two should succeed, the third should fail
        atm.Time(100.0, format='cxcsec')
        atm.Time(100.0, format='cxcsec', scale='tai')
        with pytest.raises(atm.ScaleValueError):
            atm.Time(100.0, format='cxcsec', scale='utc')

        # Check that bad scale is caught when format is specified
        with pytest.raises(atm.ScaleValueError):
            atm.Time(1950.0, format='byear')
        with pytest.raises(atm.ScaleValueError):
            atm.Time(1950.0, format='byear', scale='bad scale')

        # Check that bad scale is caught when format is auto-determined
        with pytest.raises(atm.ScaleValueError):
            atm.Time('2000:001:00:00:00')
        with pytest.raises(atm.ScaleValueError):
            atm.Time('2000:001:00:00:00', scale='bad scale')
