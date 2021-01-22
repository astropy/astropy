# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import copy
import functools
import datetime
from copy import deepcopy
from decimal import Decimal, localcontext

import numpy as np
import pytest
from numpy.testing import assert_allclose
import erfa
from erfa import ErfaWarning

from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils import isiterable, iers
from astropy.time import (Time, TimeDelta, ScaleValueError, STANDARD_TIME_SCALES,
                          TimeString, TimezoneInfo, TIME_FORMATS)
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.table import Column, Table

try:
    import pytz
    HAS_PYTZ = True
except ImportError:
    HAS_PYTZ = False


allclose_jd = functools.partial(np.allclose, rtol=np.finfo(float).eps, atol=0)
allclose_jd2 = functools.partial(np.allclose, rtol=np.finfo(float).eps,
                                 atol=np.finfo(float).eps)  # 20 ps atol
allclose_sec = functools.partial(np.allclose, rtol=np.finfo(float).eps,
                                 atol=np.finfo(float).eps * 24 * 3600)
allclose_year = functools.partial(np.allclose, rtol=np.finfo(float).eps,
                                  atol=0.)  # 14 microsec at current epoch


def setup_function(func):
    func.FORMATS_ORIG = deepcopy(Time.FORMATS)


def teardown_function(func):
    Time.FORMATS.clear()
    Time.FORMATS.update(func.FORMATS_ORIG)


class TestBasic:
    """Basic tests stemming from initial example and API reference"""

    def test_simple(self):
        times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
        t = Time(times, format='iso', scale='utc')
        assert (repr(t) == "<Time object: scale='utc' format='iso' "
                "value=['1999-01-01 00:00:00.123' '2010-01-01 00:00:00.000']>")
        assert allclose_jd(t.jd1, np.array([2451180., 2455198.]))
        assert allclose_jd2(t.jd2, np.array([-0.5 + 1.4288980208333335e-06,
                                             -0.50000000e+00]))

        # Set scale to TAI
        t = t.tai
        assert (repr(t) == "<Time object: scale='tai' format='iso' "
                "value=['1999-01-01 00:00:32.123' '2010-01-01 00:00:34.000']>")
        assert allclose_jd(t.jd1, np.array([2451180., 2455198.]))
        assert allclose_jd2(t.jd2, np.array([-0.5 + 0.00037179926839122024,
                                             -0.5 + 0.00039351851851851852]))

        # Get a new ``Time`` object which is referenced to the TT scale
        # (internal JD1 and JD1 are now with respect to TT scale)"""

        assert (repr(t.tt) == "<Time object: scale='tt' format='iso' "
                "value=['1999-01-01 00:01:04.307' '2010-01-01 00:01:06.184']>")

        # Get the representation of the ``Time`` object in a particular format
        # (in this case seconds since 1998.0).  This returns either a scalar or
        # array, depending on whether the input was a scalar or array"""

        assert allclose_sec(t.cxcsec, np.array([31536064.307456788, 378691266.18400002]))

    def test_different_dimensions(self):
        """Test scalars, vector, and higher-dimensions"""
        # scalar
        val, val1 = 2450000.0, 0.125
        t1 = Time(val, val1, format='jd')
        assert t1.isscalar is True and t1.shape == ()
        # vector
        val = np.arange(2450000., 2450010.)
        t2 = Time(val, format='jd')
        assert t2.isscalar is False and t2.shape == val.shape
        # explicitly check broadcasting for mixed vector, scalar.
        val2 = 0.
        t3 = Time(val, val2, format='jd')
        assert t3.isscalar is False and t3.shape == val.shape
        val2 = (np.arange(5.) / 10.).reshape(5, 1)
        # now see if broadcasting to two-dimensional works
        t4 = Time(val, val2, format='jd')
        assert t4.isscalar is False
        assert t4.shape == np.broadcast(val, val2).shape

    @pytest.mark.parametrize('format_', Time.FORMATS)
    def test_empty_value(self, format_):
        t = Time([], format=format_)
        assert t.size == 0
        assert t.shape == (0,)
        assert t.format == format_
        t_value = t.value
        assert t_value.size == 0
        assert t_value.shape == (0,)
        t2 = Time(t_value, format=format_)
        assert t2.size == 0
        assert t2.shape == (0,)
        assert t2.format == format_
        t3 = t2.tai
        assert t3.size == 0
        assert t3.shape == (0,)
        assert t3.format == format_
        assert t3.scale == 'tai'

    @pytest.mark.parametrize('value', [2455197.5, [2455197.5]])
    def test_copy_time(self, value):
        """Test copying the values of a Time object by passing it into the
        Time initializer.
        """
        t = Time(value, format='jd', scale='utc')

        t2 = Time(t, copy=False)
        assert np.all(t.jd - t2.jd == 0)
        assert np.all((t - t2).jd == 0)
        assert t._time.jd1 is t2._time.jd1
        assert t._time.jd2 is t2._time.jd2

        t2 = Time(t, copy=True)
        assert np.all(t.jd - t2.jd == 0)
        assert np.all((t - t2).jd == 0)
        assert t._time.jd1 is not t2._time.jd1
        assert t._time.jd2 is not t2._time.jd2

        # Include initializers
        t2 = Time(t, format='iso', scale='tai', precision=1)
        assert t2.value == '2010-01-01 00:00:34.0'
        t2 = Time(t, format='iso', scale='tai', out_subfmt='date')
        assert t2.value == '2010-01-01'

    def test_getitem(self):
        """Test that Time objects holding arrays are properly subscriptable,
        set isscalar as appropriate, and also subscript delta_ut1_utc, etc."""

        mjd = np.arange(50000, 50010)
        t = Time(mjd, format='mjd', scale='utc', location=('45d', '50d'))
        t1 = t[3]
        assert t1.isscalar is True
        assert t1._time.jd1 == t._time.jd1[3]
        assert t1.location is t.location
        t1a = Time(mjd[3], format='mjd', scale='utc')
        assert t1a.isscalar is True
        assert np.all(t1._time.jd1 == t1a._time.jd1)
        t1b = Time(t[3])
        assert t1b.isscalar is True
        assert np.all(t1._time.jd1 == t1b._time.jd1)
        t2 = t[4:6]
        assert t2.isscalar is False
        assert np.all(t2._time.jd1 == t._time.jd1[4:6])
        assert t2.location is t.location
        t2a = Time(t[4:6])
        assert t2a.isscalar is False
        assert np.all(t2a._time.jd1 == t._time.jd1[4:6])
        t2b = Time([t[4], t[5]])
        assert t2b.isscalar is False
        assert np.all(t2b._time.jd1 == t._time.jd1[4:6])
        t2c = Time((t[4], t[5]))
        assert t2c.isscalar is False
        assert np.all(t2c._time.jd1 == t._time.jd1[4:6])
        t.delta_tdb_tt = np.arange(len(t))  # Explicitly set (not testing .tdb)
        t3 = t[4:6]
        assert np.all(t3._delta_tdb_tt == t._delta_tdb_tt[4:6])
        t4 = Time(mjd, format='mjd', scale='utc',
                  location=(np.arange(len(mjd)), np.arange(len(mjd))))
        t5a = t4[3]
        assert t5a.location == t4.location[3]
        assert t5a.location.shape == ()
        t5b = t4[3:4]
        assert t5b.location.shape == (1,)
        # Check that indexing a size-1 array returns a scalar location as well;
        # see gh-10113.
        t5c = t5b[0]
        assert t5c.location.shape == ()
        t6 = t4[4:6]
        assert np.all(t6.location == t4.location[4:6])
        # check it is a view
        # (via ndarray, since quantity setter problematic for structured array)
        allzeros = np.array((0., 0., 0.), dtype=t4.location.dtype)
        assert t6.location.view(np.ndarray)[-1] != allzeros
        assert t4.location.view(np.ndarray)[5] != allzeros
        t6.location.view(np.ndarray)[-1] = allzeros
        assert t4.location.view(np.ndarray)[5] == allzeros
        # Test subscription also works for two-dimensional arrays.
        frac = np.arange(0., 0.999, 0.2)
        t7 = Time(mjd[:, np.newaxis] + frac, format='mjd', scale='utc',
                  location=('45d', '50d'))
        assert t7[0, 0]._time.jd1 == t7._time.jd1[0, 0]
        assert t7[0, 0].isscalar is True
        assert np.all(t7[5]._time.jd1 == t7._time.jd1[5])
        assert np.all(t7[5]._time.jd2 == t7._time.jd2[5])
        assert np.all(t7[:, 2]._time.jd1 == t7._time.jd1[:, 2])
        assert np.all(t7[:, 2]._time.jd2 == t7._time.jd2[:, 2])
        assert np.all(t7[:, 0]._time.jd1 == t._time.jd1)
        assert np.all(t7[:, 0]._time.jd2 == t._time.jd2)
        # Get tdb to check that delta_tdb_tt attribute is sliced properly.
        t7_tdb = t7.tdb
        assert t7_tdb[0, 0].delta_tdb_tt == t7_tdb.delta_tdb_tt[0, 0]
        assert np.all(t7_tdb[5].delta_tdb_tt == t7_tdb.delta_tdb_tt[5])
        assert np.all(t7_tdb[:, 2].delta_tdb_tt == t7_tdb.delta_tdb_tt[:, 2])
        # Explicitly set delta_tdb_tt attribute. Now it should not be sliced.
        t7.delta_tdb_tt = 0.1
        t7_tdb2 = t7.tdb
        assert t7_tdb2[0, 0].delta_tdb_tt == 0.1
        assert t7_tdb2[5].delta_tdb_tt == 0.1
        assert t7_tdb2[:, 2].delta_tdb_tt == 0.1
        # Check broadcasting of location.
        t8 = Time(mjd[:, np.newaxis] + frac, format='mjd', scale='utc',
                  location=(np.arange(len(frac)), np.arange(len(frac))))
        assert t8[0, 0].location == t8.location[0, 0]
        assert np.all(t8[5].location == t8.location[5])
        assert np.all(t8[:, 2].location == t8.location[:, 2])
        # Finally check empty array.
        t9 = t[:0]
        assert t9.isscalar is False
        assert t9.shape == (0,)
        assert t9.size == 0

    def test_properties(self):
        """Use properties to convert scales and formats.  Note that the UT1 to
        UTC transformation requires a supplementary value (``delta_ut1_utc``)
        that can be obtained by interpolating from a table supplied by IERS.
        This is tested separately."""

        t = Time('2010-01-01 00:00:00', format='iso', scale='utc')
        t.delta_ut1_utc = 0.3341  # Explicitly set one part of the xform
        assert allclose_jd(t.jd, 2455197.5)
        assert t.iso == '2010-01-01 00:00:00.000'
        assert t.tt.iso == '2010-01-01 00:01:06.184'
        assert t.tai.fits == '2010-01-01T00:00:34.000'
        assert allclose_jd(t.utc.jd, 2455197.5)
        assert allclose_jd(t.ut1.jd, 2455197.500003867)
        assert t.tcg.isot == '2010-01-01T00:01:06.910'
        assert allclose_sec(t.unix, 1262304000.0)
        assert allclose_sec(t.cxcsec, 378691266.184)
        assert allclose_sec(t.gps, 946339215.0)
        assert t.datetime == datetime.datetime(2010, 1, 1)

    def test_precision(self):
        """Set the output precision which is used for some formats.  This is
        also a test of the code that provides a dict for global and instance
        options."""

        t = Time('2010-01-01 00:00:00', format='iso', scale='utc')
        # Uses initial class-defined precision=3
        assert t.iso == '2010-01-01 00:00:00.000'

        # Set instance precision to 9
        t.precision = 9
        assert t.iso == '2010-01-01 00:00:00.000000000'
        assert t.tai.utc.iso == '2010-01-01 00:00:00.000000000'

    def test_transforms(self):
        """Transform from UTC to all supported time scales (TAI, TCB, TCG,
        TDB, TT, UT1, UTC).  This requires auxiliary information (latitude and
        longitude)."""

        lat = 19.48125
        lon = -155.933222
        t = Time('2006-01-15 21:24:37.5', format='iso', scale='utc',
                 precision=7, location=(lon, lat))
        t.delta_ut1_utc = 0.3341  # Explicitly set one part of the xform
        assert t.utc.iso == '2006-01-15 21:24:37.5000000'
        assert t.ut1.iso == '2006-01-15 21:24:37.8341000'
        assert t.tai.iso == '2006-01-15 21:25:10.5000000'
        assert t.tt.iso == '2006-01-15 21:25:42.6840000'
        assert t.tcg.iso == '2006-01-15 21:25:43.3226905'
        assert t.tdb.iso == '2006-01-15 21:25:42.6843728'
        assert t.tcb.iso == '2006-01-15 21:25:56.8939523'

    def test_transforms_no_location(self):
        """Location should default to geocenter (relevant for TDB, TCB)."""
        t = Time('2006-01-15 21:24:37.5', format='iso', scale='utc',
                 precision=7)
        t.delta_ut1_utc = 0.3341  # Explicitly set one part of the xform
        assert t.utc.iso == '2006-01-15 21:24:37.5000000'
        assert t.ut1.iso == '2006-01-15 21:24:37.8341000'
        assert t.tai.iso == '2006-01-15 21:25:10.5000000'
        assert t.tt.iso == '2006-01-15 21:25:42.6840000'
        assert t.tcg.iso == '2006-01-15 21:25:43.3226905'
        assert t.tdb.iso == '2006-01-15 21:25:42.6843725'
        assert t.tcb.iso == '2006-01-15 21:25:56.8939519'
        # Check we get the same result
        t2 = Time('2006-01-15 21:24:37.5', format='iso', scale='utc',
                  location=(0*u.m, 0*u.m, 0*u.m))
        assert t == t2
        assert t.tdb == t2.tdb

    def test_location(self):
        """Check that location creates an EarthLocation object, and that
        such objects can be used as arguments.
        """
        lat = 19.48125
        lon = -155.933222
        t = Time(['2006-01-15 21:24:37.5'], format='iso', scale='utc',
                 precision=6, location=(lon, lat))
        assert isinstance(t.location, EarthLocation)
        location = EarthLocation(lon, lat)
        t2 = Time(['2006-01-15 21:24:37.5'], format='iso', scale='utc',
                  precision=6, location=location)
        assert isinstance(t2.location, EarthLocation)
        assert t2.location == t.location
        t3 = Time(['2006-01-15 21:24:37.5'], format='iso', scale='utc',
                  precision=6, location=(location.x, location.y, location.z))
        assert isinstance(t3.location, EarthLocation)
        assert t3.location == t.location

    def test_location_array(self):
        """Check that location arrays are checked for size and used
        for the corresponding times.  Also checks that erfa
        can handle array-valued locations, and can broadcast these if needed.
        """

        lat = 19.48125
        lon = -155.933222
        t = Time(['2006-01-15 21:24:37.5'] * 2, format='iso', scale='utc',
                 precision=6, location=(lon, lat))
        assert np.all(t.utc.iso == '2006-01-15 21:24:37.500000')
        assert np.all(t.tdb.iso[0] == '2006-01-15 21:25:42.684373')
        t2 = Time(['2006-01-15 21:24:37.5'] * 2, format='iso', scale='utc',
                  precision=6, location=(np.array([lon, 0]),
                                         np.array([lat, 0])))
        assert np.all(t2.utc.iso == '2006-01-15 21:24:37.500000')
        assert t2.tdb.iso[0] == '2006-01-15 21:25:42.684373'
        assert t2.tdb.iso[1] != '2006-01-15 21:25:42.684373'
        with pytest.raises(ValueError):  # 1 time, but two locations
            Time('2006-01-15 21:24:37.5', format='iso', scale='utc',
                 precision=6, location=(np.array([lon, 0]),
                                        np.array([lat, 0])))
        with pytest.raises(ValueError):  # 3 times, but two locations
            Time(['2006-01-15 21:24:37.5'] * 3, format='iso', scale='utc',
                 precision=6, location=(np.array([lon, 0]),
                                        np.array([lat, 0])))
        # multidimensional
        mjd = np.arange(50000., 50008.).reshape(4, 2)
        t3 = Time(mjd, format='mjd', scale='utc', location=(lon, lat))
        assert t3.shape == (4, 2)
        assert t3.location.shape == ()
        assert t3.tdb.shape == t3.shape
        t4 = Time(mjd, format='mjd', scale='utc',
                  location=(np.array([lon, 0]), np.array([lat, 0])))
        assert t4.shape == (4, 2)
        assert t4.location.shape == t4.shape
        assert t4.tdb.shape == t4.shape
        t5 = Time(mjd, format='mjd', scale='utc',
                  location=(np.array([[lon], [0], [0], [0]]),
                            np.array([[lat], [0], [0], [0]])))
        assert t5.shape == (4, 2)
        assert t5.location.shape == t5.shape
        assert t5.tdb.shape == t5.shape

    def test_all_scale_transforms(self):
        """Test that standard scale transforms work.  Does not test correctness,
        except reversibility [#2074]. Also tests that standard scales can't be
        converted to local scales"""
        lat = 19.48125
        lon = -155.933222
        with iers.conf.set_temp('auto_download', False):
            for scale1 in STANDARD_TIME_SCALES:
                t1 = Time('2006-01-15 21:24:37.5', format='iso', scale=scale1,
                          location=(lon, lat))
                for scale2 in STANDARD_TIME_SCALES:
                    t2 = getattr(t1, scale2)
                    t21 = getattr(t2, scale1)
                    assert allclose_jd(t21.jd, t1.jd)

                # test for conversion to local scale
                scale3 = 'local'
                with pytest.raises(ScaleValueError):
                    t2 = getattr(t1, scale3)

    def test_creating_all_formats(self):
        """Create a time object using each defined format"""
        Time(2000.5, format='decimalyear')
        Time(100.0, format='cxcsec')
        Time(100.0, format='unix')
        Time(100.0, format='gps')
        Time(1950.0, format='byear', scale='tai')
        Time(2000.0, format='jyear', scale='tai')
        Time('B1950.0', format='byear_str', scale='tai')
        Time('J2000.0', format='jyear_str', scale='tai')
        Time('2000-01-01 12:23:34.0', format='iso', scale='tai')
        Time('2000-01-01 12:23:34.0Z', format='iso', scale='utc')
        Time('2000-01-01T12:23:34.0', format='isot', scale='tai')
        Time('2000-01-01T12:23:34.0Z', format='isot', scale='utc')
        Time('2000-01-01T12:23:34.0', format='fits')
        Time('2000-01-01T12:23:34.0', format='fits', scale='tdb')
        Time(2400000.5, 51544.0333981, format='jd', scale='tai')
        Time(0.0, 51544.0333981, format='mjd', scale='tai')
        Time('2000:001:12:23:34.0', format='yday', scale='tai')
        Time('2000:001:12:23:34.0Z', format='yday', scale='utc')
        dt = datetime.datetime(2000, 1, 2, 3, 4, 5, 123456)
        Time(dt, format='datetime', scale='tai')
        Time([dt, dt], format='datetime', scale='tai')
        dt64 = np.datetime64('2012-06-18T02:00:05.453000000')
        Time(dt64, format='datetime64', scale='tai')
        Time([dt64, dt64], format='datetime64', scale='tai')

    def test_local_format_transforms(self):
        """
        Test trasformation of local time to different formats
        Transformation to formats with reference time should give
        ScalevalueError
        """
        t = Time('2006-01-15 21:24:37.5', scale='local')
        assert_allclose(t.jd, 2453751.3921006946, atol=0.001 / 3600. / 24., rtol=0.)
        assert_allclose(t.mjd, 53750.892100694444, atol=0.001 / 3600. / 24., rtol=0.)
        assert_allclose(t.decimalyear, 2006.0408002758752, atol=0.001 / 3600. / 24. / 365., rtol=0.)
        assert t.datetime == datetime.datetime(2006, 1, 15, 21, 24, 37, 500000)
        assert t.isot == '2006-01-15T21:24:37.500'
        assert t.yday == '2006:015:21:24:37.500'
        assert t.fits == '2006-01-15T21:24:37.500'
        assert_allclose(t.byear, 2006.04217888831, atol=0.001 / 3600. / 24. / 365., rtol=0.)
        assert_allclose(t.jyear, 2006.0407723496082, atol=0.001 / 3600. / 24. / 365., rtol=0.)
        assert t.byear_str == 'B2006.042'
        assert t.jyear_str == 'J2006.041'

        # epochTimeFormats
        with pytest.raises(ScaleValueError):
            t.gps
        with pytest.raises(ScaleValueError):
            t.unix
        with pytest.raises(ScaleValueError):
            t.cxcsec
        with pytest.raises(ScaleValueError):
            t.plot_date

    def test_datetime(self):
        """
        Test datetime format, including guessing the format from the input type
        by not providing the format keyword to Time.
        """
        dt = datetime.datetime(2000, 1, 2, 3, 4, 5, 123456)
        dt2 = datetime.datetime(2001, 1, 1)
        t = Time(dt, scale='utc', precision=9)
        assert t.iso == '2000-01-02 03:04:05.123456000'
        assert t.datetime == dt
        assert t.value == dt
        t2 = Time(t.iso, scale='utc')
        assert t2.datetime == dt

        t = Time([dt, dt2], scale='utc')
        assert np.all(t.value == [dt, dt2])

        t = Time('2000-01-01 01:01:01.123456789', scale='tai')
        assert t.datetime == datetime.datetime(2000, 1, 1, 1, 1, 1, 123457)

        # broadcasting
        dt3 = (dt + (dt2-dt) * np.arange(12)).reshape(4, 3)
        t3 = Time(dt3, scale='utc')
        assert t3.shape == (4, 3)
        assert t3[2, 1].value == dt3[2, 1]
        assert t3[2, 1] == Time(dt3[2, 1])
        assert np.all(t3.value == dt3)
        assert np.all(t3[1].value == dt3[1])
        assert np.all(t3[:, 2] == Time(dt3[:, 2]))
        assert Time(t3[2, 0]) == t3[2, 0]

    def test_datetime64(self):
        dt64 = np.datetime64('2000-01-02T03:04:05.123456789')
        dt64_2 = np.datetime64('2000-01-02')
        t = Time(dt64, scale='utc', precision=9, format='datetime64')
        assert t.iso == '2000-01-02 03:04:05.123456789'
        assert t.datetime64 == dt64
        assert t.value == dt64
        t2 = Time(t.iso, scale='utc')
        assert t2.datetime64 == dt64

        t = Time(dt64_2, scale='utc', precision=3, format='datetime64')
        assert t.iso == '2000-01-02 00:00:00.000'
        assert t.datetime64 == dt64_2
        assert t.value == dt64_2
        t2 = Time(t.iso, scale='utc')
        assert t2.datetime64 == dt64_2

        t = Time([dt64, dt64_2], scale='utc', format='datetime64')
        assert np.all(t.value == [dt64, dt64_2])

        t = Time('2000-01-01 01:01:01.123456789', scale='tai')
        assert t.datetime64 == np.datetime64('2000-01-01T01:01:01.123456789')

        # broadcasting
        dt3 = (dt64 + (dt64_2-dt64) * np.arange(12)).reshape(4, 3)
        t3 = Time(dt3, scale='utc', format='datetime64')
        assert t3.shape == (4, 3)
        assert t3[2, 1].value == dt3[2, 1]
        assert t3[2, 1] == Time(dt3[2, 1], format='datetime64')
        assert np.all(t3.value == dt3)
        assert np.all(t3[1].value == dt3[1])
        assert np.all(t3[:, 2] == Time(dt3[:, 2], format='datetime64'))
        assert Time(t3[2, 0], format='datetime64') == t3[2, 0]

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
        with pytest.raises(ValueError):  # unguessable
            Time([])
        with pytest.raises(ValueError):
            Time([50000.0], ['bad'], format='mjd', scale='tai')
        with pytest.raises(ValueError):
            Time(50000.0, 'bad', format='mjd', scale='tai')
        with pytest.raises(ValueError):
            Time('2005-08-04T00:01:02.000Z', scale='tai')
        # regression test against #3396
        with pytest.raises(ValueError):
            Time(np.nan, format='jd', scale='utc')
        with pytest.raises(ValueError):
            with pytest.warns(AstropyDeprecationWarning):
                Time('2000-01-02T03:04:05(TAI)', scale='utc')
        with pytest.raises(ValueError):
            Time('2000-01-02T03:04:05(TAI')
        with pytest.raises(ValueError):
            Time('2000-01-02T03:04:05(UT(NIST)')

    def test_utc_leap_sec(self):
        """Time behaves properly near or in UTC leap second.  This
        uses the 2012-06-30 leap second for testing."""
        for year, month, day in ((2012, 6, 30), (2016, 12, 31)):
            # Start with a day without a leap second and note rollover
            yyyy_mm = f'{year:04d}-{month:02d}'
            yyyy_mm_dd = f'{year:04d}-{month:02d}-{day:02d}'
            with pytest.warns(ErfaWarning):
                t1 = Time(yyyy_mm + '-01 23:59:60.0', scale='utc')
            assert t1.iso == yyyy_mm + '-02 00:00:00.000'

            # Leap second is different
            t1 = Time(yyyy_mm_dd + ' 23:59:59.900', scale='utc')
            assert t1.iso == yyyy_mm_dd + ' 23:59:59.900'

            t1 = Time(yyyy_mm_dd + ' 23:59:60.000', scale='utc')
            assert t1.iso == yyyy_mm_dd + ' 23:59:60.000'

            t1 = Time(yyyy_mm_dd + ' 23:59:60.999', scale='utc')
            assert t1.iso == yyyy_mm_dd + ' 23:59:60.999'

            if month == 6:
                yyyy_mm_dd_plus1 = f'{year:04d}-07-01'
            else:
                yyyy_mm_dd_plus1 = '{:04d}-01-01'.format(year + 1)

            with pytest.warns(ErfaWarning):
                t1 = Time(yyyy_mm_dd + ' 23:59:61.0', scale='utc')
            assert t1.iso == yyyy_mm_dd_plus1 + ' 00:00:00.000'

            # Delta time gives 2 seconds here as expected
            t0 = Time(yyyy_mm_dd + ' 23:59:59', scale='utc')
            t1 = Time(yyyy_mm_dd_plus1 + ' 00:00:00', scale='utc')
            assert allclose_sec((t1 - t0).sec, 2.0)

    def test_init_from_time_objects(self):
        """Initialize from one or more Time objects"""
        t1 = Time('2007:001', scale='tai')
        t2 = Time(['2007-01-02', '2007-01-03'], scale='utc')
        # Init from a list of Time objects without an explicit scale
        t3 = Time([t1, t2])
        # Test that init appropriately combines a scalar (t1) and list (t2)
        # and that scale and format are same as first element.
        assert len(t3) == 3
        assert t3.scale == t1.scale
        assert t3.format == t1.format  # t1 format is yday
        assert np.all(t3.value == np.concatenate([[t1.yday], t2.tai.yday]))

        # Init from a single Time object without a scale
        t3 = Time(t1)
        assert t3.isscalar
        assert t3.scale == t1.scale
        assert t3.format == t1.format
        assert np.all(t3.value == t1.value)

        # Init from a single Time object with scale specified
        t3 = Time(t1, scale='utc')
        assert t3.scale == 'utc'
        assert np.all(t3.value == t1.utc.value)

        # Init from a list of Time object with scale specified
        t3 = Time([t1, t2], scale='tt')
        assert t3.scale == 'tt'
        assert t3.format == t1.format  # yday
        assert np.all(t3.value == np.concatenate([[t1.tt.yday], t2.tt.yday]))

        # OK, how likely is this... but might as well test.
        mjd = np.arange(50000., 50006.)
        frac = np.arange(0., 0.999, 0.2)
        t4 = Time(mjd[:, np.newaxis] + frac, format='mjd', scale='utc')
        t5 = Time([t4[:2], t4[4:5]])
        assert t5.shape == (3, 5)

        # throw error when deriving local scale time
        # from non local time scale
        with pytest.raises(ValueError):
            Time(t1, scale='local')


class TestVal2:
    """Tests related to val2"""

    @pytest.mark.parametrize("d", [
        dict(val="2001:001", val2="ignored", scale="utc"),
        dict(val={'year': 2015, 'month': 2, 'day': 3,
                  'hour': 12, 'minute': 13, 'second': 14.567},
             val2="ignored", scale="utc"),
        dict(val=np.datetime64('2005-02-25'), val2="ignored", scale="utc"),
        dict(val=datetime.datetime(2000, 1, 2, 12, 0, 0),
             val2="ignored", scale="utc"),
    ])
    def test_unused_val2_raises(self, d):
        """Test that providing val2 is for string input lets user know we won't use it"""
        with pytest.raises(ValueError):
            Time(**d)

    def test_val2(self):
        """Various tests of the val2 input"""
        t = Time([0.0, 50000.0], [50000.0, 0.0], format='mjd', scale='tai')
        assert t.mjd[0] == t.mjd[1]
        assert t.jd[0] == t.jd[1]

    def test_val_broadcasts_against_val2(self):
        mjd = np.arange(50000., 50007.)
        frac = np.arange(0., 0.999, 0.2)
        t = Time(mjd[:, np.newaxis], frac, format='mjd', scale='utc')
        assert t.shape == (7, 5)
        with pytest.raises(ValueError):
            Time([0.0, 50000.0], [0.0, 1.0, 2.0], format='mjd', scale='tai')

    def test_broadcast_not_writable(self):
        val = (2458000 + np.arange(3))[:, None]
        val2 = np.linspace(0, 1, 4, endpoint=False)
        t = Time(val=val, val2=val2, format="jd", scale="tai")
        t_b = Time(val=val + 0 * val2, val2=0 * val + val2, format="jd", scale="tai")
        t_i = Time(val=57990, val2=0.3, format="jd", scale="tai")
        t_b[1, 2] = t_i
        t[1, 2] = t_i
        assert t_b[1, 2] == t[1, 2], "writing worked"
        assert t_b[0, 2] == t[0, 2], "broadcasting didn't cause problems"
        assert t_b[1, 1] == t[1, 1], "broadcasting didn't cause problems"
        assert np.all(t_b == t), "behaved as expected"

    def test_broadcast_one_not_writable(self):
        val = (2458000 + np.arange(3))
        val2 = np.arange(1)
        t = Time(val=val, val2=val2, format="jd", scale="tai")
        t_b = Time(val=val + 0 * val2, val2=0 * val + val2, format="jd", scale="tai")
        t_i = Time(val=57990, val2=0.3, format="jd", scale="tai")
        t_b[1] = t_i
        t[1] = t_i
        assert t_b[1] == t[1], "writing worked"
        assert t_b[0] == t[0], "broadcasting didn't cause problems"
        assert np.all(t_b == t), "behaved as expected"


class TestSubFormat:
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

    def test_fits_format(self):
        """FITS format includes bigger years."""
        # Heterogeneous input formats with in_subfmt='*' (default)
        times = ['2000-01-01', '2000-01-01T01:01:01', '2000-01-01T01:01:01.123']
        t = Time(times, format='fits', scale='tai')
        assert np.all(t.fits == np.array(['2000-01-01T00:00:00.000',
                                          '2000-01-01T01:01:01.000',
                                          '2000-01-01T01:01:01.123']))
        # Explicit long format for output, default scale is UTC.
        t2 = Time(times, format='fits', out_subfmt='long*')
        assert np.all(t2.fits == np.array(['+02000-01-01T00:00:00.000',
                                           '+02000-01-01T01:01:01.000',
                                           '+02000-01-01T01:01:01.123']))
        # Implicit long format for output, because of negative year.
        times[2] = '-00594-01-01'
        t3 = Time(times, format='fits', scale='tai')
        assert np.all(t3.fits == np.array(['+02000-01-01T00:00:00.000',
                                           '+02000-01-01T01:01:01.000',
                                           '-00594-01-01T00:00:00.000']))
        # Implicit long format for output, because of large positive year.
        times[2] = '+10594-01-01'
        t4 = Time(times, format='fits', scale='tai')
        assert np.all(t4.fits == np.array(['+02000-01-01T00:00:00.000',
                                           '+02000-01-01T01:01:01.000',
                                           '+10594-01-01T00:00:00.000']))

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
        t = Time(100.0, format='gps', scale='utc')
        assert t.scale == 'utc'

        # Check that bad scale is caught when format is specified
        with pytest.raises(ScaleValueError):
            Time(1950.0, format='byear', scale='bad scale')

        # Check that bad scale is caught when format is auto-determined
        with pytest.raises(ScaleValueError):
            Time('2000:001:00:00:00', scale='bad scale')

    def test_fits_scale(self):
        """Test that the previous FITS-string formatting can still be handled
        but with a DeprecationWarning."""
        for inputs in (("2000-01-02(TAI)", "tai"),
                       ("1999-01-01T00:00:00.123(ET(NIST))", "tt"),
                       ("2014-12-12T01:00:44.1(UTC)", "utc")):
            with pytest.warns(AstropyDeprecationWarning):
                t = Time(inputs[0])
            assert t.scale == inputs[1]

            # Create Time using normal ISOT syntax and compare with FITS
            t2 = Time(inputs[0][:inputs[0].index("(")], format="isot",
                      scale=inputs[1])
            assert t == t2

        # Explicit check that conversions still work despite warning
        with pytest.warns(AstropyDeprecationWarning):
            t = Time('1999-01-01T00:00:00.123456789(UTC)')
        t = t.tai
        assert t.isot == '1999-01-01T00:00:32.123'

        with pytest.warns(AstropyDeprecationWarning):
            t = Time('1999-01-01T00:00:32.123456789(TAI)')
        t = t.utc
        assert t.isot == '1999-01-01T00:00:00.123'

        # Check scale consistency
        with pytest.warns(AstropyDeprecationWarning):
            t = Time('1999-01-01T00:00:32.123456789(TAI)', scale="tai")
        assert t.scale == "tai"
        with pytest.warns(AstropyDeprecationWarning):
            t = Time('1999-01-01T00:00:32.123456789(ET)', scale="tt")
        assert t.scale == "tt"
        with pytest.raises(ValueError), pytest.warns(AstropyDeprecationWarning):
            t = Time('1999-01-01T00:00:32.123456789(TAI)', scale="utc")

    def test_scale_default(self):
        """Test behavior when no scale is provided"""
        # These first three are TimeFromEpoch and have an intrinsic time scale
        t = Time(100.0, format='cxcsec')
        assert t.scale == 'tt'
        t = Time(100.0, format='unix')
        assert t.scale == 'utc'
        t = Time(100.0, format='gps')
        assert t.scale == 'tai'

        for date in ('2000:001', '2000-01-01T00:00:00'):
            t = Time(date)
            assert t.scale == 'utc'

        t = Time(2000.1, format='byear')
        assert t.scale == 'tt'
        t = Time('J2000')
        assert t.scale == 'tt'

    def test_epoch_times(self):
        """Test time formats derived from EpochFromTime"""
        t = Time(0.0, format='cxcsec', scale='tai')
        assert t.tt.iso == '1998-01-01 00:00:00.000'

        # Create new time object from this one and change scale, format
        t2 = Time(t, scale='tt', format='iso')
        assert t2.value == '1998-01-01 00:00:00.000'

        # Value take from Chandra.Time.DateTime('2010:001:00:00:00').secs
        t_cxcsec = 378691266.184
        t = Time(t_cxcsec, format='cxcsec', scale='utc')
        assert allclose_sec(t.value, t_cxcsec)
        assert allclose_sec(t.cxcsec, t_cxcsec)
        assert allclose_sec(t.tt.value, t_cxcsec)
        assert allclose_sec(t.tt.cxcsec, t_cxcsec)
        assert t.yday == '2010:001:00:00:00.000'
        t = Time('2010:001:00:00:00.000', scale='utc')
        assert allclose_sec(t.cxcsec, t_cxcsec)
        assert allclose_sec(t.tt.cxcsec, t_cxcsec)

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

    def test_plot_date(self):
        """Test the plot_date format.

        Depending on the situation with matplotlib, this can give different
        results because the plot date epoch time changed in matplotlib 3.3. This
        test tries to use the matplotlib date2num function to make the test
        independent of version, but if matplotlib isn't available then the code
        (and test) use the pre-3.3 epoch.
        """
        try:
            from matplotlib.dates import date2num
        except ImportError:
            # No matplotlib, in which case this uses the epoch 0000-12-31
            # as per matplotlib < 3.3.
            # Value from:
            #   matplotlib.dates.set_epoch('0000-12-31')
            #   val = matplotlib.dates.date2num('2000-01-01')
            val = 730120.0
        else:
            val = date2num(datetime.datetime(2000, 1, 1))
        t = Time('2000-01-01 00:00:00', scale='utc')
        assert np.allclose(t.plot_date, val, atol=1e-5, rtol=0)


class TestNumericalSubFormat:
    def test_explicit_example(self):
        t = Time('54321.000000000001', format='mjd')
        assert t == Time(54321, 1e-12, format='mjd')
        assert t.mjd == 54321.  # Lost precision!
        assert t.value == 54321.  # Lost precision!
        assert t.to_value('mjd') == 54321.  # Lost precision!
        assert t.to_value('mjd', subfmt='str') == '54321.000000000001'
        assert t.to_value('mjd', 'bytes') == b'54321.000000000001'
        expected_long = np.longdouble(54321.) + np.longdouble(1e-12)
        # Check we're the same to within the double holding jd2
        # (which is less precise than longdouble on arm64).
        assert np.allclose(t.to_value('mjd', subfmt='long'),
                           expected_long, rtol=0, atol=np.finfo(float).eps)
        t.out_subfmt = 'str'
        assert t.value == '54321.000000000001'
        assert t.to_value('mjd') == 54321.  # Lost precision!
        assert t.mjd == '54321.000000000001'
        assert t.to_value('mjd', subfmt='bytes') == b'54321.000000000001'
        assert t.to_value('mjd', subfmt='float') == 54321.  # Lost precision!
        t.out_subfmt = 'long'
        assert np.allclose(t.value, expected_long,
                           rtol=0., atol=np.finfo(float).eps)
        assert np.allclose(t.to_value('mjd', subfmt=None), expected_long,
                           rtol=0., atol=np.finfo(float).eps)
        assert np.allclose(t.mjd, expected_long,
                           rtol=0., atol=np.finfo(float).eps)
        assert t.to_value('mjd', subfmt='str') == '54321.000000000001'
        assert t.to_value('mjd', subfmt='float') == 54321.  # Lost precision!

    @pytest.mark.skipif(np.finfo(np.longdouble).eps >= np.finfo(float).eps,
                        reason="long double is the same as float")
    def test_explicit_longdouble(self):
        i = 54321
        # Create a different long double (which will give a different jd2
        # even when long doubles are more precise than Time, as on arm64).
        f = max(2.**(-np.finfo(np.longdouble).nmant) * 65536,
                np.finfo(float).eps)
        mjd_long = np.longdouble(i) + np.longdouble(f)
        assert mjd_long != i, "longdouble failure!"
        t = Time(mjd_long, format='mjd')
        expected = Time(i, f, format='mjd')
        assert abs(t - expected) <= 20. * u.ps
        t_float = Time(i + f, format='mjd')
        assert t_float == Time(i, format='mjd')
        assert t_float != t
        assert t.value == 54321.  # Lost precision!
        assert np.allclose(t.to_value('mjd', subfmt='long'), mjd_long,
                           rtol=0., atol=np.finfo(float).eps)
        t2 = Time(mjd_long, format='mjd', out_subfmt='long')
        assert np.allclose(t2.value, mjd_long,
                           rtol=0., atol=np.finfo(float).eps)

    @pytest.mark.skipif(np.finfo(np.longdouble).eps >= np.finfo(float).eps,
                        reason="long double is the same as float")
    def test_explicit_longdouble_one_val(self):
        """Ensure either val1 or val2 being longdouble is possible.

        Regression test for issue gh-10033.
        """
        i = 54321
        f = max(2.**(-np.finfo(np.longdouble).nmant) * 65536,
                np.finfo(float).eps)
        t1 = Time(i, f, format='mjd')
        t2 = Time(np.longdouble(i), f, format='mjd')
        t3 = Time(i, np.longdouble(f), format='mjd')
        t4 = Time(np.longdouble(i), np.longdouble(f), format='mjd')
        assert t1 == t2 == t3 == t4

    @pytest.mark.skipif(np.finfo(np.longdouble).eps >= np.finfo(float).eps,
                        reason="long double is the same as float")
    @pytest.mark.parametrize("fmt", ["mjd", "unix", "cxcsec"])
    def test_longdouble_for_other_types(self, fmt):
        t_fmt = getattr(Time(58000, format="mjd"), fmt)  # Get regular float
        t_fmt_long = np.longdouble(t_fmt)
        # Create a different long double (ensuring it will give a different jd2
        # even when long doubles are more precise than Time, as on arm64).
        atol = np.finfo(float).eps * (1. if fmt == 'mjd' else 24. * 3600.)
        t_fmt_long2 = t_fmt_long + max(
            t_fmt_long * np.finfo(np.longdouble).eps * 2, atol)
        assert t_fmt_long != t_fmt_long2, "longdouble weird!"
        tm = Time(t_fmt_long, format=fmt)
        tm2 = Time(t_fmt_long2, format=fmt)
        assert tm != tm2
        tm_long2 = tm2.to_value(fmt, subfmt='long')
        assert np.allclose(tm_long2, t_fmt_long2, rtol=0., atol=atol)

    def test_subformat_input(self):
        s = '54321.01234567890123456789'
        i, f = s.split('.')  # Note, OK only for fraction < 0.5
        t = Time(float(i), float('.' + f), format='mjd')
        t_str = Time(s, format='mjd')
        t_bytes = Time(s.encode('ascii'), format='mjd')
        t_decimal = Time(Decimal(s), format='mjd')
        assert t_str == t
        assert t_bytes == t
        assert t_decimal == t

    @pytest.mark.parametrize('out_subfmt', ('str', 'bytes'))
    def test_subformat_output(self, out_subfmt):
        i = 54321
        f = np.array([0., 1e-9, 1e-12])
        t = Time(i, f, format='mjd', out_subfmt=out_subfmt)
        t_value = t.value
        expected = np.array(['54321.0',
                             '54321.000000001',
                             '54321.000000000001'], dtype=out_subfmt)
        assert np.all(t_value == expected)
        assert np.all(Time(expected, format='mjd') == t)

        # Explicit sub-format.
        t = Time(i, f, format='mjd')
        t_mjd_subfmt = t.to_value('mjd', subfmt=out_subfmt)
        assert np.all(t_mjd_subfmt == expected)

    @pytest.mark.parametrize('fmt,string,val1,val2', [
        ('jd', '2451544.5333981', 2451544.5, .0333981),
        ('decimalyear', '2000.54321', 2000., .54321),
        ('cxcsec', '100.0123456', 100.0123456, None),
        ('unix', '100.0123456', 100.0123456, None),
        ('gps', '100.0123456', 100.0123456, None),
        ('byear', '1950.1', 1950.1, None),
        ('jyear', '2000.1', 2000.1, None)])
    def test_explicit_string_other_formats(self, fmt, string, val1, val2):
        t = Time(string, format=fmt)
        assert t == Time(val1, val2, format=fmt)
        assert t.to_value(fmt, subfmt='str') == string

    def test_basic_subformat_setting(self):
        t = Time('2001', format='jyear', scale='tai')
        t.format = "mjd"
        t.out_subfmt = "str"
        assert t.value.startswith("5")

    def test_basic_subformat_cache_does_not_crash(self):
        t = Time('2001', format='jyear', scale='tai')
        t.to_value('mjd', subfmt='str')
        assert ('mjd', 'str') in t.cache['format']
        t.to_value('mjd', 'str')

    @pytest.mark.parametrize("fmt", ["jd", "mjd", "cxcsec", "unix", "gps", "jyear"])
    def test_decimal_context_does_not_affect_string(self, fmt):
        t = Time('2001', format='jyear', scale='tai')
        t.format = fmt
        with localcontext() as ctx:
            ctx.prec = 2
            t_s_2 = t.to_value(fmt, "str")
        t2 = Time('2001', format='jyear', scale='tai')
        t2.format = fmt
        with localcontext() as ctx:
            ctx.prec = 40
            t2_s_40 = t.to_value(fmt, "str")
        assert t_s_2 == t2_s_40, "String representation should not depend on Decimal context"

    def test_decimal_context_caching(self):
        t = Time(val=58000, val2=1e-14, format='mjd', scale='tai')
        with localcontext() as ctx:
            ctx.prec = 2
            t_s_2 = t.to_value('mjd', subfmt='decimal')
        t2 = Time(val=58000, val2=1e-14, format='mjd', scale='tai')
        with localcontext() as ctx:
            ctx.prec = 40
            t_s_40 = t.to_value('mjd', subfmt='decimal')
            t2_s_40 = t2.to_value('mjd', subfmt='decimal')
        assert t_s_2 == t_s_40, "Should be the same but cache might make this automatic"
        assert t_s_2 == t2_s_40, "Different precision should produce the same results"

    @pytest.mark.parametrize("f, s, t", [("sec", "long", np.longdouble),
                                         ("sec", "decimal", Decimal),
                                         ("sec", "str", str)])
    def test_timedelta_basic(self, f, s, t):
        dt = (Time("58000", format="mjd", scale="tai")
              - Time("58001", format="mjd", scale="tai"))

        value = dt.to_value(f, s)
        assert isinstance(value, t)
        dt.format = f
        dt.out_subfmt = s
        assert isinstance(dt.value, t)
        assert isinstance(dt.to_value(f, None), t)

    def test_need_format_argument(self):
        t = Time('J2000')
        with pytest.raises(TypeError, match="missing.*required.*'format'"):
            t.to_value()
        with pytest.raises(ValueError, match='format must be one of'):
            t.to_value('julian')

    def test_wrong_in_subfmt(self):
        with pytest.raises(ValueError, match='not among selected'):
            Time("58000", format='mjd', in_subfmt='float')

        with pytest.raises(ValueError, match='not among selected'):
            Time(np.longdouble(58000), format='mjd', in_subfmt='float')

        with pytest.raises(ValueError, match='not among selected'):
            Time(58000., format='mjd', in_subfmt='str')

        with pytest.raises(ValueError, match='not among selected'):
            Time(58000., format='mjd', in_subfmt='long')

    def test_wrong_subfmt(self):
        t = Time(58000., format='mjd')
        with pytest.raises(ValueError, match='must match one'):
            t.to_value('mjd', subfmt='parrot')

        with pytest.raises(ValueError, match='must match one'):
            t.out_subfmt = 'parrot'

        with pytest.raises(ValueError, match='must match one'):
            t.in_subfmt = 'parrot'

    def test_not_allowed_subfmt(self):
        """Test case where format has no defined subfmts"""
        t = Time('J2000')
        match = 'subformat not allowed for format jyear_str'
        with pytest.raises(ValueError, match=match):
            t.to_value('jyear_str', subfmt='parrot')

        with pytest.raises(ValueError, match=match):
            t.out_subfmt = 'parrot'

        with pytest.raises(ValueError, match=match):
            Time('J2000', out_subfmt='parrot')

        with pytest.raises(ValueError, match=match):
            t.in_subfmt = 'parrot'

        with pytest.raises(ValueError, match=match):
            Time('J2000', format='jyear_str', in_subfmt='parrot')

    def test_switch_to_format_with_no_out_subfmt(self):
        t = Time('2001-01-01', out_subfmt='date_hm')
        assert t.out_subfmt == 'date_hm'

        # Now do an in-place switch to format 'jyear_str' that has no subfmts
        # where out_subfmt is changed to '*'.
        t.format = 'jyear_str'
        assert t.out_subfmt == '*'
        assert t.value == 'J2001.001'


class TestSofaErrors:
    """Test that erfa status return values are handled correctly"""

    def test_bad_time(self):
        iy = np.array([2000], dtype=np.intc)
        im = np.array([2000], dtype=np.intc)  # bad month
        id = np.array([2000], dtype=np.intc)  # bad day
        with pytest.raises(ValueError):  # bad month, fatal error
            djm0, djm = erfa.cal2jd(iy, im, id)

        iy[0] = -5000
        im[0] = 2
        with pytest.raises(ValueError):  # bad year, fatal error
            djm0, djm = erfa.cal2jd(iy, im, id)

        iy[0] = 2000
        with pytest.warns(ErfaWarning, match=r'bad day    \(JD computed\)') as w:
            djm0, djm = erfa.cal2jd(iy, im, id)
        assert len(w) == 1

        assert allclose_jd(djm0, [2400000.5])
        assert allclose_jd(djm, [53574.])


class TestCopyReplicate:
    """Test issues related to copying and replicating data"""

    def test_immutable_input(self):
        """Internals are never mutable."""
        jds = np.array([2450000.5], dtype=np.double)
        t = Time(jds, format='jd', scale='tai')
        assert allclose_jd(t.jd, jds)
        jds[0] = 2458654
        assert not allclose_jd(t.jd, jds)

        mjds = np.array([50000.0], dtype=np.double)
        t = Time(mjds, format='mjd', scale='tai')
        assert allclose_jd(t.jd, [2450000.5])
        mjds[0] = 0.0
        assert allclose_jd(t.jd, [2450000.5])

    def test_replicate(self):
        """Test replicate method"""
        t = Time(['2000:001'], format='yday', scale='tai',
                 location=('45d', '45d'))
        t_yday = t.yday
        t_loc_x = t.location.x.copy()
        t2 = t.replicate()
        assert t.yday == t2.yday
        assert t.format == t2.format
        assert t.scale == t2.scale
        assert t.location == t2.location
        # This is not allowed publicly, but here we hack the internal time
        # and location values to show that t and t2 are sharing references.
        t2._time.jd1 += 100.0

        # Need to delete the cached yday attributes (only an issue because
        # of the internal _time hack).
        del t.cache
        del t2.cache

        assert t.yday == t2.yday
        assert t.yday != t_yday  # prove that it changed
        t2_loc_x_view = t2.location.x
        t2_loc_x_view[()] = 0  # use 0 to avoid having to give units
        assert t2.location.x == t2_loc_x_view
        assert t.location.x == t2.location.x
        assert t.location.x != t_loc_x  # prove that it changed

    def test_copy(self):
        """Test copy method"""
        t = Time('2000:001', format='yday', scale='tai',
                 location=('45d', '45d'))
        t_yday = t.yday
        t_loc_x = t.location.x.copy()
        t2 = t.copy()
        assert t.yday == t2.yday
        # This is not allowed publicly, but here we hack the internal time
        # and location values to show that t and t2 are not sharing references.
        t2._time.jd1 += 100.0

        # Need to delete the cached yday attributes (only an issue because
        # of the internal _time hack).
        del t.cache
        del t2.cache

        assert t.yday != t2.yday
        assert t.yday == t_yday  # prove that it did not change
        t2_loc_x_view = t2.location.x
        t2_loc_x_view[()] = 0  # use 0 to avoid having to give units
        assert t2.location.x == t2_loc_x_view
        assert t.location.x != t2.location.x
        assert t.location.x == t_loc_x  # prove that it changed


class TestStardate:
    """Sync chronometers with Starfleet Command"""

    def test_iso_to_stardate(self):
        assert str(Time('2320-01-01', scale='tai').stardate)[:7] == '1368.99'
        assert str(Time('2330-01-01', scale='tai').stardate)[:8] == '10552.76'
        assert str(Time('2340-01-01', scale='tai').stardate)[:8] == '19734.02'

    @pytest.mark.parametrize('dates',
                             [(10000, '2329-05-26 03:02'),
                              (20000, '2340-04-15 19:05'),
                              (30000, '2351-03-07 11:08')])
    def test_stardate_to_iso(self, dates):
        stardate, iso = dates
        t_star = Time(stardate, format='stardate')
        t_iso = Time(t_star, format='iso', out_subfmt='date_hm')
        assert t_iso.value == iso


def test_python_builtin_copy():
    t = Time('2000:001', format='yday', scale='tai')
    t2 = copy.copy(t)
    t3 = copy.deepcopy(t)

    assert t.jd == t2.jd
    assert t.jd == t3.jd


def test_now():
    """
    Tests creating a Time object with the `now` class method.
    """

    now = datetime.datetime.utcnow()
    t = Time.now()

    assert t.format == 'datetime'
    assert t.scale == 'utc'

    dt = t.datetime - now  # a datetime.timedelta object

    # this gives a .1 second margin between the `utcnow` call and the `Time`
    # initializer, which is really way more generous than necessary - typical
    # times are more like microseconds.  But it seems safer in case some
    # platforms have slow clock calls or something.

    assert dt.total_seconds() < 0.1


def test_decimalyear():
    t = Time('2001:001', format='yday')
    assert t.decimalyear == 2001.0

    t = Time(2000.0, [0.5, 0.75], format='decimalyear')
    assert np.all(t.value == [2000.5, 2000.75])

    jd0 = Time('2000:001').jd
    jd1 = Time('2001:001').jd
    d_jd = jd1 - jd0
    assert np.all(t.jd == [jd0 + 0.5 * d_jd,
                           jd0 + 0.75 * d_jd])


def test_fits_year0():
    t = Time(1721425.5, format='jd', scale='tai')
    assert t.fits == '0001-01-01T00:00:00.000'
    t = Time(1721425.5 - 366., format='jd', scale='tai')
    assert t.fits == '+00000-01-01T00:00:00.000'
    t = Time(1721425.5 - 366. - 365., format='jd', scale='tai')
    assert t.fits == '-00001-01-01T00:00:00.000'


def test_fits_year10000():
    t = Time(5373484.5, format='jd', scale='tai')
    assert t.fits == '+10000-01-01T00:00:00.000'
    t = Time(5373484.5 - 365., format='jd', scale='tai')
    assert t.fits == '9999-01-01T00:00:00.000'
    t = Time(5373484.5, -1. / 24. / 3600., format='jd', scale='tai')
    assert t.fits == '9999-12-31T23:59:59.000'


def test_dir():
    t = Time('2000:001', format='yday', scale='tai')
    assert 'utc' in dir(t)


def test_time_from_epoch_jds():
    """Test that jd1/jd2 in a TimeFromEpoch format is always well-formed:
    jd1 is an integral value and abs(jd2) <= 0.5.
    """
    # From 1999:001 00:00 to 1999:002 12:00 by a non-round step. This will
    # catch jd2 == 0 and a case of abs(jd2) == 0.5.
    cxcsecs = np.linspace(0, 86400 * 1.5, 49)
    for cxcsec in cxcsecs:
        t = Time(cxcsec, format='cxcsec')
        assert np.round(t.jd1) == t.jd1
        assert np.abs(t.jd2) <= 0.5

    t = Time(cxcsecs, format='cxcsec')
    assert np.all(np.round(t.jd1) == t.jd1)
    assert np.all(np.abs(t.jd2) <= 0.5)
    assert np.any(np.abs(t.jd2) == 0.5)  # At least one exactly 0.5


def test_bool():
    """Any Time object should evaluate to True unless it is empty [#3520]."""
    t = Time(np.arange(50000, 50010), format='mjd', scale='utc')
    assert bool(t) is True
    assert bool(t[0]) is True
    assert bool(t[:0]) is False


def test_len_size():
    """Check length of Time objects and that scalar ones do not have one."""
    t = Time(np.arange(50000, 50010), format='mjd', scale='utc')
    assert len(t) == 10 and t.size == 10
    t1 = Time(np.arange(50000, 50010).reshape(2, 5), format='mjd', scale='utc')
    assert len(t1) == 2 and t1.size == 10
    # Can have length 1 or length 0 arrays.
    t2 = t[:1]
    assert len(t2) == 1 and t2.size == 1
    t3 = t[:0]
    assert len(t3) == 0 and t3.size == 0
    # But cannot get length from scalar.
    t4 = t[0]
    with pytest.raises(TypeError) as err:
        len(t4)
    # Ensure we're not just getting the old error of
    # "object of type 'float' has no len()".
    assert 'Time' in str(err.value)


def test_TimeFormat_scale():
    """guard against recurrence of #1122, where TimeFormat class looses uses
    attributes (delta_ut1_utc here), preventing conversion to unix, cxc"""
    t = Time('1900-01-01', scale='ut1')
    t.delta_ut1_utc = 0.0
    with pytest.warns(ErfaWarning):
        t.unix
        assert t.unix == t.utc.unix


@pytest.mark.remote_data
def test_scale_conversion(monkeypatch):
    # Check that if we have internet, and downloading is allowed, we
    # can get conversion to UT1 for the present, since we will download
    # IERS_A in IERS_Auto.
    monkeypatch.setattr('astropy.utils.iers.conf.auto_download', True)
    Time(Time.now().cxcsec, format='cxcsec', scale='ut1')


def test_byteorder():
    """Ensure that bigendian and little-endian both work (closes #2942)"""
    mjd = np.array([53000.00, 54000.00])
    big_endian = mjd.astype('>f8')
    little_endian = mjd.astype('<f8')
    time_mjd = Time(mjd, format='mjd')
    time_big = Time(big_endian, format='mjd')
    time_little = Time(little_endian, format='mjd')
    assert np.all(time_big == time_mjd)
    assert np.all(time_little == time_mjd)


def test_datetime_tzinfo():
    """
    Test #3160 that time zone info in datetime objects is respected.
    """
    class TZm6(datetime.tzinfo):
        def utcoffset(self, dt):
            return datetime.timedelta(hours=-6)

    d = datetime.datetime(2002, 1, 2, 10, 3, 4, tzinfo=TZm6())
    t = Time(d)
    assert t.value == datetime.datetime(2002, 1, 2, 16, 3, 4)


def test_subfmts_regex():
    """
    Test having a custom subfmts with a regular expression
    """
    class TimeLongYear(TimeString):
        name = 'longyear'
        subfmts = (('date',
                    r'(?P<year>[+-]\d{5})-%m-%d',  # hybrid
                    '{year:+06d}-{mon:02d}-{day:02d}'),)
    t = Time('+02000-02-03', format='longyear')
    assert t.value == '+02000-02-03'
    assert t.jd == Time('2000-02-03').jd


def test_set_format_basic():
    """
    Test basics of setting format attribute.
    """
    for format, value in (('jd', 2451577.5),
                          ('mjd', 51577.0),
                          ('cxcsec', 65923264.184),  # confirmed with Chandra.Time
                          ('datetime', datetime.datetime(2000, 2, 3, 0, 0)),
                          ('iso', '2000-02-03 00:00:00.000')):
        t = Time('+02000-02-03', format='fits')
        t0 = t.replicate()
        t.format = format
        assert t.value == value
        # Internal jd1 and jd2 are preserved
        assert t._time.jd1 is t0._time.jd1
        assert t._time.jd2 is t0._time.jd2


def test_unix_tai_format():
    t = Time('2020-01-01', scale='utc')
    assert allclose_sec(t.unix_tai - t.unix, 37.0)
    t = Time('1970-01-01', scale='utc')
    assert allclose_sec(t.unix_tai - t.unix, 8 + 8.2e-05)


def test_set_format_shares_subfmt():
    """
    Set format and round trip through a format that shares out_subfmt
    """
    t = Time('+02000-02-03', format='fits', out_subfmt='date_hms', precision=5)
    tc = t.copy()

    t.format = 'isot'
    assert t.precision == 5
    assert t.out_subfmt == 'date_hms'
    assert t.value == '2000-02-03T00:00:00.00000'

    t.format = 'fits'
    assert t.value == tc.value
    assert t.precision == 5


def test_set_format_does_not_share_subfmt():
    """
    Set format and round trip through a format that does not share out_subfmt
    """
    t = Time('+02000-02-03', format='fits', out_subfmt='longdate')

    t.format = 'isot'
    assert t.out_subfmt == '*'  # longdate_hms not there, goes to default
    assert t.value == '2000-02-03T00:00:00.000'

    t.format = 'fits'
    assert t.out_subfmt == '*'
    assert t.value == '2000-02-03T00:00:00.000'  # date_hms


def test_replicate_value_error():
    """
    Passing a bad format to replicate should raise ValueError, not KeyError.
    PR #3857.
    """
    t1 = Time('2007:001', scale='tai')
    with pytest.raises(ValueError) as err:
        t1.replicate(format='definitely_not_a_valid_format')
    assert 'format must be one of' in str(err.value)


def test_remove_astropy_time():
    """
    Make sure that 'astropy_time' format is really gone after #3857.  Kind of
    silly test but just to be sure.
    """
    t1 = Time('2007:001', scale='tai')
    assert 'astropy_time' not in t1.FORMATS
    with pytest.raises(ValueError) as err:
        Time(t1, format='astropy_time')
    assert 'format must be one of' in str(err.value)


def test_isiterable():
    """
    Ensure that scalar `Time` instances are not reported as iterable by the
    `isiterable` utility.

    Regression test for https://github.com/astropy/astropy/issues/4048
    """

    t1 = Time.now()
    assert not isiterable(t1)

    t2 = Time(['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00'],
              format='iso', scale='utc')
    assert isiterable(t2)


def test_to_datetime():
    tz = TimezoneInfo(utc_offset=-10 * u.hour, tzname='US/Hawaii')
    # The above lines produces a `datetime.tzinfo` object similar to:
    #     tzinfo = pytz.timezone('US/Hawaii')
    time = Time('2010-09-03 00:00:00')
    tz_aware_datetime = time.to_datetime(tz)
    assert tz_aware_datetime.time() == datetime.time(14, 0)
    forced_to_astropy_time = Time(tz_aware_datetime)
    assert tz.tzname(time.datetime) == tz_aware_datetime.tzname()
    assert time == forced_to_astropy_time

    # Test non-scalar time inputs:
    time = Time(['2010-09-03 00:00:00', '2005-09-03 06:00:00',
                 '1990-09-03 06:00:00'])
    tz_aware_datetime = time.to_datetime(tz)
    forced_to_astropy_time = Time(tz_aware_datetime)
    for dt, tz_dt in zip(time.datetime, tz_aware_datetime):
        assert tz.tzname(dt) == tz_dt.tzname()
    assert np.all(time == forced_to_astropy_time)

    with pytest.raises(ValueError, match=r'does not support leap seconds'):
        Time('2015-06-30 23:59:60.000').to_datetime()


@pytest.mark.skipif('not HAS_PYTZ')
def test_to_datetime_pytz():

    tz = pytz.timezone('US/Hawaii')
    time = Time('2010-09-03 00:00:00')
    tz_aware_datetime = time.to_datetime(tz)
    forced_to_astropy_time = Time(tz_aware_datetime)
    assert tz_aware_datetime.time() == datetime.time(14, 0)
    assert tz.tzname(time.datetime) == tz_aware_datetime.tzname()
    assert time == forced_to_astropy_time

    # Test non-scalar time inputs:
    time = Time(['2010-09-03 00:00:00', '2005-09-03 06:00:00',
                 '1990-09-03 06:00:00'])
    tz_aware_datetime = time.to_datetime(tz)
    forced_to_astropy_time = Time(tz_aware_datetime)
    for dt, tz_dt in zip(time.datetime, tz_aware_datetime):
        assert tz.tzname(dt) == tz_dt.tzname()
    assert np.all(time == forced_to_astropy_time)


def test_cache():
    t = Time('2010-09-03 00:00:00')
    t2 = Time('2010-09-03 00:00:00')

    # Time starts out without a cache
    assert 'cache' not in t._time.__dict__

    # Access the iso format and confirm that the cached version is as expected
    t.iso
    assert t.cache['format']['iso'] == t2.iso

    # Access the TAI scale and confirm that the cached version is as expected
    t.tai
    assert t.cache['scale']['tai'] == t2.tai

    # New Time object after scale transform does not have a cache yet
    assert 'cache' not in t.tt._time.__dict__

    # Clear the cache
    del t.cache
    assert 'cache' not in t._time.__dict__
    # Check accessing the cache creates an empty dictionary
    assert not t.cache
    assert 'cache' in t._time.__dict__


def test_epoch_date_jd_is_day_fraction():
    """
    Ensure that jd1 and jd2 of an epoch Time are respect the (day, fraction) convention
    (see #6638)
    """
    t0 = Time("J2000", scale="tdb")

    assert t0.jd1 == 2451545.0
    assert t0.jd2 == 0.0

    t1 = Time(datetime.datetime(2000, 1, 1, 12, 0, 0), scale="tdb")

    assert t1.jd1 == 2451545.0
    assert t1.jd2 == 0.0


def test_sum_is_equivalent():
    """
    Ensure that two equal dates defined in different ways behave equally (#6638)
    """
    t0 = Time("J2000", scale="tdb")
    t1 = Time("2000-01-01 12:00:00", scale="tdb")

    assert t0 == t1
    assert (t0 + 1 * u.second) == (t1 + 1 * u.second)


def test_string_valued_columns():
    # Columns have a nice shim that translates bytes to string as needed.
    # Ensure Time can handle these.  Use multi-d array just to be sure.
    times = [[[f'{y:04d}-{m:02d}-{d:02d}' for d in range(1, 3)]
              for m in range(5, 7)] for y in range(2012, 2014)]
    cutf32 = Column(times)
    cbytes = cutf32.astype('S')
    tutf32 = Time(cutf32)
    tbytes = Time(cbytes)
    assert np.all(tutf32 == tbytes)
    tutf32 = Time(Column(['B1950']))
    tbytes = Time(Column([b'B1950']))
    assert tutf32 == tbytes
    # Regression tests for arrays with entries with unequal length. gh-6903.
    times = Column([b'2012-01-01', b'2012-01-01T00:00:00'])
    assert np.all(Time(times) == Time(['2012-01-01', '2012-01-01T00:00:00']))


def test_bytes_input():
    tstring = '2011-01-02T03:04:05'
    tbytes = b'2011-01-02T03:04:05'
    assert tbytes.decode('ascii') == tstring
    t0 = Time(tstring)
    t1 = Time(tbytes)
    assert t1 == t0
    tarray = np.array(tbytes)
    assert tarray.dtype.kind == 'S'
    t2 = Time(tarray)
    assert t2 == t0


def test_writeable_flag():
    t = Time([1, 2, 3], format='cxcsec')
    t[1] = 5.0
    assert allclose_sec(t[1].value, 5.0)

    t.writeable = False
    with pytest.raises(ValueError) as err:
        t[1] = 5.0
    assert 'Time object is read-only. Make a copy()' in str(err.value)

    with pytest.raises(ValueError) as err:
        t[:] = 5.0
    assert 'Time object is read-only. Make a copy()' in str(err.value)

    t.writeable = True
    t[1] = 10.0
    assert allclose_sec(t[1].value, 10.0)

    # Scalar is writeable because it gets boxed into a zero-d array
    t = Time('2000:001', scale='utc')
    t[()] = '2000:002'
    assert t.value.startswith('2000:002')

    # Transformed attribute is not writeable
    t = Time(['2000:001', '2000:002'], scale='utc')
    t2 = t.tt  # t2 is read-only now because t.tt is cached
    with pytest.raises(ValueError) as err:
        t2[0] = '2005:001'
    assert 'Time object is read-only. Make a copy()' in str(err.value)


def test_setitem_location():
    loc = EarthLocation(x=[1, 2] * u.m, y=[3, 4] * u.m, z=[5, 6] * u.m)
    t = Time([[1, 2], [3, 4]], format='cxcsec', location=loc)

    # Succeeds because the right hand side makes no implication about
    # location and just inherits t.location
    t[0, 0] = 0
    assert allclose_sec(t.value, [[0, 2], [3, 4]])

    # Fails because the right hand side has location=None
    with pytest.raises(ValueError) as err:
        t[0, 0] = Time(-1, format='cxcsec')
    assert ('cannot set to Time with different location: '
            'expected location={} and '
            'got location=None'.format(loc[0])) in str(err.value)

    # Succeeds because the right hand side correctly sets location
    t[0, 0] = Time(-2, format='cxcsec', location=loc[0])
    assert allclose_sec(t.value, [[-2, 2], [3, 4]])

    # Fails because the right hand side has different location
    with pytest.raises(ValueError) as err:
        t[0, 0] = Time(-2, format='cxcsec', location=loc[1])
    assert ('cannot set to Time with different location: '
            'expected location={} and '
            'got location={}'.format(loc[0], loc[1])) in str(err.value)

    # Fails because the Time has None location and RHS has defined location
    t = Time([[1, 2], [3, 4]], format='cxcsec')
    with pytest.raises(ValueError) as err:
        t[0, 0] = Time(-2, format='cxcsec', location=loc[1])
    assert ('cannot set to Time with different location: '
            'expected location=None and '
            'got location={}'.format(loc[1])) in str(err.value)

    # Broadcasting works
    t = Time([[1, 2], [3, 4]], format='cxcsec', location=loc)
    t[0, :] = Time([-3, -4], format='cxcsec', location=loc)
    assert allclose_sec(t.value, [[-3, -4], [3, 4]])


def test_setitem_from_python_objects():
    t = Time([[1, 2], [3, 4]], format='cxcsec')
    assert t.cache == {}
    t.iso
    assert 'iso' in t.cache['format']
    assert np.all(t.iso == [['1998-01-01 00:00:01.000', '1998-01-01 00:00:02.000'],
                            ['1998-01-01 00:00:03.000', '1998-01-01 00:00:04.000']])

    # Setting item clears cache
    t[0, 1] = 100
    assert t.cache == {}
    assert allclose_sec(t.value, [[1, 100],
                                  [3, 4]])
    assert np.all(t.iso == [['1998-01-01 00:00:01.000', '1998-01-01 00:01:40.000'],
                            ['1998-01-01 00:00:03.000', '1998-01-01 00:00:04.000']])

    # Set with a float value
    t.iso
    t[1, :] = 200
    assert t.cache == {}
    assert allclose_sec(t.value, [[1, 100],
                                  [200, 200]])

    # Array of strings in yday format
    t[:, 1] = ['1998:002', '1998:003']
    assert allclose_sec(t.value, [[1, 86400 * 1],
                                  [200, 86400 * 2]])

    # Incompatible numeric value
    t = Time(['2000:001', '2000:002'])
    t[0] = '2001:001'
    with pytest.raises(ValueError) as err:
        t[0] = 100
    assert 'cannot convert value to a compatible Time object' in str(err.value)


def test_setitem_from_time_objects():
    """Set from existing Time object.
    """
    # Set from time object with different scale
    t = Time(['2000:001', '2000:002'], scale='utc')
    t2 = Time(['2000:010'], scale='tai')
    t[1] = t2[0]
    assert t.value[1] == t2.utc.value[0]

    # Time object with different scale and format
    t = Time(['2000:001', '2000:002'], scale='utc')
    t2.format = 'jyear'
    t[1] = t2[0]
    assert t.yday[1] == t2.utc.yday[0]


def test_setitem_bad_item():
    t = Time([1, 2], format='cxcsec')
    with pytest.raises(IndexError):
        t['asdf'] = 3


def test_setitem_deltas():
    """Setting invalidates any transform deltas"""
    t = Time([1, 2], format='cxcsec')
    t.delta_tdb_tt = [1, 2]
    t.delta_ut1_utc = [3, 4]
    t[1] = 3
    assert not hasattr(t, '_delta_tdb_tt')
    assert not hasattr(t, '_delta_ut1_utc')


def test_subclass():
    """Check that we can initialize subclasses with a Time instance."""
    # Ref: Issue gh-#7449 and PR gh-#7453.

    class _Time(Time):
        pass

    t1 = Time('1999-01-01T01:01:01')
    t2 = _Time(t1)

    assert t2.__class__ == _Time
    assert t1 == t2


def test_strftime_scalar():
    """Test of Time.strftime
    """
    time_string = '2010-09-03 06:00:00'
    t = Time(time_string)

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S') == time_string


def test_strftime_array():
    tstrings = ['2010-09-03 00:00:00', '2005-09-03 06:00:00',
                '1995-12-31 23:59:60']
    t = Time(tstrings)

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S').tolist() == tstrings


def test_strftime_array_2():
    tstrings = [['1998-01-01 00:00:01', '1998-01-01 00:00:02'],
                ['1998-01-01 00:00:03', '1995-12-31 23:59:60']]
    tstrings = np.array(tstrings)

    t = Time(tstrings)

    for format in t.FORMATS:
        t.format = format
        assert np.all(t.strftime('%Y-%m-%d %H:%M:%S') == tstrings)
        assert t.strftime('%Y-%m-%d %H:%M:%S').shape == tstrings.shape


def test_strftime_leapsecond():
    time_string = '1995-12-31 23:59:60'
    t = Time(time_string)

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S') == time_string


def test_strptime_scalar():
    """Test of Time.strptime
    """
    time_string = '2007-May-04 21:08:12'
    time_object = Time('2007-05-04 21:08:12')
    t = Time.strptime(time_string, '%Y-%b-%d %H:%M:%S')

    assert t == time_object


def test_strptime_array():
    """Test of Time.strptime
    """
    tstrings = [['1998-Jan-01 00:00:01', '1998-Jan-01 00:00:02'],
                ['1998-Jan-01 00:00:03', '1998-Jan-01 00:00:04']]
    tstrings = np.array(tstrings)

    time_object = Time([['1998-01-01 00:00:01', '1998-01-01 00:00:02'],
                        ['1998-01-01 00:00:03', '1998-01-01 00:00:04']])
    t = Time.strptime(tstrings, '%Y-%b-%d %H:%M:%S')

    assert np.all(t == time_object)
    assert t.shape == tstrings.shape


def test_strptime_badinput():
    tstrings = [1, 2, 3]
    with pytest.raises(TypeError):
        Time.strptime(tstrings, '%S')


def test_strptime_input_bytes_scalar():
    time_string = b'2007-May-04 21:08:12'
    time_object = Time('2007-05-04 21:08:12')
    t = Time.strptime(time_string, '%Y-%b-%d %H:%M:%S')

    assert t == time_object


def test_strptime_input_bytes_array():
    tstrings = [[b'1998-Jan-01 00:00:01', b'1998-Jan-01 00:00:02'],
                [b'1998-Jan-01 00:00:03', b'1998-Jan-01 00:00:04']]
    tstrings = np.array(tstrings)

    time_object = Time([['1998-01-01 00:00:01', '1998-01-01 00:00:02'],
                        ['1998-01-01 00:00:03', '1998-01-01 00:00:04']])
    t = Time.strptime(tstrings, '%Y-%b-%d %H:%M:%S')

    assert np.all(t == time_object)
    assert t.shape == tstrings.shape


def test_strptime_leapsecond():
    time_obj1 = Time('1995-12-31T23:59:60', format='isot')
    time_obj2 = Time.strptime('1995-Dec-31 23:59:60', '%Y-%b-%d %H:%M:%S')

    assert time_obj1 == time_obj2


def test_strptime_3_digit_year():
    time_obj1 = Time('0995-12-31T00:00:00', format='isot', scale='tai')
    time_obj2 = Time.strptime('0995-Dec-31 00:00:00', '%Y-%b-%d %H:%M:%S',
                              scale='tai')

    assert time_obj1 == time_obj2


def test_strptime_fracsec_scalar():
    time_string = '2007-May-04 21:08:12.123'
    time_object = Time('2007-05-04 21:08:12.123')
    t = Time.strptime(time_string, '%Y-%b-%d %H:%M:%S.%f')

    assert t == time_object


def test_strptime_fracsec_array():
    """Test of Time.strptime
    """
    tstrings = [['1998-Jan-01 00:00:01.123', '1998-Jan-01 00:00:02.000001'],
                ['1998-Jan-01 00:00:03.000900', '1998-Jan-01 00:00:04.123456']]
    tstrings = np.array(tstrings)

    time_object = Time([['1998-01-01 00:00:01.123', '1998-01-01 00:00:02.000001'],
                        ['1998-01-01 00:00:03.000900', '1998-01-01 00:00:04.123456']])
    t = Time.strptime(tstrings, '%Y-%b-%d %H:%M:%S.%f')

    assert np.all(t == time_object)
    assert t.shape == tstrings.shape


def test_strftime_scalar_fracsec():
    """Test of Time.strftime
    """
    time_string = '2010-09-03 06:00:00.123'
    t = Time(time_string)

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S.%f') == time_string


def test_strftime_scalar_fracsec_precision():
    time_string = '2010-09-03 06:00:00.123123123'
    t = Time(time_string)
    assert t.strftime('%Y-%m-%d %H:%M:%S.%f') == '2010-09-03 06:00:00.123'
    t.precision = 9
    assert t.strftime('%Y-%m-%d %H:%M:%S.%f') == '2010-09-03 06:00:00.123123123'


def test_strftime_array_fracsec():
    tstrings = ['2010-09-03 00:00:00.123000', '2005-09-03 06:00:00.000001',
                '1995-12-31 23:59:60.000900']
    t = Time(tstrings)
    t.precision = 6

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S.%f').tolist() == tstrings


def test_insert_time():
    tm = Time([1, 2], format='unix')

    # Insert a scalar using an auto-parsed string
    tm2 = tm.insert(1, '1970-01-01 00:01:00')
    assert np.all(tm2 == Time([1, 60, 2], format='unix'))

    # Insert scalar using a Time value
    tm2 = tm.insert(1, Time('1970-01-01 00:01:00'))
    assert np.all(tm2 == Time([1, 60, 2], format='unix'))

    # Insert length=1 array with a Time value
    tm2 = tm.insert(1, [Time('1970-01-01 00:01:00')])
    assert np.all(tm2 == Time([1, 60, 2], format='unix'))

    # Insert length=2 list with float values matching unix format.
    # Also actually provide axis=0 unlike all other tests.
    tm2 = tm.insert(1, [10, 20], axis=0)
    assert np.all(tm2 == Time([1, 10, 20, 2], format='unix'))

    # Insert length=2 np.array with float values matching unix format
    tm2 = tm.insert(1, np.array([10, 20]))
    assert np.all(tm2 == Time([1, 10, 20, 2], format='unix'))

    # Insert length=2 np.array with float values at the end
    tm2 = tm.insert(2, np.array([10, 20]))
    assert np.all(tm2 == Time([1, 2, 10, 20], format='unix'))

    # Insert length=2 np.array with float values at the beginning
    # with a negative index
    tm2 = tm.insert(-2, np.array([10, 20]))
    assert np.all(tm2 == Time([10, 20, 1, 2], format='unix'))


def test_insert_exceptions():
    tm = Time(1, format='unix')
    with pytest.raises(TypeError) as err:
        tm.insert(0, 50)
    assert 'cannot insert into scalar' in str(err.value)

    tm = Time([1, 2], format='unix')
    with pytest.raises(ValueError) as err:
        tm.insert(0, 50, axis=1)
    assert 'axis must be 0' in str(err.value)

    with pytest.raises(TypeError) as err:
        tm.insert(slice(None), 50)
    assert 'obj arg must be an integer' in str(err.value)

    with pytest.raises(IndexError) as err:
        tm.insert(-100, 50)
    assert 'index -100 is out of bounds for axis 0 with size 2' in str(err.value)


def test_datetime64_no_format():
    dt64 = np.datetime64('2000-01-02T03:04:05.123456789')
    t = Time(dt64, scale='utc', precision=9)
    assert t.iso == '2000-01-02 03:04:05.123456789'
    assert t.datetime64 == dt64
    assert t.value == dt64


def test_hash_time():
    loc1 = EarthLocation(1 * u.m, 2 * u.m, 3 * u.m)
    for loc in None, loc1:
        t = Time([1, 1, 2, 3], format='cxcsec', location=loc)
        t[3] = np.ma.masked
        h1 = hash(t[0])
        h2 = hash(t[1])
        h3 = hash(t[2])
        assert h1 == h2
        assert h1 != h3

        with pytest.raises(TypeError) as exc:
            hash(t)
        assert exc.value.args[0] == "unhashable type: 'Time' (must be scalar)"

        with pytest.raises(TypeError) as exc:
            hash(t[3])
        assert exc.value.args[0] == "unhashable type: 'Time' (value is masked)"

    t = Time(1, format='cxcsec', location=loc)
    t2 = Time(1, format='cxcsec')

    assert hash(t) != hash(t2)

    t = Time('2000:180', scale='utc')
    t2 = Time(t, scale='tai')
    assert t == t2
    assert hash(t) != hash(t2)


def test_hash_time_delta():

    t = TimeDelta([1, 1, 2, 3], format='sec')
    t[3] = np.ma.masked
    h1 = hash(t[0])
    h2 = hash(t[1])
    h3 = hash(t[2])
    assert h1 == h2
    assert h1 != h3

    with pytest.raises(TypeError) as exc:
        hash(t)
    assert exc.value.args[0] == "unhashable type: 'TimeDelta' (must be scalar)"

    with pytest.raises(TypeError) as exc:
        hash(t[3])
    assert exc.value.args[0] == "unhashable type: 'TimeDelta' (value is masked)"


def test_get_time_fmt_exception_messages():
    with pytest.raises(ValueError) as err:
        Time(10)
    assert "No time format was given, and the input is" in str(err.value)

    with pytest.raises(ValueError) as err:
        Time('2000:001', format='not-a-format')
    assert "Format 'not-a-format' is not one of the allowed" in str(err.value)

    with pytest.raises(ValueError) as err:
        Time('200')
    assert 'Input values did not match any of the formats where' in str(err.value)

    with pytest.raises(ValueError) as err:
        Time('200', format='iso')
    assert ('Input values did not match the format class iso:' + os.linesep
            + 'ValueError: Time 200 does not match iso format') == str(err.value)

    with pytest.raises(ValueError) as err:
        Time(200, format='iso')
    assert ('Input values did not match the format class iso:' + os.linesep
            + 'TypeError: Input values for iso class must be strings') == str(err.value)


def test_ymdhms_defaults():
    t1 = Time({'year': 2001}, format='ymdhms')
    assert t1 == Time('2001-01-01')


times_dict_ns = {
    'year': [2001, 2002],
    'month': [2, 3],
    'day': [4, 5],
    'hour': [6, 7],
    'minute': [8, 9],
    'second': [10, 11]
}
table_ns = Table(times_dict_ns)
struct_array_ns = table_ns.as_array()
rec_array_ns = struct_array_ns.view(np.recarray)
ymdhms_names = ('year', 'month', 'day', 'hour', 'minute', 'second')


@pytest.mark.parametrize('tm_input', [table_ns, struct_array_ns, rec_array_ns])
@pytest.mark.parametrize('kwargs', [{}, {'format': 'ymdhms'}])
@pytest.mark.parametrize('as_row', [False, True])
def test_ymdhms_init_from_table_like(tm_input, kwargs, as_row):
    time_ns = Time(['2001-02-04 06:08:10', '2002-03-05 07:09:11'])
    if as_row:
        tm_input = tm_input[0]
        time_ns = time_ns[0]

    tm = Time(tm_input, **kwargs)
    assert np.all(tm == time_ns)
    assert tm.value.dtype.names == ymdhms_names


def test_ymdhms_init_from_dict_array():
    times_dict_shape = {
        'year': [[2001, 2002],
                 [2003, 2004]],
        'month': [2, 3],
        'day': 4
    }
    time_shape = Time(
        [['2001-02-04', '2002-03-04'],
         ['2003-02-04', '2004-03-04']]
    )
    time = Time(times_dict_shape, format='ymdhms')

    assert np.all(time == time_shape)
    assert time.ymdhms.shape == time_shape.shape


@pytest.mark.parametrize('kwargs', [{}, {'format': 'ymdhms'}])
def test_ymdhms_init_from_dict_scalar(kwargs):
    """
    Test YMDHMS functionality for a dict input. This includes ensuring that
    key and attribute access work.  For extra fun use a time within a leap
    second.
    """
    time_dict = {
        'year': 2016,
        'month': 12,
        'day': 31,
        'hour': 23,
        'minute': 59,
        'second': 60.123456789}

    tm = Time(time_dict, **kwargs)

    assert tm == Time('2016-12-31T23:59:60.123456789')
    for attr in time_dict:
        for value in (tm.value[attr], getattr(tm.value, attr)):
            if attr == 'second':
                assert allclose_sec(time_dict[attr], value)
            else:
                assert time_dict[attr] == value

    # Now test initializing from a YMDHMS format time using the object
    tm_rt = Time(tm)
    assert tm_rt == tm
    assert tm_rt.format == 'ymdhms'

    # Test initializing from a YMDHMS value (np.void, i.e. recarray row)
    # without specified format.
    tm_rt = Time(tm.ymdhms)
    assert tm_rt == tm
    assert tm_rt.format == 'ymdhms'


def test_ymdhms_exceptions():
    with pytest.raises(ValueError, match='input must be dict or table-like'):
        Time(10, format='ymdhms')

    match = "'wrong' not allowed as YMDHMS key name(s)"
    # NB: for reasons unknown, using match=match in pytest.raises() fails, so we
    # fall back to old school ``match in str(err.value)``.
    with pytest.raises(ValueError) as err:
        Time({'year': 2019, 'wrong': 1}, format='ymdhms')
    assert match in str(err.value)

    match = "for 2 input key names you must supply 'year', 'month'"
    with pytest.raises(ValueError, match=match):
        Time({'year': 2019, 'minute': 1}, format='ymdhms')


def test_ymdhms_masked():
    tm = Time({'year': [2000, 2001]}, format='ymdhms')
    tm[0] = np.ma.masked
    assert isinstance(tm.value[0], np.ma.core.mvoid)
    for name in ymdhms_names:
        assert tm.value[0][name] is np.ma.masked


# Converted from doctest in astropy/test/formats.py for debugging
def test_ymdhms_output():
    t = Time({'year': 2015, 'month': 2, 'day': 3,
              'hour': 12, 'minute': 13, 'second': 14.567},
             scale='utc')
    # NOTE: actually comes back as np.void for some reason
    # NOTE: not necessarily a python int; might be an int32
    assert t.ymdhms.year == 2015


# There are two stages of validation now - one on input into a format, so that
# the format conversion code has tidy matched arrays to work with, and the
# other when object construction does not go through a format object. Or at
# least, the format object is constructed with "from_jd=True". In this case the
# normal input validation does not happen but the new input validation does,
# and can ensure that strange broadcasting anomalies can't happen.
# This form of construction uses from_jd=True.
def test_broadcasting_writeable():
    t = Time('J2015') + np.linspace(-1, 1, 10) * u.day
    t[2] = Time(58000, format="mjd")


def test_format_subformat_compatibility():
    """Test that changing format with out_subfmt defined is not a problem.
    See #9812, #9810."""
    t = Time('2019-12-20', out_subfmt='date_??')
    assert t.mjd == 58837.0
    assert t.yday == '2019:354:00:00'  # Preserves out_subfmt

    t2 = t.replicate(format='mjd')
    assert t2.out_subfmt == '*'  # Changes to default

    t2 = t.copy(format='mjd')
    assert t2.out_subfmt == '*'

    t2 = Time(t, format='mjd')
    assert t2.out_subfmt == '*'

    t2 = t.copy(format='yday')
    assert t2.out_subfmt == 'date_??'
    assert t2.value == '2019:354:00:00'

    t.format = 'yday'
    assert t.value == '2019:354:00:00'
    assert t.out_subfmt == 'date_??'

    t = Time('2019-12-20', out_subfmt='date')
    assert t.mjd == 58837.0
    assert t.yday == '2019:354'


@pytest.mark.parametrize('fmt_name,fmt_class', TIME_FORMATS.items())
def test_to_value_with_subfmt_for_every_format(fmt_name, fmt_class):
    """From a starting Time value, test that every valid combination of
    to_value(format, subfmt) works.  See #9812, #9361.
    """
    t = Time('2000-01-01')
    subfmts = list(subfmt[0] for subfmt in fmt_class.subfmts) + [None, '*']
    for subfmt in subfmts:
        t.to_value(fmt_name, subfmt)


@pytest.mark.parametrize('location', [None, (45, 45)])
def test_location_init(location):
    """Test fix in #9969 for issue #9962 where the location attribute is
    lost when initializing Time from an existing Time instance of list of
    Time instances.
    """
    tm = Time('J2010', location=location)

    # Init from a scalar Time
    tm2 = Time(tm)
    assert np.all(tm.location == tm2.location)
    assert type(tm.location) is type(tm2.location)  # noqa

    # From a list of Times
    tm2 = Time([tm, tm])
    if location is None:
        assert tm2.location is None
    else:
        for loc in tm2.location:
            assert loc == tm.location
    assert type(tm.location) is type(tm2.location)  # noqa

    # Effectively the same as a list of Times, but just to be sure that
    # Table mixin inititialization is working as expected.
    tm2 = Table([[tm, tm]])['col0']
    if location is None:
        assert tm2.location is None
    else:
        for loc in tm2.location:
            assert loc == tm.location
    assert type(tm.location) is type(tm2.location)  # noqa


def test_location_init_fail():
    """Test fix in #9969 for issue #9962 where the location attribute is
    lost when initializing Time from an existing Time instance of list of
    Time instances.  Make sure exception is correct.
    """
    tm = Time('J2010', location=(45, 45))
    tm2 = Time('J2010')

    with pytest.raises(ValueError,
                       match='cannot concatenate times unless all locations'):
        Time([tm, tm2])
