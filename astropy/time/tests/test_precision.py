import decimal
import warnings
import functools
import contextlib
from decimal import Decimal
from datetime import datetime, timedelta

import pytest
from hypothesis import assume, example, given, target
from hypothesis.extra.numpy import array_shapes, arrays
from hypothesis.strategies import (composite, datetimes, floats, integers,
                                   one_of, sampled_from, timedeltas, tuples)

import numpy as np
import erfa
from erfa import ErfaError, ErfaWarning

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import STANDARD_TIME_SCALES, Time, TimeDelta
from astropy.time.utils import day_frac, two_sum
from astropy.utils import iers


allclose_jd = functools.partial(np.allclose, rtol=np.finfo(float).eps, atol=0)
allclose_jd2 = functools.partial(np.allclose, rtol=np.finfo(float).eps,
                                 atol=np.finfo(float).eps)  # 20 ps atol
allclose_sec = functools.partial(np.allclose, rtol=np.finfo(float).eps,
                                 atol=np.finfo(float).eps * 24 * 3600)

tiny = np.finfo(float).eps
dt_tiny = TimeDelta(tiny, format='jd')


def setup_module():
    # Pre-load leap seconds table to avoid flakiness in hypothesis runs.
    # See https://github.com/astropy/astropy/issues/11030
    Time('2020-01-01').ut1


@pytest.fixture(scope='module')
def iers_b():
    """This is an expensive operation, so we share it between tests using a
    module-scoped fixture instead of using the context manager form.  This
    is particularly important for Hypothesis, which invokes the decorated
    test function many times (100 by default; see conftest.py for details).
    """
    with iers.earth_orientation_table.set(iers.IERS_B.open(iers.IERS_B_FILE)):
        yield "<using IERS-B orientation table>"


@contextlib.contextmanager
def quiet_erfa():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ErfaWarning)
        yield


def assert_almost_equal(a, b, *, rtol=None, atol=None, label=''):
    """Assert numbers are almost equal.

    This version also lets hypothesis know how far apart the inputs are, so
    that it can work towards a failure and present the worst failure ever seen
    as well as the simplest, which often just barely exceeds the threshold.
    """
    __tracebackhide__ = True
    if rtol is None or rtol == 0:
        thresh = atol
    elif atol is None:
        thresh = rtol * (abs(a) + abs(b)) / 2
    else:
        thresh = atol + rtol * (abs(a) + abs(b)) / 2

    amb = (a - b)
    if isinstance(amb, TimeDelta):
        ambv = amb.to_value(u.s)
        target(ambv, label=label + " (a-b).to_value(u.s), from TimeDelta")
        target(-ambv, label=label + " (b-a).to_value(u.s), from TimeDelta")
        if isinstance(thresh, u.Quantity):
            amb = amb.to(thresh.unit)
    else:
        try:
            target_value = float(amb)
        except TypeError:
            pass
        else:
            target(target_value, label=label + " float(a-b)")
            target(-target_value, label=label + " float(b-a)")

    assert abs(amb) < thresh


# Days that end with leap seconds
# Some time scales use a so-called "leap smear" to cope with these, others
# have times they can't represent or can represent two different ways.
# In any case these days are liable to cause trouble in time conversions.
# Note that from_erfa includes some weird non-integer steps before 1970.
leap_second_table = iers.LeapSeconds.from_iers_leap_seconds()
# Days that contain leap_seconds
leap_second_days = leap_second_table["mjd"] - 1
leap_second_deltas = list(zip(leap_second_days[1:],
                              np.diff(leap_second_table["tai_utc"])))

today = Time.now()
mjd0 = Time(0, format="mjd")


def reasonable_ordinary_jd():
    return tuples(floats(2440000, 2470000), floats(-0.5, 0.5))


@composite
def leap_second_tricky(draw):
    mjd = draw(one_of(sampled_from(leap_second_days),
                      sampled_from(leap_second_days + 1),
                      sampled_from(leap_second_days - 1)))
    return mjd + mjd0.jd1 + mjd0.jd2, draw(floats(0, 1))


def reasonable_jd():
    """Pick a reasonable JD.

    These should be not too far in the past or future (so that date conversion
    routines don't have to deal with anything too exotic), but they should
    include leap second days as a special case, and they should include several
    particularly simple cases (today, the beginning of the MJD scale, a
    reasonable date) so that hypothesis' example simplification produces
    obviously simple examples when they trigger problems.
    """
    moments = [(2455000., 0.), (mjd0.jd1, mjd0.jd2), (today.jd1, today.jd2)]
    return one_of(sampled_from(moments),
                  reasonable_ordinary_jd(),
                  leap_second_tricky())


def unreasonable_ordinary_jd():
    """JD pair that might be unordered or far away"""
    return tuples(floats(-1e7, 1e7), floats(-1e7, 1e7))


def ordered_jd():
    """JD pair that is ordered but not necessarily near now"""
    return tuples(floats(-1e7, 1e7), floats(-0.5, 0.5))


def unreasonable_jd():
    return one_of(reasonable_jd(), ordered_jd(), unreasonable_ordinary_jd())


@composite
def jd_arrays(draw, jd_values):
    s = draw(array_shapes())
    d = np.dtype([("jd1", float), ("jd2", float)])
    jdv = jd_values.map(lambda x: np.array(x, dtype=d))
    a = draw(arrays(d, s, elements=jdv))
    return a["jd1"], a["jd2"]


def unreasonable_delta():
    return tuples(floats(-1e7, 1e7), floats(-1e7, 1e7))


def reasonable_delta():
    return tuples(floats(-1e4, 1e4), floats(-0.5, 0.5))


# redundant?
def test_abs_jd2_always_less_than_half():
    """Make jd2 approach +/-0.5, and check that it doesn't go over."""
    t1 = Time(2400000.5, [-tiny, +tiny], format='jd')
    assert np.all(t1.jd1 % 1 == 0)
    assert np.all(abs(t1.jd2) < 0.5)
    t2 = Time(2400000., [[0.5 - tiny, 0.5 + tiny],
                         [-0.5 - tiny, -0.5 + tiny]], format='jd')
    assert np.all(t2.jd1 % 1 == 0)
    assert np.all(abs(t2.jd2) < 0.5)


@given(jd_arrays(unreasonable_jd()))
def test_abs_jd2_always_less_than_half_on_construction(jds):
    jd1, jd2 = jds
    t = Time(jd1, jd2, format="jd")
    target(np.amax(np.abs(t.jd2)))
    assert np.all(t.jd1 % 1 == 0)
    assert np.all(abs(t.jd2) <= 0.5)
    assert np.all((abs(t.jd2) < 0.5) | (t.jd1 % 2 == 0))


@given(integers(-10**8, 10**8), sampled_from([-0.5, 0.5]))
def test_round_to_even(jd1, jd2):
    t = Time(jd1, jd2, format="jd")
    assert (abs(t.jd2) == 0.5) and (t.jd1 % 2 == 0)


def test_addition():
    """Check that an addition at the limit of precision (2^-52) is seen"""
    t = Time(2455555., 0.5, format='jd', scale='utc')

    t_dt = t + dt_tiny
    assert t_dt.jd1 == t.jd1 and t_dt.jd2 != t.jd2

    # Check that the addition is exactly reversed by the corresponding
    # subtraction
    t2 = t_dt - dt_tiny
    assert t2.jd1 == t.jd1 and t2.jd2 == t.jd2


def test_mult_div():
    """Test precision with multiply and divide"""
    dt_small = 6 * dt_tiny
    # pick a number that will leave remainder if divided by 6.
    dt_big = TimeDelta(20000., format='jd')
    dt_big_small_by_6 = (dt_big + dt_small) / 6.
    dt_frac = dt_big_small_by_6 - TimeDelta(3333., format='jd')
    assert allclose_jd2(dt_frac.jd2, 0.33333333333333354)


def test_init_variations():
    """Check that 3 ways of specifying a time + small offset are equivalent"""
    dt_tiny_sec = dt_tiny.jd2 * 86400.
    t1 = Time(1e11, format='cxcsec') + dt_tiny
    t2 = Time(1e11, dt_tiny_sec, format='cxcsec')
    t3 = Time(dt_tiny_sec, 1e11, format='cxcsec')
    assert t1.jd1 == t2.jd1
    assert t1.jd2 == t3.jd2
    assert t1.jd1 == t2.jd1
    assert t1.jd2 == t3.jd2


def test_precision_exceeds_64bit():
    """
    Check that Time object really holds more precision than float64 by looking
    at the (naively) summed 64-bit result and asserting equality at the
    bit level.
    """
    t1 = Time(1.23456789e11, format='cxcsec')
    t2 = t1 + dt_tiny
    assert t1.jd == t2.jd


def test_through_scale_change():
    """Check that precision holds through scale change (cxcsec is TT)"""
    t0 = Time(1.0, format='cxcsec')
    t1 = Time(1.23456789e11, format='cxcsec')
    dt_tt = t1 - t0
    dt_tai = t1.tai - t0.tai
    assert allclose_jd(dt_tt.jd1, dt_tai.jd1)
    assert allclose_jd2(dt_tt.jd2, dt_tai.jd2)


def test_iso_init():
    """Check when initializing from ISO date"""
    t1 = Time('2000:001:00:00:00.00000001', scale='tai')
    t2 = Time('3000:001:13:00:00.00000002', scale='tai')
    dt = t2 - t1
    assert allclose_jd2(dt.jd2, 13. / 24. + 1e-8 / 86400. - 1.0)


def test_jd1_is_mult_of_one():
    """
    Check that jd1 is a multiple of 1.
    """
    t1 = Time('2000:001:00:00:00.00000001', scale='tai')
    assert np.round(t1.jd1) == t1.jd1
    t1 = Time(1.23456789, 12345678.90123456, format='jd', scale='tai')
    assert np.round(t1.jd1) == t1.jd1


def test_precision_neg():
    """
    Check precision when jd1 is negative.  This used to fail because ERFA
    routines use a test like jd1 > jd2 to decide which component to update.
    It was updated to abs(jd1) > abs(jd2) in erfa 1.6 (sofa 20190722).
    """
    t1 = Time(-100000.123456, format='jd', scale='tt')
    assert np.round(t1.jd1) == t1.jd1
    t1_tai = t1.tai
    assert np.round(t1_tai.jd1) == t1_tai.jd1


def test_precision_epoch():
    """
    Check that input via epoch also has full precision, i.e., against
    regression on https://github.com/astropy/astropy/pull/366
    """
    t_utc = Time(range(1980, 2001), format='jyear', scale='utc')
    t_tai = Time(range(1980, 2001), format='jyear', scale='tai')
    dt = t_utc - t_tai
    assert allclose_sec(dt.sec, np.round(dt.sec))


def test_leap_seconds_rounded_correctly():
    """Regression tests against #2083, where a leap second was rounded
    incorrectly by the underlying ERFA routine."""
    with iers.conf.set_temp('auto_download', False):
        t = Time(['2012-06-30 23:59:59.413',
                  '2012-07-01 00:00:00.413'], scale='ut1', precision=3).utc
        assert np.all(t.iso == np.array(['2012-06-30 23:59:60.000',
                                         '2012-07-01 00:00:00.000']))
    # with the bug, both yielded '2012-06-30 23:59:60.000'


@given(integers(-2**52+2, 2**52-2), floats(-1, 1))
@example(i=65536, f=3.637978807091714e-12)
def test_two_sum(i, f):
    with decimal.localcontext(decimal.Context(prec=40)):
        a = Decimal(i) + Decimal(f)
        s, r = two_sum(i, f)
        b = Decimal(s) + Decimal(r)
        assert_almost_equal(a, b, atol=Decimal(tiny), rtol=Decimal(0))


@given(floats(), floats())
def test_two_sum_symmetric(f1, f2):
    np.testing.assert_equal(two_sum(f1, f2), two_sum(f2, f1))


@given(floats(allow_nan=False, allow_infinity=False),
       floats(allow_nan=False, allow_infinity=False))
@example(f1=8.988465674311579e+307, f2=8.98846567431158e+307)
@example(f1=8.988465674311579e+307, f2=-8.98846567431158e+307)
@example(f1=-8.988465674311579e+307, f2=-8.98846567431158e+307)
def test_two_sum_size(f1, f2):
    r1, r2 = two_sum(f1, f2)
    assert (abs(r1) > abs(r2) / np.finfo(float).eps
            or r1 == r2 == 0
            or not np.isfinite(f1 + f2))


@given(integers(-2**52+2, 2**52-2), floats(-1, 1))
@example(i=65536, f=3.637978807091714e-12)
def test_day_frac_harmless(i, f):
    with decimal.localcontext(decimal.Context(prec=40)):
        a = Decimal(i) + Decimal(f)
        i_d, f_d = day_frac(i, f)
        a_d = Decimal(i_d) + Decimal(f_d)
        assert_almost_equal(a, a_d, atol=Decimal(tiny), rtol=Decimal(0))


@given(integers(-2**52+2, 2**52-2), floats(-0.5, 0.5))
@example(i=65536, f=3.637978807091714e-12)
def test_day_frac_exact(i, f):
    assume(abs(f) < 0.5 or i % 2 == 0)
    i_d, f_d = day_frac(i, f)
    assert i == i_d
    assert f == f_d


@given(integers(-2**52+2, 2**52-2), floats(-1, 1))
@example(i=65536, f=3.637978807091714e-12)
def test_day_frac_idempotent(i, f):
    i_d, f_d = day_frac(i, f)
    assert (i_d, f_d) == day_frac(i_d, f_d)


@given(integers(-2**52+2, 2**52-2), floats(-1, 1))
@example(i=65536, f=3.637978807091714e-12)
def test_mjd_initialization_precise(i, f):
    t = Time(val=i, val2=f, format="mjd", scale="tai")
    jd1, jd2 = day_frac(i + erfa.DJM0, f)
    jd1_t, jd2_t = day_frac(t.jd1, t.jd2)
    assert (abs((jd1 - jd1_t) + (jd2 - jd2_t)) * u.day).to(u.ns) < 1 * u.ns


@given(jd_arrays(unreasonable_jd()))
def test_day_frac_always_less_than_half(jds):
    jd1, jd2 = jds
    t_jd1, t_jd2 = day_frac(jd1, jd2)
    assert np.all(t_jd1 % 1 == 0)
    assert np.all(abs(t_jd2) <= 0.5)
    assert np.all((abs(t_jd2) < 0.5) | (t_jd1 % 2 == 0))


@given(integers(-10**8, 10**8), sampled_from([-0.5, 0.5]))
def test_day_frac_round_to_even(jd1, jd2):
    t_jd1, t_jd2 = day_frac(jd1, jd2)
    assert (abs(t_jd2) == 0.5) and (t_jd1 % 2 == 0)


@given(scale=sampled_from(STANDARD_TIME_SCALES), jds=unreasonable_jd())
@example(scale="tai", jds=(0.0, 0.0))
@example(scale="tai", jds=(0.0, -31738.500000000346))
def test_resolution_never_decreases(scale, jds):
    jd1, jd2 = jds
    assume(not scale == 'utc' or 2440000 < jd1 + jd2 < 2460000)
    t = Time(jd1, jd2, format="jd", scale=scale)
    with quiet_erfa():
        assert t != t + dt_tiny


@given(reasonable_jd())
def test_resolution_never_decreases_utc(jds):
    """UTC is very unhappy with unreasonable times"""
    jd1, jd2 = jds
    t = Time(jd1, jd2, format="jd", scale="utc")
    with quiet_erfa():
        assert t != t + dt_tiny


@given(scale1=sampled_from(STANDARD_TIME_SCALES),
       scale2=sampled_from(STANDARD_TIME_SCALES),
       jds=unreasonable_jd())
@example(scale1='tcg', scale2='ut1', jds=(2445149.5, 0.47187700984387526))
@example(scale1='tai', scale2='tcb', jds=(2441316.5, 0.0))
@example(scale1='tai', scale2='tcb', jds=(0.0, 0.0))
def test_conversion_preserves_jd1_jd2_invariant(iers_b, scale1, scale2, jds):
    jd1, jd2 = jds
    t = Time(jd1, jd2, scale=scale1, format="jd")
    try:
        with quiet_erfa():
            t2 = getattr(t, scale2)
    except iers.IERSRangeError:  # UT1 conversion needs IERS data
        assume(False)
    except ErfaError:
        assume(False)
    assert t2.jd1 % 1 == 0
    assert abs(t2.jd2) <= 0.5
    assert abs(t2.jd2) < 0.5 or t2.jd1 % 2 == 0


@given(scale1=sampled_from(STANDARD_TIME_SCALES),
       scale2=sampled_from(STANDARD_TIME_SCALES),
       jds=unreasonable_jd())
@example(scale1='tai', scale2='utc', jds=(0.0, 0.0))
def test_conversion_never_loses_precision(iers_b, scale1, scale2, jds):
    """Check that time ordering remains if we convert to another scale.

    Here, since scale differences can involve multiplication, we allow
    for losing one ULP, i.e., we test that two times that differ by
    two ULP will keep the same order if changed to another scale.
    """
    jd1, jd2 = jds
    t = Time(jd1, jd2, scale=scale1, format="jd")
    # Near-zero UTC JDs degrade accuracy; not clear why,
    # but also not so relevant, so ignoring.
    if (scale1 == 'utc' or scale2 == 'utc') and abs(jd1+jd2) < 1:
        tiny = 100*u.us
    else:
        tiny = 2*dt_tiny
    try:
        with quiet_erfa():
            t2 = t + tiny
            assert getattr(t, scale2) < getattr(t2, scale2)
    except iers.IERSRangeError:  # UT1 conversion needs IERS data
        assume(scale1 != 'ut1' or 2440000 < jd1 + jd2 < 2458000)
        assume(scale2 != 'ut1' or 2440000 < jd1 + jd2 < 2458000)
        raise
    except ErfaError:
        # If the generated date is too early to compute a UTC julian date,
        # and we're not converting between scales which are known to be safe,
        # tell Hypothesis that this example is invalid and to try another.
        # See https://docs.astropy.org/en/latest/time/index.html#time-scale
        barycentric = {scale1, scale2}.issubset({'tcb', 'tdb'})
        geocentric = {scale1, scale2}.issubset({'tai', 'tt', 'tcg'})
        assume(jd1 + jd2 >= -31738.5 or geocentric or barycentric)
        raise


@given(sampled_from(leap_second_deltas), floats(0.1, 0.9))
def test_leap_stretch_mjd(d, f):
    mjd, delta = d
    t0 = Time(mjd, format="mjd", scale="utc")
    th = Time(mjd + f, format="mjd", scale="utc")
    t1 = Time(mjd + 1, format="mjd", scale="utc")
    assert_quantity_allclose((t1 - t0).to(u.s), (1 * u.day + delta * u.s))
    assert_quantity_allclose((th - t0).to(u.s), f * (1 * u.day + delta * u.s))
    assert_quantity_allclose((t1 - th).to(u.s), (1 - f) * (1 * u.day + delta * u.s))


@given(scale=sampled_from(STANDARD_TIME_SCALES),
       jds=unreasonable_jd(),
       delta=floats(-10000, 10000))
@example(scale='utc',
         jds=(0.0, 2.2204460492503136e-13),
         delta=6.661338147750941e-13)
@example(scale='utc',
         jds=(2441682.5, 2.2204460492503136e-16),
         delta=7.327471962526035e-12)
@example(scale='utc', jds=(0.0, 5.787592627370942e-13), delta=0.0)
def test_jd_add_subtract_round_trip(scale, jds, delta):
    jd1, jd2 = jds
    if scale == 'utc' and abs(jd1+jd2) < 1:
        # Near-zero UTC JDs degrade accuracy; not clear why,
        # but also not so relevant, so ignoring.
        thresh = 100*u.us
    else:
        thresh = 2*dt_tiny
    t = Time(jd1, jd2, scale=scale, format="jd")
    try:
        with quiet_erfa():
            t2 = t + delta*u.day
            if abs(delta) >= np.finfo(float).eps:
                assert t2 != t
            t3 = t2 - delta*u.day
            assert_almost_equal(t3, t, atol=thresh, rtol=0)
    except ErfaError:
        assume(scale != 'utc' or 2440000 < jd1+jd2 < 2460000)
        raise


@given(scale=sampled_from(STANDARD_TIME_SCALES),
       jds=reasonable_jd(),
       delta=floats(-3*tiny, 3*tiny))
@example(scale='tai', jds=(0.0, 3.5762786865234384), delta=2.220446049250313e-16)
def test_time_argminmaxsort(scale, jds, delta):
    jd1, jd2 = jds
    t = Time(jd1, jd2+np.array([0, delta]), scale=scale, format="jd")
    imin = t.argmin()
    imax = t.argmax()
    isort = t.argsort()
    diff = (t.jd1[1]-t.jd1[0]) + (t.jd2[1]-t.jd2[0])
    if diff < 0:  # item 1 smaller
        assert delta < 0
        assert imin == 1 and imax == 0 and np.all(isort == [1, 0])
    elif diff == 0:  # identical within precision
        assert abs(delta) <= tiny
        assert imin == 0 and imax == 0 and np.all(isort == [0, 1])
    else:
        assert delta > 0
        assert imin == 0 and imax == 1 and np.all(isort == [0, 1])


@given(sampled_from(STANDARD_TIME_SCALES), unreasonable_jd(), unreasonable_jd())
@example(scale='utc',
         jds_a=(2455000.0, 0.0),
         jds_b=(2443144.5, 0.5000462962962965))
@example(scale='utc',
         jds_a=(2459003.0, 0.267502885949074),
         jds_b=(2454657.001045462, 0.49895453779026877))
def test_timedelta_full_precision(scale, jds_a, jds_b):
    jd1_a, jd2_a = jds_a
    jd1_b, jd2_b = jds_b
    assume(scale != 'utc'
           or (2440000 < jd1_a+jd2_a < 2460000
               and 2440000 < jd1_b+jd2_b < 2460000))
    if scale == 'utc':
        # UTC subtraction implies a scale change, so possible rounding errors.
        tiny = 2 * dt_tiny
    else:
        tiny = dt_tiny

    t_a = Time(jd1_a, jd2_a, scale=scale, format="jd")
    t_b = Time(jd1_b, jd2_b, scale=scale, format="jd")
    dt = t_b - t_a
    assert dt != (t_b + tiny) - t_a
    with quiet_erfa():
        assert_almost_equal(t_b-dt/2, t_a+dt/2, atol=2*dt_tiny, rtol=0,
                            label="midpoint")
        assert_almost_equal(t_b+dt, t_a+2*dt, atol=2*dt_tiny, rtol=0, label="up")
        assert_almost_equal(t_b-2*dt, t_a-dt, atol=2*dt_tiny, rtol=0, label="down")


@given(scale=sampled_from(STANDARD_TIME_SCALES),
       jds_a=unreasonable_jd(),
       jds_b=unreasonable_jd(),
       x=integers(1, 100),
       y=integers(1, 100))
def test_timedelta_full_precision_arithmetic(scale, jds_a, jds_b, x, y):
    jd1_a, jd2_a = jds_a
    jd1_b, jd2_b = jds_b
    t_a = Time(jd1_a, jd2_a, scale=scale, format="jd")
    t_b = Time(jd1_b, jd2_b, scale=scale, format="jd")
    with quiet_erfa():
        try:
            dt = t_b - t_a
            dt_x = x*dt/(x+y)
            dt_y = y*dt/(x+y)
            assert_almost_equal(dt_x + dt_y, dt, atol=(x+y)*dt_tiny, rtol=0)
        except ErfaError:
            assume(scale != 'utc'
                   or (2440000 < jd1_a+jd2_a < 2460000
                       and 2440000 < jd1_b+jd2_b < 2460000))
            raise


@given(scale1=sampled_from(STANDARD_TIME_SCALES),
       scale2=sampled_from(STANDARD_TIME_SCALES),
       jds_a=reasonable_jd(),
       jds_b=reasonable_jd())
def test_timedelta_conversion(scale1, scale2, jds_a, jds_b):
    jd1_a, jd2_a = jds_a
    jd1_b, jd2_b = jds_b
    # not translation invariant so can't convert TimeDelta
    assume('utc' not in [scale1, scale2])
    # Conversions a problem but within UT1 it should work
    assume(('ut1' not in [scale1, scale2]) or scale1 == scale2)
    t_a = Time(jd1_a, jd2_a, scale=scale1, format="jd")
    t_b = Time(jd1_b, jd2_b, scale=scale2, format="jd")
    with quiet_erfa():
        dt = t_b - t_a
        t_a_2 = getattr(t_a, scale2)
        t_b_2 = getattr(t_b, scale2)
        dt_2 = getattr(dt, scale2)
        assert_almost_equal(t_b_2 - t_a_2, dt_2, atol=dt_tiny, rtol=0,
                            label="converted")
        # Implicit conversion
        assert_almost_equal(t_b_2 - t_a_2, dt, atol=dt_tiny, rtol=0,
                            label="not converted")


# UTC disagrees when there are leap seconds
_utc_bad = [(pytest.param(s, marks=pytest.mark.xfail) if s == 'utc' else s)
            for s in STANDARD_TIME_SCALES]


@given(datetimes(), datetimes())  # datetimes have microsecond resolution
@example(dt1=datetime(1235, 1, 1, 0, 0),
         dt2=datetime(9950, 1, 1, 0, 0, 0, 890773))
@pytest.mark.parametrize("scale", _utc_bad)
def test_datetime_difference_agrees_with_timedelta(scale, dt1, dt2):
    t1 = Time(dt1, scale=scale)
    t2 = Time(dt2, scale=scale)
    assert_almost_equal(t2-t1,
                        TimeDelta(dt2-dt1,
                                  scale=None if scale == 'utc' else scale),
                        atol=2*u.us)


@given(days=integers(-3000*365, 3000*365),
       microseconds=integers(0, 24*60*60*1000000))
@pytest.mark.parametrize("scale", _utc_bad)
def test_datetime_to_timedelta(scale, days, microseconds):
    td = timedelta(days=days, microseconds=microseconds)
    assert (TimeDelta(td, scale=scale)
            == TimeDelta(days, microseconds/(86400*1e6), scale=scale, format="jd"))


@given(days=integers(-3000*365, 3000*365),
       microseconds=integers(0, 24*60*60*1000000))
@pytest.mark.parametrize("scale", _utc_bad)
def test_datetime_timedelta_roundtrip(scale, days, microseconds):
    td = timedelta(days=days, microseconds=microseconds)
    assert td == TimeDelta(td, scale=scale).value


@given(days=integers(-3000*365, 3000*365), day_frac=floats(0, 1))
@example(days=262144, day_frac=2.314815006343452e-11)
@example(days=1048576, day_frac=1.157407503171726e-10)
@pytest.mark.parametrize("scale", _utc_bad)
def test_timedelta_datetime_roundtrip(scale, days, day_frac):
    td = TimeDelta(days, day_frac, format="jd", scale=scale)
    td.format = "datetime"
    assert_almost_equal(td, TimeDelta(td.value, scale=scale), atol=2*u.us)


@given(integers(-3000*365, 3000*365), floats(0, 1))
@example(days=262144, day_frac=2.314815006343452e-11)
@pytest.mark.parametrize("scale", _utc_bad)
def test_timedelta_from_parts(scale, days, day_frac):
    kwargs = dict(format="jd", scale=scale)
    whole = TimeDelta(days, day_frac, **kwargs)
    from_parts = TimeDelta(days, **kwargs) + TimeDelta(day_frac, **kwargs)
    assert whole == from_parts


def test_datetime_difference_agrees_with_timedelta_no_hypothesis():
    scale = "tai"
    dt1 = datetime(1235, 1, 1, 0, 0)
    dt2 = datetime(9950, 1, 1, 0, 0, 0, 890773)
    t1 = Time(dt1, scale=scale)
    t2 = Time(dt2, scale=scale)
    assert(abs((t2-t1) - TimeDelta(dt2-dt1, scale=scale)) < 1*u.us)


# datetimes have microsecond resolution
@given(datetimes(), timedeltas())
@example(dt=datetime(2000, 1, 1, 0, 0),
         td=timedelta(days=-397683, microseconds=2))
@example(dt=datetime(2179, 1, 1, 0, 0),
         td=timedelta(days=-795365, microseconds=53))
@example(dt=datetime(2000, 1, 1, 0, 0),
         td=timedelta(days=1590729, microseconds=10))
@example(dt=datetime(4357, 1, 1, 0, 0),
         td=timedelta(days=-1590729, microseconds=107770))
@example(dt=datetime(4357, 1, 1, 0, 0, 0, 29),
         td=timedelta(days=-1590729, microseconds=746292))
@pytest.mark.parametrize("scale", _utc_bad)
def test_datetime_timedelta_sum(scale, dt, td):
    try:
        dt + td
    except OverflowError:
        assume(False)
    dt_a = Time(dt, scale=scale)
    td_a = TimeDelta(td, scale=None if scale == 'utc' else scale)
    assert_almost_equal(dt_a+td_a, Time(dt+td, scale=scale), atol=2*u.us)


@given(jds=reasonable_jd(),
       lat1=floats(-90, 90),
       lat2=floats(-90, 90),
       lon=floats(-180, 180))
@pytest.mark.parametrize("kind", ["apparent", "mean"])
def test_sidereal_lat_independent(iers_b, kind, jds, lat1, lat2, lon):
    jd1, jd2 = jds
    t1 = Time(jd1, jd2, scale="ut1", format="jd", location=(lon, lat1))
    t2 = Time(jd1, jd2, scale="ut1", format="jd", location=(lon, lat2))
    try:
        assert_almost_equal(t1.sidereal_time(kind),
                            t2.sidereal_time(kind),
                            atol=1*u.uas)
    except iers.IERSRangeError:
        assume(False)


@given(jds=reasonable_jd(),
       lat=floats(-90, 90),
       lon=floats(-180, 180),
       lon_delta=floats(-360, 360))
@pytest.mark.parametrize("kind", ["apparent", "mean"])
def test_sidereal_lon_independent(iers_b, kind, jds, lat, lon, lon_delta):
    jd1, jd2 = jds
    t1 = Time(jd1, jd2, scale="ut1", format="jd", location=(lon, lat))
    t2 = Time(jd1, jd2, scale="ut1", format="jd", location=(lon+lon_delta, lat))
    try:
        diff = t1.sidereal_time(kind) + lon_delta*u.degree - t2.sidereal_time(kind)
    except iers.IERSRangeError:
        assume(False)
    else:
        expected_degrees = (diff.to_value(u.degree) + 180) % 360
        assert_almost_equal(expected_degrees, 180, atol=1/(60*60*1000))
