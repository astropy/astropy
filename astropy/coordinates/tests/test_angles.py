# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Test initalization of angles not already covered by the API tests

import numpy as np
from ..angles import Angle, RA, Dec, BoundsError, Latitude, Longitude
from ...tests.helper import pytest
from ...tests.compat import assert_allclose
from ... import units as u


def test_negative_zero_dms():
    # Test for DMS parser
    a = Angle('-00:00:10', u.deg)
    assert_allclose(a.degree, -10. / 3600.)


def test_negative_zero_dm():
    # Test for DM parser
    a = Angle('-00:10', u.deg)
    assert_allclose(a.degree, -10. / 60.)


def test_negative_zero_hms():
    # Test for HMS parser
    a = Angle('-00:00:10', u.hour)
    assert_allclose(a.hour, -10. / 3600.)


def test_negative_zero_hm():
    # Test for HM parser
    a = Angle('-00:10', u.hour)
    assert_allclose(a.hour, -10. / 60.)


def test_negative_sixty_hm():
    # Test for HM parser
    a = Angle('-00:60', u.hour)
    assert_allclose(a.hour, -1.)


def test_plus_sixty_hm():
    # Test for HM parser
    a = Angle('00:60', u.hour)
    assert_allclose(a.hour, 1.)


def test_negative_fifty_nine_sixty_dms():
    # Test for DMS parser
    a = Angle('-00:59:60', u.deg)
    assert_allclose(a.degree, -1.)


def test_plus_fifty_nine_sixty_dms():
    # Test for DMS parser
    a = Angle('+00:59:60', u.deg)
    assert_allclose(a.degree, 1.)


def test_negative_sixty_dms():
    # Test for DMS parser
    a = Angle('-00:00:60', u.deg)
    assert_allclose(a.degree, -1. / 60.)


def test_plus_sixty_dms():
    # Test for DMS parser
    a = Angle('+00:00:60', u.deg)
    assert_allclose(a.degree, 1. / 60.)


def test_angle_to_is_angle():
    a = Angle('00:00:60', u.deg)
    assert isinstance(a, Angle)
    assert isinstance(a.to(u.rad), Angle)


def test_angle_to_quantity():
    a = Angle('00:00:60', u.deg)
    q = u.Quantity(a)
    assert isinstance(q, u.Quantity)
    assert q.unit is u.deg


def test_angle_string():
    a = Angle('00:00:60', u.deg)
    assert str(a) == '0d01m00.00000s'
    a = Angle('-00:00:10', u.hour)
    assert str(a) == '-0h00m10.00000s'
    a = Angle(3.2, u.radian)
    assert str(a) == '3.20000rad'
    a = Angle(4.2, u.microarcsecond)
    assert str(a) == '4.20000uarcsec'
    a = Angle('1.0uarcsec')
    assert a.value == 1.0
    assert a.unit == u.microarcsecond
    a = Angle("3d")
    assert_allclose(a.value, 3.0)
    assert a.unit == u.degree


def test_angle_repr():
    assert 'Angle' in repr(Angle(0, u.deg))
    assert 'RA' in repr(RA(0, u.deg))
    assert 'Dec' in repr(Dec(0, u.deg))

    a = Angle(0, u.deg)
    repr(a)


def test_wrap_at_in_place():
    a = Angle([-20, 150, 350, 360] * u.deg)
    out = a.wrap_at('180d', in_place=True)
    assert out is None
    assert np.all(a.degree == np.array([-20., 150., -10., 0.]))


def test_latitude():
    with pytest.raises(ValueError):
        lat = Latitude(['91d', '89d'])

    lat = Latitude(['90d', '89d'])
    with pytest.raises(ValueError):
        lat[0] = 90.001 * u.deg
    with pytest.raises(ValueError):
        lat[0] = -90.001 * u.deg

    lat = Latitude(['90d', '89d'])
    assert lat[0] == 90 * u.deg
    assert lat[1] == 89 * u.deg
    assert np.all(lat == Angle(['90d', '89d']))

    angle = lat - 190 * u.deg
    assert type(angle) is Angle
    assert angle[0] == -100 * u.deg

    lat = Latitude('80d')
    angle = lat / 2.
    assert type(angle) is Angle
    assert angle == 40 * u.deg

    angle = lat * 2.
    assert type(angle) is Angle
    assert angle == 160 * u.deg

    angle = -lat
    assert type(angle) is Angle
    assert angle == -80 * u.deg


def test_wrap_at():
    a = Angle([-20, 150, 350, 360] * u.deg)
    assert np.all(a.wrap_at(360 * u.deg).degree == np.array([340., 150., 350., 0.]))
    assert np.all(a.wrap_at(Angle(360, unit=u.deg)).degree == np.array([340., 150., 350., 0.]))
    assert np.all(a.wrap_at('360d').degree == np.array([340., 150., 350., 0.]))
    assert np.all(a.wrap_at('180d').degree == np.array([-20., 150., -10., 0.]))
    assert np.all(a.wrap_at(np.pi * u.rad).degree == np.array([-20., 150., -10., 0.]))

    a = Angle(np.arange(-1000.0, 1000.0, 0.125), unit=u.deg)
    for wrap_angle in (270, 0.2, 0.0, 360.0, 500, -2000.125):
        aw = a.wrap_at(wrap_angle * u.deg)
        assert np.all(aw.degree >= wrap_angle - 360.0)
        assert np.all(aw.degree < wrap_angle)

        aw = a.to(u.rad).wrap_at(wrap_angle * u.deg)
        assert np.all(aw.degree >= wrap_angle - 360.0)
        assert np.all(aw.degree < wrap_angle)


def test_within_bounds():
    a = Angle([-20, 150, 350] * u.deg)
    assert a.within_bounds('0d', '360d') is False
    assert a.within_bounds(None, '360d') is True
    assert a.within_bounds(-30 * u.deg, None) is True

    a = Angle('-20d')
    assert a.within_bounds('0d', '360d') is False
    assert a.within_bounds(None, '360d') is True
    assert a.within_bounds(-30 * u.deg, None) is True
