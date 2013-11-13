# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""Test initalization of angles not already covered by the API tests"""

import numpy as np
from numpy.testing.utils import assert_allclose

from ..angles import Longitude, Latitude, Angle
from ...tests.helper import pytest
from ... import units as u


def test_negative_zero_dms():
    # Test for DMS parser
    a = Angle('-00:00:10', u.deg)
    assert_allclose(a.degree, -10. / 3600.)

    # Unicode minus
    a = Angle('\u221200:00:10', u.deg)
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


def test_quantity_to_angle():
    a = Angle(1.0*u.deg)
    assert isinstance(a, Angle)
    with pytest.raises(u.UnitsError):
        Angle(1.0*u.meter)
    a = Angle(1.0*u.hour)
    assert isinstance(a, Angle)
    assert a.unit is u.hourangle
    with pytest.raises(u.UnitsError):
        Angle(1.0*u.min)


def test_angle_string():
    a = Angle('00:00:60', u.deg)
    assert str(a) == '0d01m00s'
    a = Angle('-00:00:10', u.hour)
    assert str(a) == '-0h00m10s'
    a = Angle(3.2, u.radian)
    assert str(a) == '3.2rad'
    a = Angle(4.2, u.microarcsecond)
    assert str(a) == '4.2uarcsec'
    a = Angle('1.0uarcsec')
    assert a.value == 1.0
    assert a.unit == u.microarcsecond
    a = Angle("3d")
    assert_allclose(a.value, 3.0)
    assert a.unit == u.degree
    a = Angle('10"')
    assert_allclose(a.value, 10.0)
    assert a.unit == u.arcsecond
    a = Angle("10'")
    assert_allclose(a.value, 10.0)
    assert a.unit == u.arcminute


def test_angle_repr():
    assert 'Angle' in repr(Angle(0, u.deg))
    assert 'Longitude' in repr(Longitude(0, u.deg))
    assert 'Latitude' in repr(Latitude(0, u.deg))

    a = Angle(0, u.deg)
    repr(a)


def test_large_angle_representation():
    """Test that angles above 360 degrees can be output as strings,
    in repr, str, and to_string.  (regression test for #1413)"""
    a = Angle(350, u.deg) + Angle(350, u.deg)
    a.to_string()
    a.to_string(u.hourangle)
    repr(a)
    repr(a.to(u.hourangle))
    str(a)
    str(a.to(u.hourangle))


def test_wrap_at_inplace():
    a = Angle([-20, 150, 350, 360] * u.deg)
    out = a.wrap_at('180d', inplace=True)
    assert out is None
    assert np.all(a.degree == np.array([-20., 150., -10., 0.]))


def test_latitude():
    with pytest.raises(ValueError):
        lat = Latitude(['91d', '89d'])
    with pytest.raises(ValueError):
        lat = Latitude('-91d')

    lat = Latitude(['90d', '89d'])
    with pytest.raises(ValueError):
        lat[0] = 90.001 * u.deg
    with pytest.raises(ValueError):
        lat[0] = -90.001 * u.deg

    lat = Latitude(['90d', '89d'])
    assert lat[0] == 90 * u.deg
    assert lat[1] == 89 * u.deg
    assert np.all(lat == Angle(['90d', '89d']))

    # conserve type on unit change (closes #1423)
    angle = lat.to('radian')
    assert type(angle) is Latitude
    # but not on calculations
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


def test_longitude():
    # Default wrapping at 360d with an array input
    lon = Longitude(['370d', '88d'])
    assert np.all(lon == Longitude(['10d', '88d']))
    assert np.all(lon == Angle(['10d', '88d']))

    # conserve type on unit change and keep wrap_angle (closes #1423)
    angle = lon.to('hourangle')
    assert type(angle) is Longitude
    assert angle.wrap_angle == lon.wrap_angle
    angle = lon[0]
    assert type(angle) is Longitude
    assert angle.wrap_angle == lon.wrap_angle
    angle = lon[1:]
    assert type(angle) is Longitude
    assert angle.wrap_angle == lon.wrap_angle

    # but not on calculations
    angle = lon / 2.
    assert np.all(angle == Angle(['5d', '44d']))
    assert type(angle) is Angle
    assert not hasattr(angle, 'wrap_angle')

    angle = lon * 2. + 400 * u.deg
    assert np.all(angle == Angle(['420d', '576d']))
    assert type(angle) is Angle

    # Test setting a mutable value and having it wrap
    lon[1] = -10 * u.deg
    assert np.all(lon == Angle(['10d', '350d']))

    # Test wrapping and try hitting some edge cases
    lon = Longitude(np.array([0, 0.5, 1.0, 1.5, 2.0]) * np.pi, unit=u.radian)
    assert np.all(lon.degree == np.array([0., 90, 180, 270, 0]))

    lon = Longitude(np.array([0, 0.5, 1.0, 1.5, 2.0]) * np.pi, unit=u.radian, wrap_angle='180d')
    assert np.all(lon.degree == np.array([0., 90, -180, -90, 0]))

    # Wrap on setting wrap_angle property (also test auto-conversion of wrap_angle to an Angle)
    lon = Longitude(np.array([0, 0.5, 1.0, 1.5, 2.0]) * np.pi, unit=u.radian)
    lon.wrap_angle = '180d'
    assert np.all(lon.degree == np.array([0., 90, -180, -90, 0]))

    lon = Longitude('460d')
    assert lon == Angle('100d')
    lon.wrap_angle = '90d'
    assert lon == Angle('-260d')


def test_wrap_at():
    a = Angle([-20, 150, 350, 360] * u.deg)
    assert np.all(a.wrap_at(360 * u.deg).degree == np.array([340., 150., 350., 0.]))
    assert np.all(a.wrap_at(Angle(360, unit=u.deg)).degree == np.array([340., 150., 350., 0.]))
    assert np.all(a.wrap_at('360d').degree == np.array([340., 150., 350., 0.]))
    assert np.all(a.wrap_at('180d').degree == np.array([-20., 150., -10., 0.]))
    assert np.all(a.wrap_at(np.pi * u.rad).degree == np.array([-20., 150., -10., 0.]))

    # Test wrapping a scalar Angle
    a = Angle('190d')
    assert a.wrap_at('180d') == Angle('-170d')

    a = Angle(np.arange(-1000.0, 1000.0, 0.125), unit=u.deg)
    for wrap_angle in (270, 0.2, 0.0, 360.0, 500, -2000.125):
        aw = a.wrap_at(wrap_angle * u.deg)
        assert np.all(aw.degree >= wrap_angle - 360.0)
        assert np.all(aw.degree < wrap_angle)

        aw = a.to(u.rad).wrap_at(wrap_angle * u.deg)
        assert np.all(aw.degree >= wrap_angle - 360.0)
        assert np.all(aw.degree < wrap_angle)


def test_is_within_bounds():
    a = Angle([-20, 150, 350] * u.deg)
    assert a.is_within_bounds('0d', '360d') is False
    assert a.is_within_bounds(None, '360d') is True
    assert a.is_within_bounds(-30 * u.deg, None) is True

    a = Angle('-20d')
    assert a.is_within_bounds('0d', '360d') is False
    assert a.is_within_bounds(None, '360d') is True
    assert a.is_within_bounds(-30 * u.deg, None) is True


def test_angle_mismatched_unit():
    a = Angle('+6h7m8s', unit=u.degree)
    assert_allclose(a.value, 91.78333333333332)


def test_regression_formatting_negative():
    # Regression test for a bug that caused:
    #
    # >>> Angle(-1., unit='deg').to_string()
    # '-1d00m-0s'
    assert Angle(-0., unit='deg').to_string() == '-0d00m00s'
    assert Angle(-1., unit='deg').to_string() == '-1d00m00s'
    assert Angle(-0., unit='hour').to_string() == '-0h00m00s'
    assert Angle(-1., unit='hour').to_string() == '-1h00m00s'
