# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Test initalization of angles not already covered by the API tests

from ..angles import Angle, RA, Dec, BoundsError
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
