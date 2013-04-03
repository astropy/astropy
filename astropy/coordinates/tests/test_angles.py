# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Test initalization of angles not already covered by the API tests

from ..angles import Angle
from ...tests.compat import assert_allclose
from ... import units as u


def test_negative_zero_dms():
    # Test for DMS parser
    a = Angle('-00:00:10', u.deg)
    assert_allclose(a.degrees, -10. / 3600.)


def test_negative_zero_dm():
    # Test for DM parser
    a = Angle('-00:10', u.deg)
    assert_allclose(a.degrees, -10. / 60.)


def test_negative_zero_hms():
    # Test for HMS parser
    a = Angle('-00:00:10', u.hour)
    assert_allclose(a.hours, -10. / 3600.)


def test_negative_zero_hm():
    # Test for HM parser
    a = Angle('-00:10', u.hour)
    assert_allclose(a.hours, -10. / 60.)
