# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from numpy import testing as npt
from ...tests.helper import pytest

from ... import units as u


"""
This includes tests for distances/cartesian points that are *not* in the API
tests.  Right now that's just regression tests.
"""


def test_distance_change():
    from .. import Longitude, Latitude, ICRSCoordinates, Distance

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    c = ICRSCoordinates(ra, dec)

    c.distance = Distance(1, unit=u.kpc)

    oldx = c.x
    assert (oldx - 0.35284083171901953) < 1e-10

    #now x should increase when the distance increases
    c.distance = Distance(2, unit=u.kpc)

    assert c.x == oldx * 2


def test_distance_in_coordinates():
    """
    test that distances can be created from quantities
    """
    from .. import Longitude, Latitude, ICRSCoordinates, CartesianPoints

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    c = ICRSCoordinates(ra, dec)

    c.distance = 2 * u.kpc  # auto-converts

    #make sure cartesian stuff now works
    # do the internal method first, because the properties will fail with
    # unhelpful errors because __getattr_ is overridden by coordinates
    c._make_cart()

    assert isinstance(c.cartesian, CartesianPoints)


def test_distance_is_quantity():
    """
    test that distance behaves like a proper quantity
    """
    from numpy import ndarray
    from .. import Distance

    Distance(2 * u.kpc)

    d = Distance([2, 3.1], u.kpc)

    d.shape == (2,)

    a = d.view(ndarray)
    q = d.view(u.Quantity)
    a[0] = 1.2
    q.value[1] = 5.4

    assert d[0].value == 1.2
    assert d[1].value == 5.4

    q = u.Quantity(d, copy=True)
    q.value[1] = 0
    assert q.value[1] == 0
    assert d.value[1] != 0


def test_distmod():
    from .. import Distance

    d = Distance(10, u.pc)
    assert d.distmod.value == 0

    d = Distance(distmod=20)
    assert d.distmod.value == 20
    assert d.kpc == 100

    d = Distance(distmod=-1., unit=u.au)
    npt.assert_almost_equal(d.value, 1301442.9440836983)

    with pytest.raises(ValueError):
        d = Distance(value=d, distmod=20)

    with pytest.raises(ValueError):
        d = Distance(z=.23, distmod=20)
