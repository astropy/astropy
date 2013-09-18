# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from numpy import testing as npt

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

def test_distance_from_quantity():
    from .. import Longitude, Latitude, ICRSCoordinates, Distance

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    c = ICRSCoordinates(ra, dec)

    # a Quantity object should be able to supply a distance
    q = 2 * u.kpc
    c.distance = q
