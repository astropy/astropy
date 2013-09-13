# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""
This includes tests for distances/cartesian points that are *not* in the API
tests.  Right now that's just regression tests.
"""

import numpy as np
from numpy import testing as npt

from ...tests.helper import pytest
from ... import units as u



def test_distance_change():
    from .. import Longitude, Latitude, ICRS, Distance

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    c = ICRS(ra, dec)

    c.distance = Distance(1, unit=u.kpc)

    oldx = c.x.value
    assert (oldx - 0.35284083171901953) < 1e-10

    #now x should increase when the distance increases
    c.distance = Distance(2, unit=u.kpc)

    assert c.x.value == oldx * 2


def test_distance_is_quantity():
    """
    test that distance behaves like a proper quantity
    """

    from .. import Distance

    Distance(2 * u.kpc)

    d = Distance([2, 3.1], u.kpc)

    assert d.shape == (2,)

    a = d.view(np.ndarray)
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
    npt.assert_allclose(d.value, 1301442.9440836983)

    with pytest.raises(ValueError):
        d = Distance(value=d, distmod=20)

    with pytest.raises(ValueError):
        d = Distance(z=.23, distmod=20)

def test_distance_in_coordinates():
    """
    test that distances can be created from quantities and that CartesianPoints
    can be built from them sucessfully
    """
    from .. import Longitude, Latitude, ICRS, CartesianPoints

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    c = ICRS(ra, dec)

    c.distance = 2 * u.kpc  # auto-converts

    #make sure cartesian stuff now works
    # do the internal method first, because the properties will fail with
    # unhelpful errors because __getattr_ is overridden by coordinates
    c._make_cart()

    # c.x, c.y, c.z and such are in test_api
    #make sure repr still works

    repr(c)

    assert isinstance(c.cartesian, CartesianPoints)


def test_creating_cartesian_single():
    """
    test building cartesian points with the single-argument constructor
    """

    from .. import CartesianPoints

    CartesianPoints(np.ones((3, 10)), unit=u.kpc)

    #allow dimensionless, too
    CartesianPoints(np.ones((3, 10)), unit=u.dimensionless_unscaled)

    with pytest.raises(u.UnitsError):
        CartesianPoints(np.ones((3, 10)))

    with pytest.raises(ValueError):
        CartesianPoints(np.ones((2, 10)), unit=u.kpc)

    #quantity version
    c = CartesianPoints(np.ones((3, 10)) * u.kpc)
    assert c.unit == u.kpc

    c = CartesianPoints(np.ones((3, 10)) * u.kpc, unit=u.Mpc)
    assert c.unit == u.Mpc


def test_creating_cartesian_triple():
    """
    test building cartesian points with the `x`,`y`,`z` constructor
    """

    from .. import CartesianPoints

    #make sure scalars are scalars
    c = CartesianPoints(1, 2, 3, unit=u.kpc)

    CartesianPoints(np.ones(10), np.ones(10), np.ones(10), unit=u.kpc)

    with pytest.raises(ValueError):
        #shapes must match
        CartesianPoints(np.ones(8), np.ones(12), np.ones(9), unit=u.kpc)

    #if one is a quantity, use that unit
    c = CartesianPoints(np.ones(10), np.ones(10), np.ones(10) * u.kpc)
    assert c.unit == u.kpc

    #convert when needed
    c = CartesianPoints(np.ones(10), np.ones(10), np.ones(10) * u.kpc,
                        unit=u.Mpc)
    assert c.unit == u.Mpc
    assert c[2][0].value < .1  # conversion of kpc to Mpc should give much smaller, but do this for round-off

    with pytest.raises(u.UnitsError):
        CartesianPoints(np.ones(10) * u.Mpc, np.ones(10), np.ones(10) * u.kpc)


def test_cartesian_operations():
    """
    more tests of CartesianPoints beyond those in test_api
    """

    from .. import Longitude, Latitude
    from .. import CartesianPoints

    c = CartesianPoints(np.ones(10), np.ones(10), np.ones(10), unit=u.kpc)

    c2 = c + c
    assert c2.y[2].value == 2
    assert c2.unit == c.unit

    c3 = c - c
    assert c3[1, 2].value == 0

    r, lat, lon = c.to_spherical()

    assert r.unit == c.unit
    assert isinstance(lat, Latitude)
    assert isinstance(lon, Longitude)

    c4 = c * 3
    assert c4.unit == c.unit

    #always preserve the CartesianPoint's units
    c5 = 3 * u.pc + c
    assert c5.unit == c.unit


def test_cartesian_view():
    """
    test that the cartesian subclass properly deals with new views
    """

    from .. import CartesianPoints

    c = CartesianPoints(np.ones(10), np.ones(10), np.ones(10), unit=u.kpc)

    c2 = CartesianPoints(c, copy=False)
    asarr = c.view(np.ndarray)
    assert np.all(asarr == 1)

    asarr.ravel()[0] = 2
    assert not np.all(asarr == 1)
    assert not c.x[0] == 2
