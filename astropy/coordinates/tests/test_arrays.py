# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...tests.helper import pytest

import numpy as np
from numpy import testing as npt

from ... import units as u


def test_angle_arrays():
    """
    Test arrays values with Angle objects.
    """
    from .. import Angle

    # Tests incomplete
    a1 = Angle([0, 45, 90, 180, 270, 360, 720.], unit=u.degree)
    npt.assert_almost_equal([0., 45., 90., 180., 270., 360., 720.], a1.value)

    a2 = Angle(np.array([-90, -45, 0, 45, 90, 180, 270, 360]), unit=u.degree)
    npt.assert_almost_equal([-90, -45, 0, 45, 90, 180, 270, 360],
                            a2.value)

    a3 = Angle(["12 degrees", "3 hours", "5 deg", "4rad"])
    npt.assert_almost_equal([12., 45., 5., 229.18311805],
                            a3.value)
    assert a3.unit == u.degree

    a4 = Angle(["12 degrees", "3 hours", "5 deg", "4rad"], u.radian)
    npt.assert_almost_equal(a4.degree, a3.value)
    assert a4.unit == u.radian

    a5 = Angle([0, 45, 90, 180, 270, 360], unit=u.degree)
    a6 = a5.sum()
    npt.assert_almost_equal(a6.value, 945.0)
    assert a6.unit is u.degree

    with pytest.raises(TypeError):
        # Arrays of Angle objects are not supported -- that's really
        # tricky to do correctly, if at all, due to the possibility of
        # nesting.
        a7 = Angle([a1, a2, a3], unit=u.degree)

    a8 = Angle(["04:02:02", "03:02:01", "06:02:01"], unit=u.degree)
    npt.assert_almost_equal(a8.value, [4.03388889, 3.03361111, 6.03361111])

    a9 = Angle(np.array(["04:02:02", "03:02:01", "06:02:01"]), unit=u.degree)
    npt.assert_almost_equal(a9.value, a8.value)

    with pytest.raises(u.UnitsError):
        a10 = Angle(["04:02:02", "03:02:01", "06:02:01"])


def test_dms():
    from .. import Angle
    from ..angle_utilities import dms_to_degrees

    a1 = Angle([0, 45.5, -45.5], unit=u.degree)
    d, m, s = a1.dms
    npt.assert_almost_equal(d, [0, 45, -45])
    npt.assert_almost_equal(m, [0, 30, 30])
    npt.assert_almost_equal(s, [0, 0, 0])

    dms = a1.dms
    degrees = dms_to_degrees(*dms)
    npt.assert_almost_equal(a1.degree, degrees)

    a2 = Angle(dms, unit=u.degree)

    npt.assert_almost_equal(a2.radian, a1.radian)


def test_hms():
    from .. import Angle
    from ..angle_utilities import hms_to_hours

    a1 = Angle([0, 11.5, -11.5], unit=u.hour)
    h, m, s = a1.hms
    npt.assert_almost_equal(h, [0, 11, -11])
    npt.assert_almost_equal(m, [0, 30, 30])
    npt.assert_almost_equal(s, [0, 0, 0])

    hms = a1.hms
    hours = hms_to_hours(*hms)
    npt.assert_almost_equal(a1.hour, hours)

    a2 = Angle(hms, unit=u.hour)

    npt.assert_almost_equal(a2.radian, a1.radian)


def test_array_coordinates_creation():
    """
    Test creating coordinates from arrays.
    """
    from .. import ICRSCoordinates

    c = ICRSCoordinates(np.array([1, 2]), np.array([3, 4]), unit=(u.deg, u.deg))
    assert not c.isscalar

    with pytest.raises(ValueError):
        c = ICRSCoordinates(np.array([1, 2]), np.array([3, 4, 5]), unit=(u.deg, u.deg))
    with pytest.raises(ValueError):
        c = ICRSCoordinates(np.array([1, 2]), np.array([[3, 4], [5, 6]]), unit=(u.deg, u.deg))

    #make sure cartesian initialization also works
    c = ICRSCoordinates(x=np.array([1, 2]), y=np.array([3, 4]), z=np.array([5, 6]), unit=u.kpc)

def test_array_coordinates_distances():
    """
    Test creating coordinates from arrays and distances.
    """
    from .. import ICRSCoordinates

    #correct way
    ICRSCoordinates(np.array([1, 2]), np.array([3, 4]), unit=(u.deg, u.deg), distance= [.1,.2] * u.kpc)

    with pytest.raises(ValueError):
        #scalar distance and array coordinates
        ICRSCoordinates(np.array([1, 2]), np.array([[3, 4], [5, 6]]), unit=(u.deg, u.deg), distance= 2. * u.kpc)
    with pytest.raises(ValueError):
        #scalar coordinates and array distance
        ICRSCoordinates(1., 2., unit=(u.deg, u.deg), distance= [.1,.2, 3.] * u.kpc)
    with pytest.raises(ValueError):
        #more distance values than coordinates
        ICRSCoordinates(np.array([1, 2]), np.array([[3, 4], [5, 6]]), unit=(u.deg, u.deg), distance= [.1,.2, 3.] * u.kpc)

@pytest.mark.parametrize(('arrshape'), [(2, ), (4, 2, 5)])
def test_array_coordinates_transformations(arrshape):
    """
    Test transformation on coordinates with array content (first length-2 1D, then a 3D array)
    """
    from .. import ICRSCoordinates, GalacticCoordinates
    
    #M31 coordinates from test_transformations
    raarr = np.ones(arrshape) * 10.6847929
    decarr = np.ones(arrshape) * 41.2690650
    c = ICRSCoordinates(raarr, decarr, unit=(u.deg, u.deg))
    g = c.transform_to(GalacticCoordinates)

    assert len(g.l) == 2
    assert len(g.b) == 2

    npt.assert_array_almost_equal(g.l.degree, 121.17447049007306)
    npt.assert_array_almost_equal(g.b.degree, -21.57291080408368)

    #now make sure round-tripping works through FK5 - this exercises both static and dynamic transform matricies
    c2 = c.fk5.icrs
    npt.assert_array_almost_equal(c.ra, c2.ra)
    npt.assert_array_almost_equal(c.dec, c2.dec)

    #also just make sure it's possible to get to FK4, which uses a direct trasnform function.
    #c.fk4

def test_array_coordinates_string():
    """
    tests for string representations of aarray coordinates
    """
    from .. import ICRSCoordinates

    c = ICRSCoordinates(np.array([1, 2]), np.array([3, 4]), unit=(u.deg, u.deg))
    str(c)
    unicode(c)
    repr(c)

    assert repr(c) == '<ICRSCoordinates RA=[1 2] deg, Dec=[3 4] deg>'


