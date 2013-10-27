# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six
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
    npt.assert_almost_equal(m, [0, 30, -30])
    npt.assert_almost_equal(s, [0, 0, -0])

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
    npt.assert_almost_equal(m, [0, 30, -30])
    npt.assert_almost_equal(s, [0, 0, -0])

    hms = a1.hms
    hours = hms_to_hours(*hms)
    npt.assert_almost_equal(a1.hour, hours)

    a2 = Angle(hms, unit=u.hour)

    npt.assert_almost_equal(a2.radian, a1.radian)


def test_array_coordinates_creation():
    """
    Test creating coordinates from arrays.
    """
    from .. import ICRS

    c = ICRS(np.array([1, 2]), np.array([3, 4]), unit=(u.deg, u.deg))
    assert not c.isscalar

    with pytest.raises(ValueError):
        c = ICRS(np.array([1, 2]), np.array([3, 4, 5]), unit=(u.deg, u.deg))
    with pytest.raises(ValueError):
        c = ICRS(np.array([1, 2]), np.array([[3, 4], [5, 6]]), unit=(u.deg, u.deg))

    #make sure cartesian initialization also works
    c = ICRS(x=np.array([1, 2]), y=np.array([3, 4]), z=np.array([5, 6]), unit=u.kpc)

    #also ensure strings can be arrays
    c = ICRS(np.array(['1d0m0s', '2h02m00.3s']), np.array(['3d', '4d']), unit=(u.deg, u.deg))

    #but invalid strings cannot
    with pytest.raises(ValueError):
        c = ICRS(np.array(['10m0s', '2h02m00.3s']), np.array(['3d', '4d']), unit=(u.deg, u.deg))
    with pytest.raises(ValueError):
        c = ICRS(np.array(['1d0m0s', '2h02m00.3s']), np.array(['3x', '4d']), unit=(u.deg, u.deg))


def test_array_coordinates_distances():
    """
    Test creating coordinates from arrays and distances.
    """
    from .. import ICRS

    #correct way
    ICRS(np.array([1, 2]), np.array([3, 4]), unit=(u.deg, u.deg), distance= [.1, .2] * u.kpc)

    with pytest.raises(ValueError):
        #scalar distance and array coordinates
        ICRS(np.array([1, 2]), np.array([[3, 4], [5, 6]]), unit=(u.deg, u.deg), distance= 2. * u.kpc)
    with pytest.raises(ValueError):
        #scalar coordinates and array distance
        ICRS(1., 2., unit=(u.deg, u.deg), distance= [.1, .2, 3.] * u.kpc)
    with pytest.raises(ValueError):
        #more distance values than coordinates
        ICRS(np.array([1, 2]), np.array([[3, 4], [5, 6]]), unit=(u.deg, u.deg), distance= [.1, .2, 3.] * u.kpc)


@pytest.mark.parametrize(('arrshape', 'distance'), [((2, ), None), ((4, 2, 5), None), ((4, 2, 5), 2 * u.kpc)])
def test_array_coordinates_transformations(arrshape, distance):
    """
    Test transformation on coordinates with array content (first length-2 1D, then a 3D array)
    """
    from .. import ICRS, Galactic

    #M31 coordinates from test_transformations
    raarr = np.ones(arrshape) * 10.6847929
    decarr = np.ones(arrshape) * 41.2690650
    if distance is not None:
        distance = np.ones(arrshape) * distance

    c = ICRS(raarr, decarr, unit=(u.deg, u.deg), distance=distance)
    g = c.transform_to(Galactic)

    assert g.l.shape == arrshape

    npt.assert_array_almost_equal(g.l.degree, 121.17447049007306)
    npt.assert_array_almost_equal(g.b.degree, -21.57291080408368)

    if distance is not None:
        assert g.distance.unit == c.distance.unit

    #now make sure round-tripping works through FK5 - this exercises both static and dynamic transform matricies
    c2 = c.fk5.icrs
    npt.assert_array_almost_equal(c.ra, c2.ra)
    npt.assert_array_almost_equal(c.dec, c2.dec)

    assert c2.ra.shape == arrshape

    if distance is not None:
        assert c2.distance.unit == c.distance.unit

    #also make sure it's possible to get to FK4, which uses a direct transform function.
    fk4 = c.fk4

    npt.assert_array_almost_equal(fk4.ra.degree, 10.0004, decimal=4)
    npt.assert_array_almost_equal(fk4.dec.degree, 40.9953, decimal=4)

    assert fk4.ra.shape == arrshape
    if distance is not None:
        assert fk4.distance.unit == c.distance.unit

    #now check the reverse transforms run
    cfk4 = fk4.icrs
    assert cfk4.ra.shape == arrshape


def test_array_coordinates_string():
    """
    tests for string representations of aarray coordinates
    """
    from .. import ICRS

    c = ICRS(np.array([1, 2]), np.array([3, 4]), unit=(u.deg, u.deg))
    str(c)
    six.text_type(c)
    repr(c)

    assert repr(c) == '<ICRS RA=[1 2] deg, Dec=[3 4] deg>'

    #also check with distance

    c = ICRS(np.array([1, 2]), np.array([3, 4]), unit=(u.deg, u.deg), distance= u.kpc * [0.5, 1.5])
    str(c)
    six.text_type(c)
    repr(c)

    print(repr(c))

    assert repr(c) == '<ICRS RA=[1 2] deg, Dec=[3 4] deg, Distance=[ 0.5  1.5] kpc>'


def test_array_precession():
    """
    Ensures that FK5 coordinates as arrays precess their equinoxes
    """
    from ...time import Time
    from .. import FK5

    j2000 = Time('J2000', scale='utc')
    j1975 = Time('J1975', scale='utc')

    fk5 = FK5([1, 1.1], [0.5, 0.6], unit=(u.radian, u.radian))
    assert fk5.equinox.jyear == j2000.jyear
    fk5_2 = fk5.precess_to(j1975)
    assert fk5_2.equinox.jyear == j1975.jyear

    npt.assert_array_less(0.05, np.abs(fk5.ra.degree - fk5_2.ra.degree))
    npt.assert_array_less(0.05, np.abs(fk5.dec.degree - fk5_2.dec.degree))

def test_array_separation():
    from .. import ICRS

    c1 = ICRS([0 , 0], [0, 0], unit=(u.degree, u.degree))
    c2 = ICRS([1, 2], [0, 0], unit=(u.degree, u.degree))

    npt.assert_array_almost_equal(c1.separation(c2).degree, [1, 2])

    c3 = ICRS([0 , 3.], [0., 0], unit=(u.degree, u.degree), distance=[1 ,1.] * u.kpc)
    c4 = ICRS([1, 1.], [0., 0], unit=(u.degree, u.degree), distance=[1 ,1.] * u.kpc)
    
    #the 3-1 separation should be twice the 0-1 separation, but not *exactly* the same
    sep = c3.separation_3d(c4)
    sepdiff = sep[1] - (2 * sep[0])

    assert abs(sepdiff.value) < 1e-5
    assert sepdiff != 0

def test_array_indexing():
    from .. import FK5Coordinates
    from ...time import Time

    ra = np.linspace(0, 360, 10)
    dec = np.linspace(-90, 90, 10)
    j1975 = Time(1975, format='jyear', scale='utc')

    c1 = FK5Coordinates(ra, dec, unit=(u.degree, u.degree), equinox=j1975)

    c2 = c1[4]
    assert c2.ra.degree == 160
    assert c2.dec.degree == -10

    c3 = c1[2:5]
    npt.assert_array_equal(c3.ra.degree, [80, 120, 160])
    npt.assert_array_equal(c3.dec.degree, [-50, -30, -10])

    c4 = c1[np.array([2, 5, 8])]

    npt.assert_array_equal(c4.ra.degree, [80, 200, 320])
    npt.assert_array_equal(c4.dec.degree, [-50, 10, 70])

    #now make sure the equinox is preserved
    assert c2.equinox == c1.equinox
    assert c3.equinox == c1.equinox
    assert c4.equinox == c1.equinox
