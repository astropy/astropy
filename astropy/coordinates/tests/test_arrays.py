# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
import numpy as np
from numpy import testing as npt

from ... import units as u
from ...time import Time
from ...tests.helper import assert_quantity_allclose as assert_allclose

from .. import (Angle, ICRS, FK4, FK5, Galactic, SkyCoord,
                CartesianRepresentation)
from ..angle_utilities import dms_to_degrees, hms_to_hours


def test_angle_arrays():
    """
    Test arrays values with Angle objects.
    """
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
        # Arrays where the elements are Angle objects are not supported -- it's
        # really tricky to do correctly, if at all, due to the possibility of
        # nesting.
        a7 = Angle([a1, a2, a3], unit=u.degree)

    a8 = Angle(["04:02:02", "03:02:01", "06:02:01"], unit=u.degree)
    npt.assert_almost_equal(a8.value, [4.03388889, 3.03361111, 6.03361111])

    a9 = Angle(np.array(["04:02:02", "03:02:01", "06:02:01"]), unit=u.degree)
    npt.assert_almost_equal(a9.value, a8.value)

    with pytest.raises(u.UnitsError):
        a10 = Angle(["04:02:02", "03:02:01", "06:02:01"])


def test_dms():
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
    c = ICRS(np.array([1, 2])*u.deg, np.array([3, 4])*u.deg)
    assert not c.ra.isscalar

    with pytest.raises(ValueError):
        c = ICRS(np.array([1, 2])*u.deg, np.array([3, 4, 5])*u.deg)
    with pytest.raises(ValueError):
        c = ICRS(np.array([1, 2, 4, 5])*u.deg, np.array([[3, 4], [5, 6]])*u.deg)

    # make sure cartesian initialization also works
    cart = CartesianRepresentation(x=[1., 2.]*u.kpc, y=[3., 4.]*u.kpc, z=[5., 6.]*u.kpc)
    c = ICRS(cart)

    # also ensure strings can be arrays
    c = SkyCoord(['1d0m0s', '2h02m00.3s'], ['3d', '4d'])

    # but invalid strings cannot
    with pytest.raises(ValueError):
        c = SkyCoord(Angle(['10m0s', '2h02m00.3s']), Angle(['3d', '4d']))
    with pytest.raises(ValueError):
        c = SkyCoord(Angle(['1d0m0s', '2h02m00.3s']), Angle(['3x', '4d']))


def test_array_coordinates_distances():
    """
    Test creating coordinates from arrays and distances.
    """
    # correct way
    ICRS(ra=np.array([1, 2])*u.deg, dec=np.array([3, 4])*u.deg, distance=[.1, .2] * u.kpc)

    with pytest.raises(ValueError):
        # scalar distance and mismatched array coordinates
        ICRS(ra=np.array([1, 2, 3])*u.deg, dec=np.array([[3, 4], [5, 6]])*u.deg, distance=2. * u.kpc)
    with pytest.raises(ValueError):
        # more distance values than coordinates
        ICRS(ra=np.array([1, 2])*u.deg, dec=np.array([3, 4])*u.deg, distance=[.1, .2, 3.] * u.kpc)


@pytest.mark.parametrize(('arrshape', 'distance'), [((2, ), None), ((4, 2, 5), None), ((4, 2, 5), 2 * u.kpc)])
def test_array_coordinates_transformations(arrshape, distance):
    """
    Test transformation on coordinates with array content (first length-2 1D, then a 3D array)
    """
    # M31 coordinates from test_transformations
    raarr = np.ones(arrshape) * 10.6847929
    decarr = np.ones(arrshape) * 41.2690650
    if distance is not None:
        distance = np.ones(arrshape) * distance

    print(raarr, decarr, distance)
    c = ICRS(ra=raarr*u.deg, dec=decarr*u.deg, distance=distance)
    g = c.transform_to(Galactic)

    assert g.l.shape == arrshape

    npt.assert_array_almost_equal(g.l.degree, 121.17440967)
    npt.assert_array_almost_equal(g.b.degree, -21.57299631)

    if distance is not None:
        assert g.distance.unit == c.distance.unit

    # now make sure round-tripping works through FK5
    c2 = c.transform_to(FK5).transform_to(ICRS)
    npt.assert_array_almost_equal(c.ra.radian, c2.ra.radian)
    npt.assert_array_almost_equal(c.dec.radian, c2.dec.radian)

    assert c2.ra.shape == arrshape

    if distance is not None:
        assert c2.distance.unit == c.distance.unit

    # also make sure it's possible to get to FK4, which uses a direct transform function.
    fk4 = c.transform_to(FK4)

    npt.assert_array_almost_equal(fk4.ra.degree, 10.0004, decimal=4)
    npt.assert_array_almost_equal(fk4.dec.degree, 40.9953, decimal=4)

    assert fk4.ra.shape == arrshape
    if distance is not None:
        assert fk4.distance.unit == c.distance.unit

    # now check the reverse transforms run
    cfk4 = fk4.transform_to(ICRS)
    assert cfk4.ra.shape == arrshape


def test_array_precession():
    """
    Ensures that FK5 coordinates as arrays precess their equinoxes
    """
    j2000 = Time('J2000', scale='utc')
    j1975 = Time('J1975', scale='utc')

    fk5 = FK5([1, 1.1]*u.radian, [0.5, 0.6]*u.radian)
    assert fk5.equinox.jyear == j2000.jyear
    fk5_2 = fk5.transform_to(FK5(equinox=j1975))
    assert fk5_2.equinox.jyear == j1975.jyear

    npt.assert_array_less(0.05, np.abs(fk5.ra.degree - fk5_2.ra.degree))
    npt.assert_array_less(0.05, np.abs(fk5.dec.degree - fk5_2.dec.degree))


def test_array_separation():
    c1 = ICRS([0, 0]*u.deg, [0, 0]*u.deg)
    c2 = ICRS([1, 2]*u.deg, [0, 0]*u.deg)

    npt.assert_array_almost_equal(c1.separation(c2).degree, [1, 2])

    c3 = ICRS([0, 3.]*u.deg, [0., 0]*u.deg, distance=[1, 1.] * u.kpc)
    c4 = ICRS([1, 1.]*u.deg, [0., 0]*u.deg, distance=[1, 1.] * u.kpc)

    # the 3-1 separation should be twice the 0-1 separation, but not *exactly* the same
    sep = c3.separation_3d(c4)
    sepdiff = sep[1] - (2 * sep[0])

    assert abs(sepdiff.value) < 1e-5
    assert sepdiff != 0


def test_array_indexing():
    ra = np.linspace(0, 360, 10)
    dec = np.linspace(-90, 90, 10)
    j1975 = Time(1975, format='jyear', scale='utc')

    c1 = FK5(ra*u.deg, dec*u.deg, equinox=j1975)

    c2 = c1[4]
    assert c2.ra.degree == 160
    assert c2.dec.degree == -10

    c3 = c1[2:5]
    assert_allclose(c3.ra, [80, 120, 160] * u.deg)
    assert_allclose(c3.dec, [-50, -30, -10] * u.deg)

    c4 = c1[np.array([2, 5, 8])]

    assert_allclose(c4.ra, [80, 200, 320] * u.deg)
    assert_allclose(c4.dec, [-50, 10, 70] * u.deg)

    # now make sure the equinox is preserved
    assert c2.equinox == c1.equinox
    assert c3.equinox == c1.equinox
    assert c4.equinox == c1.equinox


def test_array_len():
    input_length = [1, 5]
    for length in input_length:
        ra = np.linspace(0, 360, length)
        dec = np.linspace(0, 90, length)

        c = ICRS(ra*u.deg, dec*u.deg)

        assert len(c) == length

        assert c.shape == (length,)

    with pytest.raises(TypeError):
        c = ICRS(0*u.deg, 0*u.deg)
        len(c)

    assert c.shape == tuple()


def test_array_eq():
    c1 = ICRS([1, 2]*u.deg, [3, 4]*u.deg)
    c2 = ICRS([1, 2]*u.deg, [3, 5]*u.deg)
    c3 = ICRS([1, 3]*u.deg, [3, 4]*u.deg)
    c4 = ICRS([1, 2]*u.deg, [3, 4.2]*u.deg)

    assert c1 == c1
    assert c1 != c2
    assert c1 != c3
    assert c1 != c4
