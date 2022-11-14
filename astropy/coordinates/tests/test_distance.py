# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This includes tests for the Distance class and related calculations
"""

import numpy as np
import pytest
from numpy import testing as npt

from astropy import units as u
from astropy.coordinates import CartesianRepresentation, Distance, Latitude, Longitude
from astropy.coordinates.builtin_frames import ICRS, Galactic
from astropy.units import allclose as quantity_allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyWarning

MULTIPLE_INPUTS_ERROR_MSG = "^more than one of `.*` were given to Distance constructor$"


def test_distances():
    """
    Tests functionality for Coordinate class distances and cartesian
    transformations.
    """

    """
    Distances can also be specified, and allow for a full 3D definition of a
    coordinate.
    """

    # try all the different ways to initialize a Distance
    distance = Distance(12, u.parsec)
    Distance(40, unit=u.au)
    Distance(value=5, unit=u.kpc)

    # need to provide a unit
    with pytest.raises(u.UnitsError):
        Distance(12)

    with pytest.raises(ValueError, match="none of `value`, `z`, `distmod`,"):
        Distance(unit=u.km)

    # standard units are pre-defined
    npt.assert_allclose(distance.lyr, 39.138765325702551)
    npt.assert_allclose(distance.km, 370281309776063.0)

    # Coordinate objects can be assigned a distance object, giving them a full
    # 3D position
    c = Galactic(
        l=158.558650 * u.degree,
        b=-43.350066 * u.degree,
        distance=Distance(12, u.parsec),
    )
    assert quantity_allclose(c.distance, 12 * u.pc)

    # or initialize distances via redshifts - this is actually tested in the
    # function below that checks for scipy. This is kept here as an example
    # c.distance = Distance(z=0.2)  # uses current cosmology
    # with whatever your preferred cosmology may be
    # c.distance = Distance(z=0.2, cosmology=WMAP5)

    # Coordinate objects can be initialized with a distance using special
    # syntax
    c1 = Galactic(l=158.558650 * u.deg, b=-43.350066 * u.deg, distance=12 * u.kpc)

    # Coordinate objects can be instantiated with cartesian coordinates
    # Internally they will immediately be converted to two angles + a distance
    cart = CartesianRepresentation(x=2 * u.pc, y=4 * u.pc, z=8 * u.pc)
    c2 = Galactic(cart)

    sep12 = c1.separation_3d(c2)
    # returns a *3d* distance between the c1 and c2 coordinates
    # not that this does *not*
    assert isinstance(sep12, Distance)
    npt.assert_allclose(sep12.pc, 12005.784163916317, 10)

    """
    All spherical coordinate systems with distances can be converted to
    cartesian coordinates.
    """

    cartrep2 = c2.cartesian
    assert isinstance(cartrep2.x, u.Quantity)
    npt.assert_allclose(cartrep2.x.value, 2)
    npt.assert_allclose(cartrep2.y.value, 4)
    npt.assert_allclose(cartrep2.z.value, 8)

    # with no distance, the unit sphere is assumed when converting to cartesian
    c3 = Galactic(l=158.558650 * u.degree, b=-43.350066 * u.degree, distance=None)
    unitcart = c3.cartesian
    npt.assert_allclose(
        ((unitcart.x**2 + unitcart.y**2 + unitcart.z**2) ** 0.5).value, 1.0
    )

    # TODO: choose between these when CartesianRepresentation gets a definite
    # decision on whether or not it gets __add__
    #
    # CartesianRepresentation objects can be added and subtracted, which are
    # vector/elementwise they can also be given as arguments to a coordinate
    # system
    # csum = ICRS(c1.cartesian + c2.cartesian)
    csumrep = CartesianRepresentation(c1.cartesian.xyz + c2.cartesian.xyz)
    csum = ICRS(csumrep)

    npt.assert_allclose(csumrep.x.value, -8.12016610185)
    npt.assert_allclose(csumrep.y.value, 3.19380597435)
    npt.assert_allclose(csumrep.z.value, -8.2294483707)
    npt.assert_allclose(csum.ra.degree, 158.529401774)
    npt.assert_allclose(csum.dec.degree, -43.3235825777)
    npt.assert_allclose(csum.distance.kpc, 11.9942200501)


@pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
def test_distances_scipy():
    """
    The distance-related tests that require scipy due to the cosmology
    module needing scipy integration routines
    """
    from astropy.cosmology import WMAP5

    # try different ways to initialize a Distance
    d4 = Distance(z=0.23)  # uses default cosmology - as of writing, WMAP7
    npt.assert_allclose(d4.z, 0.23, rtol=1e-8)

    d5 = Distance(z=0.23, cosmology=WMAP5)
    npt.assert_allclose(d5.compute_z(WMAP5), 0.23, rtol=1e-8)

    d6 = Distance(z=0.23, cosmology=WMAP5, unit=u.km)
    npt.assert_allclose(d6.value, 3.5417046898762366e22)

    with pytest.raises(ValueError, match="a `cosmology` was given but `z`"):
        Distance(parallax=1 * u.mas, cosmology=WMAP5)

    # Regression test for #12531
    with pytest.raises(ValueError, match=MULTIPLE_INPUTS_ERROR_MSG):
        Distance(z=0.23, parallax=1 * u.mas)

    # vectors!  regression test for #11949
    d4 = Distance(z=[0.23, 0.45])  # as of writing, Planck18
    npt.assert_allclose(d4.z, [0.23, 0.45], rtol=1e-8)


def test_distance_change():
    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    c1 = ICRS(ra, dec, Distance(1, unit=u.kpc))

    oldx = c1.cartesian.x.value
    assert (oldx - 0.35284083171901953) < 1e-10

    # first make sure distances are immutable
    with pytest.raises(AttributeError):
        c1.distance = Distance(2, unit=u.kpc)

    # now x should increase with a bigger distance increases
    c2 = ICRS(ra, dec, Distance(2, unit=u.kpc))
    assert c2.cartesian.x.value == oldx * 2


def test_distance_is_quantity():
    """
    test that distance behaves like a proper quantity
    """

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

    # regression test against #2261
    d = Distance([2 * u.kpc, 250.0 * u.pc])
    assert d.unit is u.kpc
    assert np.all(d.value == np.array([2.0, 0.25]))


def test_distmod():
    d = Distance(10, u.pc)
    assert d.distmod.value == 0

    d = Distance(distmod=20)
    assert d.distmod.value == 20
    assert d.kpc == 100

    d = Distance(distmod=-1.0, unit=u.au)
    npt.assert_allclose(d.value, 1301442.9440836983)

    with pytest.raises(ValueError, match=MULTIPLE_INPUTS_ERROR_MSG):
        d = Distance(value=d, distmod=20)

    with pytest.raises(ValueError, match=MULTIPLE_INPUTS_ERROR_MSG):
        d = Distance(z=0.23, distmod=20)

    # check the Mpc/kpc/pc behavior
    assert Distance(distmod=1).unit == u.pc
    assert Distance(distmod=11).unit == u.kpc
    assert Distance(distmod=26).unit == u.Mpc
    assert Distance(distmod=-21).unit == u.AU

    # if an array, uses the mean of the log of the distances
    assert Distance(distmod=[1, 11, 26]).unit == u.kpc


def test_parallax():
    d = Distance(parallax=1 * u.arcsecond)
    assert d.pc == 1.0

    with pytest.raises(ValueError, match=MULTIPLE_INPUTS_ERROR_MSG):
        d = Distance(15 * u.pc, parallax=20 * u.milliarcsecond)

    with pytest.raises(ValueError, match=MULTIPLE_INPUTS_ERROR_MSG):
        d = Distance(parallax=20 * u.milliarcsecond, distmod=20)

    # array
    plx = [1, 10, 100.0] * u.mas
    d = Distance(parallax=plx)
    assert quantity_allclose(d.pc, [1000.0, 100.0, 10.0])
    assert quantity_allclose(plx, d.parallax)

    error_message = (
        r"^some parallaxes are negative, which are not interpretable as distances\. "
    )
    with pytest.raises(ValueError, match=error_message):
        Distance(parallax=-1 * u.mas)

    with pytest.raises(ValueError, match=error_message):
        Distance(parallax=[10, 1, -1] * u.mas)

    warning_message = "^negative parallaxes are converted to NaN distances even when"

    with pytest.warns(AstropyWarning, match=warning_message):
        Distance(parallax=-1 * u.mas, allow_negative=True)

    with pytest.warns(AstropyWarning, match=warning_message):
        Distance(parallax=[10, 1, -1] * u.mas, allow_negative=True)

    # Regression test for #12569; `unit` was ignored if `parallax` was given.
    d = Distance(parallax=1 * u.mas, unit=u.kpc)
    assert d.value == 1.0
    assert d.unit is u.kpc


def test_distance_in_coordinates():
    """
    test that distances can be created from quantities and that cartesian
    representations come out right
    """

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    coo = ICRS(ra, dec, distance=2 * u.kpc)

    cart = coo.cartesian

    assert isinstance(cart.xyz, u.Quantity)


def test_negative_distance():
    """Test optional kwarg allow_negative"""

    error_message = (
        r"^distance must be >= 0\. Use the argument `allow_negative=True` to allow "
        r"negative values\.$"
    )

    with pytest.raises(ValueError, match=error_message):
        Distance([-2, 3.1], u.kpc)

    with pytest.raises(ValueError, match=error_message):
        Distance([-2, -3.1], u.kpc)

    with pytest.raises(ValueError, match=error_message):
        Distance(-2, u.kpc)

    d = Distance(-2, u.kpc, allow_negative=True)
    assert d.value == -2


def test_distance_comparison():
    """Ensure comparisons of distances work (#2206, #2250)"""
    a = Distance(15 * u.kpc)
    b = Distance(15 * u.kpc)
    assert a == b
    c = Distance(1.0 * u.Mpc)
    assert a < c


def test_distance_to_quantity_when_not_units_of_length():
    """Any operation that leaves units other than those of length
    should turn a distance into a quantity (#2206, #2250)"""
    d = Distance(15 * u.kpc)
    twice = 2.0 * d
    assert isinstance(twice, Distance)
    area = 4.0 * np.pi * d**2
    assert area.unit.is_equivalent(u.m**2)
    assert not isinstance(area, Distance)
    assert type(area) is u.Quantity


def test_distance_nan():
    # Check that giving NaNs to Distance doesn't emit a warning
    Distance([0, np.nan, 1] * u.m)
