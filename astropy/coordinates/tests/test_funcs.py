# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for miscellaneous functionality in the `funcs` module
"""


import numpy as np
import pytest
from numpy import testing as npt

from astropy import units as u
from astropy.coordinates import FK5, ICRS, SkyCoord
from astropy.coordinates import representation as r
from astropy.coordinates.funcs import (
    concatenate,
    concatenate_representations,
    get_constellation,
    get_sun,
)
from astropy.time import Time

CARTESIAN_POS = r.CartesianRepresentation([1, 2, 3] * u.kpc)
CARTESIAN_VEL = r.CartesianDifferential([8, 9, 10] * u.km / u.s)
CARTESIAN_POS_AND_VEL = CARTESIAN_POS.with_differentials(CARTESIAN_VEL)

RADIAL_VEL = r.RadialDifferential(1 * u.km / u.s)
SPHERICAL_COS_LAT_VEL = r.SphericalCosLatDifferential(
    1 * u.mas / u.yr, 2 * u.mas / u.yr, 3 * u.km / u.s
)
SPHERICAL_POS = r.SphericalRepresentation(
    lon=1 * u.deg, lat=2.0 * u.deg, distance=10 * u.pc
)
UNIT_SPHERICAL_POS = r.UnitSphericalRepresentation(lon=1 * u.deg, lat=2.0 * u.deg)
CARTESIAN_POS_2D_ARR = r.CartesianRepresentation(np.ones((3, 100)) * u.kpc)
CARTESIAN_POS_3D_ARR = r.CartesianRepresentation(np.ones((3, 16, 8)) * u.kpc)
UNIT_SPHERICAL_COS_LAT_VEL = r.UnitSphericalCosLatDifferential(
    1 * u.mas / u.yr, 2 * u.mas / u.yr
)
CARTESIAN_VEL_2D_ARR = r.CartesianDifferential(*np.ones((3, 100)) * u.km / u.s)
CARTESIAN_VEL_3D_ARR = r.CartesianDifferential(*np.ones((3, 16, 8)) * u.km / u.s)


def test_sun():
    """
    Test that `get_sun` works and it behaves roughly as it should (in GCRS)
    """

    northern_summer_solstice = Time("2010-6-21")
    northern_winter_solstice = Time("2010-12-21")
    equinox_1 = Time("2010-3-21")
    equinox_2 = Time("2010-9-21")

    gcrs1 = get_sun(equinox_1)
    assert np.abs(gcrs1.dec.deg) < 1

    gcrs2 = get_sun(
        Time([northern_summer_solstice, equinox_2, northern_winter_solstice])
    )
    assert np.all(np.abs(gcrs2.dec - [23.5, 0, -23.5] * u.deg) < 1 * u.deg)


def test_constellations(recwarn):
    inuma = ICRS(9 * u.hour, 65 * u.deg)

    n_prewarn = len(recwarn)
    res = get_constellation(inuma)
    res_short = get_constellation(inuma, short_name=True)
    assert len(recwarn) == n_prewarn  # neither version should not make warnings

    assert res == "Ursa Major"
    assert res_short == "UMa"
    assert isinstance(res, str) or getattr(res, "shape", None) == tuple()

    # these are taken from the ReadMe for Roman 1987
    ras = [9, 23.5, 5.12, 9.4555, 12.8888, 15.6687, 19, 6.2222]
    decs = [65, -20, 9.12, -19.9, 22, -12.1234, -40, -81.1234]
    shortnames = ["UMa", "Aqr", "Ori", "Hya", "Com", "Lib", "CrA", "Men"]

    testcoos = FK5(ras * u.hour, decs * u.deg, equinox="B1950")
    npt.assert_equal(get_constellation(testcoos, short_name=True), shortnames)

    # test on a SkyCoord, *and* test Boötes, which is special in that it has a
    # non-ASCII character
    boores = get_constellation(SkyCoord(15 * u.hour, 30 * u.deg, frame="icrs"))
    assert boores == "Boötes"
    assert isinstance(boores, str) or getattr(boores, "shape", None) == tuple()


@pytest.mark.xfail
def test_constellation_edge_cases():
    # Test edge cases close to borders, using B1875.0 coordinates
    # Look for HMS / DMS roundoff-to-decimal issues from Roman (1987) data,
    # and misuse of PrecessedGeocentric, as documented in
    # https://github.com/astropy/astropy/issues/9855
    # Define eight test points.
    # The first four cross the boundary at 06h14m30 == 6.2416666666666... hours
    # with Monoceros on the west side of Orion at Dec +3.0.
    ras = [6.24100, 6.24160, 6.24166, 6.24171]
    # aka ['6h14m27.6s' '6h14m29.76s' '6h14m29.976s' '6h14m30.156s']

    decs = [3.0, 3.0, 3.0, 3.0]

    # Correct constellations for given RA/Dec coordinates
    shortnames = ["Ori", "Ori", "Ori", "Mon"]

    # The second four sample northward along RA 22 hours, crossing the boundary
    # at 86° 10' == 86.1666... degrees between Cepheus and Ursa Minor
    decs += [86.16, 86.1666, 86.16668, 86.1668]
    ras += [22.0, 22.0, 22.0, 22.0]
    shortnames += ["Cep", "Cep", "Umi", "Umi"]

    testcoos = FK5(ras * u.hour, decs * u.deg, equinox="B1875")
    npt.assert_equal(
        get_constellation(testcoos, short_name=True),
        shortnames,
        "get_constellation() error: misusing Roman approximations, vs IAU boundaries"
        " from Delporte?",
    )

    # TODO: When that's fixed, add other tests with coords that are in different constellations
    # depending on equinox


def test_concatenate():
    # Just positions
    fk5 = FK5(1 * u.deg, 2 * u.deg)
    sc = SkyCoord(3 * u.deg, 4 * u.deg, frame="fk5")

    res = concatenate([fk5, sc])
    np.testing.assert_allclose(res.ra, [1, 3] * u.deg)
    np.testing.assert_allclose(res.dec, [2, 4] * u.deg)

    with pytest.raises(TypeError):
        concatenate(fk5)

    with pytest.raises(TypeError):
        concatenate(1 * u.deg)

    # positions and velocities
    fr = ICRS(
        ra=10 * u.deg,
        dec=11.0 * u.deg,
        pm_ra_cosdec=12 * u.mas / u.yr,
        pm_dec=13 * u.mas / u.yr,
    )
    sc = SkyCoord(
        ra=20 * u.deg,
        dec=21.0 * u.deg,
        pm_ra_cosdec=22 * u.mas / u.yr,
        pm_dec=23 * u.mas / u.yr,
    )

    res = concatenate([fr, sc])

    with pytest.raises(ValueError):
        concatenate([fr, fk5])

    fr2 = ICRS(ra=10 * u.deg, dec=11.0 * u.deg)
    with pytest.raises(ValueError):
        concatenate([fr, fr2])


@pytest.mark.parametrize(
    "rep",
    (
        CARTESIAN_POS,
        SPHERICAL_POS,
        UNIT_SPHERICAL_POS,
        CARTESIAN_POS_2D_ARR,
        CARTESIAN_POS_3D_ARR,
        CARTESIAN_POS_AND_VEL,
        SPHERICAL_POS.with_differentials(SPHERICAL_COS_LAT_VEL),
        UNIT_SPHERICAL_POS.with_differentials(SPHERICAL_COS_LAT_VEL),
        UNIT_SPHERICAL_POS.with_differentials(UNIT_SPHERICAL_COS_LAT_VEL),
        UNIT_SPHERICAL_POS.with_differentials({"s": RADIAL_VEL}),
        CARTESIAN_POS_2D_ARR.with_differentials(CARTESIAN_VEL_2D_ARR),
        CARTESIAN_POS_3D_ARR.with_differentials(CARTESIAN_VEL_3D_ARR),
    ),
)
@pytest.mark.parametrize("n", (2, 4))
def test_concatenate_representations(rep, n):
    # Test that combining with itself succeeds
    expected_shape = (n * rep.shape[0],) + rep.shape[1:] if rep.shape else (n,)

    tmp = concatenate_representations(n * (rep,))
    assert tmp.shape == expected_shape

    if "s" in rep.differentials:
        assert tmp.differentials["s"].shape == expected_shape


def test_concatenate_representations_invalid_input():
    # Test that combining pairs fails
    with pytest.raises(TypeError):
        concatenate_representations((CARTESIAN_POS, SPHERICAL_POS))

    with pytest.raises(ValueError):
        concatenate_representations((CARTESIAN_POS, CARTESIAN_POS_AND_VEL))

    # Check that passing in a single object fails
    with pytest.raises(TypeError):
        concatenate_representations(CARTESIAN_POS)


def test_concatenate_representations_different_units():
    concat = concatenate_representations(
        [r.CartesianRepresentation([1, 2, 3] * unit) for unit in (u.pc, u.kpc)]
    )
    assert np.array_equal(concat.xyz, [[1, 1000], [2, 2000], [3, 3000]] * u.pc)
