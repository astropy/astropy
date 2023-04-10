# Licensed under a 3-clause BSD style license - see LICENSE.rst

from itertools import combinations

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.coordinates.builtin_frames import (
    FK5,
    ICRS,
    AltAz,
    Galactic,
    SkyOffsetFrame,
)
from astropy.coordinates.distances import Distance
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time

CONVERT_PRECISION = 1 * u.arcsec
ICRS_45_45 = SkyCoord(ra=45 * u.deg, dec=45 * u.deg, frame=ICRS())
M31_DISTANCE = Distance(770 * u.kpc)
POSITION_ON_SKY = {"ra": 36.4 * u.deg, "dec": -55.8 * u.deg}
DISTANCE = {"distance": 150 * u.pc}
PROPER_MOTION = {"pm_ra_cosdec": -21.2 * u.mas / u.yr, "pm_dec": 17.1 * u.mas / u.yr}


@pytest.fixture(scope="module")
def icrs_coords_with_trig_values():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges.  Explicit conversion to radians in the trig
    # functions is needed so that output would be a bare `ndarray`, not a `Quantity`.
    icrs_coord = ICRS(
        ra=np.linspace(0, 360, 12)[1:-1] * u.deg,
        dec=np.linspace(-90, 90, 12)[1:-1] * u.deg,
        distance=1.0 * u.kpc,
    )
    return (
        icrs_coord,
        np.sin(icrs_coord.dec.rad),
        np.cos(icrs_coord.dec.rad),
        np.sin(icrs_coord.ra.rad),
        np.cos(icrs_coord.ra.rad),
    )


def test_altaz_attribute_transforms():
    """Test transforms between AltAz frames with different attributes."""
    el1 = EarthLocation(0 * u.deg, 0 * u.deg, 0 * u.m)
    origin1 = AltAz(
        0 * u.deg, 0 * u.deg, obstime=Time("2000-01-01T12:00:00"), location=el1
    )
    coo1 = SkyCoord(1 * u.deg, 1 * u.deg, frame=SkyOffsetFrame(origin=origin1))

    origin2 = AltAz(
        0 * u.deg, 0 * u.deg, obstime=Time("2000-01-01T11:00:00"), location=el1
    )
    coo2 = coo1.transform_to(SkyOffsetFrame(origin=origin2))
    assert_allclose(
        [coo2.lon.wrap_at(180 * u.deg), coo2.lat],
        [1.22522446, 0.70624298] * u.deg,
        atol=CONVERT_PRECISION,
    )

    el3 = EarthLocation(0 * u.deg, 90 * u.deg, 0 * u.m)
    origin3 = AltAz(
        0 * u.deg, 90 * u.deg, obstime=Time("2000-01-01T12:00:00"), location=el3
    )
    coo3 = coo2.transform_to(SkyOffsetFrame(origin=origin3))
    assert_allclose(
        [coo3.lon.wrap_at(180 * u.deg), coo3.lat],
        [1 * u.deg, 1 * u.deg],
        atol=CONVERT_PRECISION,
    )


@pytest.mark.parametrize(
    "inradec,expectedlatlon, tolsep",
    [
        ((45, 45) * u.deg, (0, 0) * u.deg, 0.001 * u.arcsec),
        ((45, 0) * u.deg, (0, -45) * u.deg, 0.001 * u.arcsec),
        ((45, 90) * u.deg, (0, 45) * u.deg, 0.001 * u.arcsec),
        ((46, 45) * u.deg, (1 * np.cos(45 * u.deg), 0) * u.deg, 16 * u.arcsec),
    ],
)
def test_skyoffset(inradec, expectedlatlon, tolsep):
    skyoffset_frame = SkyOffsetFrame(origin=ICRS_45_45)

    skycoord = SkyCoord(*inradec, frame=ICRS)
    skycoord_inaf = skycoord.transform_to(skyoffset_frame)
    assert hasattr(skycoord_inaf, "lon")
    assert hasattr(skycoord_inaf, "lat")
    expected = SkyCoord(*expectedlatlon, frame=skyoffset_frame)

    assert skycoord_inaf.separation(expected) < tolsep
    # Check we can also transform back (regression test for gh-11254).
    roundtrip = skycoord_inaf.transform_to(ICRS())
    assert roundtrip.separation(skycoord) < 1 * u.uas


def test_skyoffset_functional_ra():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    icrs_coord = ICRS(ra=input_ra * u.deg, dec=input_dec * u.deg, distance=1.0 * u.kpc)

    for ra in np.linspace(0, 360, 24):
        # expected rotation
        expected = ICRS(
            ra=np.linspace(0 - ra, 360 - ra, 12)[1:-1] * u.deg,
            dec=np.linspace(-90, 90, 12)[1:-1] * u.deg,
            distance=1.0 * u.kpc,
        )
        expected_xyz = expected.cartesian.xyz

        # actual transformation to the frame
        skyoffset_frame = SkyOffsetFrame(origin=ICRS(ra * u.deg, 0 * u.deg))
        actual = icrs_coord.transform_to(skyoffset_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS())
        roundtrip_xyz = roundtrip.cartesian.xyz

        # Verify
        assert_allclose(actual_xyz, expected_xyz, atol=1e-5 * u.kpc)
        assert_allclose(icrs_coord.ra, roundtrip.ra, atol=1e-5 * u.deg)
        assert_allclose(icrs_coord.dec, roundtrip.dec, atol=1e-5 * u.deg)
        assert_allclose(icrs_coord.distance, roundtrip.distance, atol=1e-5 * u.kpc)


@pytest.mark.parametrize("dec", Angle(np.linspace(-90, 90, 13), u.deg))
def test_skyoffset_functional_dec(dec, icrs_coords_with_trig_values):
    icrs_coord, sin_dec_i, cos_dec_i, sin_ra_i, cos_ra_i = icrs_coords_with_trig_values
    # Dec rotations
    # Done in xyz space because dec must be [-90,90]

    # expected rotation
    sin_dec = np.sin(-dec.rad)
    cos_dec = np.cos(dec.rad)
    expected = SkyCoord(
        x=-sin_dec_i * sin_dec + cos_ra_i * cos_dec_i * cos_dec,
        y=sin_ra_i * cos_dec_i,
        z=sin_dec_i * cos_dec + sin_dec * cos_ra_i * cos_dec_i,
        unit="kpc",
        representation_type="cartesian",
    )

    # actual transformation to the frame
    actual = icrs_coord.transform_to(SkyOffsetFrame(origin=ICRS(0 * u.deg, dec)))

    # back to ICRS
    roundtrip = actual.transform_to(ICRS())

    # Verify
    assert_allclose(actual.cartesian.xyz, expected.cartesian.xyz, atol=1e-5 * u.kpc)
    assert_allclose(icrs_coord.ra, roundtrip.ra, atol=1e-5 * u.deg)
    assert_allclose(icrs_coord.dec, roundtrip.dec, atol=1e-5 * u.deg)
    assert_allclose(icrs_coord.distance, roundtrip.distance, atol=1e-5 * u.kpc)


@pytest.mark.parametrize("ra", Angle(np.linspace(0, 360, 10), u.deg))
@pytest.mark.parametrize("dec", Angle(np.linspace(-90, 90, 5), u.deg))
def test_skyoffset_functional_ra_dec(ra, dec, icrs_coords_with_trig_values):
    icrs_coord, sin_dec_i, cos_dec_i, sin_ra_i, cos_ra_i = icrs_coords_with_trig_values
    cos_dec = np.cos(dec.rad)
    sin_dec = np.sin(-dec.rad)
    cos_ra = np.cos(ra.rad)
    sin_ra = np.sin(ra.rad)
    # expected rotation
    expected = SkyCoord(
        x=(
            -sin_dec_i * sin_dec
            + cos_ra_i * cos_dec_i * cos_dec * cos_ra
            + sin_ra_i * cos_dec_i * cos_dec * sin_ra
        ),
        y=sin_ra_i * cos_dec_i * cos_ra - cos_ra_i * cos_dec_i * sin_ra,
        z=(
            sin_dec_i * cos_dec
            + sin_dec * cos_ra * cos_ra_i * cos_dec_i
            + sin_dec * sin_ra * sin_ra_i * cos_dec_i
        ),
        unit="kpc",
        representation_type="cartesian",
    )

    # actual transformation to the frame
    actual = icrs_coord.transform_to(SkyOffsetFrame(origin=ICRS(ra, dec)))

    # back to ICRS
    roundtrip = actual.transform_to(ICRS())

    # Verify
    assert_allclose(actual.cartesian.xyz, expected.cartesian.xyz, atol=1e-5 * u.kpc)
    assert_allclose(icrs_coord.ra, roundtrip.ra, atol=1e-4 * u.deg)
    assert_allclose(icrs_coord.dec, roundtrip.dec, atol=1e-5 * u.deg)
    assert_allclose(icrs_coord.distance, roundtrip.distance, atol=1e-5 * u.kpc)


def test_skycoord_skyoffset_frame():
    m31 = SkyCoord(10.6847083, 41.26875, frame="icrs", unit=u.deg)
    m33 = SkyCoord(23.4621, 30.6599417, frame="icrs", unit=u.deg)

    m31_astro = m31.skyoffset_frame()
    m31_in_m31 = m31.transform_to(m31_astro)
    m33_in_m31 = m33.transform_to(m31_astro)

    assert_allclose(
        [m31_in_m31.lon, m31_in_m31.lat], [0, 0] * u.deg, atol=1e-10 * u.deg
    )
    assert_allclose(
        [m33_in_m31.lon, m33_in_m31.lat], [11.13135175, -9.79084759] * u.deg
    )

    assert_allclose(
        m33.separation(m31), np.hypot(m33_in_m31.lon, m33_in_m31.lat), atol=0.1 * u.deg
    )


@pytest.mark.parametrize(
    "from_origin,to_origin",
    combinations(
        (
            ICRS(10.6847929 * u.deg, 41.2690650 * u.deg, M31_DISTANCE),
            FK5(10.6847929 * u.deg, 41.2690650 * u.deg, M31_DISTANCE),
            Galactic(121.1744050 * u.deg, -21.5729360 * u.deg, M31_DISTANCE),
        ),
        r=2,
    ),
)
def test_m31_coord_transforms(from_origin, to_origin):
    """
    This tests a variety of coordinate conversions for the Chandra point-source
    catalog location of M31 from NED, via SkyOffsetFrames
    """
    from_pos = SkyOffsetFrame(1 * u.deg, 1 * u.deg, origin=from_origin)

    to_astroframe = SkyOffsetFrame(origin=to_origin)
    target_pos = from_pos.transform_to(to_astroframe)

    assert_allclose(
        to_origin.separation(target_pos),
        np.hypot(from_pos.lon, from_pos.lat),
        atol=CONVERT_PRECISION,
    )
    roundtrip_pos = target_pos.transform_to(from_pos)
    assert_allclose(
        [roundtrip_pos.lon.wrap_at(180 * u.deg), roundtrip_pos.lat],
        [1.0 * u.deg, 1.0 * u.deg],
        atol=CONVERT_PRECISION,
    )


@pytest.mark.parametrize(
    "rotation, expectedlatlon",
    [
        (0 * u.deg, [0, 1] * u.deg),
        (180 * u.deg, [0, -1] * u.deg),
        (90 * u.deg, [-1, 0] * u.deg),
        (-90 * u.deg, [1, 0] * u.deg),
    ],
)
def test_rotation(rotation, expectedlatlon):
    target = ICRS(45 * u.deg, 46 * u.deg)
    trans = target.transform_to(SkyOffsetFrame(origin=ICRS_45_45, rotation=rotation))
    assert_allclose(
        [trans.lon.wrap_at(180 * u.deg), trans.lat], expectedlatlon, atol=1e-10 * u.deg
    )


@pytest.mark.parametrize(
    "rotation, expectedlatlon",
    [
        (0 * u.deg, [0, 1] * u.deg),
        (180 * u.deg, [0, -1] * u.deg),
        (90 * u.deg, [-1, 0] * u.deg),
        (-90 * u.deg, [1, 0] * u.deg),
    ],
)
def test_skycoord_skyoffset_frame_rotation(rotation, expectedlatlon):
    """Test if passing a rotation argument via SkyCoord works"""
    target = SkyCoord(45 * u.deg, 46 * u.deg)
    trans = target.transform_to(ICRS_45_45.skyoffset_frame(rotation=rotation))
    assert_allclose(
        [trans.lon.wrap_at(180 * u.deg), trans.lat], expectedlatlon, atol=1e-10 * u.deg
    )


def test_skyoffset_names():
    aframe1 = SkyOffsetFrame(origin=ICRS_45_45)
    assert type(aframe1).__name__ == "SkyOffsetICRS"

    aframe2 = SkyOffsetFrame(origin=Galactic(45 * u.deg, 45 * u.deg))
    assert type(aframe2).__name__ == "SkyOffsetGalactic"


def test_skyoffset_origindata():
    origin = ICRS()
    with pytest.raises(ValueError):
        SkyOffsetFrame(origin=origin)


@pytest.mark.parametrize("lon", (190, -10) * u.deg)
def test_skyoffset_lonwrap(lon):
    sc = SkyCoord(lon=lon, lat=-45 * u.deg, frame=SkyOffsetFrame(origin=ICRS_45_45))
    assert sc.lon < 180 * u.deg
    assert sc.realize_frame(sc.represent_as("cartesian")).lon < 180 * u.deg


def test_skyoffset_velocity():
    c = ICRS(**POSITION_ON_SKY, **PROPER_MOTION)
    c_skyoffset = c.transform_to(SkyOffsetFrame(origin=c))

    assert_allclose(c_skyoffset.pm_lon_coslat, c.pm_ra_cosdec)
    assert_allclose(c_skyoffset.pm_lat, c.pm_dec)


@pytest.mark.parametrize(
    "rotation, expectedpmlonlat",
    [
        (0 * u.deg, [1, 2] * u.mas / u.yr),
        (45 * u.deg, [-(2**-0.5), 3 * 2**-0.5] * u.mas / u.yr),
        (90 * u.deg, [-2, 1] * u.mas / u.yr),
        (180 * u.deg, [-1, -2] * u.mas / u.yr),
        (-90 * u.deg, [2, -1] * u.mas / u.yr),
    ],
)
def test_skyoffset_velocity_rotation(rotation, expectedpmlonlat):
    sc = SkyCoord(
        **POSITION_ON_SKY, pm_ra_cosdec=1 * u.mas / u.yr, pm_dec=2 * u.mas / u.yr
    )

    c_skyoffset0 = sc.transform_to(sc.skyoffset_frame(rotation=rotation))
    assert_allclose(c_skyoffset0.pm_lon_coslat, expectedpmlonlat[0])
    assert_allclose(c_skyoffset0.pm_lat, expectedpmlonlat[1])


def test_skyoffset_two_frames_interfering():
    """Regression test for gh-11277, where it turned out that the
    origin argument validation from one SkyOffsetFrame could interfere
    with that of another.

    Note that this example brought out a different bug than that at the
    top of gh-11277, viz., that an attempt was made to set origin on a SkyCoord
    when it should just be stay as part of the SkyOffsetFrame.
    """
    # Example adapted from @bmerry's minimal example at
    # https://github.com/astropy/astropy/issues/11277#issuecomment-825492335
    altaz_frame = AltAz(
        obstime=Time("2020-04-22T13:00:00Z"), location=EarthLocation(18, -30)
    )
    target = SkyCoord(alt=70 * u.deg, az=150 * u.deg, frame=altaz_frame)
    dirs_altaz_offset = SkyCoord(
        lon=[-0.02, 0.01, 0.0, 0.0, 0.0] * u.rad,
        lat=[0.0, 0.2, 0.0, -0.3, 0.1] * u.rad,
        frame=target.skyoffset_frame(),
    )
    dirs_altaz = dirs_altaz_offset.transform_to(altaz_frame)
    dirs_icrs = dirs_altaz.transform_to(ICRS())
    target_icrs = target.transform_to(ICRS())
    # The line below was almost guaranteed to fail.
    dirs_icrs.transform_to(target_icrs.skyoffset_frame())
