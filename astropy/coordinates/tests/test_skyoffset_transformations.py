# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
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


def test_altaz_attribute_transforms():
    """Test transforms between AltAz frames with different attributes."""
    el1 = EarthLocation(0 * u.deg, 0 * u.deg, 0 * u.m)
    origin1 = AltAz(
        0 * u.deg, 0 * u.deg, obstime=Time("2000-01-01T12:00:00"), location=el1
    )
    frame1 = SkyOffsetFrame(origin=origin1)
    coo1 = SkyCoord(1 * u.deg, 1 * u.deg, frame=frame1)

    el2 = EarthLocation(0 * u.deg, 0 * u.deg, 0 * u.m)
    origin2 = AltAz(
        0 * u.deg, 0 * u.deg, obstime=Time("2000-01-01T11:00:00"), location=el2
    )
    frame2 = SkyOffsetFrame(origin=origin2)
    coo2 = coo1.transform_to(frame2)
    coo2_expected = [1.22522446, 0.70624298] * u.deg
    assert_allclose(
        [coo2.lon.wrap_at(180 * u.deg), coo2.lat], coo2_expected, atol=convert_precision
    )

    el3 = EarthLocation(0 * u.deg, 90 * u.deg, 0 * u.m)
    origin3 = AltAz(
        0 * u.deg, 90 * u.deg, obstime=Time("2000-01-01T12:00:00"), location=el3
    )
    frame3 = SkyOffsetFrame(origin=origin3)
    coo3 = coo2.transform_to(frame3)
    assert_allclose(
        [coo3.lon.wrap_at(180 * u.deg), coo3.lat],
        [1 * u.deg, 1 * u.deg],
        atol=convert_precision,
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
def test_skyoffset(inradec, expectedlatlon, tolsep, originradec=(45, 45) * u.deg):
    origin = ICRS(*originradec)
    skyoffset_frame = SkyOffsetFrame(origin=origin)

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


def test_skyoffset_functional_dec():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra=input_ra * u.deg, dec=input_dec * u.deg, distance=1.0 * u.kpc)
    # Dec rotations
    # Done in xyz space because dec must be [-90,90]

    for dec in np.linspace(-90, 90, 13):
        # expected rotation
        dec_rad = -np.deg2rad(dec)
        # fmt: off
        expected_x = (-np.sin(input_dec_rad) * np.sin(dec_rad) +
                       np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad))
        expected_y = (np.sin(input_ra_rad) * np.cos(input_dec_rad))
        expected_z = (np.sin(input_dec_rad) * np.cos(dec_rad) +
                      np.sin(dec_rad) * np.cos(input_ra_rad) * np.cos(input_dec_rad))
        # fmt: on
        expected = SkyCoord(
            x=expected_x,
            y=expected_y,
            z=expected_z,
            unit="kpc",
            representation_type="cartesian",
        )
        expected_xyz = expected.cartesian.xyz

        # actual transformation to the frame
        skyoffset_frame = SkyOffsetFrame(origin=ICRS(0 * u.deg, dec * u.deg))
        actual = icrs_coord.transform_to(skyoffset_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS())

        # Verify
        assert_allclose(actual_xyz, expected_xyz, atol=1e-5 * u.kpc)
        assert_allclose(icrs_coord.ra, roundtrip.ra, atol=1e-5 * u.deg)
        assert_allclose(icrs_coord.dec, roundtrip.dec, atol=1e-5 * u.deg)
        assert_allclose(icrs_coord.distance, roundtrip.distance, atol=1e-5 * u.kpc)


def test_skyoffset_functional_ra_dec():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra=input_ra * u.deg, dec=input_dec * u.deg, distance=1.0 * u.kpc)

    for ra in np.linspace(0, 360, 10):
        for dec in np.linspace(-90, 90, 5):
            # expected rotation
            dec_rad = -np.deg2rad(dec)
            ra_rad = np.deg2rad(ra)
            # fmt: off
            expected_x = (-np.sin(input_dec_rad) * np.sin(dec_rad) +
                           np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad) * np.cos(ra_rad) +
                           np.sin(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad) * np.sin(ra_rad))
            expected_y = (np.sin(input_ra_rad) * np.cos(input_dec_rad) * np.cos(ra_rad) -
                          np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.sin(ra_rad))
            expected_z = (np.sin(input_dec_rad) * np.cos(dec_rad) +
                          np.sin(dec_rad) * np.cos(ra_rad) * np.cos(input_ra_rad) * np.cos(input_dec_rad) +
                          np.sin(dec_rad) * np.sin(ra_rad) * np.sin(input_ra_rad) * np.cos(input_dec_rad))
            # fmp: on
            expected = SkyCoord(
                x=expected_x,
                y=expected_y,
                z=expected_z,
                unit='kpc',
                representation_type='cartesian',
            )
            expected_xyz = expected.cartesian.xyz

            # actual transformation to the frame
            skyoffset_frame = SkyOffsetFrame(origin=ICRS(ra * u.deg, dec * u.deg))
            actual = icrs_coord.transform_to(skyoffset_frame)
            actual_xyz = actual.cartesian.xyz

            # back to ICRS
            roundtrip = actual.transform_to(ICRS())

            # Verify
            assert_allclose(actual_xyz, expected_xyz, atol=1e-5 * u.kpc)
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


# used below in the next parametrized test
m31_sys = [ICRS, FK5, Galactic]
m31_coo = [
    (10.6847929, 41.2690650),
    (10.6847929, 41.2690650),
    (121.1744050, -21.5729360),
]
m31_dist = Distance(770, u.kpc)
convert_precision = 1 * u.arcsec
roundtrip_precision = 1e-4 * u.degree
dist_precision = 1e-9 * u.kpc

m31_params = []
for i in range(len(m31_sys)):
    for j in range(len(m31_sys)):
        if i < j:
            m31_params.append((m31_sys[i], m31_sys[j], m31_coo[i], m31_coo[j]))


@pytest.mark.parametrize(("fromsys", "tosys", "fromcoo", "tocoo"), m31_params)
def test_m31_coord_transforms(fromsys, tosys, fromcoo, tocoo):
    """
    This tests a variety of coordinate conversions for the Chandra point-source
    catalog location of M31 from NED, via SkyOffsetFrames
    """
    from_origin = fromsys(fromcoo[0] * u.deg, fromcoo[1] * u.deg, distance=m31_dist)
    from_pos = SkyOffsetFrame(1 * u.deg, 1 * u.deg, origin=from_origin)
    to_origin = tosys(tocoo[0] * u.deg, tocoo[1] * u.deg, distance=m31_dist)

    to_astroframe = SkyOffsetFrame(origin=to_origin)
    target_pos = from_pos.transform_to(to_astroframe)

    assert_allclose(
        to_origin.separation(target_pos),
        np.hypot(from_pos.lon, from_pos.lat),
        atol=convert_precision,
    )
    roundtrip_pos = target_pos.transform_to(from_pos)
    assert_allclose(
        [roundtrip_pos.lon.wrap_at(180 * u.deg), roundtrip_pos.lat],
        [1.0 * u.deg, 1.0 * u.deg],
        atol=convert_precision,
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
    origin = ICRS(45 * u.deg, 45 * u.deg)
    target = ICRS(45 * u.deg, 46 * u.deg)

    aframe = SkyOffsetFrame(origin=origin, rotation=rotation)
    trans = target.transform_to(aframe)

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
    origin = SkyCoord(45 * u.deg, 45 * u.deg)
    target = SkyCoord(45 * u.deg, 46 * u.deg)

    aframe = origin.skyoffset_frame(rotation=rotation)
    trans = target.transform_to(aframe)

    assert_allclose(
        [trans.lon.wrap_at(180 * u.deg), trans.lat], expectedlatlon, atol=1e-10 * u.deg
    )


def test_skyoffset_names():
    origin1 = ICRS(45 * u.deg, 45 * u.deg)
    aframe1 = SkyOffsetFrame(origin=origin1)
    assert type(aframe1).__name__ == "SkyOffsetICRS"

    origin2 = Galactic(45 * u.deg, 45 * u.deg)
    aframe2 = SkyOffsetFrame(origin=origin2)
    assert type(aframe2).__name__ == "SkyOffsetGalactic"


def test_skyoffset_origindata():
    origin = ICRS()
    with pytest.raises(ValueError):
        SkyOffsetFrame(origin=origin)


def test_skyoffset_lonwrap():
    origin = ICRS(45 * u.deg, 45 * u.deg)
    sc = SkyCoord(190 * u.deg, -45 * u.deg, frame=SkyOffsetFrame(origin=origin))
    assert sc.lon < 180 * u.deg

    sc2 = SkyCoord(-10 * u.deg, -45 * u.deg, frame=SkyOffsetFrame(origin=origin))
    assert sc2.lon < 180 * u.deg

    sc3 = sc.realize_frame(sc.represent_as("cartesian"))
    assert sc3.lon < 180 * u.deg

    sc4 = sc2.realize_frame(sc2.represent_as("cartesian"))
    assert sc4.lon < 180 * u.deg


def test_skyoffset_velocity():
    c = ICRS(
        ra=170.9 * u.deg,
        dec=-78.4 * u.deg,
        pm_ra_cosdec=74.4134 * u.mas / u.yr,
        pm_dec=-93.2342 * u.mas / u.yr,
    )
    skyoffset_frame = SkyOffsetFrame(origin=c)
    c_skyoffset = c.transform_to(skyoffset_frame)

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
        ra=170.9 * u.deg,
        dec=-78.4 * u.deg,
        pm_ra_cosdec=1 * u.mas / u.yr,
        pm_dec=2 * u.mas / u.yr,
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
