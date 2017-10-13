# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
import numpy as np

from ... import units as u
from ..distances import Distance
from ..builtin_frames import ICRS, FK5, Galactic, AltAz, SkyOffsetFrame
from .. import SkyCoord, EarthLocation
from ...time import Time
from ...tests.helper import assert_quantity_allclose as assert_allclose

from ...extern.six.moves import range


@pytest.mark.parametrize("inradec,expectedlatlon, tolsep", [
    ((45, 45)*u.deg, (0, 0)*u.deg, .001*u.arcsec),
    ((45, 0)*u.deg, (0, -45)*u.deg, .001*u.arcsec),
    ((45, 90)*u.deg, (0, 45)*u.deg, .001*u.arcsec),
    ((46, 45)*u.deg, (1*np.cos(45*u.deg), 0)*u.deg, 16*u.arcsec),
    ])
def test_skyoffset(inradec, expectedlatlon, tolsep, originradec=(45, 45)*u.deg):
    origin = ICRS(*originradec)
    skyoffset_frame = SkyOffsetFrame(origin=origin)

    skycoord = SkyCoord(*inradec, frame=ICRS)
    skycoord_inaf = skycoord.transform_to(skyoffset_frame)
    assert hasattr(skycoord_inaf, 'lon')
    assert hasattr(skycoord_inaf, 'lat')
    expected = SkyCoord(*expectedlatlon, frame=skyoffset_frame)

    assert skycoord_inaf.separation(expected) < tolsep


def test_skyoffset_functional_ra():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    icrs_coord = ICRS(ra=input_ra*u.deg,
                      dec=input_dec*u.deg,
                      distance=1.*u.kpc)

    for ra in np.linspace(0, 360, 24):
        # expected rotation
        expected = ICRS(ra=np.linspace(0-ra, 360-ra, 12)[1:-1]*u.deg,
                        dec=np.linspace(-90, 90, 12)[1:-1]*u.deg,
                        distance=1.*u.kpc)
        expected_xyz = expected.cartesian.xyz

        # actual transformation to the frame
        skyoffset_frame = SkyOffsetFrame(origin=ICRS(ra*u.deg, 0*u.deg))
        actual = icrs_coord.transform_to(skyoffset_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS)
        roundtrip_xyz = roundtrip.cartesian.xyz

        # Verify
        assert_allclose(actual_xyz, expected_xyz, atol=1E-5*u.kpc)
        assert_allclose(icrs_coord.ra, roundtrip.ra, atol=1E-5*u.deg)
        assert_allclose(icrs_coord.dec, roundtrip.dec, atol=1E-5*u.deg)
        assert_allclose(icrs_coord.distance, roundtrip.distance, atol=1E-5*u.kpc)


def test_skyoffset_functional_dec():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra=input_ra*u.deg,
                      dec=input_dec*u.deg,
                      distance=1.*u.kpc)
    # Dec rotations
    # Done in xyz space because dec must be [-90,90]

    for dec in np.linspace(-90, 90, 13):
        # expected rotation
        dec_rad = -np.deg2rad(dec)
        expected_x = (-np.sin(input_dec_rad) * np.sin(dec_rad) +
                       np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad))
        expected_y = (np.sin(input_ra_rad) * np.cos(input_dec_rad))
        expected_z = (np.sin(input_dec_rad) * np.cos(dec_rad) +
                      np.sin(dec_rad) * np.cos(input_ra_rad) * np.cos(input_dec_rad))
        expected = SkyCoord(x=expected_x,
                            y=expected_y,
                            z=expected_z, unit='kpc', representation='cartesian')
        expected_xyz = expected.cartesian.xyz

        # actual transformation to the frame
        skyoffset_frame = SkyOffsetFrame(origin=ICRS(0*u.deg, dec*u.deg))
        actual = icrs_coord.transform_to(skyoffset_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS)

        # Verify
        assert_allclose(actual_xyz, expected_xyz, atol=1E-5*u.kpc)
        assert_allclose(icrs_coord.ra, roundtrip.ra, atol=1E-5*u.deg)
        assert_allclose(icrs_coord.dec, roundtrip.dec, atol=1E-5*u.deg)
        assert_allclose(icrs_coord.distance, roundtrip.distance, atol=1E-5*u.kpc)


def test_skyoffset_functional_ra_dec():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra=input_ra*u.deg,
                      dec=input_dec*u.deg,
                      distance=1.*u.kpc)

    for ra in np.linspace(0, 360, 10):
        for dec in np.linspace(-90, 90, 5):
            # expected rotation
            dec_rad = -np.deg2rad(dec)
            ra_rad = np.deg2rad(ra)
            expected_x = (-np.sin(input_dec_rad) * np.sin(dec_rad) +
                           np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad) * np.cos(ra_rad) +
                           np.sin(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad) * np.sin(ra_rad))
            expected_y = (np.sin(input_ra_rad) * np.cos(input_dec_rad) * np.cos(ra_rad) -
                          np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.sin(ra_rad))
            expected_z = (np.sin(input_dec_rad) * np.cos(dec_rad) +
                          np.sin(dec_rad) * np.cos(ra_rad) * np.cos(input_ra_rad) * np.cos(input_dec_rad) +
                          np.sin(dec_rad) * np.sin(ra_rad) * np.sin(input_ra_rad) * np.cos(input_dec_rad))
            expected = SkyCoord(x=expected_x,
                                y=expected_y,
                                z=expected_z, unit='kpc', representation='cartesian')
            expected_xyz = expected.cartesian.xyz

            # actual transformation to the frame
            skyoffset_frame = SkyOffsetFrame(origin=ICRS(ra*u.deg, dec*u.deg))
            actual = icrs_coord.transform_to(skyoffset_frame)
            actual_xyz = actual.cartesian.xyz

            # back to ICRS
            roundtrip = actual.transform_to(ICRS)

            # Verify
            assert_allclose(actual_xyz, expected_xyz, atol=1E-5*u.kpc)
            assert_allclose(icrs_coord.ra, roundtrip.ra, atol=1E-4*u.deg)
            assert_allclose(icrs_coord.dec, roundtrip.dec, atol=1E-5*u.deg)
            assert_allclose(icrs_coord.distance, roundtrip.distance, atol=1E-5*u.kpc)


def test_skycoord_skyoffset_frame():
    m31 = SkyCoord(10.6847083, 41.26875, frame='icrs', unit=u.deg)
    m33 = SkyCoord(23.4621, 30.6599417, frame='icrs', unit=u.deg)

    m31_astro = m31.skyoffset_frame()
    m31_in_m31 = m31.transform_to(m31_astro)
    m33_in_m31 = m33.transform_to(m31_astro)

    assert_allclose([m31_in_m31.lon, m31_in_m31.lat], [0, 0]*u.deg, atol=1e-10*u.deg)
    assert_allclose([m33_in_m31.lon, m33_in_m31.lat], [11.13135175, -9.79084759]*u.deg)

    assert_allclose(m33.separation(m31),
                    np.hypot(m33_in_m31.lon, m33_in_m31.lat),
                    atol=.1*u.deg)


# used below in the next parametrized test
m31_sys = [ICRS, FK5, Galactic]
m31_coo = [(10.6847929, 41.2690650), (10.6847929, 41.2690650), (121.1744050, -21.5729360)]
m31_dist = Distance(770, u.kpc)
convert_precision = 1 * u.arcsec
roundtrip_precision = 1e-4 * u.degree
dist_precision = 1e-9 * u.kpc

m31_params = []
for i in range(len(m31_sys)):
    for j in range(len(m31_sys)):
        if i < j:
            m31_params.append((m31_sys[i], m31_sys[j], m31_coo[i], m31_coo[j]))


@pytest.mark.parametrize(('fromsys', 'tosys', 'fromcoo', 'tocoo'), m31_params)
def test_m31_coord_transforms(fromsys, tosys, fromcoo, tocoo):
    """
    This tests a variety of coordinate conversions for the Chandra point-source
    catalog location of M31 from NED, via SkyOffsetFrames
    """
    from_origin = fromsys(fromcoo[0]*u.deg, fromcoo[1]*u.deg,
                          distance=m31_dist)
    from_pos = SkyOffsetFrame(1*u.deg, 1*u.deg, origin=from_origin)
    to_origin = tosys(tocoo[0]*u.deg, tocoo[1]*u.deg, distance=m31_dist)

    to_astroframe = SkyOffsetFrame(origin=to_origin)
    target_pos = from_pos.transform_to(to_astroframe)

    assert_allclose(to_origin.separation(target_pos),
                    np.hypot(from_pos.lon, from_pos.lat),
                    atol=convert_precision)
    roundtrip_pos = target_pos.transform_to(from_pos)
    assert_allclose([roundtrip_pos.lon.wrap_at(180*u.deg), roundtrip_pos.lat],
                    [1.0*u.deg, 1.0*u.deg], atol=convert_precision)


def test_altaz_attribute_transforms():
    """Test transforms between AltAz frames with different attributes."""
    el1 = EarthLocation(0*u.deg, 0*u.deg, 0*u.m)
    origin1 = AltAz(0 * u.deg, 0*u.deg, obstime=Time("2000-01-01T12:00:00"),
                    location=el1)
    frame1 = SkyOffsetFrame(origin=origin1)
    coo1 = SkyCoord(1 * u.deg, 1 * u.deg, frame=frame1)

    el2 = EarthLocation(0*u.deg, 0*u.deg, 0*u.m)
    origin2 = AltAz(0 * u.deg, 0*u.deg, obstime=Time("2000-01-01T11:00:00"),
                    location=el2)
    frame2 = SkyOffsetFrame(origin=origin2)
    coo2 = coo1.transform_to(frame2)
    coo2_expected = [1.22522446, 0.70624298] * u.deg
    assert_allclose([coo2.lon.wrap_at(180*u.deg), coo2.lat],
                    coo2_expected, atol=convert_precision)

    el3 = EarthLocation(0*u.deg, 90*u.deg, 0*u.m)
    origin3 = AltAz(0 * u.deg, 90*u.deg, obstime=Time("2000-01-01T12:00:00"),
                    location=el3)
    frame3 = SkyOffsetFrame(origin=origin3)
    coo3 = coo2.transform_to(frame3)
    assert_allclose([coo3.lon.wrap_at(180*u.deg), coo3.lat],
                    [1*u.deg, 1*u.deg], atol=convert_precision)


@pytest.mark.parametrize("rotation, expectedlatlon", [
    (0*u.deg, [0, 1]*u.deg),
    (180*u.deg, [0, -1]*u.deg),
    (90*u.deg, [-1, 0]*u.deg),
    (-90*u.deg, [1, 0]*u.deg)
    ])
def test_rotation(rotation, expectedlatlon):
    origin = ICRS(45*u.deg, 45*u.deg)
    target = ICRS(45*u.deg, 46*u.deg)

    aframe = SkyOffsetFrame(origin=origin, rotation=rotation)
    trans = target.transform_to(aframe)

    assert_allclose([trans.lon.wrap_at(180*u.deg), trans.lat],
                    expectedlatlon, atol=1e-10*u.deg)


@pytest.mark.parametrize("rotation, expectedlatlon", [
    (0*u.deg, [0, 1]*u.deg),
    (180*u.deg, [0, -1]*u.deg),
    (90*u.deg, [-1, 0]*u.deg),
    (-90*u.deg, [1, 0]*u.deg)
    ])
def test_skycoord_skyoffset_frame_rotation(rotation, expectedlatlon):
    """Test if passing a rotation argument via SkyCoord works"""
    origin = SkyCoord(45*u.deg, 45*u.deg)
    target = SkyCoord(45*u.deg, 46*u.deg)

    aframe = origin.skyoffset_frame(rotation=rotation)
    trans = target.transform_to(aframe)

    assert_allclose([trans.lon.wrap_at(180*u.deg), trans.lat],
                    expectedlatlon, atol=1e-10*u.deg)


def test_skyoffset_names():
    origin1 = ICRS(45*u.deg, 45*u.deg)
    aframe1 = SkyOffsetFrame(origin=origin1)
    assert type(aframe1).__name__ == 'SkyOffsetICRS'

    origin2 = Galactic(45*u.deg, 45*u.deg)
    aframe2 = SkyOffsetFrame(origin=origin2)
    assert type(aframe2).__name__ == 'SkyOffsetGalactic'


def test_skyoffset_origindata():
    origin = ICRS()
    with pytest.raises(ValueError):
        SkyOffsetFrame(origin=origin)


def test_skyoffset_lonwrap():
    origin = ICRS(45*u.deg, 45*u.deg)
    sc = SkyCoord(190*u.deg, -45*u.deg, frame=SkyOffsetFrame(origin=origin))
    assert sc.lon < 180 * u.deg


def test_skyoffset_velerr():
    # TODO: remove this when the SkyOffsetFrame's support velocities
    origin = ICRS(45*u.deg, 45*u.deg)
    originwvel = ICRS(45*u.deg, 45*u.deg, radial_velocity=1*u.km/u.s)

    SkyOffsetFrame(origin=origin)
    with pytest.raises(NotImplementedError):
        SkyOffsetFrame(origin=originwvel)
    SkyOffsetFrame(origin.data, origin=origin)
    with pytest.raises(NotImplementedError):
        SkyOffsetFrame(originwvel.data, origin=origin)
    with pytest.raises(NotImplementedError):
        SkyOffsetFrame(origin.data, origin=originwvel)
    with pytest.raises(NotImplementedError):
        SkyOffsetFrame(originwvel.data, origin=originwvel)
