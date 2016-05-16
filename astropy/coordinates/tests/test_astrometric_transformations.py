# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ... import units as u
from ..builtin_frames import ICRS, AstrometricICRS
from .. import SkyCoord
from ...tests.helper import (pytest, quantity_allclose as allclose,
                             assert_quantity_allclose as assert_allclose)


@pytest.mark.parametrize("inradec,expecteddradec, tolsep", [
    ((45, 45)*u.deg, (0, 0)*u.deg, .001*u.arcsec),
    ((45, 0)*u.deg, (0, -45)*u.deg, .001*u.arcsec),
    ((45, 90)*u.deg, (0, 45)*u.deg, .001*u.arcsec),
    ((46, 45)*u.deg, (1*np.cos(45*u.deg), 0)*u.deg, 16*u.arcsec),
    ])
def test_astrometric(inradec, expecteddradec, tolsep, originradec=(45, 45)*u.deg):
    origin = ICRS(*originradec)
    astrometric_frame = AstrometricICRS(origin=origin)

    skycoord = SkyCoord(*inradec, frame=ICRS)
    skycoord_inaf = skycoord.transform_to(astrometric_frame)
    assert hasattr(skycoord_inaf, 'dra')
    assert hasattr(skycoord_inaf, 'ddec')
    expected = SkyCoord(*expecteddradec, frame=astrometric_frame)

    assert skycoord_inaf.separation(expected) < tolsep


def test_astrometric_functional_ra():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra = input_ra*u.deg,
                      dec = input_dec*u.deg,
                      distance=1.*u.kpc)

    for ra in np.linspace(0,360,24):
        # expected rotation
        expected = ICRS(ra=np.linspace(0-ra, 360-ra, 12)[1:-1]*u.deg,
                        dec=np.linspace(-90, 90, 12)[1:-1]*u.deg,
                        distance=1.*u.kpc)
        expected_xyz = expected.cartesian.xyz

        # actual transformation to the frame
        astrometric_frame = AstrometricICRS(origin=ICRS(ra*u.deg, 0*u.deg))
        actual = icrs_coord.transform_to(astrometric_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS)
        roundtrip_xyz = roundtrip.cartesian.xyz

        # Verify
        assert_allclose(actual_xyz, expected_xyz, atol=1E-5*u.kpc)
        assert_allclose(icrs_coord.ra, roundtrip.ra, atol = 1E-5*u.deg)
        assert_allclose(icrs_coord.dec, roundtrip.dec, atol = 1E-5*u.deg)
        assert_allclose(icrs_coord.distance, roundtrip.distance, atol = 1E-5*u.kpc)


def test_astrometric_functional_dec():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra = input_ra*u.deg,
                      dec = input_dec*u.deg,
                      distance=1.*u.kpc)
    #Dec rotations
    #Done in xyz space because dec must be [-90,90]

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
        astrometric_frame = AstrometricICRS(origin=ICRS(0*u.deg, dec*u.deg))
        actual = icrs_coord.transform_to(astrometric_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS)
        roundtrip_xyz = roundtrip.cartesian.xyz

        # Verify
        assert_allclose(actual_xyz, expected_xyz, atol=1E-5*u.kpc)
        assert_allclose(icrs_coord.ra, roundtrip.ra, atol = 1E-5*u.deg)
        assert_allclose(icrs_coord.dec, roundtrip.dec, atol = 1E-5*u.deg)
        assert_allclose(icrs_coord.distance, roundtrip.distance, atol = 1E-5*u.kpc)


def test_astrometric_functional_ra_dec():
    # we do the 12)[1:-1] business because sometimes machine precision issues
    # lead to results that are either ~0 or ~360, which mucks up the final
    # comparison and leads to spurious failures.  So this just avoids that by
    # staying away from the edges
    input_ra = np.linspace(0, 360, 12)[1:-1]
    input_dec = np.linspace(-90, 90, 12)[1:-1]
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra = input_ra*u.deg,
                      dec = input_dec*u.deg,
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
            astrometric_frame = AstrometricICRS(origin=ICRS(ra*u.deg, dec*u.deg))
            actual = icrs_coord.transform_to(astrometric_frame)
            actual_xyz = actual.cartesian.xyz

            # back to ICRS
            roundtrip = actual.transform_to(ICRS)
            roundtrip_xyz = roundtrip.cartesian.xyz

            # Verify
            assert_allclose(actual_xyz, expected_xyz, atol=1E-5*u.kpc)
            assert_allclose(icrs_coord.ra, roundtrip.ra, atol = 1E-4*u.deg)
            assert_allclose(icrs_coord.dec, roundtrip.dec, atol = 1E-5*u.deg)
            assert_allclose(icrs_coord.distance, roundtrip.distance, atol = 1E-5*u.kpc)

def test_skycoord_astrometric_frame():
    m31 = SkyCoord(10.6847083, 41.26875, frame='icrs', unit=u.deg)
    m33 = SkyCoord(23.4621, 30.6599417, frame='icrs', unit=u.deg)

    m31_astro = m31.astrometric_frame()
    m31_in_m31 = m31.transform_to(m31_astro)
    m33_in_m31 = m33.transform_to(m31_astro)

    assert_allclose([m31_in_m31.dra, m31_in_m31.ddec], [0, 0]*u.deg, atol=1e-10*u.deg)
    assert_allclose([m33_in_m31.dra, m33_in_m31.ddec], [11.13135175, -9.79084759]*u.deg)

    assert_allclose(m33.separation(m31),
                    np.hypot(m33_in_m31.dra, m33_in_m31.ddec),
                    atol=.1*u.deg)
@pytest.mark.parametrize("rotation, expecteddradec", [
    (0*u.deg, [0, 1]*u.deg),
    (180*u.deg, [0, -1]*u.deg),
    (90*u.deg, [-1, 0]*u.deg),
    (-90*u.deg, [1, 0]*u.deg)
    ])
def test_rotation(rotation, expecteddradec):
    origin = ICRS(45*u.deg, 45*u.deg)
    target = ICRS(45*u.deg, 46*u.deg)

    aframe = AstrometricICRS(origin=origin, rotation=rotation)
    trans = target.transform_to(aframe)

    assert_allclose([trans.dra.wrap_at(180*u.deg), trans.ddec],
                    expecteddradec, atol=1e-10*u.deg)
