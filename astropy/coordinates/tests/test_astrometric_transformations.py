# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ... import units as u
from ..builtin_frames import (ICRS, AstrometricICRS)
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
    expected = SkyCoord(*expecteddradec, frame=astrometric_frame)

    assert skycoord_inaf.separation(expected) < tolsep


@pytest.mark.xfail
def test_astrometric_functional_ra():
    #Setup
    input_ra = np.linspace(0,360,10)
    input_dec = np.linspace(-90,90,10)
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra = input_ra*u.deg,
                      dec = input_dec*u.deg,
                      distance=1.*u.kpc)
    #RA rotations

    for ra in np.linspace(0,360,24):
        # expected rotation
        expected = ICRS(ra=np.linspace(0-ra,360-ra,10)*u.deg,
                        dec=np.linspace(-90,90,10)*u.deg,
                        distance=1.*u.kpc)
        expected_xyz = expected.cartesian.xyz

        # actual transformation to the frame
        astrometric_frame = AstrometricICRS(origin_ra=ra*u.deg,
                                        origin_dec=0.0*u.deg)
        actual = icrs_coord.transform_to(astrometric_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS)
        roundtrip_xyz = roundtrip.cartesian.xyz

        # Verify
        assert allclose(actual_xyz.to(u.kpc), expected_xyz.to(u.kpc), atol=1E-5*u.kpc)
        #assert allclose(actual_xyz.to(u.kpc), roundtrip_xyz.to(u.kpc), atol=1E-5*u.kpc)
        assert allclose(icrs_coord.ra.to(u.deg), roundtrip.ra.to(u.deg), atol = 1E-5*u.deg)
        assert allclose(icrs_coord.dec.to(u.deg), roundtrip.dec.to(u.deg), atol = 1E-5*u.deg)
        assert allclose(icrs_coord.distance.to(u.kpc), roundtrip.distance.to(u.kpc), atol = 1E-5*u.kpc)

@pytest.mark.xfail
def test_astrometric_functional_dec():
    #Setup
    input_ra = np.linspace(0,360,10)
    input_dec = np.linspace(-90,90,10)
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra = input_ra*u.deg,
                      dec = input_dec*u.deg,
                      distance=1.*u.kpc)
    #Dec rotations
    #Done in xyz space because dec must be [-90,90]

    for dec in np.linspace(-90,90,13):
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
        astrometric_frame = AstrometricICRS(origin_ra=0.0*u.deg,
                                        origin_dec=dec*u.deg)
        actual = icrs_coord.transform_to(astrometric_frame)
        actual_xyz = actual.cartesian.xyz

        # back to ICRS
        roundtrip = actual.transform_to(ICRS)
        roundtrip_xyz = roundtrip.cartesian.xyz

        # Verify
        assert allclose(actual_xyz.to(u.kpc), expected_xyz.to(u.kpc), atol=1E-5*u.kpc)
        assert allclose(icrs_coord.ra.to(u.deg), roundtrip.ra.to(u.deg), atol = 1E-5*u.deg)
        assert allclose(icrs_coord.dec.to(u.deg), roundtrip.dec.to(u.deg), atol = 1E-5*u.deg)
        assert allclose(icrs_coord.distance.to(u.kpc), roundtrip.distance.to(u.kpc), atol = 1E-5*u.kpc)
        #assert allclose(actual_xyz.to(u.kpc), roundtrip_xyz.to(u.kpc), atol=1E-5*u.kpc)

@pytest.mark.xfail
def test_astrometric_functional_ra_dec():
    #Setup
    input_ra = np.linspace(0,360,10)
    input_dec = np.linspace(-90,90,10)
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra = input_ra*u.deg,
                      dec = input_dec*u.deg,
                      distance=1.*u.kpc)
    #Both rotations
    for ra in np.linspace(0,360,24):
        for dec in np.linspace(-90,90,13):
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
            astrometric_frame = AstrometricICRS(origin_ra=ra*u.deg,
                                            origin_dec=dec*u.deg)
            actual = icrs_coord.transform_to(astrometric_frame)
            actual_xyz = actual.cartesian.xyz

            # back to ICRS
            roundtrip = actual.transform_to(ICRS)
            roundtrip_xyz = roundtrip.cartesian.xyz

            # Verify
            assert allclose(actual_xyz.to(u.kpc), expected_xyz.to(u.kpc), atol=1E-5*u.kpc)
            #assert allclose(actual_xyz.to(u.kpc), roundtrip_xyz.to(u.kpc), atol=1E-5*u.kpc)
            assert allclose(icrs_coord.ra.to(u.deg), roundtrip.ra.to(u.deg), atol = 1E-5*u.deg)
            assert allclose(icrs_coord.dec.to(u.deg), roundtrip.dec.to(u.deg), atol = 1E-5*u.deg)
            assert allclose(icrs_coord.distance.to(u.kpc), roundtrip.distance.to(u.kpc), atol = 1E-5*u.kpc)
