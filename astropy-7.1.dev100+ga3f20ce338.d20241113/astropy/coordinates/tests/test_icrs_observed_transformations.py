# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Accuracy tests for ICRS transformations, primarily to/from AltAz."""

import numpy as np

from astropy import units as u
from astropy.coordinates import (
    CIRS,
    ICRS,
    AltAz,
    EarthLocation,
    HADec,
    SkyCoord,
    frame_transform_graph,
    golden_spiral_grid,
)
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time


def test_icrs_altaz_consistency():
    """
    Check ICRS<->AltAz for consistency with ICRS<->CIRS<->AltAz

    The latter is extensively tested in test_intermediate_transformations.py
    """
    usph = golden_spiral_grid(200)
    dist = np.linspace(0.5, 1, len(usph)) * u.km * 1e5
    icoo = SkyCoord(ra=usph.lon, dec=usph.lat, distance=dist)

    observer = EarthLocation(28 * u.deg, 23 * u.deg, height=2000.0 * u.km)
    obstime = Time("J2010")
    aa_frame = AltAz(obstime=obstime, location=observer)

    # check we are going direct!
    trans = frame_transform_graph.get_transform(ICRS, AltAz).transforms
    assert len(trans) == 1

    # check that ICRS-AltAz and ICRS->CIRS->AltAz are consistent
    aa1 = icoo.transform_to(aa_frame)
    aa2 = icoo.transform_to(CIRS()).transform_to(aa_frame)
    assert_allclose(aa1.separation_3d(aa2), 0 * u.mm, atol=1 * u.mm)

    # check roundtrip
    roundtrip = icoo.transform_to(aa_frame).transform_to(icoo)
    assert_allclose(roundtrip.separation_3d(icoo), 0 * u.mm, atol=1 * u.mm)

    # check there and back via CIRS mish-mash
    roundtrip = icoo.transform_to(aa_frame).transform_to(CIRS()).transform_to(icoo)
    assert_allclose(roundtrip.separation_3d(icoo), 0 * u.mm, atol=1 * u.mm)


def test_icrs_hadec_consistency():
    """
    Check ICRS<->HADec for consistency with ICRS<->CIRS<->HADec
    """
    usph = golden_spiral_grid(200)
    dist = np.linspace(0.5, 1, len(usph)) * u.km * 1e5
    icoo = SkyCoord(ra=usph.lon, dec=usph.lat, distance=dist)

    observer = EarthLocation(28 * u.deg, 23 * u.deg, height=2000.0 * u.km)
    obstime = Time("J2010")
    hd_frame = HADec(obstime=obstime, location=observer)

    # check we are going direct!
    trans = frame_transform_graph.get_transform(ICRS, HADec).transforms
    assert len(trans) == 1

    # check that ICRS-HADec and ICRS->CIRS->HADec are consistent
    aa1 = icoo.transform_to(hd_frame)
    aa2 = icoo.transform_to(CIRS()).transform_to(hd_frame)
    assert_allclose(aa1.separation_3d(aa2), 0 * u.mm, atol=1 * u.mm)

    # check roundtrip
    roundtrip = icoo.transform_to(hd_frame).transform_to(icoo)
    assert_allclose(roundtrip.separation_3d(icoo), 0 * u.mm, atol=1 * u.mm)

    # check there and back via CIRS mish-mash
    roundtrip = icoo.transform_to(hd_frame).transform_to(CIRS()).transform_to(icoo)
    assert_allclose(roundtrip.separation_3d(icoo), 0 * u.mm, atol=1 * u.mm)
