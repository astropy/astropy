# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Accuracy tests for ICRS transformations, primarily to/from AltAz.

"""
import pytest
import numpy as np

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation, ICRS,  CIRS, AltAz, SkyCoord)

from astropy.coordinates.tests.utils import randomly_sample_sphere
from astropy.coordinates import frame_transform_graph

ra, dec, dist = randomly_sample_sphere(200)
icrs_coords = SkyCoord(ra=ra, dec=dec, distance=dist*u.km*1e5)


@pytest.mark.parametrize('icoo', icrs_coords)
def test_icrs_consistency(icoo):
    """
    Check ICRS<->AltAz for consistency with ICRS<->CIRS<->AltAz

    The latter is extensively tested in test_intermediate_transformations.py
    """
    observer = EarthLocation(28*u.deg, 23*u.deg, height=2000.*u.km)
    obstime = Time('J2010')
    aa_frame = AltAz(obstime=obstime, location=observer)

    # check we are going direct!
    trans = frame_transform_graph.get_transform(ICRS, AltAz).transforms
    assert(len(trans) == 1)

    # check that ICRS-AltAz and ICRS->CIRS->AltAz are consistent
    aa1 = icoo.transform_to(aa_frame)
    aa2 = icoo.transform_to(CIRS()).transform_to(aa_frame)
    assert_allclose(aa1.separation_3d(aa2), 0*u.mm, atol=1*u.mm)

    # check roundtrip
    roundtrip = icoo.transform_to(aa_frame).transform_to(icoo)
    assert_allclose(roundtrip.separation_3d(icoo),  0*u.mm, atol=1*u.mm)

    # check there and back via CIRS mish-mash
    roundtrip = icoo.transform_to(aa_frame).transform_to(
        CIRS()).transform_to(icoo)
    assert_allclose(roundtrip.separation_3d(icoo),  0*u.mm, atol=1*u.mm)
