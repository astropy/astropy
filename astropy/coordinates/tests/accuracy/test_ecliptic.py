# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Accuracy tests for Ecliptic coordinate systems.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ....tests.helper import quantity_allclose
from .... import units as u
from ...builtin_frames import FK5, ICRS, GCRS, GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic
from ....constants import R_sun, R_earth


def test_against_pytpm_doc_example():
    """
    Check that Astropy's Ecliptic systems give answers consistent with pyTPM

    Currently this is only testing against the example given in the pytpm docs
    """
    fk5_in = FK5('12h22m54.899s', '15d49m20.57s', equinox='J2000')
    pytpm_out = BarycentricTrueEcliptic(lon=178.78256462*u.deg,
                                        lat=16.7597002513*u.deg,
                                        equinox='J2000')
    astropy_out = fk5_in.transform_to(pytpm_out)

    # we check w/i 1 arcmin because there are some subtle differences in pyTPM's ecl definition
    assert pytpm_out.separation(astropy_out) < (1*u.arcmin)


def test_ecliptic_heliobary():
    """
    Check that the ecliptic transformations for heliocentric and barycentric
    at least more or less make sense
    """
    icrs = ICRS(1*u.deg, 2*u.deg, distance=1.5*R_sun)

    bary = icrs.transform_to(BarycentricTrueEcliptic)
    helio = icrs.transform_to(HeliocentricTrueEcliptic)

    # make sure there's a sizable distance shift - in 3d hundreds of km, but
    # this is 1D so we allow it to be somewhat smaller
    assert np.abs(bary.distance - helio.distance) > 1*u.km

    # now make something that's got the location of helio but in bary's frame.
    # this is a convenience to allow `separation` to work as expected
    helio_in_bary_frame = bary.realize_frame(helio.cartesian)
    assert bary.separation(helio_in_bary_frame) > 1*u.arcmin


def test_ecl_geo():
    """
    Check that the geocentric version at least gets well away from GCRS.  For a
    true "accuracy" test we need a comparison dataset that is similar to the
    geocentric/GCRS comparison we want to do here.  Contributions welcome!
    """
    gcrs = GCRS(10*u.deg, 20*u.deg, distance=1.5*R_earth)
    gecl = gcrs.transform_to(GeocentricTrueEcliptic)

    assert quantity_allclose(gecl.distance, gcrs.distance)


def test_arraytransforms():
    """
    Test that transforms to/from ecliptic coordinates work on array coordinates
    (not testing for accuracy.)
    """
    ra = np.ones((3,), dtype=float) * u.deg
    dec = np.zeros((3,), dtype=float) * u.deg
    distance = np.ones((3, ), dtype=float) * u.au

    test_icrs = ICRS(ra=ra, dec=dec, distance=distance)

    helio_arr = test_icrs.transform_to(HeliocentricTrueEcliptic)
    assert len(helio_arr) == 3
    #really, if hasn't thrown an exception, it passes

    bary_arr = test_icrs.transform_to(BarycentricTrueEcliptic)
    assert len(bary_arr) == 3
    #really, if hasn't thrown an exception, it passes
