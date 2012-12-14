# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains tests for the name resolve convenience module.
"""

# Standard library
import time
import urllib,urllib2

# Third party
import numpy as np
import pytest

# Astropy
from ..name_resolve import get_icrs_coordinates
from ..builtin_systems import ICRSCoordinates, FK5Coordinates, FK4Coordinates, GalacticCoordinates

def test_names():

    # First check that sesame is up
    if urllib.urlopen("http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame").getcode() != 200:
    	pytest.skip("SESAME appears to be down, skipping test_name_resolve.py:test_names()...")

    with pytest.raises(ValueError):
        get_icrs_coordinates("m87h34hhh")

    #for name in ["ngc 3642", "m42", "castor", "pollux"]:
    icrs = get_icrs_coordinates("ngc 3642")
    icrs_true = ICRSCoordinates("11h 22m 18.014s", "59d 04m 27.27s")
    np.testing.assert_almost_equal(icrs.ra.degrees, icrs_true.ra.degrees, 3)
    np.testing.assert_almost_equal(icrs.dec.degrees, icrs_true.dec.degrees, 3)

    icrs = get_icrs_coordinates("castor")
    icrs_true = ICRSCoordinates("07h 34m 35.87s", "+31d 53m 17.8s")
    np.testing.assert_almost_equal(icrs.ra.degrees, icrs_true.ra.degrees, 3)
    np.testing.assert_almost_equal(icrs.dec.degrees, icrs_true.dec.degrees, 3)

def test_database_specify():

    # First check that sesame is up
    if urllib.urlopen("http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame").getcode() != 200:
    	pytest.skip("SESAME appears to be down, skipping test_database_specify.py:test_names()...")

    name = "ngc 3642"
    for db in ["simbad", "ned", "vizier", "all"]:
        icrs = ICRSCoordinates.resolve_name(name, database=db)
        gal = GalacticCoordinates.resolve_name(name).transform_to(ICRSCoordinates)
        np.testing.assert_almost_equal(icrs.ra.degrees, gal.ra.degrees, 1)
        np.testing.assert_almost_equal(icrs.dec.degrees, gal.dec.degrees, 1)

        time.sleep(1)

    name = "castor"
    # Don't search ned or vizier since castor doesn't seem to be in either
    for db in ["simbad",  "all"]:
        icrs = ICRSCoordinates.resolve_name(name, database=db)
        gal = GalacticCoordinates.resolve_name(name).transform_to(ICRSCoordinates)
        np.testing.assert_almost_equal(icrs.ra.degrees, gal.ra.degrees, 1)
        np.testing.assert_almost_equal(icrs.dec.degrees, gal.dec.degrees, 1)

        time.sleep(1)
