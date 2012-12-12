# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains tests for the name resolve convenience module.
"""

# Standard library
import time
import urllib2

# Third party
import numpy as np
import pytest

# Astropy
from ..name_resolve import get_icrs_coordinates
from ..builtin_systems import ICRSCoordinates, FK5Coordinates, FK4Coordinates, GalacticCoordinates

def test_names():

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

def test_transforms():

    for name in ["ngc 3642", "m42", "castor", "pollux"]:
        icrs = ICRSCoordinates.resolve_name(name)
        gal = GalacticCoordinates.resolve_name(name).transform_to(ICRSCoordinates)
        np.testing.assert_almost_equal(icrs.ra.degrees, gal.ra.degrees, 7)
        np.testing.assert_almost_equal(icrs.dec.degrees, gal.dec.degrees, 7)

def test_database_specify():

    for db in ["simbad", "ned", "vizier", "all"]:
        for name in ["ngc 3642", "m42"]:
            icrs = ICRSCoordinates.resolve_name(name, database=db)
            gal = GalacticCoordinates.resolve_name(name).transform_to(ICRSCoordinates)
            np.testing.assert_almost_equal(icrs.ra.degrees, gal.ra.degrees, 1)
            np.testing.assert_almost_equal(icrs.dec.degrees, gal.dec.degrees, 1)

            time.sleep(1)

def test_database_castor():
    name = "castor"
    for db in ["simbad", "ned", "vizier", "all"]:
        if db == "ned":
            with pytest.raises(urllib2.URLError):
                icrs = ICRSCoordinates.resolve_name(name, database=db)
            continue

        icrs = ICRSCoordinates.resolve_name(name, database=db)
        gal = GalacticCoordinates.resolve_name(name).transform_to(ICRSCoordinates)
        np.testing.assert_almost_equal(icrs.ra.degrees, gal.ra.degrees, 1)
        np.testing.assert_almost_equal(icrs.dec.degrees, gal.dec.degrees, 1)

        time.sleep(1)