# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import urllib.request
import warnings

import pytest
import numpy as np
from numpy.testing import assert_equal

from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import download_file, clear_download_cache
from astropy.utils.iers import (
    IERS_A, IERS_Auto, IERS_B, IERS_B_URL, IERS_B_FILE
)
from astropy import units as u
from astropy.table import QTable
from astropy.time import Time, TimeDelta

try:
    IERS_A.open("finals2000A.all")  # check if IERS_A is available
except OSError:
    HAS_IERS_A = False
else:
    HAS_IERS_A = True


test_columns = [
    # These are very stringent thresholds
    # Can mark individual parameter sets with xfail, see docs
    ("UT1_UTC", 10 * u.us),
    ("PM_x", 1e-3 * u.marcsec),
    ("PM_y", 1e-3 * u.marcsec),
    ("dX_2000A", 1e-4 * u.marcsec),
    ("dY_2000A", 1e-4 * u.marcsec),
    # Maybe repeat tests with more practical thresholds
    # but not xfail?
]


@pytest.mark.parametrize("keyword,atol", test_columns)
@pytest.mark.xfail(
        reason="IERS B has changed old values since the version in astropy")
@pytest.mark.remote_data
def test_IERS_B_download_agrees_with_IERS_B_local(keyword, atol):
    B_local = IERS_B.open(IERS_B_FILE)
    B = IERS_B.open(download_file(IERS_B_URL, cache=True))
    if B["MJD"][-1] < B_local["MJD"][-1]:
        clear_download_cache(IERS_B_URL)
        B = IERS_B.open(download_file(IERS_B_URL, cache=True))
    assert (
        B["MJD"][-1] >= B_local["MJD"][-1]
    ), "Local IERS B data should never be newer than the internet version!"
    if B["MJD"][-1] == B_local["MJD"][-1]:
        # Unable to run test because current version is in astropy already
        return

    mjd = B_local["MJD"][-1]

    ok_B = B["MJD"] <= mjd
    assert_quantity_allclose(
        B["MJD"][ok_B],
        B_local["MJD"],
        atol=0.01 * u.day,
        rtol=0,
        err_msg="MJDs don't match",
    )
    assert_quantity_allclose(
        B[keyword][ok_B],
        B_local[keyword],
        atol=atol,
        rtol=0,
        err_msg="Built-in IERS B {} values don't match downloaded "
        "IERS B values".format(keyword),
    )


copy_columns = [
    ("UT1_UTC", "UT1_UTC_B"),
    ("PM_x", "PM_X_B"),
    ("PM_y", "PM_Y_B"),
    ("dX_2000A", "dX_2000A_B"),
    ("dY_2000A", "dY_2000A_B"),
]


@pytest.mark.parametrize("b_name,a_name", copy_columns)
def test_IERS_B_parameters_loaded_into_IERS_Auto(b_name, a_name):
    B = IERS_B.open(IERS_B_FILE)
    B[b_name]
    A = IERS_Auto.open()
    try:
        A[a_name]
    except KeyError:
        if A["MJD"][-1] < 59100 * u.d:
            pytest.xfail(
                "Bug #9205 IERS B data is available but not merged into "
                "IERS_Auto object unless new IERS A data is available."
            )
        else:
            raise


copy_columns_atol = [
    ("UT1_UTC", "UT1_UTC_B", 10 * u.us),
    ("PM_x", "PM_X_B", 1e-3 * u.marcsec),
    ("PM_y", "PM_Y_B", 1e-3 * u.marcsec),
    pytest.param("dX_2000A", "dX_2000A_B", 1e-4 * u.marcsec,
                 marks=pytest.mark.xfail),
    pytest.param("dY_2000A", "dY_2000A_B", 1e-4 * u.marcsec,
                 marks=pytest.mark.xfail),
]


@pytest.mark.remote_data
@pytest.mark.parametrize("b_name,a_name,atol", copy_columns_atol)
def test_IERS_B_parameters_loaded_into_IERS_Auto_correctly(
        b_name,
        a_name,
        atol):
    A = IERS_Auto.open()
    A[a_name]
    B = IERS_B.open(IERS_B_FILE)

    ok_A = A["MJD"] < B["MJD"][-1]

    mjds_A = A["MJD"][ok_A].to(u.day).value
    i_B = np.searchsorted(B["MJD"].to(u.day).value, mjds_A)

    assert_equal(np.diff(i_B), 1, err_msg="Valid region not contiguous")
    assert_equal(A["MJD"][ok_A], B["MJD"][i_B],
                 err_msg="MJDs don't make sense")
    assert_quantity_allclose(
        A[a_name][ok_A],
        B[b_name][i_B],
        atol=atol,
        rtol=0,
        err_msg="Bug #9206 IERS B parameter {} not copied over "
                "IERS A parameter {}".format(b_name, a_name),
    )
