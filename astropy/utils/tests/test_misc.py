# Licensed under a 3-clause BSD style license - see LICENSE.rst

import json
import locale
import os
import urllib.error
from datetime import datetime

import numpy as np
import pytest

from astropy.io import fits
from astropy.tests.helper import CI
from astropy.utils import data, misc
from astropy.utils.exceptions import AstropyDeprecationWarning


def test_isiterable():
    assert misc.isiterable(2) is False
    assert misc.isiterable([2]) is True
    assert misc.isiterable([1, 2, 3]) is True
    assert misc.isiterable(np.array(2)) is False
    assert misc.isiterable(np.array([1, 2, 3])) is True


def test_signal_number_to_name_no_failure():
    # Regression test for #5340: ensure signal_number_to_name throws no
    # AttributeError (it used ".iteritems()" which was removed in Python3).
    misc.signal_number_to_name(0)


@pytest.mark.remote_data
def test_api_lookup():
    try:
        strurl = misc.find_api_page("astropy.utils.misc", "dev", False, timeout=5)
        objurl = misc.find_api_page(misc, "dev", False, timeout=5)
    except (urllib.error.URLError, TimeoutError):
        if CI:
            pytest.xfail("Timed out in CI")
        else:
            raise

    assert strurl == objurl
    assert (
        strurl
        == "http://devdocs.astropy.org/utils/ref_api.html#module-astropy.utils.misc"
    )

    # Try a non-dev version
    objurl = misc.find_api_page(misc, "stable", False, timeout=3)
    assert (
        objurl
        == "https://docs.astropy.org/en/stable/utils/ref_api.html#module-astropy.utils.misc"
    )


def test_is_path_hidden_deprecation():
    with pytest.warns(
        AstropyDeprecationWarning, match="^The is_path_hidden function is deprecated"
    ):
        misc.is_path_hidden("data")


# This is the only test that uses astropy/utils/tests/data/.hidden_file.txt
def test_skip_hidden():
    path = data.get_pkg_data_path("data")
    for root, dirs, files in os.walk(path):
        assert ".hidden_file.txt" in files
        assert "local.dat" in files
        # break after the first level since the data dir contains some other
        # subdirectories that don't have these files
        break
    with pytest.warns(
        AstropyDeprecationWarning, match="^The .*_hidden function is deprecated"
    ):
        for root, dirs, files in misc.walk_skip_hidden(path):
            assert ".hidden_file.txt" not in files
            assert "local.dat" in files
            break


def test_JsonCustomEncoder():
    from astropy import units as u

    assert json.dumps(np.arange(3), cls=misc.JsonCustomEncoder) == "[0, 1, 2]"
    assert json.dumps(1 + 2j, cls=misc.JsonCustomEncoder) == "[1.0, 2.0]"
    assert json.dumps({1, 2}, cls=misc.JsonCustomEncoder) == "[1, 2]"
    assert (
        json.dumps(b"hello world \xc3\x85", cls=misc.JsonCustomEncoder)
        == '"hello world \\u00c5"'
    )
    assert json.dumps({1: 2}, cls=misc.JsonCustomEncoder) == '{"1": 2}'  # default
    assert json.dumps({1: u.m}, cls=misc.JsonCustomEncoder) == '{"1": "m"}'
    # Quantities
    tmp = json.dumps({"a": 5 * u.cm}, cls=misc.JsonCustomEncoder)
    newd = json.loads(tmp)
    tmpd = {"a": {"unit": "cm", "value": 5.0}}
    assert newd == tmpd
    tmp2 = json.dumps({"a": np.arange(2) * u.cm}, cls=misc.JsonCustomEncoder)
    newd = json.loads(tmp2)
    tmpd = {"a": {"unit": "cm", "value": [0.0, 1.0]}}
    assert newd == tmpd
    tmp3 = json.dumps({"a": np.arange(2) * u.erg / u.s}, cls=misc.JsonCustomEncoder)
    newd = json.loads(tmp3)
    tmpd = {"a": {"unit": "erg / s", "value": [0.0, 1.0]}}
    assert newd == tmpd


def test_JsonCustomEncoder_FITS_rec_from_files():
    with fits.open(
        fits.util.get_testdata_filepath("variable_length_table.fits")
    ) as hdul:
        assert (
            json.dumps(hdul[1].data, cls=misc.JsonCustomEncoder)
            == "[[[45, 56], [11, 3]], [[11, 12, 13], [12, 4]]]"
        )

    with fits.open(fits.util.get_testdata_filepath("btable.fits")) as hdul:
        assert (
            json.dumps(hdul[1].data, cls=misc.JsonCustomEncoder)
            == '[[1, "Sirius", -1.4500000476837158, "A1V"], '
            '[2, "Canopus", -0.7300000190734863, "F0Ib"], '
            '[3, "Rigil Kent", -0.10000000149011612, "G2V"]]'
        )

    with fits.open(fits.util.get_testdata_filepath("table.fits")) as hdul:
        assert (
            json.dumps(hdul[1].data, cls=misc.JsonCustomEncoder)
            == '[["NGC1001", 11.100000381469727], '
            '["NGC1002", 12.300000190734863], '
            '["NGC1003", 15.199999809265137]]'
        )


def test_set_locale():
    # First, test if the required locales are available
    current = locale.setlocale(locale.LC_ALL)
    try:
        locale.setlocale(locale.LC_ALL, "en_US.utf8")
        locale.setlocale(locale.LC_ALL, "fr_FR.utf8")
    except locale.Error as e:
        pytest.skip(f"Locale error: {e}")
    finally:
        locale.setlocale(locale.LC_ALL, current)

    date = datetime(2000, 10, 1, 0, 0, 0)
    day_mon = date.strftime("%a, %b")

    with misc._set_locale("en_US.utf8"):
        assert date.strftime("%a, %b") == "Sun, Oct"

    with misc._set_locale("fr_FR.utf8"):
        assert date.strftime("%a, %b") == "dim., oct."

    # Back to original
    assert date.strftime("%a, %b") == day_mon

    with misc._set_locale(current):
        assert date.strftime("%a, %b") == day_mon


def test_dtype_bytes_or_chars():
    assert misc.dtype_bytes_or_chars(np.dtype(np.float64)) == 8
    assert misc.dtype_bytes_or_chars(np.dtype(object)) is None
    assert misc.dtype_bytes_or_chars(np.dtype(np.int32)) == 4
    assert misc.dtype_bytes_or_chars(np.array(b"12345").dtype) == 5
    assert misc.dtype_bytes_or_chars(np.array("12345").dtype) == 5


def test_indent_deprecation():
    with pytest.warns(AstropyDeprecationWarning, match=r"Use textwrap\.indent"):
        misc.indent("Obsolete since Python 3.3")
