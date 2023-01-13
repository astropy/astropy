# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A set of tests for the util.py module.
"""
import pytest

from astropy.io.votable import util


def test_range_list():
    assert util.coerce_range_list_param((5,)) == ("5.0", 1)


def test_range_list2():
    assert util.coerce_range_list_param((5e-7, 8e-7)) == ("5e-07,8e-07", 2)


def test_range_list3():
    assert util.coerce_range_list_param((5e-7, 8e-7, "FOO")) == ("5e-07,8e-07;FOO", 3)


def test_range_list4a():
    with pytest.raises(ValueError):
        util.coerce_range_list_param(
            (5e-7, (None, 8e-7), (4, None), (4, 5), "J", "FOO")
        )


def test_range_list4():
    assert util.coerce_range_list_param(
        (5e-7, (None, 8e-7), (4, None), (4, 5), "J", "FOO"), numeric=False
    ) == ("5e-07,/8e-07,4/,4/5,J;FOO", 6)


def test_range_list5():
    with pytest.raises(ValueError):
        util.coerce_range_list_param(("FOO",))


def test_range_list6():
    with pytest.raises(ValueError):
        print(util.coerce_range_list_param((5, "FOO"), util.stc_reference_frames))


def test_range_list7():
    assert util.coerce_range_list_param(("J",), numeric=False) == ("J", 1)


def test_range_list8():
    for s in [
        "5.0",
        "5e-07,8e-07",
        "5e-07,8e-07;FOO",
        "5e-07,/8e-07,4.0/,4.0/5.0;FOO",
        "J",
    ]:
        assert util.coerce_range_list_param(s, numeric=False)[0] == s


def test_range_list9a():
    with pytest.raises(ValueError):
        util.coerce_range_list_param("52,-27.8;FOO", util.stc_reference_frames)


def test_range_list9():
    assert util.coerce_range_list_param("52,-27.8;GALACTIC", util.stc_reference_frames)
