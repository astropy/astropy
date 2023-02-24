# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import pytest

asdf = pytest.importorskip("asdf")

from asdf.exceptions import AsdfDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        category=AsdfDeprecationWarning,
        message=r"asdf.tests.helpers is deprecated.*",
    )
    from asdf.tests.helpers import assert_roundtrip_tree

import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude


def test_angle(tmpdir):
    tree = {"angle": Angle(100, u.deg)}
    assert_roundtrip_tree(tree, tmpdir)


def test_latitude(tmpdir):
    tree = {"angle": Latitude(10, u.deg)}
    assert_roundtrip_tree(tree, tmpdir)


def test_longitude(tmpdir):
    tree = {"angle": Longitude(-100, u.deg, wrap_angle=180 * u.deg)}
    assert_roundtrip_tree(tree, tmpdir)
