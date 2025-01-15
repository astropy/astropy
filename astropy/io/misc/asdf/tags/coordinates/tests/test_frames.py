# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

asdf = pytest.importorskip("asdf")
from asdf.tests.helpers import assert_roundtrip_tree  # noqa: E402

from astropy import units  # noqa: E402
from astropy.coordinates import FK5, ICRS, Angle, Latitude, Longitude  # noqa: E402


def test_hcrs_basic(tmpdir):
    ra = Longitude(25, unit=units.deg)
    dec = Latitude(45, unit=units.deg)

    tree = {"coord": ICRS(ra=ra, dec=dec)}

    assert_roundtrip_tree(tree, tmpdir)


def test_icrs_basic(tmpdir):
    wrap_angle = Angle(1.5, unit=units.rad)
    ra = Longitude(25, unit=units.deg, wrap_angle=wrap_angle)
    dec = Latitude(45, unit=units.deg)

    tree = {"coord": ICRS(ra=ra, dec=dec)}

    assert_roundtrip_tree(tree, tmpdir)


def test_icrs_nodata(tmpdir):
    tree = {"coord": ICRS()}

    assert_roundtrip_tree(tree, tmpdir)


def test_icrs_compound(tmpdir):
    icrs = ICRS(ra=[0, 1, 2] * units.deg, dec=[3, 4, 5] * units.deg)

    tree = {"coord": icrs}

    assert_roundtrip_tree(tree, tmpdir)


def test_fk5_time(tmpdir):
    tree = {"coord": FK5(equinox="2011-01-01T00:00:00")}

    assert_roundtrip_tree(tree, tmpdir)
