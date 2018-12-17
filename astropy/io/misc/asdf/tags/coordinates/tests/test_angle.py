# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

from astropy import __minimum_asdf_version__

asdf = pytest.importorskip('asdf', minversion=__minimum_asdf_version__)

import astropy.units as u

from asdf.tests.helpers import assert_roundtrip_tree

from astropy.coordinates import Longitude, Latitude, Angle

from astropy.io.misc.asdf.extension import AstropyExtension


def test_angle(tmpdir):
    tree = {'angle': Angle(100, u.deg)}
    assert_roundtrip_tree(tree, tmpdir)


def test_latitude(tmpdir):
    tree = {'angle': Latitude(10, u.deg)}
    assert_roundtrip_tree(tree, tmpdir)


def test_longitude(tmpdir):
    tree = {'angle': Longitude(-100, u.deg, wrap_angle=180*u.deg)}
    assert_roundtrip_tree(tree, tmpdir)
