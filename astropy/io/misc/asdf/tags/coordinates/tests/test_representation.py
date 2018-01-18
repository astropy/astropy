# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from collections import OrderedDict
import pytest


import astropy.units as u

from astropy.coordinates import Longitude, Latitude
import astropy.coordinates.representation as r

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree


class UnitSphericalWrap180Representationasdf(r.UnitSphericalRepresentation):
    attr_classes = OrderedDict([('lon', Longitude), ('lat', Latitude)])


def test_unitspherical(tmpdir):
    tree = {'representation': r.UnitSphericalRepresentation(10*u.deg, 10*u.deg)}
    assert_roundtrip_tree(tree, tmpdir)


def test_spherical(tmpdir):
    tree = {'representation': r.SphericalRepresentation(10*u.deg, 10*u.deg, 2.2*u.AU)}
    assert_roundtrip_tree(tree, tmpdir)


def test_unitspherical180(tmpdir):
    x = UnitSphericalWrap180Representationasdf(-10*u.deg, 10*u.deg)
    tree = {'representation': x}
    with pytest.raises(ValueError):
        assert_roundtrip_tree(tree, tmpdir)


def test_cart_differential(tmpdir):
    d = r.CartesianDifferential([11.1, 12.24, 7.25]*u.km/u.s)

    tree = {'representation': d}
    assert_roundtrip_tree(tree, tmpdir)


def test_sunpy180(tmpdir):
    scoord = pytest.importorskip('sunpy.coordinates')

    x = scoord.representation.SphericalWrap180Representation(-10*u.deg, 10*u.deg, 0.98*u.AU)
    tree = {'representation': x}
    assert_roundtrip_tree(tree, tmpdir)
