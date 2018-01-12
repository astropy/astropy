# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from collections import OrderedDict
import pytest


import astropy.units as u

from astropy.coordinates import SphericalRepresentation, UnitSphericalRepresentation, Longitude, Latitude

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree


class UnitSphericalWrap180Representationasdf(UnitSphericalRepresentation):
    attr_classes = OrderedDict([('lon', Longitude), ('lat', Latitude)])


def test_unitspherical(tmpdir):
    tree = {'representation': UnitSphericalRepresentation(10*u.deg, 10*u.deg)}
    assert_roundtrip_tree(tree, tmpdir)


def test_spherical(tmpdir):
    tree = {'representation': SphericalRepresentation(10*u.deg, 10*u.deg, 2.2*u.AU)}
    assert_roundtrip_tree(tree, tmpdir)


def test_unitspherical180(tmpdir):
    x = UnitSphericalWrap180Representationasdf(-10*u.deg, 10*u.deg)
    tree = {'representation': x}
    with pytest.raises(ValueError):
        assert_roundtrip_tree(tree, tmpdir)


def test_sunpy180(tmpdir):
    scoord = pytest.importorskip('sunpy.coordinates')

    x = scoord.representation.SphericalWrap180Representation(-10*u.deg, 10*u.deg, 0.98*u.AU)
    tree = {'representation': x}
    assert_roundtrip_tree(tree, tmpdir)
