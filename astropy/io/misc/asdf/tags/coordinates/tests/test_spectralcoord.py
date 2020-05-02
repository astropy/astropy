# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

from astropy import units as u
from astropy.coordinates import SpectralCoord, ICRS, Galactic
from astropy.tests.helper import assert_quantity_allclose

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree  # noqa


def test_scalar_spectralcoord(tmpdir):

    sc = SpectralCoord(565 * u.nm)
    tree = dict(spectralcoord=sc)

    def check(asdffile):
        assert isinstance(asdffile['spectralcoord'], SpectralCoord)
        assert_quantity_allclose(asdffile['spectralcoord'].quantity, 565 * u.nm)

    assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check)


def test_vector_spectralcoord(tmpdir):

    sc = SpectralCoord([100, 200, 300] * u.GHz)
    tree = dict(spectralcoord=sc)

    def check(asdffile):
        assert isinstance(asdffile['spectralcoord'], SpectralCoord)
        assert_quantity_allclose(asdffile['spectralcoord'].quantity, [100, 200, 300] * u.GHz)

    assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check, tree_match_func=assert_quantity_allclose)


@pytest.mark.filterwarnings("ignore:No velocity")
def test_spectralcoord_with_obstarget(tmpdir):

    sc = SpectralCoord(10 * u.GHz,
                       observer=ICRS(1 * u.km, 2 * u.km, 3 * u.km, representation_type='cartesian'),
                       target=Galactic(10 * u.deg, 20 * u.deg, distance=30 * u.pc))
    tree = dict(spectralcoord=sc)

    def check(asdffile):
        assert isinstance(asdffile['spectralcoord'], SpectralCoord)
        assert_quantity_allclose(asdffile['spectralcoord'].quantity, 10 * u.GHz)
        assert isinstance(asdffile['spectralcoord'].observer, ICRS)
        assert isinstance(asdffile['spectralcoord'].target, Galactic)

    assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check)
