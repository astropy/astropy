# -*- coding: utf-8 -*-

"""Testing :mod:`astropy.cosmology.units`."""

##############################################################################
# IMPORTS

import pytest

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import default_cosmology
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_ASDF

##############################################################################
# TESTS
##############################################################################


def test_has_expected_units():
    """
    Test that this module has the expected set of units. Some of the units are
    imported from :mod:`astropy.units`, or vice versa. Here we test presence,
    not usage. Units from :mod:`astropy.units` are tested in that module. Units
    defined in :mod:`astropy.cosmology` will be tested subsequently.
    """
    assert u.astrophys.littleh is cu.littleh


def test_has_expected_equivalencies():
    """
    Test that this module has the expected set of equivalencies. Many of the
    equivalencies are imported from :mod:`astropy.units`, so here we test
    presence, not usage. Equivalencies from :mod:`astropy.units` are tested in
    that module. Equivalencies defined in :mod:`astropy.cosmology` will be
    tested subsequently.
    """
    assert u.equivalencies.with_H0 is cu.with_H0


def test_littleh():
    H0_70 = 70 * u.km / u.s / u.Mpc
    h70dist = 70 * u.Mpc / cu.littleh

    assert_quantity_allclose(h70dist.to(u.Mpc, cu.with_H0(H0_70)), 100 * u.Mpc)

    # make sure using the default cosmology works
    cosmodist = default_cosmology.get().H0.value * u.Mpc / cu.littleh
    assert_quantity_allclose(cosmodist.to(u.Mpc, cu.with_H0()), 100 * u.Mpc)

    # Now try a luminosity scaling
    h1lum = 0.49 * u.Lsun * cu.littleh ** -2
    assert_quantity_allclose(h1lum.to(u.Lsun, cu.with_H0(H0_70)), 1 * u.Lsun)

    # And the trickiest one: magnitudes.  Using H0=10 here for the round numbers
    H0_10 = 10 * u.km / u.s / u.Mpc
    # assume the "true" magnitude M = 12.
    # Then M - 5*log_10(h)  = M + 5 = 17
    withlittlehmag = 17 * (u.mag - u.MagUnit(cu.littleh ** 2))
    assert_quantity_allclose(withlittlehmag.to(u.mag, cu.with_H0(H0_10)), 12 * u.mag)


@pytest.mark.skipif(not HAS_ASDF, reason="requires ASDF")
@pytest.mark.parametrize('equiv', [cu.with_H0()])
def test_equivalencies(tmpdir, equiv):
    from asdf.tests import helpers

    tree = {'equiv': equiv}
    helpers.assert_roundtrip_tree(tree, tmpdir)
