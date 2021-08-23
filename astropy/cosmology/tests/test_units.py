# -*- coding: utf-8 -*-

"""Testing :mod:`astropy.cosmology.units`."""

##############################################################################
# IMPORTS

import astropy.cosmology.units as cu
import astropy.units as u

##############################################################################
# TESTS
##############################################################################


def test_has_expected_units():
    """
    Test that this module has the expected set of units. Many of the units are
    imported from :mod:`astropy.units`, so here we test presence, not usage.
    Units from :mod:`astropy.units` are tested in that module. Units defined in
    :mod:`astropy.cosmology` will be tested subsequently.
    """
    assert cu.littleh is u.astrophys.littleh


def test_has_expected_equivalencies():
    """
    Test that this module has the expected set of equivalencies. Many of the
    equivalencies are imported from :mod:`astropy.units`, so here we test
    presence, not usage. Equivalencies from :mod:`astropy.units` are tested in
    that module. Equivalencies defined in :mod:`astropy.cosmology` will be
    tested subsequently.
    """
    assert cu.with_H0 is u.equivalencies.with_H0
