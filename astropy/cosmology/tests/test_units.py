# -*- coding: utf-8 -*-

"""Testing :mod:`astropy.cosmology.units`."""

##############################################################################
# IMPORTS

import contextlib

import pytest

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import default_cosmology, Planck13
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_ASDF, HAS_SCIPY
from astropy.utils.exceptions import AstropyDeprecationWarning

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
    with pytest.warns(AstropyDeprecationWarning, match="`littleh`"):
        assert u.astrophys.littleh is cu.littleh


def test_has_expected_equivalencies():
    """
    Test that this module has the expected set of equivalencies. Many of the
    equivalencies are imported from :mod:`astropy.units`, so here we test
    presence, not usage. Equivalencies from :mod:`astropy.units` are tested in
    that module. Equivalencies defined in :mod:`astropy.cosmology` will be
    tested subsequently.
    """
    with pytest.warns(AstropyDeprecationWarning, match="`with_H0`"):
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


@pytest.mark.skipif(not HAS_SCIPY, reason="Cosmology needs scipy")
def test_dimensionless_redshift():
    """Test the equivalency  ``dimensionless_redshift``.
    """
    z = 3 * cu.redshift
    val = 3 * u.one

    # show units not equal
    assert z.unit == cu.redshift
    assert z.unit != u.one

    # test equivalency enabled by default
    assert z.to(u.one) == val

    # also test that it works for powers
    assert (3 * cu.redshift ** 3).to(u.one) == val

    # and in composite units
    assert (3 * u.km / cu.redshift**3).to(u.km) == 3 * u.km

    # test it also works as an equivalency
    with u.set_enabled_equivalencies([]):  # turn off default equivalencies
        assert z.to(u.one, equivalencies=cu.dimensionless_redshift()) == val

    # if this fails, something is really wrong
    with u.add_enabled_equivalencies(cu.dimensionless_redshift()):
        assert z.to(u.one) == val


@pytest.mark.skipif(not HAS_SCIPY, reason="Cosmology needs scipy")
def test_with_redshift():
    """Test the equivalency  ``with_redshift``.
    """
    cosmo = Planck13.clone(Tcmb0=3*u.K)  # make non-default cosmology
    z = 15 * cu.redshift
    Tcmb = cosmo.Tcmb(z)

    default_cosmo = default_cosmology.get()
    assert default_cosmo != cosmo  # shows changing default

    # 1) Default (without specifying the cosmology)
    with default_cosmology.set(cosmo):
        # A) Default (Tcmb=True)
        equivalency = cu.with_redshift()
        assert_quantity_allclose(z.to(u.K, equivalency), Tcmb)
        assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)

        # B) Tcmb=False
        equivalency = cu.with_redshift(Tcmb=False)
        with pytest.raises(u.UnitConversionError, match="'redshift' and 'K'"):
            z.to(u.K, equivalency)

    # showing the answer changes if the cosmology changes
    # this test uses the default cosmology
    equivalency = cu.with_redshift()
    assert_quantity_allclose(z.to(u.K, equivalency),
                             default_cosmo.Tcmb(z))
    assert default_cosmo.Tcmb(z) != Tcmb

    # 2) Specifying the cosmology
    # A) Default (Tcmb=True)
    equivalency = cu.with_redshift(cosmo)
    assert_quantity_allclose(z.to(u.K, equivalency), Tcmb)
    assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)

    # B) Tcmb=False
    equivalency = cu.with_redshift(cosmo, Tcmb=False)
    with pytest.raises(u.UnitConversionError, match="'redshift' and 'K'"):
        z.to(u.K, equivalency)

    # Test atzkw
    # this is really just a test that 'atzkw' doesn't fail
    equivalency = cu.with_redshift(cosmo, atzkw={"ztol": 1e-10})
    assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)


# FIXME! get "dimensionless_redshift", "with_redshift" to work in this
# they are not in ``astropy.units.equivalencies``, so the following fails
@pytest.mark.skipif(not HAS_ASDF, reason="requires ASDF")
@pytest.mark.parametrize('equiv',[cu.with_H0])
def test_equivalencies_asdf(tmpdir, equiv):
    from asdf.tests import helpers

    tree = {'equiv': equiv()}
    with (pytest.warns(AstropyDeprecationWarning, match="`with_H0`")
          if equiv.__name__ == "with_H0" else contextlib.nullcontext()):

        helpers.assert_roundtrip_tree(tree, tmpdir)
