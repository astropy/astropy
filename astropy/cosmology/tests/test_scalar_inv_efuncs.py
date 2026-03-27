# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for the scalar_inv_efuncs Cython extension in astropy.cosmology."""

from numpy.testing import assert_allclose

from astropy.cosmology._src.flrw.scalar_inv_efuncs import (
    flcdm_inv_efunc_nomnu,
    flcdm_inv_efunc_norel,
    fwcdm_inv_efunc_norel,
    lcdm_inv_efunc_nomnu,
    lcdm_inv_efunc_norel,
    wcdm_inv_efunc_norel,
)

# ---------------------------------------------------------------------------
# LambdaCDM (no relativistic species)
# At z=0: E(0)=1 => inv_efunc=1 for flat universe (Om0+Ode0=1, Ok0=0)
# ---------------------------------------------------------------------------


def test_lcdm_inv_efunc_norel_z0():
    # Flat universe at z=0: Om0+Ode0=1, Ok0=0 => E(0)=1
    assert_allclose(lcdm_inv_efunc_norel(0.0, 0.3, 0.7, 0.0), 1.0)


def test_lcdm_inv_efunc_norel_matter_dominated():
    # Matter-only universe (Ode0=0, Ok0=0, Om0=1): E(z)=(1+z)^1.5
    # => inv_efunc(z) = (1+z)^-1.5
    z = 1.0
    expected = (1.0 + z) ** -1.5
    assert_allclose(lcdm_inv_efunc_norel(z, 1.0, 0.0, 0.0), expected)


def test_lcdm_inv_efunc_norel_positive():
    """inv_efunc must always be positive for physical parameters."""
    assert lcdm_inv_efunc_norel(0.5, 0.3, 0.7, 0.0) > 0.0


# ---------------------------------------------------------------------------
# LambdaCDM (massless neutrinos)
# ---------------------------------------------------------------------------


def test_lcdm_inv_efunc_nomnu_z0():
    # Flat universe at z=0: Om0+Ode0+Or0=1, Ok0=0 => E(0)=1
    assert_allclose(lcdm_inv_efunc_nomnu(0.0, 0.3, 0.6999, 0.0, 0.0001), 1.0)


def test_lcdm_inv_efunc_nomnu_positive():
    """inv_efunc must always be positive for physical parameters."""
    assert lcdm_inv_efunc_nomnu(1.0, 0.3, 0.7, 0.0, 0.0) > 0.0


# ---------------------------------------------------------------------------
# FlatLambdaCDM (no relativistic species)
# ---------------------------------------------------------------------------


def test_flcdm_inv_efunc_norel_z0():
    # Flat: Om0+Ode0=1 => E(0)=1
    assert_allclose(flcdm_inv_efunc_norel(0.0, 0.3, 0.7), 1.0)


def test_flcdm_inv_efunc_norel_matter_dominated():
    # Matter-only flat: Ode0=0, Om0=1 => inv_efunc(z)=(1+z)^-1.5
    z = 2.0
    expected = (1.0 + z) ** -1.5
    assert_allclose(flcdm_inv_efunc_norel(z, 1.0, 0.0), expected)


def test_flcdm_inv_efunc_norel_positive():
    """inv_efunc must always be positive for physical parameters."""
    assert flcdm_inv_efunc_norel(0.5, 0.3, 0.7) > 0.0


# ---------------------------------------------------------------------------
# FlatLambdaCDM (massless neutrinos)
# ---------------------------------------------------------------------------


def test_flcdm_inv_efunc_nomnu_z0():
    # Flat: Om0+Ode0+Or0=1 => E(0)=1
    assert_allclose(flcdm_inv_efunc_nomnu(0.0, 0.3, 0.6999, 0.0001), 1.0)


# ---------------------------------------------------------------------------
# wCDM (no relativistic species)
# ---------------------------------------------------------------------------


def test_wcdm_inv_efunc_norel_z0():
    # Flat wCDM at z=0 with Om0+Ode0=1, Ok0=0 => E(0)=1
    assert_allclose(wcdm_inv_efunc_norel(0.0, 0.3, 0.7, 0.0, -1.0), 1.0)


def test_wcdm_inv_efunc_norel_lcdm_limit():
    # w0=-1 => wCDM reduces to LambdaCDM
    z = 0.5
    Om0, Ode0, Ok0, w0 = 0.3, 0.7, 0.0, -1.0
    assert_allclose(
        wcdm_inv_efunc_norel(z, Om0, Ode0, Ok0, w0),
        lcdm_inv_efunc_norel(z, Om0, Ode0, Ok0),
    )


def test_wcdm_inv_efunc_norel_positive():
    """inv_efunc must always be positive for physical parameters."""
    assert wcdm_inv_efunc_norel(1.0, 0.3, 0.7, 0.0, -0.8) > 0.0


# ---------------------------------------------------------------------------
# FlatwCDM (no relativistic species)
# ---------------------------------------------------------------------------


def test_fwcdm_inv_efunc_norel_z0():
    # Flat wCDM at z=0: Om0+Ode0=1 => E(0)=1
    assert_allclose(fwcdm_inv_efunc_norel(0.0, 0.3, 0.7, -1.0), 1.0)


def test_fwcdm_inv_efunc_norel_lcdm_limit():
    # w0=-1 => flat wCDM reduces to flat LambdaCDM
    z = 1.0
    Om0, Ode0, w0 = 0.3, 0.7, -1.0
    assert_allclose(
        fwcdm_inv_efunc_norel(z, Om0, Ode0, w0),
        flcdm_inv_efunc_norel(z, Om0, Ode0),
    )
