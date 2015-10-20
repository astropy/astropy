# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for blackbody functions."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import warnings
# THIRD-PARTY
import numpy as np

# LOCAL
from ..blackbody import blackbody_nu, blackbody_lambda, FNU
from ... import units as u
from ... import constants as const
from ...tests.helper import pytest, catch_warnings
from ...utils.exceptions import AstropyUserWarning

try:
    from scipy import integrate  # pylint: disable=W0611
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


__doctest_skip__ = ['*']


@pytest.mark.skipif("os.environ.get('APPVEYOR')",  reason="fails on AppVeyor")
@pytest.mark.skipif('not HAS_SCIPY')
def test_blackbody_scipy():
    """Test Planck function.

    .. note:: Needs ``scipy`` to work.

    """
    flux_unit = u.Watt / (u.m ** 2 * u.um)
    wave = np.logspace(0, 8, 1e6) * u.AA
    temp = 100. * u.K
    with np.errstate(all='ignore'):
        bb_nu = blackbody_nu(wave, temp) * u.sr
    flux = bb_nu.to(flux_unit, u.spectral_density(wave)) / u.sr

    lum = wave.to(u.um)
    intflux = integrate.trapz(flux.value, x=lum.value)
    ans = const.sigma_sb * temp**4 / np.pi
    np.testing.assert_allclose(intflux, ans.value, rtol=0.01)  # 1% accuracy


@pytest.mark.skipif("os.environ.get('APPVEYOR')",  reason="fails on AppVeyor")
def test_blackbody_overflow():
    """Test Planck function with overflow."""
    photlam = u.photon / (u.cm**2 * u.s * u.AA)
    wave = [0, 1000.0, 100000.0, 1e55]  # Angstrom
    temp = 10000.0  # Kelvin
    with np.errstate(all='ignore'):
        bb_lam = blackbody_lambda(wave, temp) * u.sr
    flux = bb_lam.to(photlam, u.spectral_density(wave * u.AA)) / u.sr

    # First element is NaN, last element is very small, others normal
    assert np.isnan(flux[0])
    assert np.log10(flux[-1].value) < -134
    np.testing.assert_allclose(
        flux.value[1:-1], [3.38131732e+16, 3.87451317e+15],
        rtol=1e-3)  # 0.1% accuracy in PHOTLAM/sr

    with np.errstate(all='ignore'):
        flux = blackbody_lambda(1, 1e4)
    assert flux.value == 0


@pytest.mark.skipif("os.environ.get('APPVEYOR')",  reason="fails on AppVeyor")
def test_blackbody_synphot():
    """Test that it is consistent with IRAF SYNPHOT BBFUNC."""
    # Solid angle of solar radius at 1 kpc
    fac = np.pi * (const.R_sun / const.kpc) ** 2 * u.sr

    with np.errstate(all='ignore'):
        flux = blackbody_nu([100, 1, 1000, 1e4, 1e5] * u.AA, 5000) * fac
    assert flux.unit == FNU

    # Special check for overflow value (SYNPHOT gives 0)
    assert np.log10(flux[0].value) < -143

    np.testing.assert_allclose(
        flux.value[1:], [0, 2.01950807e-34, 3.78584515e-26, 1.90431881e-27],
        rtol=0.01)  # 1% accuracy


def test_blackbody_exceptions_and_warnings():
    """Test exceptions."""

    # Negative temperature
    with pytest.raises(ValueError) as exc:
        blackbody_nu(1000 * u.AA, -100)
    assert exc.value.args[0] == 'Temperature should be positive: -100.0 K'

    # Zero wavelength given for conversion to Hz
    with catch_warnings(AstropyUserWarning) as w:
        blackbody_nu(0 * u.AA, 5000)
    assert len(w) == 1
    assert 'invalid' in w[0].message.args[0]

    # Negative wavelength given for conversion to Hz
    with catch_warnings(AstropyUserWarning) as w:
        blackbody_nu(-1. * u.AA, 5000)
    assert len(w) == 1
    assert 'invalid' in w[0].message.args[0]


def test_blackbody_array_temperature():

    # Regression test to make sure that the temperature can be an array

    flux = blackbody_nu(1.2 * u.mm, [100, 200, 300] * u.K)

    np.testing.assert_allclose(flux.value,
                               [1.804908e-12,   3.721328e-12,   5.638513e-12],
                                rtol=1e-5)

    flux = blackbody_nu([2, 4, 6] * u.mm, [100, 200, 300] * u.K)

    np.testing.assert_allclose(flux.value,
                               [6.657915e-13,   3.420677e-13,   2.291897e-13],
                                rtol=1e-5)

    flux = blackbody_nu(np.ones((3, 4)) * u.mm, np.ones(4) * u.K)
    assert flux.shape == (3, 4)

