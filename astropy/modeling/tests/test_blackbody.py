# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for blackbody model and functions."""
# pylint: disable=invalid-name, no-member

import pytest
import numpy as np

from astropy import constants as const
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose, catch_warnings
from astropy.utils.exceptions import AstropyUserWarning

from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.blackbody import (
    BlackBody1D, blackbody_nu, blackbody_lambda, FNU)

try:
    from scipy import integrate  # noqa
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

__doctest_skip__ = ['*']


@pytest.mark.filterwarnings("ignore:BlackBody provides the same capabilities")
class TestBlackbody1D:

    # Make sure the temperature equivalency automatically applies by trying
    # to pass temperatures in celsius

    @pytest.mark.parametrize('temperature', (3000 * u.K, 2726.85 * u.deg_C))
    def test_evaluate(self, temperature):

        bolometric_flux = 1000 * u.L_sun / (4 * np.pi * (1.5 * u.pc) ** 2)

        b = BlackBody1D(temperature=temperature,
                        bolometric_flux=bolometric_flux)

        assert_quantity_allclose(b(1.4 * u.micron), 4734463.93280426 * u.Jy)
        assert_quantity_allclose(b(214.13747 * u.THz), 4734463.93280426 * u.Jy)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_fit(self):

        fitter = LevMarLSQFitter()

        b = BlackBody1D(3000 * u.K)

        wav = np.array([0.5, 5, 10]) * u.micron
        fnu = np.array([1, 10, 5]) * u.Jy

        b_fit = fitter(b, wav, fnu)

        assert_quantity_allclose(b_fit.temperature, 2840.7438339457754 * u.K)
        assert_quantity_allclose(b_fit.bolometric_flux, 6.821837075583734e-08 * u.erg / u.cm**2 / u.s)


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.filterwarnings("ignore:BlackBody provides the same capabilities")
def test_blackbody_scipy():
    """Test Planck function.

    .. note:: Needs ``scipy`` to work.

    """
    flux_unit = u.Watt / (u.m ** 2 * u.um)
    wave = np.logspace(0, 8, 100000) * u.AA
    temp = 100. * u.K
    with np.errstate(all='ignore'):
        bb_nu = blackbody_nu(wave, temp) * u.sr
    flux = bb_nu.to(flux_unit, u.spectral_density(wave)) / u.sr

    lum = wave.to(u.um)
    intflux = integrate.trapz(flux.value, x=lum.value)
    ans = const.sigma_sb * temp ** 4 / np.pi
    np.testing.assert_allclose(intflux, ans.value, rtol=0.01)  # 1% accuracy


@pytest.mark.filterwarnings("ignore:BlackBody provides the same capabilities")
def test_blackbody_overflow():
    """Test Planck function with overflow."""
    photlam = u.photon / (u.cm**2 * u.s * u.AA)
    wave = [0, 1000.0, 100000.0, 1e55]  # Angstrom
    temp = 10000.0  # Kelvin
    with pytest.warns(
            AstropyUserWarning,
            match=r'Input contains invalid wavelength/frequency value\(s\)'):
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


@pytest.mark.filterwarnings("ignore:BlackBody provides the same capabilities")
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


@pytest.mark.filterwarnings("ignore:BlackBody provides the same capabilities")
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


@pytest.mark.filterwarnings("ignore:BlackBody provides the same capabilities")
def test_blackbody_array_temperature():
    """Regression test to make sure that the temperature can be an array."""
    flux = blackbody_nu(1.2 * u.mm, [100, 200, 300] * u.K)
    np.testing.assert_allclose(
        flux.value, [1.804908e-12, 3.721328e-12, 5.638513e-12], rtol=1e-5)

    flux = blackbody_nu([2, 4, 6] * u.mm, [100, 200, 300] * u.K)
    np.testing.assert_allclose(
        flux.value, [6.657915e-13, 3.420677e-13, 2.291897e-13], rtol=1e-5)

    flux = blackbody_nu(np.ones((3, 4)) * u.mm, np.ones(4) * u.K)
    assert flux.shape == (3, 4)
