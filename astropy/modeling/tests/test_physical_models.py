# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for physical functions."""
# pylint: disable=no-member, invalid-name
import pytest
import numpy as np

from astropy.modeling.physical_models import BlackBody
from astropy.modeling.fitting import LevMarLSQFitter

from astropy.tests.helper import assert_quantity_allclose, catch_warnings
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning

try:
    from scipy import optimize, integrate  # noqa

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

__doctest_skip__ = ["*"]


# BlackBody tests


@pytest.mark.parametrize("temperature", (3000 * u.K, 2726.85 * u.deg_C))
def test_blackbody_evaluate(temperature):

    b = BlackBody(temperature=temperature, scale=1.0)

    assert_quantity_allclose(b(1.4 * u.micron), 486787299458.15656 * u.MJy / u.sr)
    assert_quantity_allclose(b(214.13747 * u.THz), 486787299458.15656 * u.MJy / u.sr)


def test_blackbody_weins_law():
    b = BlackBody(293.0 * u.K)
    assert_quantity_allclose(b.lambda_max, 9.890006672986939 * u.micron)
    assert_quantity_allclose(b.nu_max, 17.22525080856469 * u.THz)


def test_blackbody_sefanboltzman_law():
    b = BlackBody(293.0 * u.K)
    assert_quantity_allclose(b.bolometric_flux, 133.02471751812573 * u.W / (u.m * u.m))


def test_blackbody_return_units():
    # return of evaluate has no units when temperature has no units
    b = BlackBody(1000.0 * u.K, scale=1.0)
    assert not isinstance(b.evaluate(1.0 * u.micron, 1000.0, 1.0), u.Quantity)

    # return has "standard" units when scale has no units
    b = BlackBody(1000.0 * u.K, scale=1.0)
    assert isinstance(b(1.0 * u.micron), u.Quantity)
    assert b(1.0 * u.micron).unit == u.erg / (u.cm ** 2 * u.s * u.Hz * u.sr)

    # return has scale units when scale has units
    b = BlackBody(1000.0 * u.K, scale=1.0 * u.MJy / u.sr)
    assert isinstance(b(1.0 * u.micron), u.Quantity)
    assert b(1.0 * u.micron).unit == u.MJy / u.sr


@pytest.mark.skipif("not HAS_SCIPY")
def test_blackbody_fit():

    fitter = LevMarLSQFitter()

    b = BlackBody(3000 * u.K, scale=5e-17 * u.Jy / u.sr)

    wav = np.array([0.5, 5, 10]) * u.micron
    fnu = np.array([1, 10, 5]) * u.Jy / u.sr

    b_fit = fitter(b, wav, fnu, maxiter=1000)

    assert_quantity_allclose(b_fit.temperature, 2840.7438355865065 * u.K)
    assert_quantity_allclose(b_fit.scale, 5.803783292762381e-17 * u.Jy / u.sr)


def test_blackbody_overflow():
    """Test Planck function with overflow."""
    photlam = u.photon / (u.cm ** 2 * u.s * u.AA)
    wave = [0.0, 1000.0, 100000.0, 1e55]  # Angstrom
    temp = 10000.0  # Kelvin
    bb = BlackBody(temperature=temp * u.K, scale=1.0)
    with pytest.warns(
            AstropyUserWarning,
            match=r'Input contains invalid wavelength/frequency value\(s\)'):
        with np.errstate(all="ignore"):
            bb_lam = bb(wave) * u.sr
    flux = bb_lam.to(photlam, u.spectral_density(wave * u.AA)) / u.sr

    # First element is NaN, last element is very small, others normal
    assert np.isnan(flux[0])
    with np.errstate(all="ignore"):
        assert np.log10(flux[-1].value) < -134
    np.testing.assert_allclose(
        flux.value[1:-1], [0.00046368, 0.04636773], rtol=1e-3
    )  # 0.1% accuracy in PHOTLAM/sr
    with np.errstate(all="ignore"):
        flux = bb(1.0 * u.AA)
    assert flux.value == 0


def test_blackbody_exceptions_and_warnings():
    """Test exceptions."""

    # Negative temperature
    with pytest.raises(ValueError) as exc:
        bb = BlackBody(-100 * u.K)
        bb(1.0 * u.micron)
    assert exc.value.args[0] == "Temperature should be positive: [-100.] K"

    bb = BlackBody(5000 * u.K)

    # Zero wavelength given for conversion to Hz
    with catch_warnings(AstropyUserWarning) as w:
        bb(0 * u.AA)
    assert len(w) == 1
    assert "invalid" in w[0].message.args[0]

    # Negative wavelength given for conversion to Hz
    with catch_warnings(AstropyUserWarning) as w:
        bb(-1.0 * u.AA)
    assert len(w) == 1
    assert "invalid" in w[0].message.args[0]

    # Test that a non surface brightness converatable scale unit
    with pytest.raises(ValueError) as exc:
        bb = BlackBody(5000 * u.K, scale=1.0 * u.Jy)
        bb(1.0 * u.micron)
    assert exc.value.args[0] == "scale units not surface brightness: Jy"


def test_blackbody_array_temperature():
    """Regression test to make sure that the temperature can be an array."""
    multibb = BlackBody([100, 200, 300] * u.K)
    flux = multibb(1.2 * u.mm)
    np.testing.assert_allclose(
        flux.value, [1.804908e-12, 3.721328e-12, 5.638513e-12], rtol=1e-5
    )

    flux = multibb([2, 4, 6] * u.mm)
    np.testing.assert_allclose(
        flux.value, [6.657915e-13, 3.420677e-13, 2.291897e-13], rtol=1e-5
    )

    multibb = BlackBody(np.ones(4) * u.K)
    flux = multibb(np.ones((3, 4)) * u.mm)
    assert flux.shape == (3, 4)
