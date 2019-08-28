# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for physical functions."""

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


class TestBlackbody:

    # Make sure the temperature equivalency automatically applies by trying
    # to pass temperatures in celsius

    @pytest.mark.parametrize("temperature", (3000 * u.K, 2726.85 * u.deg_C))
    def test_evaluate(self, temperature):

        b = BlackBody(temperature=temperature, scale=1.0)

        assert_quantity_allclose(b(1.4 * u.micron), 486787299458.15656 * u.MJy / u.sr)
        assert_quantity_allclose(
            b(214.13747 * u.THz), 486787299458.15656 * u.MJy / u.sr
        )

    @pytest.mark.skipif("not HAS_SCIPY")
    def test_fit(self):

        fitter = LevMarLSQFitter()

        b = BlackBody(3000 * u.K)

        wav = np.array([0.5, 5, 10]) * u.micron
        fnu = np.array([1, 10, 5]) * u.Jy / u.sr

        b_fit = fitter(b, wav, fnu, maxiter=1000)

        assert_quantity_allclose(b_fit.temperature, 2840.74382317 * u.K)
        assert_quantity_allclose(b_fit.scale, 5803783.33328)


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
