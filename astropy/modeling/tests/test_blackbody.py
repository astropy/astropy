# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import pytest
import numpy as np

from ..blackbody import BlackBody1D
from ..fitting import LevMarLSQFitter

from ...tests.helper import assert_quantity_allclose
from ... import units as u

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestBlackbody1D():

    # Make sure the temperature equivalency automatically applies by trying
    # to pass temperatures in celsius

    @pytest.mark.parametrize('temperature', (3000 * u.K, 2726.85 * u.deg_C))
    def test_evaluate(self, temperature):

        bolometric_flux = 1000 * u.L_sun / (4 * np.pi * (1.5 * u.pc)**2)

        b = BlackBody1D(temperature=temperature,
                        bolometric_flux=bolometric_flux)

        assert_quantity_allclose(b(1.4 * u.micron), 4734464.498937388 * u.Jy)
        assert_quantity_allclose(b(214.13747 * u.THz), 4734464.498937388 * u.Jy)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_fit(self):

        fitter = LevMarLSQFitter()

        b = BlackBody1D(3000 * u.K)

        wav = np.array([0.5, 5, 10]) * u.micron
        fnu = np.array([1, 10, 5]) * u.Jy

        b_fit = fitter(b, wav, fnu)

        assert_quantity_allclose(b_fit.temperature, 2840.744774408546 * u.K)
        assert_quantity_allclose(b_fit.bolometric_flux, 6.821837296857152e-08 * u.erg / u.cm**2 / u.s)
