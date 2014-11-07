# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for blackbody functions."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# THIRD-PARTY
import numpy as np

# LOCAL
from ..blackbody import planck_lambda
from ... import units as u
from ... import constants as const
from ...tests.helper import pytest


__doctest_skip__ = ['*']


def test_planck_lambda():
    """Test Planck function.

    .. note:: Needs ``scipy`` to work.

    """
    try:
        from scipy import integrate
    except ImportError:
        pytest.xfail('No scipy')

    wave = np.logspace(0, 8, 1e6) * u.AA
    temp = 100. * u.K
    flux = planck_lambda(wave, temp)

    lum = wave.to(u.um)
    intflux = integrate.trapz(flux.value, x=lum.value)
    ans = const.sigma_sb * temp**4 / np.pi

    # 1% accuracy
    np.testing.assert_allclose(intflux, ans.value, rtol=0.01)
