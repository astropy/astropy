# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for blackbody functions."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# THIRD-PARTY
import numpy as np

# LOCAL
from ..blackbody import planck, FNU
from ... import units as u
from ... import constants as const
from ...tests.helper import pytest


__doctest_skip__ = ['*']


def test_planck_scipy():
    """Test Planck function.

    .. note:: Needs ``scipy`` to work.

    """
    try:
        from scipy import integrate
    except ImportError:
        pytest.xfail('No scipy')

    flux_unit = u.Watt / (u.m ** 2 * u.um)
    wave = np.logspace(0, 8, 1e6) * u.AA
    temp = 100. * u.K
    flux = planck(wave, temp, flux_unit=flux_unit)
    assert flux.unit == flux_unit / u.sr

    lum = wave.to(u.um)
    intflux = integrate.trapz(flux.value, x=lum.value)
    ans = const.sigma_sb * temp**4 / np.pi
    np.testing.assert_allclose(intflux, ans.value, rtol=0.01)  # 1% accuracy


def test_planck_overflow():
    """Test Planck function with overflow catching."""
    photlam = u.photon / (u.cm**2 * u.s * u.AA)
    wave = [0, 1000.0, 100000.0, 1e55]  # Angstrom
    temp = 10000.0  # Kelvin
    flux = planck(wave, temp, flux_unit=photlam)
    ans = [0, 3.38131732e+16, 3.87451317e+15, 0]  # PHOTLAM/sr
    assert flux.unit == photlam / u.sr
    np.testing.assert_allclose(flux.value, ans, rtol=1e-3)  # 0.1% accuracy

    flux = planck(1, 1e4)
    assert flux.value == 0

    flux = planck(10 * u.m, 1.e10)
    assert flux.value == 0


def test_planck_synphot():
    """Test that it is consistent with IRAF SYNPHOT BBFUNC."""
    # Solid angle of solar radius at 1 kpc
    fac = np.pi * (const.R_sun / const.kpc) ** 2 * u.sr

    flux = planck([1, 100, 1000, 1e4, 1e5], 5000) * fac
    assert flux.unit == FNU
    np.testing.assert_allclose(
        flux.value, [0, 0, 2.01950807e-34, 3.78584515e-26, 1.90431881e-27],
        rtol=0.01)
