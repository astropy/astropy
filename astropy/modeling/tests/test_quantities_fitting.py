# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to fitting models with quantity parameters
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
import pytest

from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError
from ...tests.helper import assert_quantity_allclose
from ...utils import NumpyRNGContext
from .. import fitting

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# Fitting should be as intuitive as possible to the user. Essentially, models
# and fitting should work without units, but if one has units, the other should
# have units too, and the resulting fitted parameters will also have units.


def _fake_gaussian_data():

    # Generate fake data
    with NumpyRNGContext(12345):
        x = np.linspace(-5., 5., 2000)
        y = 3 * np.exp(-0.5 * (x - 1.3)**2 / 0.8**2)
        y += np.random.normal(0., 0.2, x.shape)

    # Attach units to data
    x = x * u.m
    y = y * u.Jy

    return x, y


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitting_simple():

    x, y = _fake_gaussian_data()

    # Fit the data using a Gaussian with units
    g_init = Gaussian1D()
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)

    # TODO: update actual numerical results once implemented, but these should
    # be close to the values below.
    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitting_with_initial_values():

    x, y = _fake_gaussian_data()

    # Fit the data using a Gaussian with units
    g_init = Gaussian1D(amplitude=1. * u.mJy, mean=3 * u.cm, stddev=2 * u.mm)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)

    # TODO: update actual numerical results once implemented, but these should
    # be close to the values below.
    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitting_missing_data_units():
    """
    Raise an error if the model has units but the data doesn't
    """
    g_init = Gaussian1D(amplitude=1. * u.mJy, mean=3 * u.cm, stddev=2 * u.mm)
    fit_g = fitting.LevMarLSQFitter()

    with pytest.raises(UnitsError) as exc:
        fit_g(g_init, [1, 2, 3], [4, 5, 6])
    assert exc.value.args[0] == ("'cm' (length) and '' (dimensionless) are not "
                                 "convertible")

    with pytest.raises(UnitsError) as exc:
        fit_g(g_init, [1, 2, 3] * u.m, [4, 5, 6])
    assert exc.value.args[0] == ("'mJy' (spectral flux density) and '' "
                                 "(dimensionless) are not convertible")


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitting_missing_model_units():
    """
    Proceed if the data has units but the model doesn't
    """

    x, y = _fake_gaussian_data()

    g_init = Gaussian1D(amplitude=1., mean=3, stddev=2)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)

    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)

    g_init = Gaussian1D(amplitude=1., mean=3 * u.m, stddev=2 * u.m)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)

    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitting_incompatible_units():
    """
    Raise an error if the data and model have incompatible units
    """

    g_init = Gaussian1D(amplitude=1. * u.Jy, mean=3 * u.m, stddev=2 * u.cm)
    fit_g = fitting.LevMarLSQFitter()

    with pytest.raises(UnitsError) as exc:
        fit_g(g_init, [1, 2, 3] * u.Hz, [4, 5, 6] * u.Jy)
    assert exc.value.args[0] == ("'Hz' (frequency) and 'm' (length) are not convertible")
