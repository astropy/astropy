# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to fitting models with quantity parameters
"""
import numpy as np
import pytest


from astropy.modeling import models
from astropy import units as u
from astropy.units import UnitsError
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils import NumpyRNGContext
from astropy.modeling import fitting


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


bad_compound_models_no_units = [
    models.Linear1D() + models.Gaussian1D() | models.Scale(),
    models.Linear1D() + models.Gaussian1D() | models.Shift()
]

compound_models_no_units = [
    models.Linear1D() + models.Gaussian1D() + models.Gaussian1D()
]


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitting_simple():

    x, y = _fake_gaussian_data()

    # Fit the data using a Gaussian with units
    g_init = models.Gaussian1D()
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
    g_init = models.Gaussian1D(amplitude=1. * u.mJy,
                               mean=3 * u.cm,
                               stddev=2 * u.mm)
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
    g_init = models.Gaussian1D(amplitude=1. * u.mJy,
                               mean=3 * u.cm,
                               stddev=2 * u.mm)
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

    g_init = models.Gaussian1D(amplitude=1., mean=3, stddev=2)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)

    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)

    g_init = models.Gaussian1D(amplitude=1., mean=3 * u.m, stddev=2 * u.m)
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

    g_init = models.Gaussian1D(amplitude=1. * u.Jy,
                               mean=3 * u.m,
                               stddev=2 * u.cm)
    fit_g = fitting.LevMarLSQFitter()

    with pytest.raises(UnitsError) as exc:
        fit_g(g_init, [1, 2, 3] * u.Hz, [4, 5, 6] * u.Jy)
    assert exc.value.args[0] == ("'Hz' (frequency) and 'm' (length) are not convertible")


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.parametrize('model', compound_models_no_units)
def test_compound_without_units(model):
    x = np.linspace(-5, 5, 10) * u.Angstrom
    with NumpyRNGContext(12345):
        y = np.random.sample(10)

    fitter = fitting.LevMarLSQFitter()
    res_fit = fitter(model, x, y * u.Hz)
    for param_name in res_fit.param_names:
        print(getattr(res_fit, param_name))
    assert all([res_fit[i]._has_units for i in range(3)])
    z = res_fit(x)
    assert isinstance(z, u.Quantity)

    res_fit = fitter(model, np.arange(10) * u.Unit('Angstrom'), y)
    assert all([res_fit[i]._has_units for i in range(3)])
    z = res_fit(x)
    assert isinstance(z, np.ndarray)


@pytest.mark.skipif('not HAS_SCIPY')
def test_compound_fitting_with_units():
    x = np.linspace(-5, 5, 15) * u.Angstrom
    y = np.linspace(-5, 5, 15) * u.Angstrom

    fitter = fitting.LevMarLSQFitter()
    m = models.Gaussian2D(10*u.Hz,
                          3*u.Angstrom, 4*u.Angstrom,
                          1*u.Angstrom, 2*u.Angstrom)
    p = models.Planar2D(3*u.Hz/u.Angstrom, 4*u.Hz/u.Angstrom, 1*u.Hz)
    model = m + p

    z = model(x, y)
    res = fitter(model, x, y, z)
    assert isinstance(res(x, y), np.ndarray)
    assert all([res[i]._has_units for i in range(2)])

    model = models.Gaussian2D() + models.Planar2D()
    res = fitter(model, x, y, z)
    assert isinstance(res(x, y), np.ndarray)
    assert all([res[i]._has_units for i in range(2)])


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.parametrize('model', bad_compound_models_no_units)
def test_bad_compound_without_units(model):
    with pytest.raises(ValueError):
        x = np.linspace(-5, 5, 10) * u.Angstrom
        with NumpyRNGContext(12345):
            y = np.random.sample(10)

        fitter = fitting.LevMarLSQFitter()
        res_fit = fitter(model, x, y * u.Hz)
