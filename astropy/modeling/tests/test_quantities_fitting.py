# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests that relate to fitting models with quantity parameters
"""
import numpy as np
import pytest

from astropy import units as u
from astropy.modeling import fitting, models
from astropy.modeling.core import Fittable1DModel
from astropy.modeling.parameters import Parameter
from astropy.tests.helper import assert_quantity_allclose
from astropy.units import UnitsError
from astropy.utils import NumpyRNGContext
from astropy.utils.compat.optional_deps import HAS_SCIPY

# Fitting should be as intuitive as possible to the user. Essentially, models
# and fitting should work without units, but if one has units, the other should
# have units too, and the resulting fitted parameters will also have units.

fitters = [
    fitting.LevMarLSQFitter,
    fitting.TRFLSQFitter,
    fitting.LMLSQFitter,
    fitting.DogBoxLSQFitter,
]


def _fake_gaussian_data():
    # Generate fake data
    with NumpyRNGContext(12345):
        x = np.linspace(-5.0, 5.0, 2000)
        y = 3 * np.exp(-0.5 * (x - 1.3) ** 2 / 0.8**2)
        y += np.random.normal(0.0, 0.2, x.shape)

    # Attach units to data
    x = x * u.m
    y = y * u.Jy

    return x, y


compound_models_no_units = [
    models.Linear1D() + models.Gaussian1D() + models.Gaussian1D(),
    models.Linear1D() + models.Gaussian1D() | models.Scale(),
    models.Linear1D() + models.Gaussian1D() | models.Shift(),
]


class CustomInputNamesModel(Fittable1DModel):
    n_inputs = 1
    n_outputs = 1

    a = Parameter(default=1.0)
    b = Parameter(default=1.0)

    def __init__(self, a=a, b=b):
        super().__init__(a=a, b=b)
        self.inputs = ("inn",)
        self.outputs = ("out",)

    @staticmethod
    def evaluate(inn, a, b):
        return a * inn + b

    @property
    def input_units(self):
        if self.a.unit is None and self.b.unit is None:
            return None
        else:
            return {"inn": self.b.unit / self.a.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {"a": outputs_unit["out"] / inputs_unit["inn"], "b": outputs_unit["out"]}


def models_with_custom_names():
    line = models.Linear1D(1 * u.m / u.s, 2 * u.m)
    line.inputs = ("inn",)
    line.outputs = ("out",)

    custom_names_model = CustomInputNamesModel(1 * u.m / u.s, 2 * u.m)
    return [line, custom_names_model]


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_fitting_simple(fitter):
    fitter = fitter()

    x, y = _fake_gaussian_data()

    # Fit the data using a Gaussian with units
    g_init = models.Gaussian1D()
    g = fitter(g_init, x, y)

    # TODO: update actual numerical results once implemented, but these should
    # be close to the values below.
    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_fitting_with_initial_values(fitter):
    fitter = fitter()

    x, y = _fake_gaussian_data()

    # Fit the data using a Gaussian with units
    g_init = models.Gaussian1D(amplitude=1.0 * u.mJy, mean=3 * u.cm, stddev=2 * u.mm)
    g = fitter(g_init, x, y)

    # TODO: update actual numerical results once implemented, but these should
    # be close to the values below.
    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_fitting_missing_data_units(fitter):
    """
    Raise an error if the model has units but the data doesn't
    """
    fitter = fitter()

    class UnorderedGaussian1D(models.Gaussian1D):
        # Parameters are ordered differently here from Gaussian1D
        # to ensure the order does not break functionality.
        def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
            return {
                "amplitude": outputs_unit["y"],
                "mean": inputs_unit["x"],
                "stddev": inputs_unit["x"],
            }

    g_init = UnorderedGaussian1D(amplitude=1.0 * u.mJy, mean=3 * u.cm, stddev=2 * u.mm)
    # We define flux unit so that conversion fails at wavelength unit.
    # This is because the order of parameter unit conversion seems to
    # follow the order defined in _parameter_units_for_data_units method.
    MESSAGE = r"'cm' .* and '' .* are not convertible"
    with pytest.raises(UnitsError, match=MESSAGE):
        fitter(g_init, [1, 2, 3], [4, 5, 6] * (u.erg / (u.s * u.cm * u.cm * u.Hz)))

    MESSAGE = r"'mJy' .* and '' .* are not convertible"
    with pytest.raises(UnitsError, match=MESSAGE):
        fitter(g_init, [1, 2, 3] * u.m, [4, 5, 6])


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_fitting_missing_model_units(fitter):
    """
    Proceed if the data has units but the model doesn't
    """
    fitter = fitter()

    x, y = _fake_gaussian_data()

    g_init = models.Gaussian1D(amplitude=1.0, mean=3, stddev=2)
    g = fitter(g_init, x, y)

    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)

    g_init = models.Gaussian1D(amplitude=1.0, mean=3 * u.m, stddev=2 * u.m)
    g = fitter(g_init, x, y)

    assert_quantity_allclose(g.amplitude, 3 * u.Jy, rtol=0.05)
    assert_quantity_allclose(g.mean, 1.3 * u.m, rtol=0.05)
    assert_quantity_allclose(g.stddev, 0.8 * u.m, rtol=0.05)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_fitting_incompatible_units(fitter):
    """
    Raise an error if the data and model have incompatible units
    """
    fitter = fitter()

    g_init = models.Gaussian1D(amplitude=1.0 * u.Jy, mean=3 * u.m, stddev=2 * u.cm)
    MESSAGE = r"'Hz' .* and 'm' .* are not convertible"
    with pytest.raises(UnitsError, match=MESSAGE):
        fitter(g_init, [1, 2, 3] * u.Hz, [4, 5, 6] * u.Jy)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.filterwarnings(r"ignore:The fit may be unsuccessful.*")
@pytest.mark.filterwarnings(r"ignore:divide by zero encountered.*")
@pytest.mark.parametrize("model", compound_models_no_units)
@pytest.mark.parametrize("fitter", fitters)
def test_compound_without_units(model, fitter):
    fitter = fitter()

    x = np.linspace(-5, 5, 10) * u.Angstrom
    with NumpyRNGContext(12345):
        y = np.random.sample(10)

    res_fit = fitter(model, x, y * u.Hz)
    for param_name in res_fit.param_names:
        print(getattr(res_fit, param_name))
    assert all([res_fit[i]._has_units for i in range(3)])
    z = res_fit(x)
    assert isinstance(z, u.Quantity)

    res_fit = fitter(model, np.arange(10) * u.Unit("Angstrom"), y)
    assert all([res_fit[i]._has_units for i in range(3)])
    z = res_fit(x)
    assert isinstance(z, np.ndarray)


# FIXME: See https://github.com/astropy/astropy/issues/10675
# @pytest.mark.skipif(not HAS_SCIPY, reason='requires scipy')
@pytest.mark.skip(reason="Flaky and ill-conditioned")
@pytest.mark.parametrize("fitter", fitters)
def test_compound_fitting_with_units(fitter):
    fitter = fitter()

    x = np.linspace(-5, 5, 15) * u.Angstrom
    y = np.linspace(-5, 5, 15) * u.Angstrom

    fitter = fitter()
    m = models.Gaussian2D(
        10 * u.Hz, 3 * u.Angstrom, 4 * u.Angstrom, 1 * u.Angstrom, 2 * u.Angstrom
    )
    p = models.Planar2D(3 * u.Hz / u.Angstrom, 4 * u.Hz / u.Angstrom, 1 * u.Hz)
    model = m + p

    z = model(x, y)
    res = fitter(model, x, y, z)
    assert isinstance(res(x, y), np.ndarray)
    assert all([res[i]._has_units for i in range(2)])

    model = models.Gaussian2D() + models.Planar2D()
    res = fitter(model, x, y, z)
    assert isinstance(res(x, y), np.ndarray)
    assert all([res[i]._has_units for i in range(2)])

    # A case of a mixture of models with and without units
    model = models.BlackBody(temperature=3000 * u.K) * models.Const1D(amplitude=1.0)
    x = np.linspace(1, 3, 10000) * u.micron

    with NumpyRNGContext(12345):
        n = np.random.normal(3)

    y = model(x)
    res = fitter(model, x, y * (1 + n))
    # The large rtol here is due to different results on linux and macosx, likely
    # the model is ill-conditioned.
    np.testing.assert_allclose(
        res.parameters, [3000, 2.1433621e00, 2.647347e00], rtol=0.4
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.filterwarnings(r"ignore:Model is linear in parameters*")
@pytest.mark.parametrize("model", models_with_custom_names())
@pytest.mark.parametrize("fitter", fitters)
def test_fitting_custom_names(model, fitter):
    """Tests fitting of models with custom inputs and outsputs names."""
    fitter = fitter()

    x = np.linspace(0, 10, 100) * u.s
    y = model(x)
    new_model = fitter(model, x, y)
    for param_name in model.param_names:
        assert_quantity_allclose(
            getattr(new_model, param_name).quantity, getattr(model, param_name).quantity
        )
