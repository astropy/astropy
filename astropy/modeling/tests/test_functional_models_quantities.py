import pytest
import numpy as np

from ... import units as u
from ...tests.helper import assert_quantity_allclose

from ..functional_models import (Gaussian1D, GaussianAbsorption1D,
                                 Sersic1D, Sine1D, Linear1D,
                                 Lorentz1D, Voigt1D, Const1D,
                                 Box1D, Trapezoid1D, MexicanHat1D,
                                 Moffat1D, Gaussian2D, Const2D, Ellipse2D,
                                 Disk2D, Ring2D, Box2D, TrapezoidDisk2D,
                                 MexicanHat2D, AiryDisk2D, Moffat2D, Sersic2D)
from ..fitting import LevMarLSQFitter

# TODO: GaussianAbsorption1D doesn't work with units because the 1- part doesn't
# have units. How do we want to deal with that?

MODELS_1D = [
{'class': Gaussian1D,
 'parameters': {'amplitude': 3 * u.Jy, 'mean': 2 * u.m, 'stddev': 30 * u.cm},
 'evaluation':[(2600 * u.mm, 3 * u.Jy * np.exp(-2))],
 'bounding_box': True},
{'class': Sersic1D,
 'parameters': {'amplitude': 3 * u.MJy / u.sr, 'r_eff': 2 * u.arcsec, 'n': 4},
 'evaluation':[(3 * u.arcsec, 1.3237148119468918 * u.MJy/u.sr)],
 'bounding_box': False},
{'class': Sine1D,
 'parameters': {'amplitude': 3 * u.km / u.s, 'frequency': 0.25 * u.Hz, 'phase': 0.5},
 'evaluation':[(1 * u.s, -3 * u.km / u.s)],
 'bounding_box': False},
{'class': Linear1D,
 'parameters': {'slope': 3 * u.km / u.s, 'intercept': 5000 * u.m},
 'evaluation':[(6000 * u.ms, 23 * u.km)],
 'bounding_box': False},
{'class': Lorentz1D,
 'parameters': {'amplitude': 2 * u.Jy, 'x_0': 505 * u.nm, 'fwhm': 100 * u.AA},
 'evaluation':[(0.51 * u.micron, 1 * u.Jy)],
 'bounding_box': True},
{'class': Voigt1D,
 'parameters': {'amplitude_L': 2 * u.Jy, 'x_0': 505 * u.nm,
                'fwhm_L': 100 * u.AA, 'fwhm_G': 50 * u.AA},
 'evaluation':[(0.51 * u.micron, 1.06264568 * u.Jy)],
 'bounding_box': False},
{'class': Const1D,
 'parameters': {'amplitude': 3 * u.Jy},
 'evaluation':[(0.6 * u.micron, 3 * u.Jy)],
 'bounding_box': False},
{'class': Box1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'width': 1 * u.um},
 'evaluation':[(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)],
 'bounding_box': True},
{'class': Trapezoid1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'width': 1 * u.um, 'slope': 5 * u.Jy / u.um},
 'evaluation':[(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)],
 'bounding_box': True},
{'class': MexicanHat1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'sigma': 1e-3 * u.mm},
 'evaluation':[(1000 * u.nm, -0.09785050 * u.Jy)],
 'bounding_box': True},
{'class': Moffat1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'gamma': 1e-3 * u.mm, 'alpha':1},
 'evaluation':[(1000 * u.nm, 0.238853503 * u.Jy)],
 'bounding_box': False},
 ]

MODELS_2D = [
 {'class': Gaussian2D,
  'parameters': {'amplitude': 3 * u.Jy, 'x_mean': 2 * u.m, 'y_mean': 1 * u.m,
                 'x_stddev': 3 * u.m, 'y_stddev': 2 * u.m, 'theta': 45 * u.deg},
  'evaluation':[(412.1320343 * u.cm, 3.121320343 * u.m, 3 * u.Jy * np.exp(-0.5))],
  'bounding_box': True},
{'class': Const2D,
 'parameters': {'amplitude': 3 * u.Jy},
 'evaluation':[(0.6 * u.micron, 0.2 * u.m, 3 * u.Jy)],
 'bounding_box': False},
{'class': Disk2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0':3 * u.m, 'y_0': 2 * u.m,
                'R_0': 300 * u.cm},
 'evaluation':[(5.8 * u.m, 201 * u.cm, 3 * u.Jy)],
 'bounding_box': True},
{'class': TrapezoidDisk2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0':1 * u.m, 'y_0': 2 * u.m,
                'R_0': 100 * u.cm, 'slope':1 * u.Jy / u.m},
 'evaluation':[(3.5 * u.m, 2 * u.m, 1.5 * u.Jy)],
 'bounding_box': True},
{'class': Ellipse2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0':3 * u.m, 'y_0': 2 * u.m,
                'a': 300 * u.cm, 'b': 200 * u.cm, 'theta': 45 * u.deg},
 'evaluation':[(4 * u.m, 300 * u.cm, 3 * u.Jy)],
 'bounding_box': True},
{'class': Ring2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0':3 * u.m, 'y_0': 2 * u.m,
                'r_in': 2 * u.cm, 'r_out': 2.1 * u.cm, 'width': None},
 'evaluation':[(302.05 * u.cm, 2 * u.m + 10 * u.um, 3 * u.Jy)],
 'bounding_box': True},
{'class': Box2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0':3 * u.m, 'y_0': 2 * u.s,
                'x_width': 2 * u.cm, 'y_width': 3 * u.s},
 'evaluation':[(301 * u.cm, 3 * u.s, 3 * u.Jy)],
 'bounding_box': True},
{'class': MexicanHat2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0':3 * u.m, 'y_0': 2 * u.m,
                'sigma': 1 * u.m},
 'evaluation':[(4 * u.m, 2.5 * u.m, 0.602169107 * u.Jy)],
 'bounding_box': False},
{'class': AiryDisk2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0':3 * u.m, 'y_0': 2 * u.m,
                'radius': 1 * u.m},
 'evaluation':[(4 * u.m, 2.1 * u.m, 4.76998480e-05 * u.Jy)],
 'bounding_box': False},
{'class': Moffat2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'y_0': 3.5 * u.um,
                'gamma': 1e-3 * u.mm, 'alpha':1},
 'evaluation':[(1000 * u.nm, 2 * u.um, 0.202565833 * u.Jy)],
 'bounding_box': False},
{'class': Sersic2D,
 'parameters': {'amplitude': 3 * u.MJy / u.sr, 'x_0': 1 * u.arcsec,
                'y_0': 2 * u.arcsec, 'r_eff': 2 * u.arcsec, 'n': 4},
 'evaluation':[(3 * u.arcsec, 2.5 * u.arcsec, 2.829990489 * u.MJy/u.sr)],
 'bounding_box': False},
]

MODELS = MODELS_1D + MODELS_2D

@pytest.mark.parametrize('model', MODELS)
def test_models_evaluatate_with_units(model):
    m = model['class'](**model['parameters'])
    for args in model['evaluation']:
        assert_quantity_allclose(m(*args[:-1]), args[-1])


@pytest.mark.parametrize('model', MODELS)
def test_models_evaluatate_with_units_x_array(model):

    m = model['class'](**model['parameters'])
    for args in model['evaluation']:
        if len(args) == 2:
            x, y = args
            x_arr = u.Quantity([x, x])
            result = m(x_arr)
            assert_quantity_allclose(result, u.Quantity([y, y]))
        else:
            x, y, z = args
            x_arr = u.Quantity([x, x])
            y_arr = u.Quantity([y, y])
            result = m(x_arr, y_arr)
            assert_quantity_allclose(result, u.Quantity([z, z]))


@pytest.mark.parametrize('model', MODELS)
def test_models_evaluatate_with_units_param_array(model):

    params = {}
    for key, value in model['parameters'].items():
        if value is None:
            params[key] = None
        else:
            params[key] = np.repeat(value, 2)

    m = model['class'](**params)
    for args in model['evaluation']:
        if len(args) == 2:
            x, y = args
            x_arr = u.Quantity([x, x])
            result = m(x_arr)
            assert_quantity_allclose(result, u.Quantity([y, y]))
        else:
            x, y, z = args
            x_arr = u.Quantity([x, x])
            y_arr = u.Quantity([y, y])
            result = m(x_arr, y_arr)
            assert_quantity_allclose(result, u.Quantity([z, z]))


@pytest.mark.parametrize('model', MODELS)
def test_models_bounding_box(model):
    m = model['class'](**model['parameters'])
    if model['bounding_box']:
        m.bounding_box
    else:
        # Check that NotImplementedError is raised, so that if bounding_box is
        # implemented we remember to set bounding_box=True in the list of models
        # above
        with pytest.raises(NotImplementedError):
            m.bounding_box


@pytest.mark.parametrize('model', MODELS_1D)
def test_1d_models_fitting(model):
    m = model['class'](**model['parameters'])
    x = np.linspace(-3, 3, 100) * model['evaluation'][0][0].unit
    y = np.ones(100) * model['evaluation'][0][1].unit
    fitter = LevMarLSQFitter()
    m_new = fitter(m, x, y)
