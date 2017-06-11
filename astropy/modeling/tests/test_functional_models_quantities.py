import pytest
import numpy as np

from ... import units as u
from ...tests.helper import assert_quantity_allclose

from ..functional_models import (Gaussian1D, GaussianAbsorption1D,
                                                Sersic1D, Sine1D, Linear1D,
                                                Lorentz1D, Voigt1D, Const1D,
                                                Box1D, Trapezoid1D, MexicanHat1D,
                                                Moffat1D)

# TODO: GaussianAbsorption1D doesn't work with units because the 1- part doesn't
# have units. How do we want to deal with that?

MODELS_1D = [
{'class': Gaussian1D,
 'parameters': {'amplitude': 3 * u.Jy, 'mean': 2 * u.m, 'stddev': 30 * u.cm},
 'evaluation':[(2600 * u.mm, 3 * u.Jy * np.exp(-2))]},
{'class': Sersic1D,
 'parameters': {'amplitude': 3 * u.MJy / u.sr, 'r_eff': 2 * u.arcsec, 'n': 4},
 'evaluation':[(3 * u.arcsec, 1.3237148119468918 * u.MJy/u.sr)]},
{'class': Sine1D,
 'parameters': {'amplitude': 3 * u.km / u.s, 'frequency': 0.25 * u.Hz, 'phase': 0.5},
 'evaluation':[(1 * u.s, -3 * u.km / u.s)]},
{'class': Linear1D,
 'parameters': {'slope': 3 * u.km / u.s, 'intercept': 5000 * u.m},
 'evaluation':[(6000 * u.ms, 23 * u.km)]},
{'class': Lorentz1D,
 'parameters': {'amplitude': 2 * u.Jy, 'x_0': 505 * u.nm, 'fwhm': 100 * u.AA},
 'evaluation':[(0.51 * u.micron, 1 * u.Jy)]},
{'class': Voigt1D,
 'parameters': {'amplitude_L': 2 * u.Jy, 'x_0': 505 * u.nm,
                'fwhm_L': 100 * u.AA, 'fwhm_G': 50 * u.AA},
 'evaluation':[(0.51 * u.micron, 1.06264568 * u.Jy)]},
{'class': Const1D,
 'parameters': {'amplitude': 3 * u.Jy},
 'evaluation':[(0.6 * u.micron, 3 * u.Jy)]},
{'class': Box1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'width': 1 * u.um},
 'evaluation':[(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)]},
{'class': Trapezoid1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'width': 1 * u.um, 'slope': 5 * u.Jy / u.um},
 'evaluation':[(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)]},
{'class': MexicanHat1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'sigma': 1e-3 * u.mm},
 'evaluation':[(1000 * u.nm, -0.09785050 * u.Jy)]},
{'class': Moffat1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'gamma': 1e-3 * u.mm, 'alpha':1},
 'evaluation':[(1000 * u.nm, 0.238853503 * u.Jy)]},
 ]

@pytest.mark.parametrize('model', MODELS_1D)
def test_1d_models_evaluatate_with_units(model):
    m = model['class'](**model['parameters'])
    for x, y in model['evaluation']:
        assert_quantity_allclose(m(x), y)


@pytest.mark.parametrize('model', MODELS_1D)
def test_1d_models_evaluatate_with_units_x_array(model):

    m = model['class'](**model['parameters'])
    for x, y in model['evaluation']:
        x_arr = u.Quantity([x, x])
        result = m(x_arr)
        assert_quantity_allclose(result, u.Quantity([y, y]))

@pytest.mark.parametrize('model', MODELS_1D)
def test_1d_models_evaluatate_with_units_param_array(model):

    params = {}
    for key, value in model['parameters'].items():
        params[key] = np.repeat(value, 2)

    m = model['class'](**params)
    for x, y in model['evaluation']:
        x_arr = u.Quantity([x, x])
        result = m(x_arr)
        assert_quantity_allclose(result, u.Quantity([y, y]))


# class Const2D(Fittable2DModel):
# class Ellipse2D(Fittable2DModel):
# class Disk2D(Fittable2DModel):
# class Ring2D(Fittable2DModel):
# class Delta1D(Fittable1DModel):
# class Delta2D(Fittable2DModel):
# class Box1D(Fittable1DModel):
#     @classmethod
# class Box2D(Fittable2DModel):
# class Trapezoid1D(Fittable1DModel):
# class TrapezoidDisk2D(Fittable2DModel):
# class MexicanHat1D(Fittable1DModel):
# class MexicanHat2D(Fittable2DModel):
# class AiryDisk2D(Fittable2DModel):
#     @classmethod
# class Moffat1D(Fittable1DModel):
# class Moffat2D(Fittable2DModel):
# class Sersic2D(Fittable2DModel):
#     @classmethod
