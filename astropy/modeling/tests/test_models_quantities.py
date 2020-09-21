# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name, no-member
from collections import OrderedDict

import pytest
import numpy as np

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from astropy.modeling.functional_models import (
    Gaussian1D,
    Sersic1D, Sine1D, Linear1D,
    Lorentz1D, Voigt1D, Const1D,
    Box1D, Trapezoid1D, RickerWavelet1D,
    Moffat1D, Gaussian2D, Const2D, Ellipse2D,
    Disk2D, Ring2D, Box2D, TrapezoidDisk2D,
    RickerWavelet2D, AiryDisk2D, Moffat2D, Sersic2D,
    KingProjectedAnalytic1D)

from astropy.modeling.physical_models import Plummer1D

from astropy.modeling.powerlaws import (
    PowerLaw1D, BrokenPowerLaw1D, SmoothlyBrokenPowerLaw1D,
    ExponentialCutoffPowerLaw1D, LogParabola1D)

from astropy.modeling.polynomial import Polynomial1D, Polynomial2D

from astropy.modeling.fitting import LevMarLSQFitter

try:
    from scipy import optimize  # noqa
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

FUNC_MODELS_1D = [
{'class': Gaussian1D,
 'parameters': {'amplitude': 3 * u.Jy, 'mean': 2 * u.m, 'stddev': 30 * u.cm},
 'evaluation': [(2600 * u.mm, 3 * u.Jy * np.exp(-2))],
 'bounding_box': [0.35, 3.65] * u.m},
{'class': Sersic1D,
 'parameters': {'amplitude': 3 * u.MJy / u.sr, 'r_eff': 2 * u.arcsec, 'n': 4},
 'evaluation': [(3 * u.arcsec, 1.3237148119468918 * u.MJy/u.sr)],
 'bounding_box': False},
{'class': Sine1D,
 'parameters': {'amplitude': 3 * u.km / u.s, 'frequency': 0.25 * u.Hz, 'phase': 0.5},
 'evaluation': [(1 * u.s, -3 * u.km / u.s)],
 'bounding_box': False},
{'class': Linear1D,
 'parameters': {'slope': 3 * u.km / u.s, 'intercept': 5000 * u.m},
 'evaluation': [(6000 * u.ms, 23 * u.km)],
 'bounding_box': False},
{'class': Lorentz1D,
 'parameters': {'amplitude': 2 * u.Jy, 'x_0': 505 * u.nm, 'fwhm': 100 * u.AA},
 'evaluation': [(0.51 * u.micron, 1 * u.Jy)],
 'bounding_box': [255, 755] * u.nm},
{'class': Voigt1D,
 'parameters': {'amplitude_L': 2 * u.Jy, 'x_0': 505 * u.nm,
                'fwhm_L': 100 * u.AA, 'fwhm_G': 50 * u.AA},
 'evaluation': [(0.51 * u.micron, 1.06264568 * u.Jy)],
 'bounding_box': False},
{'class': Const1D,
 'parameters': {'amplitude': 3 * u.Jy},
 'evaluation': [(0.6 * u.micron, 3 * u.Jy)],
 'bounding_box': False},
{'class': Box1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'width': 1 * u.um},
 'evaluation': [(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)],
 'bounding_box': [3.9, 4.9] * u.um},
{'class': Trapezoid1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'width': 1 * u.um, 'slope': 5 * u.Jy / u.um},
 'evaluation': [(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)],
 'bounding_box': [3.3, 5.5] * u.um},
{'class': RickerWavelet1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'sigma': 1e-3 * u.mm},
 'evaluation': [(1000 * u.nm, -0.09785050 * u.Jy)],
 'bounding_box': [-5.6, 14.4] * u.um},
{'class': Moffat1D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'gamma': 1e-3 * u.mm, 'alpha': 1},
 'evaluation': [(1000 * u.nm, 0.238853503 * u.Jy)],
 'bounding_box': False},
{'class': KingProjectedAnalytic1D,
 'parameters': {'amplitude': 1. * u.Msun/u.pc**2, 'r_core': 1. * u.pc, 'r_tide': 2. * u.pc},
 'evaluation': [(0.5 * u.pc, 0.2 * u.Msun/u.pc**2)],
 'bounding_box': [0. * u.pc, 2. * u.pc]}
 ]

PHYS_MODELS_1D = [
{'class': Plummer1D,
 'parameters': {'mass': 3 * u.kg, 'r_plum': 0.5 * u.m},
 'evaluation': [(1* u.m, 0.10249381 * u.kg / (u.m **3))],
 'bounding_box': False}
 ]

FUNC_MODELS_2D = [
{'class': Gaussian2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_mean': 2 * u.m, 'y_mean': 1 * u.m,
                'x_stddev': 3 * u.m, 'y_stddev': 2 * u.m, 'theta': 45 * u.deg},
 'evaluation': [(412.1320343 * u.cm, 3.121320343 * u.m, 3 * u.Jy * np.exp(-0.5))],
 'bounding_box': [[-14.18257445, 16.18257445], [-10.75693665, 14.75693665]] * u.m},
{'class': Const2D,
 'parameters': {'amplitude': 3 * u.Jy},
 'evaluation': [(0.6 * u.micron, 0.2 * u.m, 3 * u.Jy)],
 'bounding_box': False},
{'class': Disk2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 3 * u.m, 'y_0': 2 * u.m,
                'R_0': 300 * u.cm},
 'evaluation': [(5.8 * u.m, 201 * u.cm, 3 * u.Jy)],
 'bounding_box': [[-1, 5], [0, 6]] * u.m},
{'class': TrapezoidDisk2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 1 * u.m, 'y_0': 2 * u.m,
                'R_0': 100 * u.cm, 'slope': 1 * u.Jy / u.m},
 'evaluation': [(3.5 * u.m, 2 * u.m, 1.5 * u.Jy)],
 'bounding_box': [[-2, 6], [-3, 5]] * u.m},
{'class': Ellipse2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 3 * u.m, 'y_0': 2 * u.m,
                'a': 300 * u.cm, 'b': 200 * u.cm, 'theta': 45 * u.deg},
 'evaluation': [(4 * u.m, 300 * u.cm, 3 * u.Jy)],
 'bounding_box': [[-0.76046808, 4.76046808], [0.68055697, 5.31944302]] * u.m},
{'class': Ring2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 3 * u.m, 'y_0': 2 * u.m,
                'r_in': 2 * u.cm, 'r_out': 2.1 * u.cm},
 'evaluation': [(302.05 * u.cm, 2 * u.m + 10 * u.um, 3 * u.Jy)],
 'bounding_box': [[1.979, 2.021], [2.979, 3.021]] * u.m},
{'class': Box2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 3 * u.m, 'y_0': 2 * u.s,
                'x_width': 4 * u.cm, 'y_width': 3 * u.s},
 'evaluation': [(301 * u.cm, 3 * u.s, 3 * u.Jy)],
 'bounding_box': [[0.5 * u.s, 3.5 * u.s], [2.98 * u.m, 3.02 * u.m]]},
{'class': RickerWavelet2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 3 * u.m, 'y_0': 2 * u.m,
                'sigma': 1 * u.m},
 'evaluation': [(4 * u.m, 2.5 * u.m, 0.602169107 * u.Jy)],
 'bounding_box': False},
{'class': AiryDisk2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 3 * u.m, 'y_0': 2 * u.m,
                'radius': 1 * u.m},
 'evaluation': [(4 * u.m, 2.1 * u.m, 4.76998480e-05 * u.Jy)],
 'bounding_box': False},
{'class': Moffat2D,
 'parameters': {'amplitude': 3 * u.Jy, 'x_0': 4.4 * u.um, 'y_0': 3.5 * u.um,
                'gamma': 1e-3 * u.mm, 'alpha': 1},
 'evaluation': [(1000 * u.nm, 2 * u.um, 0.202565833 * u.Jy)],
 'bounding_box': False},
{'class': Sersic2D,
 'parameters': {'amplitude': 3 * u.MJy / u.sr, 'x_0': 1 * u.arcsec,
                'y_0': 2 * u.arcsec, 'r_eff': 2 * u.arcsec, 'n': 4,
                'ellip': 0, 'theta': 0},
 'evaluation': [(3 * u.arcsec, 2.5 * u.arcsec, 2.829990489 * u.MJy/u.sr)],
 'bounding_box': False},
]

POWERLAW_MODELS = [
{'class': PowerLaw1D,
 'parameters': {'amplitude': 5 * u.kg, 'x_0': 10 * u.cm, 'alpha': 1},
 'evaluation': [(1 * u.m, 500 * u.g)],
 'bounding_box': False},
{'class': BrokenPowerLaw1D,
 'parameters': {'amplitude': 5 * u.kg, 'x_break': 10 * u.cm, 'alpha_1': 1, 'alpha_2': -1},
 'evaluation': [(1 * u.m, 50 * u.kg), (1 * u.cm, 50 * u.kg)],
 'bounding_box': False},
{'class': SmoothlyBrokenPowerLaw1D,
 'parameters': {'amplitude': 5 * u.kg, 'x_break': 10 * u.cm, 'alpha_1': 1, 'alpha_2': -1, 'delta': 1},
 'evaluation': [(1 * u.cm, 15.125 * u.kg), (1 * u.m, 15.125 * u.kg)],
 'bounding_box': False},
{'class': ExponentialCutoffPowerLaw1D,
 'parameters': {'amplitude': 5 * u.kg, 'x_0': 10 * u.cm, 'alpha': 1, 'x_cutoff': 1 * u.m},
 'evaluation': [(1 * u.um, 499999.5 * u.kg), (10 * u.m, 50 * np.exp(-10) * u.g)],
 'bounding_box': False},
{'class': LogParabola1D,
 'parameters': {'amplitude': 5 * u.kg, 'x_0': 10 * u.cm, 'alpha': 1, 'beta': 2},
 'evaluation': [(1 * u.cm, 5 * 0.1 ** (-1 - 2 * np.log(0.1)) * u.kg)],
 'bounding_box': False}
]

POLY_MODELS = [
{'class': Polynomial1D,
 'parameters': {'degree': 2, 'c0': 3 * u.one, 'c1': 2 / u.m, 'c2': 3 / u.m**2},
 'evaluation': [(3 * u.m, 36 * u.one)],
 'bounding_box': False},
{'class': Polynomial1D,
 'parameters': {'degree': 2, 'c0': 3 * u.kg, 'c1': 2 * u.kg / u.m, 'c2': 3 * u.kg / u.m**2},
 'evaluation': [(3 * u.m, 36 * u.kg)],
 'bounding_box': False},
{'class': Polynomial1D,
 'parameters': {'degree': 2, 'c0': 3 * u.kg, 'c1': 2 * u.kg, 'c2': 3 * u.kg},
 'evaluation': [(3 * u.one, 36 * u.kg)],
 'bounding_box': False},
{'class': Polynomial2D,
 'parameters': {'degree': 2, 'c0_0': 3 * u.one, 'c1_0': 2 / u.m, 'c2_0': 3 / u.m**2,
                'c0_1': 3 / u.s, 'c0_2': -2 / u.s**2, 'c1_1': 5 / u.m / u.s},
 'evaluation': [(3 * u.m, 2 * u.s, 64 * u.one)],
 'bounding_box': False},
{'class': Polynomial2D,
 'parameters': {'degree': 2, 'c0_0': 3 * u.kg, 'c1_0': 2 * u.kg / u.m, 'c2_0': 3 * u.kg / u.m**2,
                'c0_1': 3 * u.kg / u.s, 'c0_2': -2 * u.kg / u.s**2, 'c1_1': 5 * u.kg / u.m / u.s},
 'evaluation': [(3 * u.m, 2 * u.s, 64 * u.kg)],
 'bounding_box': False},
{'class': Polynomial2D,
 'parameters': {'degree': 2, 'c0_0': 3 * u.kg, 'c1_0': 2 * u.kg, 'c2_0': 3 * u.kg,
                'c0_1': 3 * u.kg, 'c0_2': -2 * u.kg, 'c1_1': 5 * u.kg},
 'evaluation': [(3 * u.one, 2 * u.one, 64 * u.kg)],
 'bounding_box': False},
 ]


MODELS = FUNC_MODELS_1D + FUNC_MODELS_2D + POWERLAW_MODELS + PHYS_MODELS_1D

SCIPY_MODELS = set([Sersic1D, Sersic2D, AiryDisk2D])


@pytest.mark.parametrize('model', MODELS)
def test_models_evaluate_without_units(model):
    if not HAS_SCIPY and model['class'] in SCIPY_MODELS:
        pytest.skip()
    m = model['class'](**model['parameters'])
    for args in model['evaluation']:
        if len(args) == 2:
            kwargs = OrderedDict(zip(('x', 'y'), args))
        else:
            kwargs = OrderedDict(zip(('x', 'y', 'z'), args))
            if kwargs['x'].unit.is_equivalent(kwargs['y'].unit):
                kwargs['x'] = kwargs['x'].to(kwargs['y'].unit)
        mnu = m.without_units_for_data(**kwargs)
        args = [x.value for x in kwargs.values()]
        assert_quantity_allclose(mnu(*args[:-1]), args[-1])


@pytest.mark.parametrize('model', MODELS)
def test_models_evaluate_with_units(model):
    if not HAS_SCIPY and model['class'] in SCIPY_MODELS:
        pytest.skip()
    m = model['class'](**model['parameters'])
    for args in model['evaluation']:
        assert_quantity_allclose(m(*args[:-1]), args[-1])


@pytest.mark.parametrize('model', MODELS)
def test_models_evaluate_with_units_x_array(model):

    if not HAS_SCIPY and model['class'] in SCIPY_MODELS:
        pytest.skip()

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
def test_models_evaluate_with_units_param_array(model):

    if not HAS_SCIPY and model['class'] in SCIPY_MODELS:
        pytest.skip()

    params = {}
    for key, value in model['parameters'].items():
        if value is None or key == 'degree':
            params[key] = value
        else:
            params[key] = np.repeat(value, 2)

    params['n_models'] = 2

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

    # In some cases, having units in parameters caused bounding_box to break,
    # so this is to ensure that it works correctly.

    if not HAS_SCIPY and model['class'] in SCIPY_MODELS:
        pytest.skip()

    m = model['class'](**model['parameters'])

    # In the following we need to explicitly test that the value is False
    # since Quantities no longer evaluate as as True
    if model['bounding_box'] is False:
        # Check that NotImplementedError is raised, so that if bounding_box is
        # implemented we remember to set bounding_box=True in the list of models
        # above
        with pytest.raises(NotImplementedError):
            m.bounding_box
    else:
        # A bounding box may have inhomogeneous units so we need to check the
        # values one by one.
        for i in range(len(model['bounding_box'])):
            bbox = m.bounding_box
            assert_quantity_allclose(bbox[i], model['bounding_box'][i])


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.filterwarnings(r'ignore:.*:RuntimeWarning')
@pytest.mark.filterwarnings(r'ignore:Model is linear in parameters.*')
@pytest.mark.filterwarnings(r'ignore:The fit may be unsuccessful.*')
@pytest.mark.parametrize('model', MODELS)
def test_models_fitting(model):

    m = model['class'](**model['parameters'])
    if len(model['evaluation'][0]) == 2:
        x = np.linspace(1, 3, 100) * model['evaluation'][0][0].unit
        y = np.exp(-x.value ** 2) * model['evaluation'][0][1].unit
        args = [x, y]
    else:
        x = np.linspace(1, 3, 100) * model['evaluation'][0][0].unit
        y = np.linspace(1, 3, 100) * model['evaluation'][0][1].unit
        z = np.exp(-x.value**2 - y.value**2) * model['evaluation'][0][2].unit
        args = [x, y, z]

    # Test that the model fits even if it has units on parameters
    fitter = LevMarLSQFitter()
    m_new = fitter(m, *args)

    # Check that units have been put back correctly
    for param_name in m.param_names:
        par_bef = getattr(m, param_name)
        par_aft = getattr(m_new, param_name)
        if par_bef.unit is None:
            # If the parameter used to not have a unit then had a radian unit
            # for example, then we should allow that
            assert par_aft.unit is None or par_aft.unit is u.rad
        else:
            assert par_aft.unit.is_equivalent(par_bef.unit)
