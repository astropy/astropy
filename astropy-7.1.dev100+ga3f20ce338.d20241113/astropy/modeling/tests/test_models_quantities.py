# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name, no-member

import numpy as np
import pytest

from astropy import units as u
from astropy.modeling.bounding_box import ModelBoundingBox
from astropy.modeling.core import fix_inputs
from astropy.modeling.fitting import (
    DogBoxLSQFitter,
    LevMarLSQFitter,
    LMLSQFitter,
    TRFLSQFitter,
)
from astropy.modeling.functional_models import (
    AiryDisk2D,
    ArcCosine1D,
    ArcSine1D,
    ArcTangent1D,
    Box1D,
    Box2D,
    Const1D,
    Const2D,
    Cosine1D,
    Disk2D,
    Ellipse2D,
    Exponential1D,
    Gaussian1D,
    Gaussian2D,
    GeneralSersic2D,
    KingProjectedAnalytic1D,
    Linear1D,
    Logarithmic1D,
    Lorentz1D,
    Lorentz2D,
    Moffat1D,
    Moffat2D,
    Multiply,
    Planar2D,
    RickerWavelet1D,
    RickerWavelet2D,
    Ring2D,
    Scale,
    Sersic1D,
    Sersic2D,
    Sine1D,
    Tangent1D,
    Trapezoid1D,
    TrapezoidDisk2D,
    Voigt1D,
)
from astropy.modeling.parameters import InputParameterError
from astropy.modeling.physical_models import Drude1D, Plummer1D
from astropy.modeling.polynomial import Polynomial1D, Polynomial2D
from astropy.modeling.powerlaws import (
    BrokenPowerLaw1D,
    ExponentialCutoffPowerLaw1D,
    LogParabola1D,
    PowerLaw1D,
    Schechter1D,
    SmoothlyBrokenPowerLaw1D,
)
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY

fitters = [LevMarLSQFitter, TRFLSQFitter, LMLSQFitter, DogBoxLSQFitter]

FUNC_MODELS_1D = [
    {
        "class": Gaussian1D,
        "parameters": {"amplitude": 3 * u.Jy, "mean": 2 * u.m, "stddev": 30 * u.cm},
        "evaluation": [(2600 * u.mm, 3 * u.Jy * np.exp(-2))],
        "bounding_box": [0.35, 3.65] * u.m,
    },
    {
        "class": Sersic1D,
        "parameters": {"amplitude": 3 * u.MJy / u.sr, "r_eff": 2 * u.arcsec, "n": 4},
        "evaluation": [(3 * u.arcsec, 1.3237148119468918 * u.MJy / u.sr)],
        "bounding_box": False,
    },
    {
        "class": Sine1D,
        "parameters": {
            "amplitude": 3 * u.km / u.s,
            "frequency": 0.25 * u.Hz,
            "phase": 0.5,
        },
        "evaluation": [(1 * u.s, -3 * u.km / u.s)],
        "bounding_box": False,
    },
    {
        "class": Cosine1D,
        "parameters": {
            "amplitude": 3 * u.km / u.s,
            "frequency": 0.25 * u.Hz,
            "phase": 0.25,
        },
        "evaluation": [(1 * u.s, -3 * u.km / u.s)],
        "bounding_box": False,
    },
    {
        "class": Tangent1D,
        "parameters": {
            "amplitude": 3 * u.km / u.s,
            "frequency": 0.125 * u.Hz,
            "phase": 0.25,
        },
        "evaluation": [(1 * u.s, -3 * u.km / u.s)],
        "bounding_box": [-4, 0] / u.Hz,
    },
    {
        "class": ArcSine1D,
        "parameters": {
            "amplitude": 3 * u.km / u.s,
            "frequency": 0.25 * u.Hz,
            "phase": 0.5,
        },
        "evaluation": [(0 * u.km / u.s, -2 * u.s)],
        "bounding_box": [-3, 3] * u.km / u.s,
    },
    {
        "class": ArcCosine1D,
        "parameters": {
            "amplitude": 3 * u.km / u.s,
            "frequency": 0.25 * u.Hz,
            "phase": 0.5,
        },
        "evaluation": [(0 * u.km / u.s, -1 * u.s)],
        "bounding_box": [-3, 3] * u.km / u.s,
    },
    {
        "class": ArcTangent1D,
        "parameters": {
            "amplitude": 3 * u.km / u.s,
            "frequency": 0.125 * u.Hz,
            "phase": 0.25,
        },
        "evaluation": [(0 * u.km / u.s, -2 * u.s)],
        "bounding_box": False,
    },
    {
        "class": Linear1D,
        "parameters": {"slope": 3 * u.km / u.s, "intercept": 5000 * u.m},
        "evaluation": [(6000 * u.ms, 23 * u.km)],
        "bounding_box": False,
    },
    {
        "class": Lorentz1D,
        "parameters": {"amplitude": 2 * u.Jy, "x_0": 505 * u.nm, "fwhm": 100 * u.AA},
        "evaluation": [(0.51 * u.micron, 1 * u.Jy)],
        "bounding_box": [255, 755] * u.nm,
    },
    {
        "class": Voigt1D,
        "parameters": {
            "amplitude_L": 2 * u.Jy,
            "x_0": 505 * u.nm,
            "fwhm_L": 100 * u.AA,
            "fwhm_G": 50 * u.AA,
        },
        "evaluation": [(0.51 * u.micron, 1.0621795524 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": Voigt1D,
        "parameters": {
            "amplitude_L": 2 * u.Jy,
            "x_0": 505 * u.nm,
            "fwhm_L": 100 * u.AA,
            "fwhm_G": 50 * u.AA,
            "method": "humlicek2",
        },
        "evaluation": [(0.51 * u.micron, 1.0621795524 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": Const1D,
        "parameters": {"amplitude": 3 * u.Jy},
        "evaluation": [(0.6 * u.micron, 3 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": Box1D,
        "parameters": {"amplitude": 3 * u.Jy, "x_0": 4.4 * u.um, "width": 1 * u.um},
        "evaluation": [(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)],
        "bounding_box": [3.9, 4.9] * u.um,
    },
    {
        "class": Trapezoid1D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 4.4 * u.um,
            "width": 1 * u.um,
            "slope": 5 * u.Jy / u.um,
        },
        "evaluation": [(4200 * u.nm, 3 * u.Jy), (1 * u.m, 0 * u.Jy)],
        "bounding_box": [3.3, 5.5] * u.um,
    },
    {
        "class": RickerWavelet1D,
        "parameters": {"amplitude": 3 * u.Jy, "x_0": 4.4 * u.um, "sigma": 1e-3 * u.mm},
        "evaluation": [(1000 * u.nm, -0.09785050 * u.Jy)],
        "bounding_box": [-5.6, 14.4] * u.um,
    },
    {
        "class": Moffat1D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 4.4 * u.um,
            "gamma": 1e-3 * u.mm,
            "alpha": 1,
        },
        "evaluation": [(1000 * u.nm, 0.238853503 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": KingProjectedAnalytic1D,
        "parameters": {
            "amplitude": 1.0 * u.Msun / u.pc**2,
            "r_core": 1.0 * u.pc,
            "r_tide": 2.0 * u.pc,
        },
        "evaluation": [(0.5 * u.pc, 0.2 * u.Msun / u.pc**2)],
        "bounding_box": [0.0 * u.pc, 2.0 * u.pc],
    },
    {
        "class": Logarithmic1D,
        "parameters": {"amplitude": 5 * u.m, "tau": 2 * u.m},
        "evaluation": [(4 * u.m, 3.4657359027997265 * u.m)],
        "bounding_box": False,
    },
    {
        "class": Exponential1D,
        "parameters": {"amplitude": 5 * u.m, "tau": 2 * u.m},
        "evaluation": [(4 * u.m, 36.945280494653254 * u.m)],
        "bounding_box": False,
    },
]

SCALE_MODELS = [
    {
        "class": Scale,
        "parameters": {"factor": 2 * u.m},
        "evaluation": [(1 * u.m, 2 * u.m)],
        "bounding_box": False,
    },
    {
        "class": Multiply,
        "parameters": {"factor": 2 * u.m},
        "evaluation": [(1 * u.m / u.m, 2 * u.m)],
        "bounding_box": False,
    },
]

PHYS_MODELS_1D = [
    {
        "class": Plummer1D,
        "parameters": {"mass": 3 * u.kg, "r_plum": 0.5 * u.m},
        "evaluation": [(1 * u.m, 0.10249381 * u.kg / (u.m**3))],
        "bounding_box": False,
    },
    {
        "class": Drude1D,
        "parameters": {
            "amplitude": 1.0 * u.m,
            "x_0": 2175.0 * u.AA,
            "fwhm": 400.0 * u.AA,
        },
        "evaluation": [(2000 * u.AA, 0.5452317018423869 * u.m)],
        "bounding_box": [-17825, 22175] * u.AA,
    },
]

FUNC_MODELS_2D = [
    {
        "class": Gaussian2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_mean": 2 * u.m,
            "y_mean": 1 * u.m,
            "x_stddev": 3 * u.m,
            "y_stddev": 2 * u.m,
            "theta": 45 * u.deg,
        },
        "evaluation": [
            (412.1320343 * u.cm, 3.121320343 * u.m, 3 * u.Jy * np.exp(-0.5))
        ],
        "bounding_box": [[-13.02230366, 15.02230366], [-12.02230366, 16.02230366]]
        * u.m,
    },
    {
        "class": Lorentz2D,
        "parameters": {
            "amplitude": 2 * u.Jy,
            "x_0": 505 * u.nm,
            "y_0": 507 * u.nm,
            "fwhm": 100 * u.AA,
        },
        "evaluation": [(0.51 * u.micron, 0.53 * u.micron, 0.08635579 * u.Jy)],
        "bounding_box": [[255, 755], [257, 757]] * u.nm,
    },
    {
        "class": Const2D,
        "parameters": {"amplitude": 3 * u.Jy},
        "evaluation": [(0.6 * u.micron, 0.2 * u.m, 3 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": Disk2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "R_0": 300 * u.cm,
        },
        "evaluation": [(5.8 * u.m, 201 * u.cm, 3 * u.Jy)],
        "bounding_box": [[-1, 5], [0, 6]] * u.m,
    },
    {
        "class": TrapezoidDisk2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 1 * u.m,
            "y_0": 2 * u.m,
            "R_0": 100 * u.cm,
            "slope": 1 * u.Jy / u.m,
        },
        "evaluation": [(3.5 * u.m, 2 * u.m, 1.5 * u.Jy)],
        "bounding_box": [[-2, 6], [-3, 5]] * u.m,
    },
    {
        "class": Ellipse2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "a": 300 * u.cm,
            "b": 200 * u.cm,
            "theta": 45 * u.deg,
        },
        "evaluation": [(4 * u.m, 300 * u.cm, 3 * u.Jy)],
        "bounding_box": [
            [-0.5495097567963922, 4.549509756796392],
            [0.4504902432036073, 5.549509756796393],
        ]
        * u.m,
    },
    {
        "class": Ring2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "r_in": 2 * u.cm,
            "r_out": 2.1 * u.cm,
        },
        "evaluation": [(302.05 * u.cm, 2 * u.m + 10 * u.um, 3 * u.Jy)],
        "bounding_box": [[1.979, 2.021], [2.979, 3.021]] * u.m,
    },
    {
        "class": Box2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.s,
            "x_width": 4 * u.cm,
            "y_width": 3 * u.s,
        },
        "evaluation": [(301 * u.cm, 3 * u.s, 3 * u.Jy)],
        "bounding_box": [[0.5 * u.s, 3.5 * u.s], [2.98 * u.m, 3.02 * u.m]],
    },
    {
        "class": RickerWavelet2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "sigma": 1 * u.m,
        },
        "evaluation": [(4 * u.m, 2.5 * u.m, 0.602169107 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": AiryDisk2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "radius": 1 * u.m,
        },
        "evaluation": [(4 * u.m, 2.1 * u.m, 4.76998480e-05 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": Moffat2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 4.4 * u.um,
            "y_0": 3.5 * u.um,
            "gamma": 1e-3 * u.mm,
            "alpha": 1,
        },
        "evaluation": [(1000 * u.nm, 2 * u.um, 0.202565833 * u.Jy)],
        "bounding_box": False,
    },
    {
        "class": Sersic2D,
        "parameters": {
            "amplitude": 3 * u.MJy / u.sr,
            "x_0": 1 * u.arcsec,
            "y_0": 2 * u.arcsec,
            "r_eff": 2 * u.arcsec,
            "n": 4,
            "ellip": 0,
            "theta": 0,
        },
        "evaluation": [(3 * u.arcsec, 2.5 * u.arcsec, 2.829990489 * u.MJy / u.sr)],
        "bounding_box": False,
    },
    {
        "class": GeneralSersic2D,
        "parameters": {
            "amplitude": 3 * u.MJy / u.sr,
            "x_0": 1 * u.arcsec,
            "y_0": 2 * u.arcsec,
            "r_eff": 2 * u.arcsec,
            "n": 4,
            "c": 1,
            "ellip": 0,
            "theta": 0,
        },
        "evaluation": [
            (3 * u.arcsec, 2.5 * u.arcsec, 2.9704014001846475 * u.MJy / u.sr)
        ],
        "bounding_box": False,
    },
    {
        "class": Planar2D,
        "parameters": {"slope_x": 2 * u.m, "slope_y": 3 * u.m, "intercept": 4 * u.m},
        "evaluation": [(5 * u.m / u.m, 6 * u.m / u.m, 32 * u.m)],
        "bounding_box": False,
    },
]

POWERLAW_MODELS = [
    {
        "class": PowerLaw1D,
        "parameters": {"amplitude": 5 * u.kg, "x_0": 10 * u.cm, "alpha": 1},
        "evaluation": [(1 * u.m, 500 * u.g)],
        "bounding_box": False,
    },
    {
        "class": BrokenPowerLaw1D,
        "parameters": {
            "amplitude": 5 * u.kg,
            "x_break": 10 * u.cm,
            "alpha_1": 1,
            "alpha_2": -1,
        },
        "evaluation": [(1 * u.m, 50 * u.kg), (1 * u.cm, 50 * u.kg)],
        "bounding_box": False,
    },
    {
        "class": SmoothlyBrokenPowerLaw1D,
        "parameters": {
            "amplitude": 5 * u.kg,
            "x_break": 10 * u.cm,
            "alpha_1": 1,
            "alpha_2": -1,
            "delta": 1,
        },
        "evaluation": [(1 * u.cm, 15.125 * u.kg), (1 * u.m, 15.125 * u.kg)],
        "bounding_box": False,
    },
    {
        "class": ExponentialCutoffPowerLaw1D,
        "parameters": {
            "amplitude": 5 * u.kg,
            "x_0": 10 * u.cm,
            "alpha": 1,
            "x_cutoff": 1 * u.m,
        },
        "evaluation": [(1 * u.um, 499999.5 * u.kg), (10 * u.m, 50 * np.exp(-10) * u.g)],
        "bounding_box": False,
    },
    {
        "class": LogParabola1D,
        "parameters": {"amplitude": 5 * u.kg, "x_0": 10 * u.cm, "alpha": 1, "beta": 2},
        "evaluation": [(1 * u.cm, 5 * 0.1 ** (-1 - 2 * np.log(0.1)) * u.kg)],
        "bounding_box": False,
    },
    {
        "class": Schechter1D,
        "parameters": {
            "phi_star": 1.0e-4 * (u.Mpc**-3),
            "m_star": -20.0 * u.ABmag,
            "alpha": -1.9,
        },
        "evaluation": [(-23 * u.ABmag, 1.002702276867279e-12 * (u.Mpc**-3))],
        "bounding_box": False,
    },
]

POLY_MODELS = [
    {
        "class": Polynomial1D,
        "parameters": {"degree": 2, "c0": 3 * u.one, "c1": 2 / u.m, "c2": 3 / u.m**2},
        "evaluation": [(3 * u.m, 36 * u.one)],
        "bounding_box": False,
    },
    {
        "class": Polynomial1D,
        "parameters": {
            "degree": 2,
            "c0": 3 * u.kg,
            "c1": 2 * u.kg / u.m,
            "c2": 3 * u.kg / u.m**2,
        },
        "evaluation": [(3 * u.m, 36 * u.kg)],
        "bounding_box": False,
    },
    {
        "class": Polynomial1D,
        "parameters": {"degree": 2, "c0": 3 * u.kg, "c1": 2 * u.kg, "c2": 3 * u.kg},
        "evaluation": [(3 * u.one, 36 * u.kg)],
        "bounding_box": False,
    },
    {
        "class": Polynomial2D,
        "parameters": {
            "degree": 2,
            "c0_0": 3 * u.one,
            "c1_0": 2 / u.m,
            "c2_0": 3 / u.m**2,
            "c0_1": 3 / u.s,
            "c0_2": -2 / u.s**2,
            "c1_1": 5 / u.m / u.s,
        },
        "evaluation": [(3 * u.m, 2 * u.s, 64 * u.one)],
        "bounding_box": False,
    },
    {
        "class": Polynomial2D,
        "parameters": {
            "degree": 2,
            "c0_0": 3 * u.kg,
            "c1_0": 2 * u.kg / u.m,
            "c2_0": 3 * u.kg / u.m**2,
            "c0_1": 3 * u.kg / u.s,
            "c0_2": -2 * u.kg / u.s**2,
            "c1_1": 5 * u.kg / u.m / u.s,
        },
        "evaluation": [(3 * u.m, 2 * u.s, 64 * u.kg)],
        "bounding_box": False,
    },
    {
        "class": Polynomial2D,
        "parameters": {
            "degree": 2,
            "c0_0": 3 * u.kg,
            "c1_0": 2 * u.kg,
            "c2_0": 3 * u.kg,
            "c0_1": 3 * u.kg,
            "c0_2": -2 * u.kg,
            "c1_1": 5 * u.kg,
        },
        "evaluation": [(3 * u.one, 2 * u.one, 64 * u.kg)],
        "bounding_box": False,
    },
]


MODELS = (
    FUNC_MODELS_1D
    + SCALE_MODELS
    + FUNC_MODELS_2D
    + POWERLAW_MODELS
    + PHYS_MODELS_1D
    + POLY_MODELS
)

SCIPY_MODELS = {Sersic1D, Sersic2D, GeneralSersic2D, AiryDisk2D}

# These models will fail fitting test, because built in fitting data
#   will produce non-finite values
NON_FINITE_LevMar_MODELS = [
    Sersic1D,
    Sersic2D,
    GeneralSersic2D,
    ArcSine1D,
    ArcCosine1D,
    PowerLaw1D,
    ExponentialCutoffPowerLaw1D,
    BrokenPowerLaw1D,
    LogParabola1D,
    Schechter1D,
]

# These models will fail the TRFLSQFitter fitting test due to non-finite
NON_FINITE_TRF_MODELS = [
    ArcSine1D,
    ArcCosine1D,
    Sersic1D,
    Sersic2D,
    GeneralSersic2D,
    PowerLaw1D,
    ExponentialCutoffPowerLaw1D,
    BrokenPowerLaw1D,
]

# These models will fail the LMLSQFitter fitting test due to non-finite
NON_FINITE_LM_MODELS = [
    Sersic1D,
    ArcSine1D,
    ArcCosine1D,
    PowerLaw1D,
    LogParabola1D,
    Schechter1D,
    ExponentialCutoffPowerLaw1D,
    BrokenPowerLaw1D,
    GeneralSersic2D,
]

# These models will fail the DogBoxLSQFitter fitting test due to non-finite
NON_FINITE_DogBox_MODELS = [
    Sersic1D,
    Sersic2D,
    GeneralSersic2D,
    ArcSine1D,
    ArcCosine1D,
    SmoothlyBrokenPowerLaw1D,
    ExponentialCutoffPowerLaw1D,
    LogParabola1D,
    Gaussian2D,
]


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.filterwarnings(r"ignore:humlicek2 has been deprecated since .*")
def test_models_evaluate_without_units(model):
    if not HAS_SCIPY and model["class"] in SCIPY_MODELS:
        pytest.skip()
    m = model["class"](**model["parameters"])
    for args in model["evaluation"]:
        if len(args) == 2:
            kwargs = dict(zip(("x", "y"), args))
        else:
            kwargs = dict(zip(("x", "y", "z"), args))
            if kwargs["x"].unit.is_equivalent(kwargs["y"].unit):
                kwargs["x"] = kwargs["x"].to(kwargs["y"].unit)
        mnu = m.without_units_for_data(**kwargs)
        args = [x.value for x in kwargs.values()]
        assert_quantity_allclose(mnu(*args[:-1]), args[-1])


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.filterwarnings(r"ignore:humlicek2 has been deprecated since .*")
def test_models_evaluate_with_units(model):
    if not HAS_SCIPY and model["class"] in SCIPY_MODELS:
        pytest.skip()
    m = model["class"](**model["parameters"])
    for args in model["evaluation"]:
        assert_quantity_allclose(m(*args[:-1]), args[-1])


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.filterwarnings(r"ignore:humlicek2 has been deprecated since .*")
def test_models_evaluate_with_units_x_array(model):
    if not HAS_SCIPY and model["class"] in SCIPY_MODELS:
        pytest.skip()

    m = model["class"](**model["parameters"])

    for args in model["evaluation"]:
        if len(args) == 2:
            x, y = args
            x_arr = u.Quantity([x, x], subok=True)
            result = m(x_arr)
            assert_quantity_allclose(result, u.Quantity([y, y], subok=True))
        else:
            x, y, z = args
            x_arr = u.Quantity([x, x])
            y_arr = u.Quantity([y, y])
            result = m(x_arr, y_arr)
            assert_quantity_allclose(result, u.Quantity([z, z]))


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.filterwarnings(r"ignore:humlicek2 has been deprecated since .*")
def test_models_evaluate_with_units_param_array(model):
    if not HAS_SCIPY and model["class"] in SCIPY_MODELS:
        pytest.skip()

    params = {}
    for key, value in model["parameters"].items():
        if value is None or key in ("degree", "method"):
            params[key] = value
        else:
            params[key] = np.repeat(value, 2)

    params["n_models"] = 2

    m = model["class"](**params)

    for args in model["evaluation"]:
        if len(args) == 2:
            x, y = args
            x_arr = u.Quantity([x, x], subok=True)
            result = m(x_arr)
            assert_quantity_allclose(result, u.Quantity([y, y], subok=True))
        else:
            x, y, z = args
            x_arr = u.Quantity([x, x])
            y_arr = u.Quantity([y, y])
            result = m(x_arr, y_arr)
            assert_quantity_allclose(result, u.Quantity([z, z]))

    if model["class"] == Drude1D:
        params["x_0"][-1] = 0 * u.AA
        MESSAGE = r"0 is not an allowed value for x_0"
        with pytest.raises(InputParameterError, match=MESSAGE):
            model["class"](**params)


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.filterwarnings(r"ignore:humlicek2 has been deprecated since .*")
def test_models_bounding_box(model):
    # In some cases, having units in parameters caused bounding_box to break,
    # so this is to ensure that it works correctly.

    if not HAS_SCIPY and model["class"] in SCIPY_MODELS:
        pytest.skip()

    m = model["class"](**model["parameters"])

    # In the following we need to explicitly test that the value is False
    # since Quantities no longer evaluate as as True
    if model["bounding_box"] is False:
        # Check that NotImplementedError is raised, so that if bounding_box is
        # implemented we remember to set bounding_box=True in the list of models
        # above
        MESSAGE = r"No bounding box is defined for this model"
        with pytest.raises(NotImplementedError, match=MESSAGE):
            m.bounding_box
    else:
        # A bounding box may have inhomogeneous units so we need to check the
        # values one by one.
        for i in range(len(model["bounding_box"])):
            bbox = m.bounding_box
            if isinstance(bbox, ModelBoundingBox):
                bbox = bbox.bounding_box()
            assert_quantity_allclose(bbox[i], model["bounding_box"][i])


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.filterwarnings(r"ignore:humlicek2 has been deprecated since .*")
def test_compound_model_input_units_equivalencies_defaults(model):
    m = model["class"](**model["parameters"])

    assert m.input_units_equivalencies is None

    compound_model = m + m
    assert compound_model.inputs_map()["x"][0].input_units_equivalencies is None
    fixed_input_model = fix_inputs(compound_model, {"x": 1})

    assert fixed_input_model.input_units_equivalencies is None

    compound_model = m - m
    assert compound_model.inputs_map()["x"][0].input_units_equivalencies is None
    fixed_input_model = fix_inputs(compound_model, {"x": 1})

    assert fixed_input_model.input_units_equivalencies is None

    compound_model = m & m
    assert compound_model.inputs_map()["x1"][0].input_units_equivalencies is None
    fixed_input_model = fix_inputs(compound_model, {"x0": 1})
    assert fixed_input_model.inputs_map()["x1"][0].input_units_equivalencies is None

    assert fixed_input_model.input_units_equivalencies is None

    if m.n_outputs == m.n_inputs:
        compound_model = m | m
        assert compound_model.inputs_map()["x"][0].input_units_equivalencies is None
        fixed_input_model = fix_inputs(compound_model, {"x": 1})

        assert fixed_input_model.input_units_equivalencies is None


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.filterwarnings(r"ignore:.*:RuntimeWarning")
@pytest.mark.filterwarnings(r"ignore:Model is linear in parameters.*")
@pytest.mark.filterwarnings(r"ignore:The fit may be unsuccessful.*")
@pytest.mark.filterwarnings(r"ignore:humlicek2 has been deprecated since .*")
@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("fitter", fitters)
def test_models_fitting(model, fitter):
    fitter = fitter()

    bad_voigt = model["class"] == Voigt1D and ("method" not in model["parameters"])
    if (
        (
            isinstance(fitter, LevMarLSQFitter)
            and model["class"] in NON_FINITE_LevMar_MODELS
        )
        or (
            isinstance(fitter, TRFLSQFitter)
            and (model["class"] in NON_FINITE_TRF_MODELS or bad_voigt)
        )
        or (
            isinstance(fitter, LMLSQFitter)
            and (model["class"] in NON_FINITE_LM_MODELS or bad_voigt)
        )
        or (
            isinstance(fitter, DogBoxLSQFitter)
            and model["class"] in NON_FINITE_DogBox_MODELS
        )
    ):
        return

    m = model["class"](**model["parameters"])

    if m.has_bounds and isinstance(fitter, LMLSQFitter):
        pytest.skip("The LMLSQFitter fitter does not support models with bounds")

    if len(model["evaluation"][0]) == 2:
        x = np.linspace(1, 3, 100) * model["evaluation"][0][0].unit
        y = np.exp(-(x.value**2)) * model["evaluation"][0][1].unit
        args = [x, y]
    else:
        x = np.linspace(1, 3, 100) * model["evaluation"][0][0].unit
        y = np.linspace(1, 3, 100) * model["evaluation"][0][1].unit
        z = np.exp(-(x.value**2 + y.value**2)) * model["evaluation"][0][2].unit
        args = [x, y, z]

    # Test that the model fits even if it has units on parameters
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


unit_mismatch_models = [
    {
        "class": Gaussian2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_mean": 2 * u.m,
            "y_mean": 1 * u.m,
            "x_stddev": 3 * u.m,
            "y_stddev": 2 * u.m,
            "theta": 45 * u.deg,
        },
        "evaluation": [
            (412.1320343 * u.cm, 3.121320343 * u.K, 3 * u.Jy * np.exp(-0.5)),
            (412.1320343 * u.K, 3.121320343 * u.m, 3 * u.Jy * np.exp(-0.5)),
        ],
        "bounding_box": [[-14.18257445, 16.18257445], [-10.75693665, 14.75693665]]
        * u.m,
    },
    {
        "class": Ellipse2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "a": 300 * u.cm,
            "b": 200 * u.cm,
            "theta": 45 * u.deg,
        },
        "evaluation": [(4 * u.m, 300 * u.K, 3 * u.Jy), (4 * u.K, 300 * u.cm, 3 * u.Jy)],
        "bounding_box": [[-0.76046808, 4.76046808], [0.68055697, 5.31944302]] * u.m,
    },
    {
        "class": Disk2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "R_0": 300 * u.cm,
        },
        "evaluation": [
            (5.8 * u.m, 201 * u.K, 3 * u.Jy),
            (5.8 * u.K, 201 * u.cm, 3 * u.Jy),
        ],
        "bounding_box": [[-1, 5], [0, 6]] * u.m,
    },
    {
        "class": Ring2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "r_in": 2 * u.cm,
            "r_out": 2.1 * u.cm,
        },
        "evaluation": [
            (302.05 * u.cm, 2 * u.K + 10 * u.K, 3 * u.Jy),
            (302.05 * u.K, 2 * u.m + 10 * u.um, 3 * u.Jy),
        ],
        "bounding_box": [[1.979, 2.021], [2.979, 3.021]] * u.m,
    },
    {
        "class": TrapezoidDisk2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 1 * u.m,
            "y_0": 2 * u.m,
            "R_0": 100 * u.cm,
            "slope": 1 * u.Jy / u.m,
        },
        "evaluation": [
            (3.5 * u.m, 2 * u.K, 1.5 * u.Jy),
            (3.5 * u.K, 2 * u.m, 1.5 * u.Jy),
        ],
        "bounding_box": [[-2, 6], [-3, 5]] * u.m,
    },
    {
        "class": RickerWavelet2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "sigma": 1 * u.m,
        },
        "evaluation": [
            (4 * u.m, 2.5 * u.K, 0.602169107 * u.Jy),
            (4 * u.K, 2.5 * u.m, 0.602169107 * u.Jy),
        ],
        "bounding_box": False,
    },
    {
        "class": AiryDisk2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "radius": 1 * u.m,
        },
        "evaluation": [
            (4 * u.m, 2.1 * u.K, 4.76998480e-05 * u.Jy),
            (4 * u.K, 2.1 * u.m, 4.76998480e-05 * u.Jy),
        ],
        "bounding_box": False,
    },
    {
        "class": Moffat2D,
        "parameters": {
            "amplitude": 3 * u.Jy,
            "x_0": 4.4 * u.um,
            "y_0": 3.5 * u.um,
            "gamma": 1e-3 * u.mm,
            "alpha": 1,
        },
        "evaluation": [
            (1000 * u.nm, 2 * u.K, 0.202565833 * u.Jy),
            (1000 * u.K, 2 * u.um, 0.202565833 * u.Jy),
        ],
        "bounding_box": False,
    },
    {
        "class": Sersic2D,
        "parameters": {
            "amplitude": 3 * u.MJy / u.sr,
            "x_0": 1 * u.arcsec,
            "y_0": 2 * u.arcsec,
            "r_eff": 2 * u.arcsec,
            "n": 4,
            "ellip": 0,
            "theta": 0,
        },
        "evaluation": [
            (3 * u.arcsec, 2.5 * u.m, 2.829990489 * u.MJy / u.sr),
            (3 * u.m, 2.5 * u.arcsec, 2.829990489 * u.MJy / u.sr),
        ],
        "bounding_box": False,
    },
]


@pytest.mark.parametrize("model", unit_mismatch_models)
def test_input_unit_mismatch_error(model):
    if not HAS_SCIPY and model["class"] in SCIPY_MODELS:
        pytest.skip()

    MESSAGE = "Units of 'x' and 'y' inputs should match"

    m = model["class"](**model["parameters"])

    for args in model["evaluation"]:
        if len(args) == 2:
            kwargs = dict(zip(("x", "y"), args))
        else:
            kwargs = dict(zip(("x", "y", "z"), args))
            if kwargs["x"].unit.is_equivalent(kwargs["y"].unit):
                kwargs["x"] = kwargs["x"].to(kwargs["y"].unit)
        with pytest.raises(u.UnitsError, match=MESSAGE):
            m.without_units_for_data(**kwargs)


mag_models = [
    {
        "class": Const1D,
        "parameters": {"amplitude": 3 * u.ABmag},
        "evaluation": [(0.6 * u.ABmag, 3 * u.ABmag)],
    },
    {
        "class": Const1D,
        "parameters": {"amplitude": 3 * u.ABmag},
        "evaluation": [(0.6 * u.mag, 3 * u.ABmag)],
    },
    {
        "class": Const1D,
        "parameters": {"amplitude": 3 * u.mag},
        "evaluation": [(0.6 * u.ABmag, 3 * u.mag)],
    },
    {
        "class": Const1D,
        "parameters": {"amplitude": 3 * u.mag},
        "evaluation": [(0.6 * u.mag, 3 * u.mag)],
    },
    {
        "class": Const2D,
        "parameters": {"amplitude": 3 * u.ABmag},
        "evaluation": [(0.6 * u.micron, 0.2 * u.m, 3 * u.ABmag)],
    },
    {
        "class": Ellipse2D,
        "parameters": {
            "amplitude": 3 * u.ABmag,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "a": 300 * u.cm,
            "b": 200 * u.cm,
            "theta": 45 * u.deg,
        },
        "evaluation": [(4 * u.m, 300 * u.cm, 3 * u.ABmag)],
    },
    {
        "class": Disk2D,
        "parameters": {
            "amplitude": 3 * u.ABmag,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "R_0": 300 * u.cm,
        },
        "evaluation": [(5.8 * u.m, 201 * u.cm, 3 * u.ABmag)],
    },
    {
        "class": Ring2D,
        "parameters": {
            "amplitude": 3 * u.ABmag,
            "x_0": 3 * u.m,
            "y_0": 2 * u.m,
            "r_in": 2 * u.cm,
            "r_out": 2.1 * u.cm,
        },
        "evaluation": [(302.05 * u.cm, 2 * u.m + 10 * u.um, 3 * u.ABmag)],
    },
    {
        "class": Box2D,
        "parameters": {
            "amplitude": 3 * u.ABmag,
            "x_0": 3 * u.m,
            "y_0": 2 * u.s,
            "x_width": 4 * u.cm,
            "y_width": 3 * u.s,
        },
        "evaluation": [(301 * u.cm, 3 * u.s, 3 * u.ABmag)],
    },
    {
        "class": SmoothlyBrokenPowerLaw1D,
        "parameters": {
            "amplitude": 5 * u.ABmag,
            "x_break": 10 * u.cm,
            "alpha_1": 1,
            "alpha_2": -1,
            "delta": 1,
        },
        "evaluation": [(1 * u.cm, 15.125 * u.ABmag), (1 * u.m, 15.125 * u.ABmag)],
    },
    {
        "class": Box1D,
        "parameters": {"amplitude": 3 * u.ABmag, "x_0": 4.4 * u.um, "width": 1 * u.um},
        "evaluation": [(4200 * u.nm, 3 * u.ABmag), (1 * u.m, 0 * u.ABmag)],
        "bounding_box": [3.9, 4.9] * u.um,
    },
    {
        "class": Schechter1D,
        "parameters": {
            "phi_star": 1.0e-4 * (u.Mpc**-3),
            "m_star": -20.0 * u.ABmag,
            "alpha": -1.9,
        },
        "evaluation": [(-23 * u.ABmag, 1.002702276867279e-12 * (u.Mpc**-3))],
    },
    {
        "class": Schechter1D,
        "parameters": {
            "phi_star": 1.0e-4 * (u.Mpc**-3),
            "m_star": -20.0 * u.mag,
            "alpha": -1.9,
        },
        "evaluation": [(-23 * u.mag, 1.002702276867279e-12 * (u.Mpc**-3))],
    },
]


@pytest.mark.parametrize("model", mag_models)
def test_models_evaluate_magunits(model):
    if not HAS_SCIPY and model["class"] in SCIPY_MODELS:
        pytest.skip()

    m = model["class"](**model["parameters"])
    for args in model["evaluation"]:
        assert_quantity_allclose(m(*args[:-1]), args[-1])


def test_Schechter1D_errors():
    # Non magnitude units are bad
    model = Schechter1D(phi_star=1.0e-4 * (u.Mpc**-3), m_star=-20.0 * u.km, alpha=-1.9)
    MESSAGE = r"The units of magnitude and m_star must be a magnitude"
    with pytest.raises(u.UnitsError, match=MESSAGE):
        model(-23 * u.km)

    # Differing magnitude systems are bad
    model = Schechter1D(
        phi_star=1.0e-4 * (u.Mpc**-3), m_star=-20.0 * u.ABmag, alpha=-1.9
    )
    MESSAGE = (
        r".*: Units of input 'x', .*, could not be converted to required input units"
        r" of .*"
    )
    with pytest.raises(u.UnitsError, match=MESSAGE):
        model(-23 * u.STmag)

    # Differing magnitude systems are bad
    model = Schechter1D(
        phi_star=1.0e-4 * (u.Mpc**-3), m_star=-20.0 * u.ABmag, alpha=-1.9
    )
    with pytest.raises(u.UnitsError, match=MESSAGE):
        model(-23 * u.mag)


def test_compound_without_units_for_data_parameters():
    # Regression test for a bug that caused models returned by
    # CompoundModel.without_units_for_data to return a model that has top-level
    # parameters decoupled from the parameters on the individual models.

    g1 = Gaussian1D(amplitude=2 * u.Jy, stddev=4 * u.nm, mean=1000 * u.nm)
    g2 = Gaussian1D(amplitude=1 * u.Jy, stddev=2 * u.nm, mean=500 * u.nm)

    gg = g1 * g2

    gg_nounit = gg.without_units_for_data(x=1 * u.nm, y=2 * u.Jy**2)[0]

    gg_nounit.amplitude_0 = 5
    assert gg_nounit.left.amplitude == 5

    gg_nounit.amplitude_1 = 6
    assert gg_nounit.right.amplitude == 6
