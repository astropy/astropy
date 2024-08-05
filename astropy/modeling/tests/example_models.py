# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Here are all the test parameters and values for the each
`~astropy.modeling.FittableModel` defined. There is a dictionary for 1D and a
dictionary for 2D models.

Explanation of keywords of the dictionaries:

"parameters" : list or dict
    Model parameters, the model is tested with. Make sure you keep the right
    order.  For polynomials you can also use a dict to specify the
    coefficients. See examples below.

"x_values" : list
    x values where the model is evaluated.

"y_values" : list
    Reference y values for the in x_values given positions.

"z_values" : list
    Reference z values for the in x_values and y_values given positions.
    (2D model option)

"x_lim" : list
    x test range for the model fitter. Depending on the model this can differ
    e.g. the PowerLaw model should be tested over a few magnitudes.

"y_lim" : list
    y test range for the model fitter. Depending on the model this can differ
    e.g. the PowerLaw model should be tested over a few magnitudes.  (2D model
    option)

"log_fit" : bool
    PowerLaw models should be tested over a few magnitudes. So log_fit should
    be true.

"requires_scipy" : bool
    If a model requires scipy (Bessel functions etc.) set this flag.

"integral" : float
    Approximate value of the integral in the range x_lim (and y_lim).

"deriv_parameters" : list
    If given the test of the derivative will use these parameters to create a
    model (optional)

"deriv_initial" : list
    If given the test of the derivative will use these parameters as initial
    values for the fit (optional)
"""

import numpy as np

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
    Exponential1D,
    Gaussian1D,
    Gaussian2D,
    GeneralSersic2D,
    KingProjectedAnalytic1D,
    Linear1D,
    Logarithmic1D,
    Lorentz1D,
    Moffat1D,
    Moffat2D,
    Planar2D,
    RickerWavelet1D,
    RickerWavelet2D,
    Ring2D,
    Sersic1D,
    Sersic2D,
    Sine1D,
    Tangent1D,
    Trapezoid1D,
    TrapezoidDisk2D,
    Voigt1D,
)
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

# 1D Models
models_1D = {
    Gaussian1D: {
        "parameters": [1, 0, 1],
        "x_values": [0, np.sqrt(2), -np.sqrt(2)],
        "y_values": [1.0, 0.367879, 0.367879],
        "x_lim": [-10, 10],
        "integral": np.sqrt(2 * np.pi),
        "bbox_peak": True,
    },
    Sine1D: {
        "parameters": [1, 0.1, 0],
        "x_values": [0, 2.5],
        "y_values": [0, 1],
        "x_lim": [-10, 10],
        "integral": 0,
    },
    Cosine1D: {
        "parameters": [1, 0.1, 0],
        "x_values": [0, 2.5],
        "y_values": [1, 0],
        "x_lim": [-10, 10],
        "integral": 0,
    },
    Tangent1D: {
        "parameters": [1, 0.1, 0],
        "x_values": [0, 1.25],
        "y_values": [0, 1],
        "x_lim": [-10, 10],
        "integral": 0,
    },
    ArcSine1D: {
        "parameters": [1, 0.1, 0],
        "x_values": [0, 1],
        "y_values": [0, 2.5],
        "x_lim": [-0.5, 0.5],
        "integral": 0,
    },
    ArcCosine1D: {
        "parameters": [1, 0.1, 0],
        "x_values": [1, 0],
        "y_values": [0, 2.5],
        "x_lim": [-0.5, 0.5],
        "integral": 0,
    },
    ArcTangent1D: {
        "parameters": [1, 0.1, 0],
        "x_values": [0, 1],
        "y_values": [0, 1.25],
        "x_lim": [-10, 10],
        "integral": 0,
    },
    Box1D: {
        "parameters": [1, 0, 10],
        "x_values": [-5, 5, 0, -10, 10],
        "y_values": [1, 1, 1, 0, 0],
        "x_lim": [-10, 10],
        "integral": 10,
        "bbox_peak": True,
    },
    Linear1D: {
        "parameters": [1, 0],
        "x_values": [0, np.pi, 42, -1],
        "y_values": [0, np.pi, 42, -1],
        "x_lim": [-10, 10],
        "integral": 0,
    },
    Lorentz1D: {
        "parameters": [1, 0, 1],
        "x_values": [0, -1, 1, 0.5, -0.5],
        "y_values": [1.0, 0.2, 0.2, 0.5, 0.5],
        "x_lim": [-10, 10],
        "integral": 1,
        "bbox_peak": True,
        "deriv_initial": [10, 0.5, 4],
    },
    RickerWavelet1D: {
        "parameters": [1, 0, 1],
        "x_values": [0, 1, -1, 3, -3],
        "y_values": [1.0, 0.0, 0.0, -0.088872, -0.088872],
        "x_lim": [-20, 20],
        "integral": 0,
        "bbox_peak": True,
    },
    Trapezoid1D: {
        "parameters": [1, 0, 2, 1],
        "x_values": [0, 1, -1, 1.5, -1.5, 2, 2],
        "y_values": [1, 1, 1, 0.5, 0.5, 0, 0],
        "x_lim": [-10, 10],
        "integral": 3,
        "bbox_peak": True,
    },
    Const1D: {
        "parameters": [1],
        "x_values": [-1, 1, np.pi, -42.0, 0],
        "y_values": [1, 1, 1, 1, 1],
        "x_lim": [-10, 10],
        "integral": 20,
    },
    Moffat1D: {
        "parameters": [1, 0, 1, 2],
        "x_values": [0, 1, -1, 3, -3],
        "y_values": [1.0, 0.25, 0.25, 0.01, 0.01],
        "x_lim": [-10, 10],
        "integral": 1,
        "deriv_parameters": [23.4, 1.2, 2.1, 2.3],
        "deriv_initial": [10, 1, 1, 1],
    },
    PowerLaw1D: {
        "parameters": [1, 1, 2],
        "constraints": {"fixed": {"x_0": True}},
        "x_values": [1, 10, 100],
        "y_values": [1.0, 0.01, 0.0001],
        "x_lim": [1, 10],
        "log_fit": True,
        "integral": 0.99,
    },
    BrokenPowerLaw1D: {
        "parameters": [1, 1, 2, 3],
        "constraints": {"fixed": {"x_break": True}},
        "x_values": [0.1, 1, 10, 100],
        "y_values": [1e2, 1.0, 1e-3, 1e-6],
        "x_lim": [0.1, 100],
        "log_fit": True,
    },
    SmoothlyBrokenPowerLaw1D: {
        "parameters": [1, 1, -2, 2, 0.5],
        "constraints": {"fixed": {"x_break": True, "delta": True}},
        "x_values": [0.01, 1, 100],
        "y_values": [3.99920012e-04, 1.0, 3.99920012e-04],
        "x_lim": [0.01, 100],
        "log_fit": True,
    },
    ExponentialCutoffPowerLaw1D: {
        "parameters": [1, 1, 2, 3],
        "constraints": {"fixed": {"x_0": True}},
        "x_values": [0.1, 1, 10, 100],
        "y_values": [9.67216100e01, 7.16531311e-01, 3.56739933e-04, 3.33823780e-19],
        "x_lim": [0.01, 100],
        "log_fit": True,
    },
    Schechter1D: {
        "parameters": [1.0, -20.0, -1.0],
        "x_values": [-25.0, -23.0, -21.0, -19.0, -17.0],
        "y_values": [
            3.42631659e-44,
            1.20551329e-07,
            7.47097466e-02,
            6.18557294e-01,
            8.64716111e-01,
        ],
        "x_lim": [-25, -17.0],
    },
    LogParabola1D: {
        "parameters": [1, 2, 3, 0.1],
        "constraints": {"fixed": {"x_0": True}},
        "x_values": [0.1, 1, 10, 100],
        "y_values": [3.26089063e03, 7.62472488e00, 6.17440488e-03, 1.73160572e-06],
        "x_lim": [0.1, 100],
        "log_fit": True,
    },
    Polynomial1D: {
        "parameters": {"degree": 2, "c0": 1.0, "c1": 1.0, "c2": 1.0},
        "x_values": [1, 10, 100],
        "y_values": [3, 111, 10101],
        "x_lim": [-3, 3],
    },
    Sersic1D: {
        "parameters": [1, 20, 4],
        "x_values": [0.1, 1, 10, 100],
        "y_values": [2.78629391e02, 5.69791430e01, 3.38788244e00, 2.23941982e-02],
        "requires_scipy": True,
        "x_lim": [0, 10],
        "log_fit": True,
    },
    Voigt1D: {
        "parameters": [0, 1, 0.5, 0.9],
        "x_values": [0, 0.2, 0.5, 1, 2, 4, 8, 20],
        "y_values": [
            0.52092360,
            0.479697445,
            0.317550374,
            0.0988079347,
            1.73876624e-2,
            4.00173216e-3,
            9.82351731e-4,
            1.56396993e-4,
        ],
        "x_lim": [-3, 3],
    },
    KingProjectedAnalytic1D: {
        "parameters": [1, 1, 2],
        "x_values": [0, 0.1, 0.5, 0.8],
        "y_values": [0.30557281, 0.30011069, 0.2, 0.1113258],
        "x_lim": [0, 10],
        "y_lim": [0, 10],
        "bbox_peak": True,
    },
    Drude1D: {
        "parameters": [1.0, 8.0, 1.0],
        "x_values": [7.0, 8.0, 9.0, 10.0],
        "y_values": [0.17883212, 1.0, 0.21891892, 0.07163324],
        "x_lim": [1.0, 20.0],
        "y_lim": [0.0, 10.0],
        "bbox_peak": True,
    },
    Plummer1D: {
        "parameters": [10.0, 0.5],
        "x_values": [1.0000e-03, 2.5005e00, 5.0000e00],
        "y_values": [1.90984022e01, 5.53541843e-03, 1.86293603e-04],
        "x_lim": [0.001, 100],
    },
    Exponential1D: {
        "parameters": [1, 1],
        "x_values": [0, 0.5, 1],
        "y_values": [1, np.sqrt(np.e), np.e],
        "x_lim": [0, 2],
        "integral": (np.e**2 - 1.0),
    },
    Logarithmic1D: {
        "parameters": [1, 1],
        "x_values": [1, np.e, np.e**2],
        "y_values": [0, 1, 2],
        "x_lim": [1, np.e**2],
        "integral": (np.e**2 + 1),
    },
}


# 2D Models
models_2D = {
    Gaussian2D: {
        "parameters": [1, 0, 0, 1, 1],
        "constraints": {"fixed": {"theta": True}},
        "x_values": [0, np.sqrt(2), -np.sqrt(2)],
        "y_values": [0, np.sqrt(2), -np.sqrt(2)],
        "z_values": [1, 1.0 / np.exp(1) ** 2, 1.0 / np.exp(1) ** 2],
        "x_lim": [-10, 10],
        "y_lim": [-10, 10],
        "integral": 2 * np.pi,
        "deriv_parameters": [137.0, 5.1, 5.4, 1.5, 2.0, np.pi / 4],
        "deriv_initial": [10, 5, 5, 4, 4, 0.5],
        "bbox_peak": True,
    },
    Const2D: {
        "parameters": [1],
        "x_values": [-1, 1, np.pi, -42.0, 0],
        "y_values": [0, 1, 42, np.pi, -1],
        "z_values": [1, 1, 1, 1, 1],
        "x_lim": [-10, 10],
        "y_lim": [-10, 10],
        "integral": 400,
    },
    Box2D: {
        "parameters": [1, 0, 0, 10, 10],
        "x_values": [-5, 5, -5, 5, 0, -10, 10],
        "y_values": [-5, 5, 0, 0, 0, -10, 10],
        "z_values": [1, 1, 1, 1, 1, 0, 0],
        "x_lim": [-10, 10],
        "y_lim": [-10, 10],
        "integral": 100,
        "bbox_peak": True,
    },
    RickerWavelet2D: {
        "parameters": [1, 0, 0, 1],
        "x_values": [0, 0, 0, 0, 0, 1, -1, 3, -3],
        "y_values": [0, 1, -1, 3, -3, 0, 0, 0, 0],
        "z_values": [
            1.0,
            0.303265,
            0.303265,
            -0.038881,
            -0.038881,
            0.303265,
            0.303265,
            -0.038881,
            -0.038881,
        ],
        "x_lim": [-10, 11],
        "y_lim": [-10, 11],
        "integral": 0,
    },
    TrapezoidDisk2D: {
        "parameters": [1, 0, 0, 1, 1],
        "x_values": [0, 0.5, 0, 1.5],
        "y_values": [0, 0.5, 1.5, 0],
        "z_values": [1, 1, 0.5, 0.5],
        "x_lim": [-3, 3],
        "y_lim": [-3, 3],
        "bbox_peak": True,
    },
    AiryDisk2D: {
        "parameters": [7, 0, 0, 10],
        "x_values": [0, 1, -1, -0.5, -0.5],
        "y_values": [0, -1, 0.5, 0.5, -0.5],
        "z_values": [7.0, 6.50158267, 6.68490643, 6.87251093, 6.87251093],
        "x_lim": [-10, 10],
        "y_lim": [-10, 10],
        "requires_scipy": True,
    },
    Moffat2D: {
        "parameters": [1, 0, 0, 1, 2],
        "x_values": [0, 1, -1, 3, -3],
        "y_values": [0, -1, 3, 1, -3],
        "z_values": [1.0, 0.111111, 0.008264, 0.008264, 0.00277],
        "x_lim": [-3, 3],
        "y_lim": [-3, 3],
    },
    Polynomial2D: {
        "parameters": {"degree": 1, "c0_0": 1.0, "c1_0": 1.0, "c0_1": 1.0},
        "x_values": [1, 2, 3],
        "y_values": [1, 3, 2],
        "z_values": [3, 6, 6],
        "x_lim": [1, 100],
        "y_lim": [1, 100],
    },
    Disk2D: {
        "parameters": [1, 0, 0, 5],
        "x_values": [-5, 5, -5, 5, 0, -10, 10],
        "y_values": [-5, 5, 0, 0, 0, -10, 10],
        "z_values": [0, 0, 1, 1, 1, 0, 0],
        "x_lim": [-10, 10],
        "y_lim": [-10, 10],
        "integral": np.pi * 5**2,
        "bbox_peak": True,
    },
    Ring2D: {
        "parameters": [1, 0, 0, 5, 5],
        "x_values": [-5, 5, -5, 5, 0, -10, 10],
        "y_values": [-5, 5, 0, 0, 0, -10, 10],
        "z_values": [1, 1, 1, 1, 0, 0, 0],
        "x_lim": [-10, 10],
        "y_lim": [-10, 10],
        "integral": np.pi * (10**2 - 5**2),
        "bbox_peak": True,
    },
    Sersic2D: {
        "parameters": [1, 25, 4, 50, 50, 0.5, -1],
        "x_values": [0.0, 1, 10, 100],
        "y_values": [1, 100, 0.0, 10],
        "z_values": [1.686398e-02, 9.095221e-02, 2.341879e-02, 9.419231e-02],
        "requires_scipy": True,
        "x_lim": [1, 1e10],
        "y_lim": [1, 1e10],
    },
    GeneralSersic2D: {
        "parameters": [1, 25, 4, 50, 50, 0.5, -1, 1.0],
        "x_values": [0.0, 1, 10, 100],
        "y_values": [1, 100, 0.0, 10],
        "z_values": [0.01711221, 0.10471922, 0.02426327, 0.1201565],
        "requires_scipy": True,
        "x_lim": [1, 1e10],
        "y_lim": [1, 1e10],
    },
    Planar2D: {
        "parameters": [1, 1, 0],
        "x_values": [0, np.pi, 42, -1],
        "y_values": [np.pi, 0, -1, 42],
        "z_values": [np.pi, np.pi, 41, 41],
        "x_lim": [-10, 10],
        "y_lim": [-10, 10],
        "integral": 0,
    },
}
