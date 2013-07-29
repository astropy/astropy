# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Here are all the test parameters and values for the ParametricModels
defined. There is a dictionary for 1D and a dictionary for 2D models.

Explanation of keywords of the dictionaries:

"parameters" : list or dict
    Model parameters, the model is tested with. Make sure you keep the right order.
    For polynomials you can also use a dict to specify the coefficients. See examples
    below.

"x_values" : list
    x values where the model is evaluated.

"y_values" : list
    Reference y values for the in x_values given positions.

"z_values" : list
    Reference z values for the in x_values and y_values given positions.
    (2D model option)

"x_lim" : list
    x test range for the model fitter. Depending on the model this can
    differ e.g. the PowerLaw model should be tested over a few magnitudes.

"y_lim" : list
    y test range for the model fitter. Depending on the model this can
    differ e.g. the PowerLaw model should be tested over a few magnitudes.
    (2D model option)

"log_fit" : bool
    PowerLaw models should be tested over a few magnitudes. So log_fit
    should be true.

"requires_scipy" : bool
    If a model requires scipy (Bessel functions etc.) set this flag.
"""

from ..functional_models import *
from ..models import *
import numpy as np

#1D Models
models_1D = {}
models_1D[Gaussian1DModel] = {'parameters': [1, 0, 1],
                           'x_values': [0, np.sqrt(2), -np.sqrt(2)],
                           'y_values': [1.0, 0.367879, 0.367879],
                           'x_lim': [-10, 10]}

models_1D[Sine1DModel] = {'parameters': [1, 1],
                           'x_values': [0, 0.25],
                           'y_values': [0, 1],
                           'x_lim': [-10, 10]}

models_1D[Box1DModel] = {'parameters': [1, 0, 1],
                           'x_values': [-0.5, 0.5, 0, -1, 1],
                           'y_values': [1, 1, 1, 0, 0],
                           'x_lim': [-2, 2]}

models_1D[Linear1DModel] = {'parameters': [1, 0],
                           'x_values': [0, np.pi, 42, -1],
                           'y_values': [0, np.pi, 42, -1],
                           'x_lim': [-10, 10]}

models_1D[Lorentz1DModel] = {'parameters': [1, 0, 1],
                           'x_values': [0, -1, 1, 0.5, -0.5],
                           'y_values': [1., 0.2, 0.2, 0.5, 0.5],
                           'x_lim': [-10, 10]}

models_1D[MexicanHat1DModel] = {'parameters': [1, 0, 1],
                           'x_values': [0, 1, -1, 3, -3],
                           'y_values': [1.0, 0.303265, 0.303265, -0.038881, -0.038881],
                           'x_lim': [-10, 10]}

models_1D[Trapezoid1DModel] = {'parameters': [1, 0, 2, 1],
                           'x_values': [0, 1, -1, 1.5, -1.5, 2, 2],
                           'y_values': [1, 1, 1, 0.5, 0.5, 0, 0],
                           'x_lim': [-10, 10]}

models_1D[Const1DModel] = {'parameters': [1],
                           'x_values': [-1, 1, np.pi, -42., 0],
                           'y_values': [1, 1, 1, 1, 1],
                           'x_lim': [-10, 10]}

models_1D[Beta1DModel] = {'parameters': [1, 0, 1, 2],
                           'x_values': [0, 1, -1, 3, -3],
                           'y_values': [1.0, 0.25, 0.25, 0.01, 0.01],
                           'x_lim': [-10, 10]}

models_1D[PowerLaw1DModel] = {'parameters': [1, 2],
                              'x_values': [1, 10, 100],
                              'y_values': [1.0, 0.01, 0.0001],
                              'x_lim': [1, 100],
                              'log_fit': True}

models_1D[Poly1DModel] = {'parameters': {'degree': 2, 'c0': 1., 'c1': 1., 'c2': 1.},
                            'x_values': [1, 10, 100],
                            'y_values': [3, 111, 10101],
                            'x_lim': [1, 100]}

#2D Models
models_2D = {}
models_2D[Gaussian2DModel] = {'parameters': [1, 0, 0, 1, 1],
                              'constraints': {'fixed': {'theta': True}},
                              'x_values': [0, np.sqrt(2), -np.sqrt(2)],
                              'y_values': [0, np.sqrt(2), -np.sqrt(2)],
                              'z_values': [1, 1. / np.exp(1) ** 2, 1. / np.exp(1) ** 2],
                              'x_lim': [-10, 10],
                              'y_lim': [-10, 10]}

models_2D[Const2DModel] = {'parameters': [1],
                           'x_values': [-1, 1, np.pi, -42., 0],
                           'y_values': [0, 1, 42, np.pi, -1],
                           'z_values': [1, 1, 1, 1, 1],
                           'x_lim': [-10, 10],
                           'y_lim': [-10, 10]}

models_2D[Box2DModel] = {'parameters': [1, 0, 0, 1, 1],
                         'x_values': [-0.5, 0.5, 0, -1, 1],
                         'y_values': [-0.5, 0.5, 0, -1, 1],
                         'z_values': [1, 1, 1, 0, 0],
                         'x_lim': [-2, 2],
                         'y_lim': [-2, 2]}

models_2D[MexicanHat2DModel] = {'parameters': [1, 0, 0, 1],
                                'x_values': [0, 0, 0, 0, 0, 1, -1, 3, -3],
                                'y_values': [0, 1, -1, 3, -3, 0, 0, 0, 0],
                                'z_values': [1.0, 0.303265, 0.303265, -0.038881, -0.038881,
                                             0.303265, 0.303265, -0.038881, -0.038881],
                                'x_lim': [-10, 10],
                                'y_lim': [-10, 10]}

models_2D[TrapezoidDisk2DModel] = {'parameters': [1, 0, 0, 1, 1],
                                   'x_values': [0, 0.5, 0, 1.5],
                                   'y_values': [0, 0.5, 1.5, 0],
                                   'z_values': [1, 1, 0.5, 0.5],
                                   'x_lim': [-3, 3],
                                   'y_lim': [-3, 3]}

models_2D[AiryDisk2DModel] = {'parameters': [1, 0, 0, 1],
                          'x_values': [0, 1, -1, -0.5, -0.5],
                          'y_values': [0, -1, 0.5, 0.5, -0.5],
                          'z_values': [1, 0.057894, 0.000788, -0.096890, -0.096890],
                          'x_lim': [-10, 10],
                          'y_lim': [-10, 10],
                          'requires_scipy': True}

models_2D[Beta2DModel] = {'parameters': [1, 0, 0, 1, 2],
                          'x_values': [0, 1, -1, 3, -3],
                          'y_values': [0, -1, 3, 1, -3],
                          'z_values': [1.0, 0.111111, 0.008264, 0.008264, 0.00277],
                          'x_lim': [-3, 3],
                          'y_lim': [-3, 3]}

models_2D[Poly2DModel] = {'parameters': {'degree': 1, 'c0_0': 1., 'c1_0': 1., 'c0_1': 1.},
                            'x_values': [1, 2, 3],
                            'y_values': [1, 3, 2],
                            'z_values': [3, 6, 6],
                            'x_lim': [1, 100],
                            'y_lim': [1, 100]}
