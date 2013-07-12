# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Here are all the test parameters and values for the ParametricModels
defined. There is a dictionary for 1D and a dictionary for 2D models.
"""
from ..functional_models import *
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

models_1D[Trapezoid1DModel] = {'parameters': [1, 0, 1, 1],
                           'x_values': [0],
                           'y_values': [1],
                           'x_lim': [-10, 10]}

models_1D[Const1DModel] = {'parameters': [1],
                           'x_values': [-1, 1, np.pi, -42., 0],
                           'y_values': [1, 1, 1, 1, 1],
                           'x_lim': [-10, 10]}

models_1D[PowerLaw1DModel] = {'parameters': [1, 2],
                           'x_values': [1, 10, 100],
                           'y_values': [1.0, 0.01, 0.0001],
                           'x_lim': [1, 100],
                           'log_fit': True}

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
                                'x_values': [0],
                                'y_values': [0],
                                'z_values': [1],
                                'x_lim': [-10, 10],
                                'y_lim': [-10, 10]}

models_2D[TrapezoidDisk2DModel] = {'parameters': [1, 0, 0, 1, 1],
                                   'x_values': [0],
                                   'y_values': [0],
                                   'z_values': [1],
                                   'x_lim': [-3, 3],
                                   'y_lim': [-3, 3]}

models_2D[Airy2DModel] = {'parameters': [1, 0, 0, 1],
                          'x_values': [0],
                          'y_values': [0],
                          'z_values': [1],
                          'x_lim': [-10, 10],
                          'y_lim': [-10, 10],
                          'requires_scipy': True}

 