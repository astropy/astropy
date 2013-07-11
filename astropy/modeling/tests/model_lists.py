from ..functional_models import *
import numpy as np

models_1D = {}
models_1D[Gaussian1DModel] = {'parameters': [1, 0, 1],
                           'x_values': [0, np.sqrt(2), -np.sqrt(2)],
                           'y_values': [1.0, 0.135335, 0.135335],
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

models_1D[Lorentz1DModel] = {'parameters': [1, 1],
                           'x_values': [0, -1, 1, 0.5, -0.5],
                           'y_values': [1., 0.2, 0.2, 0.5, 0.5],
                           'x_lim': [-10, 10]}

models_1D[MexicanHat1DModel] = {'parameters': [1, 0, 1],
                           'x_values': [0],
                           'y_values': [1],
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
                           'x_values': [2, 1, 10.],
                           'y_values': [0.25, 1, 0.01],
                           'x_lim': [1, 100],
                           'log_fit': True}


models2D = {}

