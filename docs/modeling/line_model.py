import numpy as np
from astropy.modeling import models, Parameter, format_input

class LineModel(models.PolynomialModel):
    slope = Parameter('slope')
    intercept = Parameter('intercept')
    linear = True

    def __init__(self, slope, intercept, param_dim=1, **constraints):
        super(LineModel, self).__init__(slope=slope, intercept=intercept,
                                        param_dim=param_dim, **constraints)
        self.domain = [-1, 1]
        self.window = [-1, 1]
        self._order = 2

    @staticmethod
    def eval(x, slope, intercept):
        return slope * x + intercept

    @staticmethod
    def fit_deriv(x, slope, intercept):
        d_slope = x
        d_intercept = np.ones_like(x)
        return [d_slope, d_intercept]

    @format_input
    def __call__(self, x):
        return self.eval(x, *self.param_sets)
