"""
An example of how to define a new user-defined linear model
Note that this will probably not be needed in practice, you can
simply use the built-in polynomial model.
"""
import numpy as np
from astropy import models
from astropy.models.parameters import Parameter

class LineModel(models.ParametricModel):
    """Define a line model: y ~ slope * x + intercept"""
    parnames = ['slope', 'intercept']

    def __init__(self, slope, intercept, paramdim=1):
        self.linear = True
        self.deriv = None
        self._slope = Parameter(name='slope', val=slope, mclass=self, paramdim=paramdim)
        self._intercept = Parameter(name='intercept', val=intercept, mclass=self, paramdim=paramdim)
        models.ParametricModel.__init__(self, self.parnames, ndim=1, outdim=1, paramdim=paramdim)
        self.domain = [-1, 1]
        self.window = [-1, 1]
        self._order = 2

    def eval(self, x, params):
        return params[0] * x + params[1]

    def __call__(self, x):
        x, format = models._convert_input(x, self.paramdim)
        result = self.eval(x, self.psets)
        return models._convert_output(result, format)


if __name__ == '__main__':
    # Run an example to show that it works
    # The new user-defined model is used in the same way as built-in models are used
    from astropy.models.fitting import LinearLSQFitter
    line_model = LineModel(slope=1, intercept=2)

    # Generate some example data
    x = np.arange(1, 10, .1)
    y = line_model(x) + np.random.randn(len(x))
    line_chi2_fit = LinearLSQFitter(line_model)
    line_chi2_fit(x, y)  # This call of the fitter actually performs the fit
    print(line_model)
