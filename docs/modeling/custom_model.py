__doctest_requires__ = {('.'): ['matplotlib', 'scipy']}
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.models import custom_model_1d
from astropy.modeling.fitting import NonLinearLSQFitter

# Define model
@custom_model_1d
def sum_of_gaussians(x, amplitude1=1., mean1=-1., sigma1=1.,
                     amplitude2=1., mean2=1., sigma2=1.):
    return (amplitude1 * np.exp(-0.5 * ((x - mean1) / sigma1) ** 2) +
            amplitude2 * np.exp(-0.5 * ((x - mean2) / sigma2) ** 2))

# Generate fake data
np.random.seed(0)
x = np.linspace(-5., 5., 200)
m_ref = sum_of_gaussians(amplitude1=2., mean1=-0.5, sigma1=0.4,
                         amplitude2=0.5, mean2=2., sigma2=1.0)
y = m_ref(x) + np.random.normal(0., 0.1, x.shape)

# Fit model to data
m_init = sum_of_gaussians()
fit = NonLinearLSQFitter()
m = fit(m_init, x, y)

# Plot the data and the best fit
plt.plot(x, y, 'o', color='k')
plt.plot(x, m(x), color='r', lw=2)
