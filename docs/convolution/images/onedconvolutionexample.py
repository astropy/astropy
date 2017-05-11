import numpy as np
import matplotlib.pyplot as plt

from astropy.convolution import Gaussian1DKernel, convolve

# Generate fake data
x = np.arange(1000).astype(float)
y = np.sin(x / 100.) + np.random.normal(0., 1., x.shape)
y[::10] = np.nan

# Create kernel
g = Gaussian1DKernel(stddev=50)

# Convolve data
z = convolve(y, g)

# Plot data before and after convolution
plt.plot(x, y, 'k.', label='Before')
plt.plot(x, z, 'b-', label='After')
plt.legend(loc='best')
plt.show()
