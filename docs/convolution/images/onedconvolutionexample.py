import numpy as np
import matplotlib.pyplot as plt

from astropy.convolution import Gaussian1DKernel, convolve

plt.figure(3).clf()

# Generate fake data
x = np.arange(1000).astype(float)
y = np.sin(x / 100.) + np.random.normal(0., 1., x.shape)
y[::3] = np.nan

# Create kernel
g = Gaussian1DKernel(stddev=50)

# Convolve data
z = convolve(y, g)

# Plot data before and after convolution
plt.plot(x, y, 'k-', label='Before')
plt.plot(x, z, 'b-', label='After', alpha=0.5, linewidth=2)
plt.legend(loc='best')
plt.show()
