import numpy as np
import matplotlib.pyplot as plt

from astropy.nddata.convolution import convolve, Tophat2DKernel
from astropy.nddata.convolution import Box2DKernel, Gaussian2DKernel
from astropy.nddata.convolution import MexicanHat2DKernel, AiryDisk2DKernel

from astropy.modeling.models import Gaussian2DModel


# Small Gaussian source in the middle of the image
gauss = Gaussian2DModel(1, 0, 0, 3, 3)

# Fake data including noise
x = np.arange(-100, 101)
y = np.arange(-100, 101)
x, y = np.meshgrid(x, y)
data = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)

# Setup kernels, including unity kernel for original image
kernels = [[[1]],
           Tophat2DKernel(11),
           Gaussian2DKernel(11),
           Box2DKernel(11),
           MexicanHat2DKernel(11),
           AiryDisk2DKernel(11)]

# Plot kernels
axisNum = 0
for row in range(3):
    for col in range(2):
        axisNum += 1
        ax = plt.subplot(2, 3, axisNum)
        smoothed = convolve(data, kernels[axisNum - 1])
        plt.imshow(smoothed)
        title = kernels[axisNum - 1].__class__.__name__
        if axisNum == 1:
            title = 'Original'
        plt.title(title)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
plt.show()
