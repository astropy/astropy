**********************************
Creating color RGB images
**********************************

`Lupton et al. (2004)`_ describe an "optimal" algorithm for producing red-green-
blue composite images from three separate high-dynamic range arrays. This method
is implemented in `astropy.visualization.lupton_rgb` as a set of classes
providing different scalings and a convenience wraper function. To generate a
color PNG file with the default (arcsinh) scaling:

.. plot::
    :include-source:

>>> import numpy as np
>>> from astropy.visualization import lupton_rgb
>>> imageR = np.random.random((100,100))
>>> imageG = np.random.random((100,100))
>>> imageB = np.random.random((100,100))
>>> image = lupton_rgb.makeRGB(imageR, imageG, imageB, fileName='randoms.png')

This method requires that the three images be aligned and have the same pixel
scale and size.

Changing `minimum` and `dataRange` will change the black (`minimum`) and white
(`minimum+dataRange`) levels of the resulting image, while changing `Q` will
change how the values between black and white are scaled.

The SDSS SkyServer color images were made with this technique.

.. _Lupton et al. (2004): http://adsabs.harvard.edu/abs/2004PASP..116..133L
