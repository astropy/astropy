.. _astropy-visualization-rgb:

*************************
Creating color RGB images
*************************

`Lupton et al. (2004)`_ describe an "optimal" algorithm for producing red-green-
blue composite images from three separate high-dynamic range arrays. This method
is implemented in `~astropy.visualization.make_lupton_rgb` as a convenience wraper function and an associated set of classes to provide alternate scalings. To
generate a color PNG file with the default (arcsinh) scaling:

.. doctest-skip::

>>> import numpy as np
>>> from astropy.visualization import make_lupton_rgb
>>> imageR = np.random.random((100,100))
>>> imageG = np.random.random((100,100))
>>> imageB = np.random.random((100,100))
>>> image = lupton_rgb.make_lupton_rgb(image_r, image_g, image_b, filename='randoms.png')

This method requires that the three images be aligned and have the same pixel
scale and size.

Changing ``minimum`` and ``data_range`` will change the black (``minimum``) and white
(``minimum+data_range``) levels of the resulting image, while changing ``Q`` will
change how the values between black and white are scaled.

The SDSS SkyServer color images were made with this technique.

.. _Lupton et al. (2004): http://adsabs.harvard.edu/abs/2004PASP..116..133L
