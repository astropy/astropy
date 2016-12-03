**********************************
Color RGB image creation
**********************************

`Lupton et al. (2004)`_ describe an "optimal" algorithm for producing red-green-
blue composite images from three separate high-dynamic range images. This method
is implemented in `astropy.visualization.lupton_rgb`. To generate a color PNG
file  with the default (arcsinh) scaling:

>>> from astropy.visualization import lupton_rgb
>>> imageR = np.random.random((100,100))
>>> imageG = np.random.random((100,100))
>>> imageB = np.random.random((100,100))
>>> image = lupton_rgb.makeRGB(imageR, imageG, imageB, fileName='randoms.png')

.. _Lupton et al. (2004): http://adsabs.harvard.edu/abs/2004PASP..116..133L

This method requires that the three images be aligned and have the same pixel
scale and size.
