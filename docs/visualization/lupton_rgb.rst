.. _astropy-visualization-rgb:

*************************
Creating color RGB images
*************************

`Lupton et al. (2004)`_ describe an "optimal" algorithm for producing red-green-
blue composite images from three separate high-dynamic range arrays. This method
is implemented in `~astropy.visualization.make_lupton_rgb` as a convenience
wraper function and an associated set of classes to provide alternate scalings.
The SDSS SkyServer color images were made using a variation on this technique.
To generate a color PNG file with the default (arcsinh) scaling:

.. _Lupton et al. (2004): http://adsabs.harvard.edu/abs/2004PASP..116..133L

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.visualization import make_lupton_rgb
    image_r = np.random.random((100,100))
    image_g = np.random.random((100,100))
    image_b = np.random.random((100,100))
    image = make_lupton_rgb(image_r, image_g, image_b, stretch=0.5)
    plt.imshow(image)

This method requires that the three images be aligned and have the same pixel
scale and size. Changing ``minimum`` will change the black level, while
``stretch`` and ``Q`` will change how the values between black and white are
scaled.

For a more in-depth example, download the `g`_, `r`_, `i`_ SDSS frames of the
area around the Hickson 88 group and try the example below (requires the
`reproject package`_, installable via the astropy conda channel, or via pip),
and compare it with Figure 1 of Lupton et al. (2004):

.. _reproject package: https://reproject.readthedocs.io/

.. _g: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-g-001737-5-0039.fits.bz2
.. _r: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-r-001737-5-0039.fits.bz2
.. _i: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-i-001737-5-0039.fits.bz2

.. doctest-skip::

   import numpy as np
   from astropy.visualization import make_lupton_rgb
   from astropy.io import fits
   from reproject import reproject_interp

   # Read in the three images downloaded from here:
   # g: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-g-001737-5-0039.fits.bz2
   # r: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-r-001737-5-0039.fits.bz2
   # i: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-i-001737-5-0039.fits.bz2
   g = fits.open('frame-g-001737-5-0039.fits.bz2')[0]
   r = fits.open('frame-r-001737-5-0039.fits.bz2')[0]
   i = fits.open('frame-i-001737-5-0039.fits.bz2')[0]

   # remap r and i onto g
   r_new, r_mask = reproject_interp(r, g.header)
   i_new, i_mask = reproject_interp(i, g.header)

   # zero out the unmapped values
   i_new[np.logical_not(i_mask)] = 0
   r_new[np.logical_not(r_mask)] = 0

   # red=i, green=r, blue=g
   # make a file with the default scaling
   rgb_default = make_lupton_rgb(i_new, r_new, g.data, filename="ngc6976-default.jpeg")
   # this scaling is very similar to the one used in Lupton et al. (2004)
   rgb = make_lupton_rgb(i_new, r_new, g.data, Q=10, stretch=0.5, filename="ngc6976.jpeg")

This will produce the following two images. The first is the image generated
with the default parameters.

.. raw:: html

    <a class="reference internal image-reference" href="http://data.astropy.org/visualization/ngc6976-default.jpeg"><img alt="default rgb image" src="http://data.astropy.org/visualization/ngc6976-default-small.jpeg" /></a>

The second is the image generated with Q=10, stretch=0.5, showing faint features
of the galaxies. Compare with Fig. 1 of `Lupton et al. (2004)`_ or the
`SDSS Skyserver image`_.

.. raw:: html

    <a class="reference internal image-reference" href="http://data.astropy.org/visualization/ngc6976.jpeg"><img alt="wider stretch image" src="http://data.astropy.org/visualization/ngc6976-small.jpeg" /></a>

.. _SDSS Skyserver image: http://skyserver.sdss.org/dr13/en/tools/chart/navi.aspx?ra=179.68929&dec=-0.45438&opt=
