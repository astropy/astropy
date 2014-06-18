=============================
Slicing Multidimensional Data
=============================

WCSAxes can utlimately only plot two-dimensional data. If we have an n-dimensional dataset, we have to select which dimensions to use for the x and y axis of the image. This example will show how to slice a FITS data cube and plot an image from it.

Slicing the WCS object
======================

Like the example introduced in :doc:`getting_started`, we will read in the data using `astropy.io.fits
<http://docs.astropy.org/en/stable/io/fits/index.html>`_ and parse the WCS information. 

.. plot::
   :context:
   :nofigs:
   :include-source:
   :align: center

    from astropy.wcs import WCS
    from wcsaxes import datasets
    hdu = datasets.l1448_co_hdu()
    wcs = WCS(hdu.header)
    image_data = hdu.data

If you have the original FITS file, which can be download from `here
<http://astrofrog.github.io/wcsaxes-datasets/L1448_13CO.fits>`_, this is equivalent to::

    from astropy.io import fits
    from astropy.wcs import WCS
    hdu = fits.open('L1448_13CO.fits')[0]
    wcs = WCS(hdu.header)

We then instantiate the :class:`~wcsaxes.wcsaxes.WCSAxes` using the :class:`~astropy.wcs.WCS` object and select the slices we want to plot:

.. plot::
   :context:
   :include-source:
   :align: center

    import matplotlib.pyplot as plt
    fig = plt.figure()
    from wcsaxes import WCSAxes
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs, slices=(50, 'y', 'x'))

This is a three-dimensional dataset which you can check by looking at the header information by::
    
    hdu.header

By setting ``slices=(50, 'x', 'y')``, we have chosen to plot the second dimension on the y-axis and the third dimension on the x-axis. Even though we are not plotting the all the dimensions, we have to specify which slices to select for the dimensions that are not shown. In this example, we are not plotting the first dimension so we have selected the slice 50 to display. You can experiment with this by changing the selected slice and looking at how the plotted image changes. 

Plotting the image
==================

We then add the axes to the image and plot it using the matplotlib method :meth:`~wcsaxes.wcsaxes.WCSAxes.imshow`.

.. plot::
   :context:
   :include-source:
   :align: center

    fig.add_axes(ax)
    ax.imshow(image_data[:, :, 50].transpose(), cmap=plt.cm.gist_heat)

Here, ``image_data`` is a `Numpy ndarry
<http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html>`_. In Numpy, the order of the axes is reversed so the first dimension in the FITS file appears last, the last dimension appears first and so on. Therefore the index passes to :meth:`~wcsaxes.wcsaxes.WCSAxes.imshow` should be the same as passed to ``slices`` but in reversed order. We also need to :meth:`~numpy.numpy.ndarray.transpose` ``image_data`` as we have reversed the dimensions plotted on the x and y axes in the slice.

If we don't want to reverse the dimensions plotted, we can simply do:

.. plot::
   :context:
   :include-source:
   :align: center

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs, slices=(50, 'x', 'y'))
    fig.add_axes(ax)
    ax.imshow(image_data[:, :, 50], cmap=plt.cm.gist_heat)
