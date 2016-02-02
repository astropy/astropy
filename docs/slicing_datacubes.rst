=============================
Slicing Multidimensional Data
=============================

WCSAxes can ultimately only plot two-dimensional data. If we have an
n-dimensional dataset, we have to select which dimensions to use for
the x and y axis of the image. This example will show how to slice a FITS
data cube and plot an image from it.

Slicing the WCS object
======================

Like the example introduced in :doc:`getting_started`, we will read in the
data using `astropy.io.fits
<http://docs.astropy.org/en/stable/io/fits/index.html>`_ and parse the WCS
information. The original FITS file can be downloaded from `here
<http://astrofrog.github.io/wcsaxes-datasets/L1448_13CO.fits>`_.

.. plot::
   :context: reset
   :include-source:
   :align: center
   :nofigs:

    from astropy.wcs import WCS
    from wcsaxes import datasets
    hdu = datasets.fetch_l1448_co_hdu()
    wcs = WCS(hdu.header)
    image_data = hdu.data

This is a three-dimensional dataset which you can check by looking at the
header information by::

    >>> hdu.header  # doctest: +SKIP
    ...
    NAXIS = 3 /number of axes
    CTYPE1  = 'RA---SFL'           /
    CTYPE2  = 'DEC--SFL'           /
    CTYPE3  = 'VELO-LSR'           /
    ...

The header keyword 'NAXIS' gives the number of dimensions of the dataset. The keywords 'CTYPE1', 'CTYPE2' and 'CTYPE3' give the data type of these dimensions to be right ascension, declination and velocity respectively.

We then instantiate the `~wcsaxes.WCSAxes` using the
:class:`~astropy.wcs.WCS` object and select the slices we want to plot:

.. plot::
   :context:
   :include-source:
   :align: center
   :nofigs:

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs, slices=(50, 'y', 'x'))

By setting ``slices=(50, 'y', 'x')``, we have chosen to plot the second
dimension on the y-axis and the third dimension on the x-axis. Even though we
are not plotting the all the dimensions, we have to specify which slices to
select for the dimensions that are not shown. In this example, we are not
plotting the first dimension so we have selected the slice 50 to display. You
can experiment with this by changing the selected slice and looking at how the
plotted image changes.

Plotting the image
==================

We then add the axes to the image and plot it using the method
:meth:`~matplotlib.axes.Axes.imshow`.

.. plot::
   :context:
   :include-source:
   :align: center

    ax.coords[2].set_ticks(exclude_overlapping=True)
    ax.imshow(image_data[:, :, 50].transpose(), cmap=plt.cm.gist_heat)

Here, ``image_data`` is an :class:`~numpy.ndarray` object. In Numpy, the order
of the axes is reversed so the first dimension in the FITS file appears last,
the last dimension appears first and so on. Therefore the index passed to
:meth:`~matplotlib.axes.Axes.imshow` should be the same as passed to
``slices`` but in reversed order. We also need to
:meth:`~numpy.ndarray.transpose` ``image_data`` as we have reversed the
dimensions plotted on the x and y axes in the slice.

If we don't want to reverse the dimensions plotted, we can simply do:

.. plot::
   :context: reset
   :align: center
   :nofigs:

    from astropy.wcs import WCS
    from wcsaxes import datasets
    hdu = datasets.fetch_l1448_co_hdu()
    wcs = WCS(hdu.header)
    image_data = hdu.data

.. plot::
   :context:
   :include-source:
   :align: center

    fig = plt.figure(figsize=(6,3))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs, slices=(50, 'x', 'y'))
    ax.imshow(image_data[:, :, 50], cmap=plt.cm.gist_heat)
