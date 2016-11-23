===============
Getting started
===============

Initialization
==============

To make a plot using `~astropy.visualization.wcsaxes.WCSAxes`, we first read in
the data using `astropy.io.fits
<http://docs.astropy.org/en/stable/io/fits/index.html>`_ and parse the WCS
information. In this example, we will use an example FITS file from the
http://data.astropy.org server (the
:func:`~astropy.utils.data.get_pkg_data_filename` function downloads the file
and returns a filename):

.. plot::
   :context: reset
   :nofigs:
   :include-source:
   :align: center

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

We then create a figure using Matplotlib and create the axes using the
:class:`~astropy.wcs.WCS` object created above:

.. plot::
   :context:
   :include-source:
   :align: center

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)

The ``ax`` object created is an instance of the
:class:`~astropy.visualization.wcsaxes.WCSAxes`
class. For more information about the different ways of initializing axes,
see :doc:`initializing_axes`.

The field of view shown is, as for standard matplotlib axes, 0 to
1 in both directions, in pixel coordinates. The
:meth:`~matplotlib.axes.Axes.set_xlim` and
:meth:`~matplotlib.axes.Axes.set_ylim` methods can be used to re-set the
pixel coordinates. For example, we can set the limits to the edge of the FITS
image in pixel coordinates:

.. plot::
   :context:
   :include-source:
   :align: center

    ax.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
    ax.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

If no WCS transformation is specified, the transformation will default to
identity, meaning that the world coordinates will match the pixel coordinates.

Plotting images and contours
============================

Plotting images as bitmaps or contours should be done via the usual matplotlib
methods such as :meth:`~matplotlib.axes.Axes.imshow` or
:meth:`~matplotlib.axes.Axes.contour`. For example, continuing from the
example in `Initialization`_, you can do:

.. plot::
   :context:
   :include-source:
   :align: center

    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, cmap=plt.cm.gist_heat,
              origin='lower')

and we can also add contours corresponding to the same image using:

.. plot::
   :context:
   :include-source:
   :align: center

    import numpy as np
    ax.contour(hdu.data, levels=np.logspace(-4.7, -3., 10), colors='white', alpha=0.5)

To show contours for an image in a different coordinate system, see
:doc:`overlays`.
