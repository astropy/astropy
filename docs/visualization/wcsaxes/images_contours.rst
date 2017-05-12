****************************
Plotting images and contours
****************************

For the example in the following page we start from the example introduced in
:ref:`initialization`.

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    import matplotlib.pyplot as plt

    ax = plt.subplot(projection=wcs)
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

Plotting images as bitmaps or contours should be done via the usual matplotlib
methods such as :meth:`~matplotlib.axes.Axes.imshow` or
:meth:`~matplotlib.axes.Axes.contour`. For example, continuing from the
example in :ref:`initialization`, you can do:

.. plot::
   :context:
   :include-source:
   :align: center

    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

and we can also add contours corresponding to the same image using:

.. plot::
   :context:
   :include-source:
   :align: center

    import numpy as np
    ax.contour(hdu.data, levels=np.logspace(-4.7, -3., 10), colors='white', alpha=0.5)

To show contours for an image in a different coordinate system, see
:doc:`overlays`.

.. note:: If you like using the pyplot interface, you can also call
          ``plt.imshow`` and ``plt.contour`` instead of ``ax.imshow`` and
          ``ax.contour``.
