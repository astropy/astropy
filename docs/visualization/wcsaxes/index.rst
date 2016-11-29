.. _wcsaxes:

*********************************************
Making plots with world coordinates (WCSAxes)
*********************************************

WCSAxes is a framework for making plots of Astronomical data in Matplotlib. It
was previously distributed as a standalone package, but is now included in
:ref:`astropy.visualization <astropy-visualization>`.

Compared to packages such as `APLpy <http://aplpy.readthedocs.io>`__, WCSAxes is
a framework that is closer to the Matplotlib API and allows more advanced plots
to be made, but on the other hand requires more code for simple plots. APLpy
will in future be updated to make use of WCSAxes and provide a simpler
interface for making common plots.

Getting started
===============

The following is a simple example of plotting an image with the WCSAxes
package, and overlaying the coordinate grid from the image as well as an
equatorial coordinate grid.

.. plot::
   :context: reset
   :include-source:
   :align: center

    import matplotlib.pyplot as plt
    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    fig = plt.figure()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)

    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

    ax.coords.grid(color='white', ls='solid')
    ax.coords[0].set_axislabel('Galactic Longitude')
    ax.coords[1].set_axislabel('Galactic Latitude')

    overlay = ax.get_coords_overlay('fk5')
    overlay.grid(color='white', ls='dotted')
    overlay[0].set_axislabel('Right Ascension (J2000)')
    overlay[1].set_axislabel('Declination (J2000)')

Using WCSAxes
=============

.. toctree::
   :maxdepth: 1

   initializing_axes
   images_contours
   ticks_labels_grid
   overlays
   overlaying_coordinate_systems
   slicing_datacubes
   controlling_axes
   custom_frames

Reference/API
=============

.. automodapi:: astropy.visualization.wcsaxes
   :no-inheritance-diagram:

.. automodapi:: astropy.visualization.wcsaxes.frame
   :no-inheritance-diagram:
