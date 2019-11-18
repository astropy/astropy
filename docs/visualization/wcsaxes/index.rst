.. _wcsaxes:

*********************************************
Making plots with world coordinates (WCSAxes)
*********************************************

WCSAxes is a framework for making plots of Astronomical data in
`Matplotlib <https://matplotlib.org/>`_. It was previously distributed
as a standalone package, but is now included in
:ref:`astropy.visualization <astropy-visualization>`.

.. _wcsaxes-getting-started:

Getting started
===============

The following is a very simple example of plotting an image with the WCSAxes
package:

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

    plt.subplot(projection=wcs)
    plt.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
    plt.grid(color='white', ls='solid')
    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')

This example uses the :mod:`matplotlib.pyplot` interface to Matplotlib, but WCSAxes
can be used with any of the other ways of using Matplotlib (some examples of which
are given in :ref:`initialization`). For example, using the partially object-oriented
interface, you can do::

    ax = plt.subplot(projection=wcs)
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
    ax.grid(color='white', ls='solid')
    ax.set_xlabel('Galactic Longitude')
    ax.set_ylabel('Galactic Latitude')

However, the axes object is needed to access some of the more advanced functionality
of WCSAxes.  An example of this usage is:

.. plot::
   :context:
   :include-source:
   :align: center

    ax = plt.subplot(projection=wcs, label='overlays')

    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

    ax.coords.grid(True, color='white', ls='solid')
    ax.coords[0].set_axislabel('Galactic Longitude')
    ax.coords[1].set_axislabel('Galactic Latitude')

    overlay = ax.get_coords_overlay('fk5')
    overlay.grid(color='white', ls='dotted')
    overlay[0].set_axislabel('Right Ascension (J2000)')
    overlay[1].set_axislabel('Declination (J2000)')

In the rest of this documentation we will assume that you have kept a reference
to the axes object, which we will refer to as ``ax``. However, we also note
when something can be done directly with the pyplot interface.

WCSAxes supports a number of advanced plotting options, including the ability to
control which axes to show labels on for which coordinates, overlaying contours
from data with different coordinate systems, overlaying grids for different
coordinate systems, dealing with plotting slices from data with more dimensions
than the plot, and defining custom (non-rectangular) frames.

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
   generic_transforms
   custom_frames

Reference/API
=============

.. automodapi:: astropy.visualization.wcsaxes
   :no-inheritance-diagram:

.. automodapi:: astropy.visualization.wcsaxes.frame
   :no-inheritance-diagram:
