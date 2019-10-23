******************
Controlling Axes
******************

Changing Axis Units
*******************

WCSAxes also allows users to change the units of the axes of an image. In the
example in :doc:`slicing_datacubes`, the x axis represents velocity in m/s. We
can change the unit to an equivalent one by:


.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename

    filename = get_pkg_data_filename('l1448/l1448_13co.fits')
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    import matplotlib.pyplot as plt

    ax = plt.subplot(projection=wcs, slices=(50, 'y', 'x'))
    ax.imshow(hdu.data[:, :, 50].transpose())

.. plot::
   :context:
   :include-source:
   :align: center

    import astropy.units as u
    ax.coords[2].set_major_formatter('x.x') # Otherwise values round to the nearest whole number
    ax.coords[2].set_format_unit(u.km / u.s)


Disabling Automatic Labelling
*****************************

By default WCSAxes adds labels to the axes to indicate what world coordinate is
being represented on that axis, and what unit is being used to display it. If
you want to disable this behavior you can either set an explicit label for that
axis with `~astropy.visualization.wcsaxes.CoordinateHelper.set_axislabel` or you
can disable the feature per coordinate with::

  ax = plt.subplot(projection=wcs)  # doctest: +SKIP
  ax.coords[0].set_auto_axislabel(False)  # doctest: +SKIP


Changing Axis Directions
************************

Sometimes astronomy FITS files don't follow the convention of having the longitude increase to the left,
so we want to flip an axis so that it goes in the opposite direction. To do this on our example image:

.. plot::
   :context:
   :include-source:
   :align: center

    ax.invert_xaxis()
