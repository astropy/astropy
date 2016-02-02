=============================
Overlaying coordinate systems
=============================

For the example in the following page we start from the example introduced in
:doc:`getting_started`.

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from wcsaxes import datasets

    hdu = datasets.fetch_msx_hdu()
    wcs = WCS(hdu.header)

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_axes([0.25, 0.25, 0.6, 0.6], projection=wcs)

    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, cmap=plt.cm.gist_heat,
              origin='lower')

The coordinates shown by default in a plot will be those derived from the WCS
or transformation passed to the `wcsaxes.WCSAxes` class.
However, it is possible to overlay different coordinate systems using the
:meth:`wcsaxes.WCSAxes.get_coords_overlay` method:

.. plot::
   :context:
   :include-source:
   :align: center

    overlay = ax.get_coords_overlay('fk5')

The object returned is a :class:`~wcsaxes.coordinates_map.CoordinatesMap`, the
same type of object as ``ax.coord``. It can therefore be used in the same way
as ``ax.coord`` to set the ticks, tick labels, and axis labels properties:

.. plot::
   :context:
   :include-source:
   :align: center

    ax.coords['glon'].set_ticks(color='white')
    ax.coords['glat'].set_ticks(color='white')

    ax.coords['glon'].set_axislabel('Galactic Longitude')
    ax.coords['glat'].set_axislabel('Galactic Latitude')

    ax.coords.grid(color='yellow', linestyle='solid', alpha=0.5)

    overlay['ra'].set_ticks(color='white')
    overlay['dec'].set_ticks(color='white')

    overlay['ra'].set_axislabel('Right Ascension')
    overlay['dec'].set_axislabel('Declination')

    overlay.grid(color='white', linestyle='solid', alpha=0.5)
