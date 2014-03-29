=================
Plotting overlays
=================

For the example in the following page we start from the example introduced in
:doc:`getting_started`.

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    hdu = fits.open('msx.fits')[0]
    wcs = WCS(hdu.header)

    import matplotlib.pyplot as plt
    fig = plt.figure()

    from wcsaxes import WCSAxes

    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
    fig.add_axes(ax)  # note that the axes have to be added to the figure

    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, cmap=plt.cm.gist_heat,
              origin='lower')


Transforms
==========

Apart from the handling of the ticks, tick labels, and grid lines, the
:class:`wcsaxes.wcsaxes.WCSAxes` class behaves like a normal Matplotlib ``Axes``
instance, and methods such as :meth:`~wcsaxes.wcsaxes.WCSAxes.imshow`, :meth:`~wcsaxes.wcsaxes.WCSAxes.contour`, :meth:`~wcsaxes.wcsaxes.WCSAxes.plot`,
:meth:`~wcsaxes.wcsaxes.WCSAxes.scatter`, and so on will work and plot the data in
pixel coordinates. However, all such Matplotlib commands allow a
``transform=`` argument to be passed, and the
:meth:`~wcsaxes.wcsaxes.WCSAxes.get_transform` method can be used to get the
appropriate transformation object.

The following example shows how to get the transformation object for the FK5
coordinate system::

    tr_fk5 = ax.get_transform("fk5")

To plot in the FK5 system one would then do::

    ax.method(..., transform=tr_fk5)

where ``method`` is whatever Matplotlib command you are running.

To plot in the default world coordinates system, you can use::

    ax.get_transform("world")

By specifying a WCS object, you can also define a transformation from the
current WCS system to another file's WCS system, allowing e.g. overplotting of
contours in a different system::

    ax.get_transform(<WCS instance>)

If the world coordinate system of the plot is a celestial coordinate system,
the following built-in sky coordinate systems would be available from the
``get_transform`` method:

* ``'fk4'`` or ``'b1950'``: B1950 equatorial coordinates
* ``'fk5'`` or ``'j2000'``: J2000 equatorial coordinates

and more coordinate systems will be added in future.

Patches/shapes/lines
====================

As mentioned above, matplotlib methods will by default work in pixel
coordinates:

.. plot::
   :context:
   :include-source:
   :align: center
    
    from matplotlib.patches import Rectangle
    r = Rectangle((60., 20.), 10., 12., edgecolor='yellow', facecolor='none')
    ax.add_patch(r)

but we can use the :meth:`~wcsaxes.wcsaxes.WCSAxes.get_transform` method above to plot for example in FK5 equatorial coordinates:

.. plot::
   :context:
   :include-source:
   :align: center
   
    r = Rectangle((266.0, -28.9), 0.3, 0.15, edgecolor='green', facecolor='none',
                  transform=ax.get_transform('fk5'))
    ax.add_patch(r)

Many Matplotlib methods accept the ``transform=`` option, so
:meth:`~wcsaxes.wcsaxes.WCSAxes.get_transform` can be used in many cases to
plot overlays in various coordinate systems.

..     ax.add_collection(c, transform=ax.get_transform('gal'))
..     ax.add_line(l, transform=ax.get_transform('fk4'))
..     ax.scatter(l, b, transform=ax.get_transform('gal'))