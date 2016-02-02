====================
Overplotting artists
====================

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
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)

    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, cmap=plt.cm.gist_heat,
              origin='lower')


Transforms
==========

Apart from the handling of the ticks, tick labels, and grid lines, the
`wcsaxes.WCSAxes` class behaves like a normal Matplotlib
``Axes`` instance, and methods such as
:meth:`~matplotlib.axes.Axes.imshow`,
:meth:`~matplotlib.axes.Axes.contour`,
:meth:`~matplotlib.axes.Axes.plot`,
:meth:`~matplotlib.axes.Axes.scatter`, and so on will work and plot the
data in pixel coordinates. However, all such Matplotlib commands allow a
``transform=`` argument to be passed, and the
:meth:`~wcsaxes.WCSAxes.get_transform` method can be used to get the
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

* ``'fk4'``: B1950 FK4 equatorial coordinates
* ``'fk5'``: J2000 FK5 equatorial coordinates
* ``'icrs'``: ICRS equatorial coordinates
* ``'galactic'``: Galactic coordinates

It is also possible to directly pass a frame object from
:mod:`astropy.coordinates`.

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

but we can use the :meth:`~wcsaxes.WCSAxes.get_transform` method above
to plot for example in FK5 equatorial coordinates:

.. plot::
   :context:
   :include-source:
   :align: center

    r = Rectangle((266.0, -28.9), 0.3, 0.15, edgecolor='green', facecolor='none',
                  transform=ax.get_transform('fk5'))
    ax.add_patch(r)

Many Matplotlib methods accept the ``transform=`` option, so
:meth:`~wcsaxes.WCSAxes.get_transform` can be used in many cases to
plot overlays in various coordinate systems. A few examples are shown below.

Contours
========

Overplotting contours is also simple using the
:meth:`~wcsaxes.WCSAxes.get_transform` method. For contours,
:meth:`~wcsaxes.WCSAxes.get_transform` should be given the WCS of the
image to plot the contours for:

.. plot::
   :context:
   :include-source:
   :align: center

    hdu = datasets.fetch_bolocam_hdu()
    ax.contour(hdu.data, transform=ax.get_transform(WCS(hdu.header)),
               levels=[1,2,3,4,5,6], colors='white')

    ax.set_xlim(-0.5, 148.5)
    ax.set_ylim(-0.5, 148.5)

The calls to ``set_xlim`` and ``set_ylim`` are included here as the contours
cover a larger region than the image, so we want to make sure we focus just on
the image.

Scatter plots
=============

Since the ``ax.scatter`` Matplotlib routine can take the ``transform`` option,
it can also be used to plot objects in various coordinate systems:

.. plot::
   :context:
   :include-source:
   :align: center

    l = [0.25, 0.20, 0.30, 0.27]
    b = [0.20, 0.23, 0.27, 0.30]

    ax.scatter(l, b, transform=ax.get_transform('galactic'), s=100,
               edgecolor='white', facecolor='yellow', alpha=0.5)
