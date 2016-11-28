.. _initialization:

========================================
Initializing axes with world coordinates
========================================

Basic initialization
====================

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

Alternative methods
===================

As in Matplotlib, there are in fact several ways you can initialize the
:class:`~astropy.visualization.wcsaxes.WCSAxes`.

As shown above, the simplest way is to make use of the :class:`~astropy.wcs.WCS`
class and pass this to the :meth:`~matplotlib.figure.Figure.add_subplot`
method::

    from astropy.wcs import WCS
    import matplotlib.pyplot as plt

    wcs = WCS(...)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    ax.imshow(...)

However, if you normally make plots directly with pyplot directly instead of
creating axes and figure instances, you can also do::


    plt.subplot(1, 1, 1, projection=wcs)
    plt.imshow(...)

Note that this also works with :meth:`~matplotlib.figure.Figure.add_axes` and
:func:`~matplotlib.pyplot.axes`, e.g.::

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)

or::

    plt.axes([0.1, 0.1, 0.8, 0.8], projection=wcs)

Any additional arguments passed to
:meth:`~matplotlib.figure.Figure.add_subplot`,
:meth:`~matplotlib.figure.Figure.add_axes`,
:func:`~matplotlib.pyplot.subplot`, or :func:`~matplotlib.pyplot.axes`, such
as ``slices`` or ``frame_class``, will be passed on to the
:class:`~astropy.visualization.wcsaxes.WCSAxes` class.

.. _initialize_alternative:

Directly initializing WCSAxes
=============================

As an alternative to the above methods of initializing
:class:`~astropy.visualization.wcsaxes.WCSAxes`, you can also instantiate
:class:`~astropy.visualization.wcsaxes.WCSAxes` directly and add it to the
figure::

    from astropy.wcs import WCS
    from astropy.visualization.wcsaxes import WCSAxes
    import matplotlib.pyplot as plt

    wcs = WCS(...)

    fig = plt.figure()

    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
    fig.add_axes(ax)  # note that the axes have to be explicitly added to the figure
