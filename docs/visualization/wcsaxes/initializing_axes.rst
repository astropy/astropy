.. _initialization:

****************************************
Initializing axes with world coordinates
****************************************

Basic initialization
********************

To make a plot using `~astropy.visualization.wcsaxes.WCSAxes`, we first read in
the data using `astropy.io.fits
<https://docs.astropy.org/en/stable/io/fits/index.html>`_ and parse the WCS
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
:class:`~astropy.wcs.WCS` object created above. The following example shows how
to do this with the Matplotlib 'pyplot' interface, keeping a reference to the
axes object:

.. plot::
   :context:
   :include-source:
   :align: center

    import matplotlib.pyplot as plt
    ax = plt.subplot(projection=wcs)

The ``ax`` object created is an instance of the
:class:`~astropy.visualization.wcsaxes.WCSAxes` class. Note that if no WCS
transformation is specified, the transformation will default to identity,
meaning that the world coordinates will match the pixel coordinates.

The field of view shown is, as for standard matplotlib axes, 0 to 1 in both
directions, in pixel coordinates. As soon as you show an image (see
:doc:`images_contours`), the limits will be adjusted, but if you want you can
also adjust the limits manually. Adjusting the limits is done using the
same functions/methods as for a normal Matplotlib plot:

.. plot::
   :context:
   :include-source:
   :align: center

    ax.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
    ax.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

.. note:: If you use the pyplot interface, you can also replace ``ax.set_xlim`` and
          ``ax.set_ylim`` by ``plt.xlim`` and ``plt.ylim``.

Alternative methods
*******************

As in Matplotlib, there are in fact several ways you can initialize the
:class:`~astropy.visualization.wcsaxes.WCSAxes`.

As shown above, the simplest way is to make use of the :class:`~astropy.wcs.WCS`
class and pass this to ``plt.subplot``. If you normally use the (partially)
object-oriented interface of Matplotlib, you can also do::

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

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
*****************************

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
