=================
Initializing axes
=================

As in Matplotlib, there are several ways you can initialize the
:class:`~wcsaxes.WCSAxes`.

As shown in the rest of the documentation, the
simplest way is to make use of the :class:`wcsaxes.WCS` class (instead of
:class:`astropy.wcs.WCS`) and pass this to the
:meth:`~matplotlib.figure.Figure.add_subplot` method::

    from astropy.wcs import WCS
    import matplotlib.pyplot as plt
    
    wcs = WCS(...)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    ax.imshow(...)

If you normally make plots directly with pyplot directly instead of creating
axes and figure instances, you can do::


    plt.subplot(1, 1, 1, projection=wcs)
    plt.imshow(...)

Note that this also works with :meth:`~matplotlib.figure.Figure.add_axes` and :func:`~matplotlib.pyplot.axes`, e.g.::

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)
 
or::

    plt.axes([0.1, 0.1, 0.8, 0.8], projection=wcs)

Note that in the above example, there is no explicit call to the ``wcsaxes``
package. This is because the :class:`astropy.wcs.WCS` class knows about
``wcsaxes``, and is able to automatically instantiate a
:class:`wcsaxes.WCSAxes` object.

Note that any additional arguments passed to
:meth:`~matplotlib.figure.Figure.add_subplot`,
:meth:`~matplotlib.figure.Figure.add_axes`,
:func:`~matplotlib.pyplot.subplot`, or :func:`~matplotlib.pyplot.axes`, such
as ``slices`` or ``frame_class``, will be passed on to the
:class:`~wcsaxes.WCSAxes` class.

.. _initialize_alternative:

Alternative
===========

As an alternative to the above methods of initializing
:class:`~wcsaxes.WCSAxes`, you can also instantiate :class:`~wcsaxes.WCSAxes`
directly and add it to the figure::

    from astropy.wcs import WCS
    from wcsaxes import WCSAxes
    import matplotlib.pyplot as plt
    
    wcs = WCS(...)

    fig = plt.figure()

    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
    fig.add_axes(ax)  # note that the axes have to be explicitly added to the figure

Note that in this example, we can use :class:`astropy.wcs.WCS` (but
:class:`wcsaxes.WCS` will also work).
