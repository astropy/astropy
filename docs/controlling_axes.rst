==================
Controlling Axes
==================

Changing Axis Units
===================

WCSAxes also allows users to change the units of the axes of an image. In the
example in :doc:`slicing_datacubes`, the x axis represents velocity in m/s. We
can change the unit to an equivalent one by:


.. plot::
   :context: reset
   :nofigs:


    from wcsaxes import datasets, WCS

    hdu = datasets.fetch_l1448_co_hdu()
    wcs = WCS(hdu.header)

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs,
                      slices=(50, 'y', 'x'))

    ax.imshow(hdu.data[:, :, 50].transpose(), cmap=plt.cm.gist_heat)

.. plot::
   :context:
   :include-source:
   :align: center

    import astropy.units as u
    ax.coords[2].set_major_formatter('x.x') # Otherwise values round to the nearest whole number
    ax.coords[2].set_format_unit(u.km / u.s)


This feature is only for non-angular coordinate axes. To change the format of
angles, refer to :ref:`tick_label_format`.

Flipping axis direction
=======================
