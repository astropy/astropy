====================
Using a custom frame
====================

By default, ``WCSAxes`` will make use of a rectangular frame for a plot, but
this can be changed to provide any custom frame:

.. plot::
   :context: reset
   :include-source:
   :align: center

    from astropy.wcs import WCS
    from wcsaxes import datasets
    from wcsaxes.frame import EllipticalFrame
    
    hdu = datasets.msx_hdu()
    wcs = WCS(hdu.header)

    import matplotlib.pyplot as plt
    fig = plt.figure()

    from wcsaxes import WCSAxes

    ax = WCSAxes(fig, [0.15, 0.15, 0.7, 0.7], wcs=wcs, 
                 frame_class=EllipticalFrame)
    fig.add_axes(ax)  # note that the axes have to be added to the figure

    ax.coords.grid(color='white')

    im = ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, cmap=plt.cm.gist_heat,
              origin='lower')

    # Clip the image to the frame
    im.set_clip_path(ax.coords.frame.patch)



