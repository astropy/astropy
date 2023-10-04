*****************************
Overlaying coordinate systems
*****************************

For the example in the following page we start from the example introduced in
:ref:`initialization`.

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    import matplotlib.pyplot as plt

    ax = plt.subplot(projection=wcs)
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

The coordinates shown by default in a plot will be those derived from the WCS
or transformation passed to the :class:`~astropy.visualization.wcsaxes.WCSAxes` class.
However, it is possible to overlay different coordinate systems using the
:meth:`~astropy.visualization.wcsaxes.WCSAxes.get_coords_overlay` method:

.. plot::
   :context:
   :include-source:
   :align: center

    overlay = ax.get_coords_overlay('fk5')

The object returned is a :class:`~astropy.visualization.wcsaxes.coordinates_map.CoordinatesMap`, the
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


Interior ticks and tick labels
******************************

The tick labels for an overlay grid can be difficult to associate correctly with
gridlines because the default locations at the edges of the rectangular frame
may result in multiple gridlines intersecting an edge near the same tick label
or too few gridlines intersecting an edge.  As with the base grid, it is
possible to add interior ticks or tick labels for the overlay grid.  Here we add
a "tickable" gridline at constant RA (``const-ra``) and one at constant
declination (``const-dec``).  Note that when you use a multi-character string as
the name for one of these gridlines, you need to specify that name as a part of
a tuple to other methods.

.. plot::
   :context:
   :include-source:

    from astropy.coordinates import Angle

    overlay['ra'].grid(color='red')
    overlay['dec'].grid(color='magenta')

    overlay['ra'].add_tickable_gridline('const-ra', Angle('266d20m'))
    overlay['dec'].add_tickable_gridline('const-dec', Angle('-29d00m'))

    overlay['ra'].set_ticks_position(('const-dec', 't'))
    overlay['ra'].set_ticks(color='red')
    overlay['ra'].set_ticklabel_position(('const-dec',))
    overlay['ra'].set_ticklabel(color='red', size=6)
    overlay['ra'].set_axislabel_position('r')
    overlay['ra'].set_axislabel('Right Ascension', color='red')

    overlay['dec'].set_ticks_position(('const-ra', 'r'))
    overlay['dec'].set_ticks(color='magenta')
    overlay['dec'].set_ticklabel_position(('const-ra',))
    overlay['dec'].set_ticklabel(color='magenta', size=6)
    overlay['dec'].set_axislabel_position('t')
    overlay['dec'].set_axislabel('Declination', color='magenta')

## Overlaying a Compass Arrow
*****************************

It is often useful to add compass arrows to your images, denoting which directions correspond to North and East on the sky.
The function :meth:`~astropy.wcs.utils.north_pole_angle` calculates the correct angle for this compass, which can easily be displayed using a matplotlib :class:`~matplotlib.mpl_toolkits.axes_grid1.anchored_artistsAnchoredDirectionArrows()` artist.


.. plot::
   :context:
   :include-source:

   from astropy.wcs import utils
   from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDirectionArrows

   # Find the pixel coordinates for a point in the image
   # This is the point where we'll calculate the North direction from
   image_point = (300, 300)
   north_angle = utils.north_pole_angle(image_point, wcs)

   # Make the matplotlib artist
   # For East to point to the left, we require aspect_ratio=-1
   # Note that the angle calculated from utils.north_pole_angle
   # is from the positive x-axis to North, but maptlotlib assumes
   # the angle is from the positive x-axis to East. We therefore
   # subtract 90 degrees.
   arrow = AnchoredDirectionArrows(
            ax.transAxes,
            xlabel='E',
            ylabel='N',
            loc = 'upper right',
            length = -0.15,
            aspect_ratio = -1,
            sep_y = -0.1,
            sep_x = 0.04,
            color='white',
            angle=north_angle - 90,
            back_length=0
            )
   ax.add_artist(arrow)
