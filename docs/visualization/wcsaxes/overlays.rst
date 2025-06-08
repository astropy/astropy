********************************
Overplotting markers and artists
********************************

The class :class:`~astropy.visualization.wcsaxes.WCSAxes` provides two handy methods:
:meth:`~astropy.visualization.wcsaxes.WCSAxes.plot_coord`,
:meth:`~astropy.visualization.wcsaxes.WCSAxes.scatter_coord`

Used to plots and scatter respectively :class:`~astropy.coordinates.SkyCoord` or :class:`~astropy.coordinates.BaseCoordinateFrame` coordinates on the axes. The ``transform`` keyword argument will be created based on the coordinate, specifying it here will throw a :class:`~TypeError`.

For the example in the following page we start from the example introduced in
:ref:`initialization`.

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    import matplotlib.pyplot as plt

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')


Pixel coordinates
*****************

Apart from the handling of the ticks, tick labels, and grid lines, the
`~astropy.visualization.wcsaxes.WCSAxes` class behaves like a normal Matplotlib
``Axes`` instance, and methods such as
:meth:`~matplotlib.axes.Axes.imshow`,
:meth:`~matplotlib.axes.Axes.contour`,
:meth:`~matplotlib.axes.Axes.plot`,
:meth:`~matplotlib.axes.Axes.scatter`, and so on will work and plot the
data in **pixel coordinates** by default.

In the following example, the scatter markers and the rectangle will be plotted
in pixel coordinates:

.. plot::
   :context:
   :include-source:
   :align: center

    # The following line makes it so that the zoom level no longer changes,
    # otherwise Matplotlib has a tendency to zoom out when adding overlays.
    ax.set_autoscale_on(False)

    # Add a rectangle with bottom left corner at pixel position (30, 50) with a
    # width and height of 60 and 50 pixels respectively.
    from matplotlib.patches import Rectangle
    r = Rectangle((30., 50.), 60., 50., edgecolor='yellow', facecolor='none')
    ax.add_patch(r)

    # Add three markers at (40, 30), (100, 130), and (130, 60). The facecolor is
    # a transparent white (0.5 is the alpha value).
    ax.scatter([40, 100, 130], [30, 130, 60], s=100, edgecolor='white', facecolor=(1, 1, 1, 0.5))

World coordinates
*****************

All such Matplotlib commands allow a ``transform=`` argument to be passed,
which will transform the input from world to pixel coordinates before it is
passed to Matplotlib and plotted. For instance::

    ax.scatter(..., transform=...)

will take the values passed to :meth:`~matplotlib.axes.Axes.scatter` and will
transform them using the transformation passed to ``transform=``, in order to
end up with the final pixel coordinates.

The `~astropy.visualization.wcsaxes.WCSAxes` class includes a :meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform`
method that can be used to get the appropriate transformation object to convert
from various world coordinate systems to the final display coordinate system
required by Matplotlib. The :meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform` method can
take a number of different inputs, which are described in this and subsequent
sections. The two simplest inputs to this method are ``'world'`` and
``'pixel'``.

For example, if your WCS defines an image where the coordinate system consists of an angle in degrees and a wavelength in nanometers, you can do::

    ax.scatter([34], [3.2], transform=ax.get_transform('world'))

to plot a marker at (34deg, 3.2nm).

Using ``ax.get_transform('pixel')`` is equivalent to not using any
transformation at all (and things then behave as described in the `Pixel
coordinates`_ section).

Celestial coordinates
*********************

For the special case where the WCS represents celestial coordinates, a number
of other inputs can be passed to :meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform`. These
are:

* ``'fk4'``: B1950 FK4 equatorial coordinates
* ``'fk5'``: J2000 FK5 equatorial coordinates
* ``'icrs'``: ICRS equatorial coordinates
* ``'galactic'``: Galactic coordinates

In addition, any valid `astropy.coordinates` coordinate frame can be passed.

For example, you can add markers with positions defined in the FK5 system using:

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    import matplotlib.pyplot as plt

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

    ax.set_autoscale_on(False)

.. plot::
   :context:
   :include-source:
   :align: center

    ax.scatter(266.78238, -28.769255, transform=ax.get_transform('fk5'), s=300,
               edgecolor='white', facecolor='none')

In the case of :meth:`~matplotlib.axes.Axes.scatter` and :meth:`~matplotlib.axes.Axes.plot`, the positions of the center of the markers is transformed, but the markers themselves are drawn in the frame of reference of the image, which means that they will not look distorted.

Patches/shapes/lines
********************

Transformations can also be passed to Astropy or Matplotlib patches. For example, we can
use the :meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform` method above to plot a quadrangle
in FK5 equatorial coordinates:

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    import matplotlib.pyplot as plt

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

    ax.set_autoscale_on(False)

.. plot::
   :context:
   :include-source:
   :align: center

    from astropy import units as u
    from astropy.visualization.wcsaxes import Quadrangle

    r = Quadrangle((266.0, -28.9)*u.deg, 0.3*u.deg, 0.15*u.deg,
                   edgecolor='green', facecolor='none',
                   transform=ax.get_transform('fk5'))
    ax.add_patch(r)

In this case, the quadrangle will be plotted at FK5 J2000 coordinates (266deg, -28.9deg).
See the `Quadrangles`_ section for more information on `~astropy.visualization.wcsaxes.Quadrangle`.

However, it is **very important** to note that while the height will indeed be 0.15 degrees, the width will not strictly represent 0.3 degrees on the sky, but an interval of 0.3 degrees in longitude (which, depending on the latitude, will represent a different angle on the sky).
In other words, if the width and height are set to the same value, the resulting polygon will not be a square.
The same applies to the `~matplotlib.patches.Circle` patch, which will not actually produce a circle:

.. plot::
   :context:
   :include-source:
   :align: center

    from matplotlib.patches import Circle

    r = Quadrangle((266.4, -28.9)*u.deg, 0.3*u.deg, 0.3*u.deg,
                   edgecolor='cyan', facecolor='none',
                   transform=ax.get_transform('fk5'))
    ax.add_patch(r)

    c = Circle((266.4, -29.1), 0.15, edgecolor='yellow', facecolor='none',
               transform=ax.get_transform('fk5'))
    ax.add_patch(c)



.. important:: If what you are interested is simply plotting circles around
               sources to highlight them, then we recommend using
               :meth:`~matplotlib.axes.Axes.scatter`, since for the circular
               marker (the default), the circles will be guaranteed to be
               circles in the plot, and only the position of the center is
               transformed.

               To plot 'true' spherical circles, see the `Spherical patches`_
               section.

Quadrangles
***********

`~astropy.visualization.wcsaxes.Quadrangle` is the recommended patch for plotting a quadrangle, as opposed to Matplotlib's `~matplotlib.patches.Rectangle`.
The edges of a quadrangle lie on two lines of constant longitude and two lines of constant latitude (or the equivalent component names in the coordinate frame of interest, such as right ascension and declination).
The edges of `~astropy.visualization.wcsaxes.Quadrangle` will render as curved lines if appropriate for the WCS transformation.
In contrast, `~matplotlib.patches.Rectangle` will always have straight edges.
Here's a comparison of the two types of patches for plotting a quadrangle in `~astropy.coordinates.ICRS` coordinates on `~astropy.coordinates.Galactic` axes:

.. plot::
   :context: reset
   :nofigs:

    from astropy import units as u
    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    from astropy.visualization.wcsaxes import Quadrangle
    import matplotlib.pyplot as plt

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

.. plot::
   :context:
   :include-source:
   :align: center

    from matplotlib.patches import Rectangle

    # Set the Galactic axes such that the plot includes the ICRS south pole
    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.set_xlim(0, 10000)
    ax.set_ylim(-10000, 0)

    # Overlay the ICRS coordinate grid
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='black', ls='dotted')

    # Add a quadrangle patch (100 degrees by 20 degrees)
    q = Quadrangle((255, -70)*u.deg, 100*u.deg, 20*u.deg,
                   label='Quadrangle', edgecolor='blue', facecolor='none',
                   transform=ax.get_transform('icrs'))
    ax.add_patch(q)

    # Add a rectangle patch (100 degrees by 20 degrees)
    r = Rectangle((255, -70), 100, 20,
                  label='Rectangle', edgecolor='red', facecolor='none', linestyle='--',
                  transform=ax.get_transform('icrs'))
    ax.add_patch(r)

    ax.legend(loc='upper right')

Contours
********

Overplotting contours is also simple using the
:meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform` method. For contours,
:meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform` should be given the WCS of the
image to plot the contours for:

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    from matplotlib.patches import Rectangle
    import matplotlib.pyplot as plt

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

    ax.set_autoscale_on(False)

.. plot::
   :context:
   :include-source:
   :align: center

    filename = get_pkg_data_filename('galactic_center/gc_bolocam_gps.fits')
    hdu = fits.open(filename)[0]
    ax.contour(hdu.data, transform=ax.get_transform(WCS(hdu.header)),
               levels=[1,2,3,4,5,6], colors='white')

Spherical patches
*****************

In the case where you are making a plot of a celestial image, and want to plot a circle that represents the area within a certain angle of a longitude/latitude,
the `~matplotlib.patches.Circle` patch is not appropriate, since it will result in a distorted shape (because longitude is not the same as the angle on the sky).
For this use case, you can instead use `~astropy.visualization.wcsaxes.SphericalCircle`, which takes a tuple of |Quantity| or a |SkyCoord| object as the input,
and a |Quantity| as the radius:

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    import matplotlib.pyplot as plt

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

    ax.set_autoscale_on(False)

.. plot::
   :context:
   :include-source:
   :align: center

    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.visualization.wcsaxes import SphericalCircle


    r = SphericalCircle((266.4 * u.deg, -29.1 * u.deg), 0.15 * u.degree,
                         edgecolor='yellow', facecolor='none',
                         transform=ax.get_transform('fk5'))

    ax.add_patch(r)

    #The following lines show the usage of a SkyCoord object as the input.
    skycoord_object = SkyCoord(266.4 * u.deg, -28.7 * u.deg)
    s = SphericalCircle(skycoord_object, 0.15 * u.degree,
                        edgecolor='white', facecolor='none',
                        transform=ax.get_transform('fk5'))

    ax.add_patch(s)

Beam shape and scale bar
************************

Adding an ellipse that represents the shape of the beam on a celestial
image can be done with the
:func:`~astropy.visualization.wcsaxes.add_beam` function:

.. plot::
   :context: reset
   :nofigs:

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    import matplotlib.pyplot as plt

    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

    ax.set_autoscale_on(False)

.. plot::
   :context:
   :include-source:
   :align: center

    from astropy import units as u
    from astropy.visualization.wcsaxes import add_beam, add_scalebar

    add_beam(ax, major=1.2 * u.arcmin, minor=1.2 * u.arcmin, angle=0, frame=True)

To add a segment that shows a physical scale, you can use the
:func:`~astropy.visualization.wcsaxes.add_scalebar` function:

.. plot::
   :context:
   :include-source:
   :align: center

    # Compute the angle corresponding to 10 pc at the distance of the galactic center
    gc_distance = 8.2 * u.kpc
    scalebar_length = 10 * u.pc
    scalebar_angle = (scalebar_length / gc_distance).to(
        u.deg, equivalencies=u.dimensionless_angles()
    )

    # Add a scale bar
    add_scalebar(ax, scalebar_angle, label="10 pc", color="white")
