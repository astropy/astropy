.. _astropy-coordinates-example-gallery:

Example Gallery
***************

This gallery of examples shows a variety of relatively small snippets or
examples of tasks that can be done with the ``astropy.coordinates`` sub-package.


.. _sphx_glr_generated_examples_coordinates_plot_galactocentric-frame.py:

Transforming positions and velocities to and from a Galactocentric frame
========================================================================

..
  EXAMPLE START
  Transforming positions and velocities to and from a Galactocentric frames

This example shows a few examples of how to use and customize the
`~astropy.coordinates.Galactocentric` frame to transform Heliocentric sky
positions, distance, proper motions, and radial velocities to a Galactocentric,
Cartesian frame, and the same in reverse.

The main configurable parameters of the `~astropy.coordinates.Galactocentric`
frame control the position and velocity of the solar system barycenter within
the Galaxy. These are specified by setting the ICRS coordinates of the
Galactic center, the distance to the Galactic center (the sun-galactic center
line is always assumed to be the x-axis of the Galactocentric frame), and the
Cartesian 3-velocity of the sun in the Galactocentric frame. We will first
demonstrate how to customize these values, then show how to set the solar motion
instead by inputting the proper motion of Sgr A*.

Note that, for brevity, we may refer to the solar system barycenter as just "the
sun" in the examples below.

Let's first define a barycentric coordinate and velocity in the ICRS frame.
We will use the data for the star HD 39881 from the
`Simbad <https://simbad.unistra.fr/simbad/>`_ database:


>>> import astropy.coordinates as coord
>>> from astropy import units as u
>>> c1 = coord.SkyCoord(
...     ra=89.014303 * u.degree,
...     dec=13.924912 * u.degree,
...     distance=(37.59 * u.mas).to(u.pc, u.parallax()),
...     pm_ra_cosdec=372.72 * (u.mas / u.yr),
...     pm_dec=-483.69 * (u.mas / u.yr),
...     radial_velocity=0.37 * (u.km / u.s),
...     frame="icrs",
... )

This is a high proper-motion star; suppose we'd like to transform its position
and velocity to a Galactocentric frame to see if it has a large 3D velocity
as well. To use the Astropy default solar position and motion parameters, we
can do the following:

>>> gc1 = c1.transform_to(coord.Galactocentric)

From here, we can access the components of the resulting
`~astropy.coordinates.Galactocentric` instance to see the 3D Cartesian
velocity components:

>>> print(gc1.v_x, gc1.v_y, gc1.v_z)  # doctest: +FLOAT_CMP
30.254684717897074 km / s 171.29916086104885 km / s 18.19390627095307 km / s

The default parameters for the `~astropy.coordinates.Galactocentric` frame
are detailed in the linked documentation, but we can modify the most commonly
changed values using the keywords ``galcen_distance``, ``galcen_v_sun``, and
``z_sun`` which set the sun-Galactic center distance, the 3D velocity vector
of the sun, and the height of the sun above the Galactic midplane,
respectively. The velocity of the sun can be specified as an
`~astropy.units.Quantity` object with velocity units and is interpreted as a
Cartesian velocity, as in the example below. Note that, as with the positions,
the Galactocentric frame is a right-handed system (i.e., the Sun is at negative
x values) so ``v_x`` is opposite of the Galactocentric radial velocity:

>>> v_sun = [11.1, 244, 7.25] * (u.km / u.s)  # [vx, vy, vz]
>>> gc_frame = coord.Galactocentric(
...     galcen_distance=8 * u.kpc, galcen_v_sun=v_sun, z_sun=0 * u.pc
... )

We can then transform to this frame instead, with our custom parameters:

>>> gc2 = c1.transform_to(gc_frame)
>>> print(gc2.v_x, gc2.v_y, gc2.v_z)  # doctest: +FLOAT_CMP
28.427958360720748 km / s 169.69916086104888 km / s 17.70831652451455 km / s

It is sometimes useful to specify the solar motion using the
`proper motion of Sgr A* <https://arxiv.org/abs/astro-ph/0408107>`_
instead of Cartesian velocity components. With an assumed distance, we can convert
proper motion components to Cartesian velocity components using `astropy.units`:

>>> galcen_distance = 8 * u.kpc
>>> pm_gal_sgrA = [-6.379, -0.202] * (u.mas / u.yr)  # from Reid & Brunthaler 2004
>>> vy, vz = -(galcen_distance * pm_gal_sgrA).to(u.km / u.s, u.dimensionless_angles())

We still have to assume a line-of-sight velocity for the Galactic center,
which we will again take to be 11 km/s:

>>> vx = 11.1 * (u.km / u.s)
>>> v_sun2 = u.Quantity([vx, vy, vz])  # List of Quantity -> a single Quantity
>>> gc_frame2 = coord.Galactocentric(
...     galcen_distance=galcen_distance, galcen_v_sun=v_sun2, z_sun=0 * u.pc
... )
>>> gc3 = c1.transform_to(gc_frame2)
>>> print(gc3.v_x, gc3.v_y, gc3.v_z)  # doctest: +FLOAT_CMP
28.427958360720748 km / s 167.61484955608267 km / s 18.118916793584443 km / s

The transformations also work in the opposite direction. This can be useful
for transforming simulated or theoretical data to observable quantities. As
an example, we will generate 4 theoretical circular orbits at different
Galactocentric radii with the same circular velocity, and transform them to
Heliocentric coordinates:

.. plot::
   :include-source:

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import astropy.coordinates as coord
    >>> from astropy import units as u
    >>> ring_distances = np.arange(10, 26, 5) * u.kpc
    >>> circ_velocity = 220 * (u.km / u.s)
    >>> phi_grid = np.linspace(90, 270, 512) * u.degree  # grid of azimuths
    >>> ring_rep = coord.CylindricalRepresentation(
    ...     rho=ring_distances[:, np.newaxis],
    ...     phi=phi_grid[np.newaxis],
    ...     z=np.zeros_like(ring_distances)[:, np.newaxis],
    ... )
    >>> angular_velocity = (-circ_velocity / ring_distances).to(
    ...     u.mas / u.yr, u.dimensionless_angles()
    ... )
    >>> ring_dif = coord.CylindricalDifferential(
    ...     d_rho=np.zeros(phi_grid.shape)[np.newaxis] * (u.km / u.s),
    ...     d_phi=angular_velocity[:, np.newaxis],
    ...     d_z=np.zeros(phi_grid.shape)[np.newaxis] * (u.km / u.s),
    ... )
    >>> ring_rep = ring_rep.with_differentials(ring_dif)
    >>> gc_rings = coord.SkyCoord(ring_rep, frame=coord.Galactocentric)

    First, let's visualize the geometry in Galactocentric coordinates. Here are
    the positions and velocities of the rings; note that in the velocity plot,
    the velocities of the 4 rings are identical and thus overlaid under the same
    curve:

    >>> fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    >>> _ = axes[0].plot(gc_rings.x.T, gc_rings.y.T, marker="None", linewidth=3)
    >>> _ = axes[0].text(-8.0, 0, r"$\odot$", fontsize=20)
    >>> _ = axes[0].set_xlim(-30, 30)
    >>> _ = axes[0].set_ylim(-30, 30)
    >>> _ = axes[0].set_xlabel("$x$ [kpc]")
    >>> _ = axes[0].set_ylabel("$y$ [kpc]")
    >>> _ = axes[0].set_title("Positions")
    >>> _ = axes[1].plot(gc_rings.v_x.T, gc_rings.v_y.T, marker="None", linewidth=3)
    >>> _ = axes[1].set_xlim(-250, 250)
    >>> _ = axes[1].set_ylim(-250, 250)
    >>> _ = axes[1].set_xlabel(f"$v_x$ [{(u.km / u.s).to_string('latex_inline')}]")
    >>> _ = axes[1].set_ylabel(f"$v_y$ [{(u.km / u.s).to_string('latex_inline')}]")
    >>> _ = axes[1].set_title("Velocities")
    >>> fig.tight_layout()

    Now we can transform to Galactic coordinates and visualize the rings in
    observable coordinates:

    >>> gal_rings = gc_rings.transform_to(coord.Galactic)
    >>> fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    >>> for i in range(len(ring_distances)):
    ...     _ = ax.plot(
    ...         gal_rings[i].l.degree,
    ...         gal_rings[i].pm_l_cosb.value,
    ...         label=str(ring_distances[i]),
    ...         marker="None",
    ...         linewidth=3,
    ...     )
    >>> _ = ax.set_xlim(360, 0)
    >>> _ = ax.set_xlabel("$l$ [deg]")
    >>> _ = ax.set_ylabel(rf'$\mu_l \, \cos b$ [{(u.mas/u.yr).to_string("latex_inline")}]')
    >>> _ = ax.legend()
    >>> plt.draw()

..
  EXAMPLE END


.. _sphx_glr_generated_examples_coordinates_plot_mars-coordinate-frame.py:

Create a new coordinate frame class for Mars
============================================

..
  EXAMPLE START
  Create a new coordinate frame class for Mars

This example describes how to subclass and define a custom coordinate frame for a
planetary body which can be described by a geodetic or bodycentric representation,
as discussed in :ref:`astropy:astropy-coordinates-design` and
:ref:`astropy-coordinates-create-geodetic`.

Note that we use the frame here only to store coordinates. To use it to determine, e.g.,
where to point a telescope on Earth to observe Olympus Mons, one would need to add the
frame to the transfer graph, which is beyond the scope of this example.

To do this, first we need to define a subclass of a
`~astropy.coordinates.BaseGeodeticRepresentation` and
`~astropy.coordinates.BaseBodycentricRepresentation`, then a subclass of
`~astropy.coordinates.BaseCoordinateFrame` using the previous defined
representations.

.. plot::
   :include-source:

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.coordinates.baseframe import BaseCoordinateFrame
    >>> from astropy.coordinates.representation import CartesianRepresentation
    >>> from astropy.coordinates.representation.geodetic import (
    ...     BaseBodycentricRepresentation,
    ...     BaseGeodeticRepresentation,
    ... )
    >>> from astropy.visualization import quantity_support

    The first step is to create a new class, and make it a subclass of
    `~astropy.coordinates.BaseGeodeticRepresentation`.
    Geodetic latitudes are used and longitudes span from 0 to 360 degrees east positive
    It represent a best fit of the Mars spheroid to the martian geoid (areoid):

    >>> class MarsBestFitAeroid(BaseGeodeticRepresentation):
    ...     """A Spheroidal representation of Mars that minimized deviations with respect to the
    ...     areoid following
    ...         Ardalan A. A, R. Karimi, and E. W. Grafarend (2010)
    ...         https://doi.org/10.1007/s11038-009-9342-7
    ...     """
    ...     _equatorial_radius = 3395.4280 * u.km
    ...     _flattening = 0.5227617843759314 * u.percent

    Now let's define a new geodetic representation obtained from MarsBestFitAeroid but
    described by planetocentric latitudes:

    >>> class MarsBestFitOcentricAeroid(BaseBodycentricRepresentation):
    ...     """A Spheroidal planetocentric representation of Mars that minimized deviations with
    ...     respect to the areoid following
    ...         Ardalan A. A, R. Karimi, and E. W. Grafarend (2010)
    ...         https://doi.org/10.1007/s11038-009-9342-7
    ...     """
    ...     _equatorial_radius = 3395.4280 * u.km
    ...     _flattening = 0.5227617843759314 * u.percent

    As a comparison we define a new spherical frame representation, we could
    have based it on `~astropy.coordinates.BaseBodycentricRepresentation` too:

    >>> class MarsSphere(BaseGeodeticRepresentation):
    ...     """A Spherical representation of Mars."""
    ...     _equatorial_radius = 3395.4280 * u.km
    ...     _flattening = 0.0 * u.percent

    The new planetary body-fixed reference system will be described using the
    previous defined representations:

    >>> class MarsCoordinateFrame(BaseCoordinateFrame):
    ...     """A reference system for Mars."""
    ...     name = "Mars"

    Now we plot the differences between each component of the cartesian
    representation with respect to the spherical model, assuming the point on the
    surface of the body (``height = 0``):

    >>> mars_sphere = MarsCoordinateFrame(
    ...     lon=np.linspace(0, 360, 128) * u.deg,
    ...     lat=np.linspace(-90, 90, 128) * u.deg,
    ...     representation_type=MarsSphere,
    ... )
    >>> mars = MarsCoordinateFrame(
    ...     lon=np.linspace(0, 360, 128) * u.deg,
    ...     lat=np.linspace(-90, 90, 128) * u.deg,
    ...     representation_type=MarsBestFitAeroid,
    ... )
    >>> mars_ocentric = MarsCoordinateFrame(
    ...     lon=np.linspace(0, 360, 128) * u.deg,
    ...     lat=np.linspace(-90, 90, 128) * u.deg,
    ...     representation_type=MarsBestFitOcentricAeroid,
    ... )
    >>> xyz_sphere = mars_sphere.represent_as(CartesianRepresentation)
    >>> xyz = mars.represent_as(CartesianRepresentation)
    >>> xyz_ocentric = mars_ocentric.represent_as(CartesianRepresentation)
    >>> with quantity_support():
    ...     fig, ax = plt.subplots(2, subplot_kw={"projection": "3d"})
    ...     _ = ax[0].scatter(*((xyz - xyz_sphere).xyz << u.km))
    ...     _ = ax[0].tick_params(labelsize=8)
    ...     _ = ax[0].set(xlabel="x [km]", ylabel="y [km]", zlabel="z [km]")
    ...     _ = ax[0].set_title("Mars-odetic spheroid difference from sphere")
    ...     _ = ax[1].scatter(*((xyz_ocentric - xyz_sphere).xyz << u.km))
    ...     _ = ax[1].tick_params(labelsize=8)
    ...     _ = ax[1].set(xlabel="x [km]", ylabel="y [km]", zlabel="z [km]")
    ...     _ = ax[1].set_title("Mars-ocentric spheroid difference from sphere")
    ...     plt.draw()

..
  EXAMPLE END


.. _sphx_glr_generated_examples_coordinates_plot_obs-planning.py:

Determining and plotting the altitude/azimuth of a celestial object
===================================================================

..
  EXAMPLE START
  Determining and plotting the altitude/azimuth of a celestial object

This example demonstrates coordinate transformations and the creation of
visibility curves to assist with observing run planning.

In this example, we make a `~astropy.coordinates.SkyCoord` instance for M33.
The altitude-azimuth coordinates are then found using
`astropy.coordinates.EarthLocation` and `astropy.time.Time` objects.

This example is meant to demonstrate the capabilities of the
`astropy.coordinates` package. For more convenient and/or complex observation
planning, consider the `astroplan <https://astroplan.readthedocs.io/>`_
package.

Let's suppose you are planning to visit picturesque Bear Mountain State Park
in New York, USA. You are bringing your telescope with you (of course), and
someone told you M33 is a great target to observe there. You happen to know
you are free at 11:00 PM local time, and you want to know if it will be up.
Astropy can answer that.

.. plot::
   :include-source:

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
    >>> from astropy.time import Time
    >>> from astropy.visualization import quantity_support

    :meth:`astropy.coordinates.SkyCoord.from_name` uses Simbad to resolve object
    names and retrieve coordinates.

    Get the coordinates of M33:

    >>> m33 = SkyCoord.from_name("M33")  # doctest: +SKIP

    .. testsetup:: sphx_glr_generated_examples_coordinates_plot_obs-planning

        >>> m33 = SkyCoord(23.46206906, 30.66017511, unit="deg")

    Use `astropy.coordinates.EarthLocation` to provide the location of Bear
    Mountain and set the time to 11pm Eastern Daylight Time (EDT) on 2012 July 12:

    >>> bear_mountain = EarthLocation(lat=41.3 * u.deg, lon=-74 * u.deg, height=390 * u.m)
    >>> utcoffset = -4 * u.hour  # EDT
    >>> time = Time("2012-7-12 23:00:00") - utcoffset

    :meth:`astropy.coordinates.EarthLocation.get_site_names` can be used to get
    locations of major observatories.

    Use `astropy.coordinates` to find the Alt, Az coordinates of M33 at as
    observed from Bear Mountain at 11pm on 2012 July 12:

    >>> m33altaz = m33.transform_to(AltAz(obstime=time, location=bear_mountain))
    >>> print(f"M33's Altitude = {m33altaz.alt:.2}")
    M33's Altitude = 0.13 deg

    This is helpful since it turns out M33 is barely above the horizon at this
    time. It is more informative to find M33's airmass over the course of
    the night.

    Find the Alt, Az coordinates of M33 at 100 times evenly spaced between 10 PM
    and 7 AM EDT:

    >>> midnight = Time("2012-7-13 00:00:00") - utcoffset
    >>> delta_midnight = np.linspace(-2, 10, 100) * u.hour
    >>> frame_July13night = AltAz(obstime=midnight + delta_midnight, location=bear_mountain)
    >>> m33altazs_July13night = m33.transform_to(frame_July13night)

    Convert Alt, Az to airmass with `~astropy.coordinates.AltAz.secz` attribute:

    >>> m33airmasss_July13night = m33altazs_July13night.secz

    Plot the airmass as a function of time:

    >>> with quantity_support():
    ...     fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ...     _ = ax.plot(delta_midnight, m33airmasss_July13night)
    ...     _ = ax.set_xlim(-2, 10)
    ...     _ = ax.set_ylim(1, 4)
    ...     _ = ax.set_xlabel("Hours from EDT Midnight")
    ...     _ = ax.set_ylabel("Airmass [Sec(z)]")
    ...     plt.draw()

    Use  :func:`~astropy.coordinates.get_sun` to find the location of the Sun at 1000
    evenly spaced times between noon on July 12 and noon on July 13:

    >>> delta_midnight = np.linspace(-12, 12, 1000) * u.hour
    >>> times_July12_to_13 = midnight + delta_midnight
    >>> frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=bear_mountain)
    >>> sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)

    Do the same with :func:`~astropy.coordinates.get_body` to find when the moon is
    up. Be aware that this will need to download a 10 MB file from the internet
    to get a precise location of the moon.

    >>> moon_July12_to_13 = get_body("moon", times_July12_to_13)
    >>> moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)

    Find the Alt, Az coordinates of M33 at those same times:

    >>> m33altazs_July12_to_13 = m33.transform_to(frame_July12_to_13)

    Make a figure illustrating nighttime and the altitudes of M33 and
    the Sun over that time:

    >>> with quantity_support():
    ...     fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ...     _ = ax.plot(delta_midnight, sunaltazs_July12_to_13.alt, color="r", label="Sun")
    ...     _ = ax.plot(
    ...         delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75] * 3, ls="--", label="Moon"
    ...     )
    ...     mappable = ax.scatter(
    ...         delta_midnight,
    ...         m33altazs_July12_to_13.alt,
    ...         c=m33altazs_July12_to_13.az.value,
    ...         label="M33",
    ...         lw=0,
    ...         s=8,
    ...         cmap="viridis",
    ...     )
    ...     _ = ax.fill_between(
    ...         delta_midnight,
    ...         0 * u.deg,
    ...         90 * u.deg,
    ...         sunaltazs_July12_to_13.alt < (-0 * u.deg),
    ...         color="0.5",
    ...         zorder=0,
    ...     )
    ...     _ = ax.fill_between(
    ...         delta_midnight,
    ...         0 * u.deg,
    ...         90 * u.deg,
    ...         sunaltazs_July12_to_13.alt < (-18 * u.deg),
    ...         color="k",
    ...         zorder=0,
    ...     )
    ...     _ = fig.colorbar(mappable).set_label("Azimuth [deg]")
    ...     _ = ax.legend(loc="upper left")
    ...     _ = ax.set_xlim(-12 * u.hour, 12 * u.hour)
    ...     _ = ax.set_xticks((np.arange(13) * 2 - 12) * u.hour)
    ...     _ = ax.set_ylim(0 * u.deg, 90 * u.deg)
    ...     _ = ax.set_xlabel("Hours from EDT Midnight")
    ...     _ = ax.set_ylabel("Altitude [deg]")
    ...     _ = ax.grid(visible=True)
    ...     plt.draw()

..
  EXAMPLE END


.. _sphx_glr_generated_examples_coordinates_plot_sgr-coordinate-frame.py:

Create a new coordinate class (for the Sagittarius stream)
==========================================================

..
  EXAMPLE START
  Create a new coordinate class (for the Sagittarius stream)

This document describes in detail how to subclass and define a custom spherical
coordinate frame, as discussed in :ref:`astropy:astropy-coordinates-design` and
the docstring for `~astropy.coordinates.BaseCoordinateFrame`. In this example,
we will define a coordinate system defined by the plane of orbit of the
Sagittarius Dwarf Galaxy (hereafter Sgr; as defined in Majewski et al. 2003).
The Sgr coordinate system is often referred to in terms of two angular
coordinates, :math:`\Lambda,B`.

To do this, we need to define a subclass of
`~astropy.coordinates.BaseCoordinateFrame` that knows the names and units of the
coordinate system angles in each of the supported representations. In this case
we support `~astropy.coordinates.SphericalRepresentation` with "Lambda" and
"Beta". Then we have to define the transformation from this coordinate system to
some other built-in system. Here we will use Galactic coordinates, represented
by the `~astropy.coordinates.Galactic` class.

.. seealso::

    The `gala package <http://gala.adrian.pw/>`_
        Defines a number of Astropy coordinate frames for
        stellar stream coordinate systems.

    Majewski et al. 2003
        "A Two Micron All Sky Survey View of the Sagittarius
        Dwarf Galaxy. I. Morphology of the Sagittarius Core and Tidal Arms",
        https://arxiv.org/abs/astro-ph/0304198

    Law & Majewski 2010
        "The Sagittarius Dwarf Galaxy: A Model for Evolution in a
        Triaxial Milky Way Halo", https://arxiv.org/abs/1003.1132

    David Law's Sgr info page
        https://www.stsci.edu/~dlaw/Sgr/

.. plot::
   :include-source:

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import astropy.coordinates as coord
    >>> from astropy import units as u
    >>> from astropy.coordinates import frame_transform_graph
    >>> from astropy.coordinates.matrix_utilities import matrix_transpose, rotation_matrix

    The first step is to create a new class, which we'll call
    ``Sagittarius`` and make it a subclass of
    `~astropy.coordinates.BaseCoordinateFrame`:

    >>> class Sagittarius(coord.BaseCoordinateFrame):
    ...     """A Heliocentric spherical coordinate system defined by the orbit
    ...     of the Sagittarius dwarf galaxy, as described in
    ...         https://ui.adsabs.harvard.edu/abs/2003ApJ...599.1082M
    ...     and further explained in
    ...         https://www.stsci.edu/~dlaw/Sgr/.
    ...
    ...     Parameters
    ...     ----------
    ...     representation : `~astropy.coordinates.BaseRepresentation` or None
    ...         A representation object or None to have no data (or use the other keywords)
    ...     Lambda : `~astropy.coordinates.Angle`, optional, must be keyword
    ...         The longitude-like angle corresponding to Sagittarius' orbit.
    ...     Beta : `~astropy.coordinates.Angle`, optional, must be keyword
    ...         The latitude-like angle corresponding to Sagittarius' orbit.
    ...     distance : `~astropy.units.Quantity`, optional, must be keyword
    ...         The Distance for this object along the line-of-sight.
    ...     pm_Lambda_cosBeta : `~astropy.units.Quantity`, optional, must be keyword
    ...         The proper motion along the stream in ``Lambda`` (including the
    ...         ``cos(Beta)`` factor) for this object (``pm_Beta`` must also be given).
    ...     pm_Beta : `~astropy.units.Quantity`, optional, must be keyword
    ...         The proper motion in Declination for this object (``pm_ra_cosdec`` must
    ...         also be given).
    ...     radial_velocity : `~astropy.units.Quantity`, optional, keyword-only
    ...         The radial velocity of this object.
    ...     """
    ...     default_representation = coord.SphericalRepresentation
    ...     default_differential = coord.SphericalCosLatDifferential
    ...     frame_specific_representation_info = {
    ...         coord.SphericalRepresentation: [
    ...             coord.RepresentationMapping("lon", "Lambda"),
    ...             coord.RepresentationMapping("lat", "Beta"),
    ...             coord.RepresentationMapping("distance", "distance"),
    ...         ]
    ...     }

    Breaking this down line-by-line, we define the class as a subclass of
    `~astropy.coordinates.BaseCoordinateFrame`. Then we include a descriptive
    docstring. The final lines are class-level attributes that specify the
    default representation for the data, default differential for the velocity
    information, and mappings from the attribute names used by representation
    objects to the names that are to be used by the ``Sagittarius`` frame. In this
    case we override the names in the spherical representations but do not do
    anything with other representations like cartesian or cylindrical.

    Next we have to define the transformation from this coordinate system to some
    other built-in coordinate system; we will use Galactic coordinates. We can do
    this by defining functions that return transformation matrices, or by simply
    defining a function that accepts a coordinate and returns a new coordinate in
    the new system. Because the transformation to the Sagittarius coordinate
    stem is just a spherical rotation from Galactic coordinates, we will
    define a function that returns this matrix. We will start by constructing the
    transformation matrix using pre-determined Euler angles and the
    ``rotation_matrix`` helper function:

    >>> SGR_PHI = (180 + 3.75) * u.degree  # Euler angles (from Law & Majewski 2010)
    >>> SGR_THETA = (90 - 13.46) * u.degree
    >>> SGR_PSI = (180 + 14.111534) * u.degree

    Generate the rotation matrix using the x-convention (see Goldstein):

    >>> SGR_MATRIX = (
    ...     np.diag([1.0, 1.0, -1.0])
    ...     @ rotation_matrix(SGR_PSI, "z")
    ...     @ rotation_matrix(SGR_THETA, "x")
    ...     @ rotation_matrix(SGR_PHI, "z")
    ... )

    Since we already constructed the transformation (rotation) matrix above, and
    the inverse of a rotation matrix is just its transpose, the required
    transformation functions are very simple:

    >>> @frame_transform_graph.transform(
    ...     coord.StaticMatrixTransform, coord.Galactic, Sagittarius
    ... )
    ... def galactic_to_sgr():
    ...     """Compute the Galactic spherical to heliocentric Sgr transformation matrix."""
    ...     return SGR_MATRIX

    The decorator ``@frame_transform_graph.transform(coord.StaticMatrixTransform, coord.Galactic, Sagittarius)``
    registers this function on the
    ``frame_transform_graph`` as a coordinate transformation. Inside the function,
    we return the previously defined rotation matrix.

    We then register the inverse transformation by using the transpose of the
    rotation matrix (which is faster to compute than the inverse):

    >>> @frame_transform_graph.transform(
    ...     coord.StaticMatrixTransform, Sagittarius, coord.Galactic
    ... )
    ... def sgr_to_galactic():
    ...     """Compute the heliocentric Sgr to spherical Galactic transformation matrix."""
    ...     return matrix_transpose(SGR_MATRIX)

    Now that we have registered these transformations between ``Sagittarius`` and
    `~astropy.coordinates.Galactic`, we can transform between *any* coordinate
    system and ``Sagittarius`` (as long as the other system has a path to
    transform to `~astropy.coordinates.Galactic`). For example, to transform from
    ICRS coordinates to ``Sagittarius``, we would do:

    >>> icrs = coord.SkyCoord(280.161732 * u.degree, 11.91934 * u.degree, frame="icrs")
    >>> sgr = icrs.transform_to(Sagittarius)
    >>> print(sgr)
    <SkyCoord (Sagittarius): (Lambda, Beta) in deg
        (346.81830652, -39.28360407)>

    Or, to transform from the ``Sagittarius`` frame to ICRS coordinates (in this
    case, a line along the ``Sagittarius`` x-y plane):

    >>> sgr = coord.SkyCoord(
    ...     Lambda=np.linspace(0, 2 * np.pi, 128) * u.radian,
    ...     Beta=np.zeros(128) * u.radian,
    ...     frame="sagittarius",
    ... )
    >>> icrs = sgr.transform_to(coord.ICRS)
    >>> print(icrs)  # doctest: +ELLIPSIS
    <SkyCoord (ICRS): (ra, dec) in deg
        [(284.03876751, -29.00408353), (287.24685769, -29.44848352),
         (290.48068369, -29.81535572), (293.7357366 , -30.1029631 ),
         ...]>

    As an example, we will now plot the points in both coordinate systems:

    >>> fig, axes = plt.subplots(2, 1, figsize=(8, 10), subplot_kw={"projection": "aitoff"})
    >>> _ = axes[0].set_title("Sagittarius")
    >>> _ = axes[0].plot(
    ...     sgr.Lambda.wrap_at(180 * u.deg).radian,
    ...     sgr.Beta.radian,
    ...     linestyle="none",
    ...     marker=".",
    ... )
    >>> _ = axes[0].grid(visible=True)
    >>> _ = axes[1].set_title("ICRS")
    >>> _ = axes[1].plot(
    ...     icrs.ra.wrap_at(180 * u.deg).radian, icrs.dec.radian, linestyle="none", marker="."
    ... )
    >>> _ = axes[1].grid(visible=True)

    This particular transformation is just a spherical rotation, which is a
    special case of an Affine transformation with no vector offset. The
    transformation of velocity components is therefore natively supported as
    well:

    >>> sgr = coord.SkyCoord(
    ...     Lambda=np.linspace(0, 2 * np.pi, 128) * u.radian,
    ...     Beta=np.zeros(128) * u.radian,
    ...     pm_Lambda_cosBeta=np.random.uniform(-5, 5, 128) * (u.mas / u.yr),
    ...     pm_Beta=np.zeros(128) * (u.mas / u.yr),
    ...     frame="sagittarius",
    ... )
    >>> icrs = sgr.transform_to(coord.ICRS)
    >>> print(icrs)  # doctest: +ELLIPSIS
    <SkyCoord (ICRS): (ra, dec) in deg
        [(284.03876751, -29.00408353), (287.24685769, -29.44848352),
         ...,
         ...]>
    >>> fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    >>> _ = axes[0].set_title("Sagittarius")
    >>> _ = axes[0].plot(
    ...     sgr.Lambda.degree, sgr.pm_Lambda_cosBeta.value, linestyle="none", marker="."
    ... )
    >>> _ = axes[0].set_xlabel(r"$\Lambda$ [deg]")
    >>> _ = axes[0].set_ylabel(
    ...     rf"$\mu_\Lambda \, \cos B$ [{sgr.pm_Lambda_cosBeta.unit.to_string('latex_inline')}]"
    ... )
    >>> _ = axes[0].grid(visible=True)
    >>> _ = axes[1].set_title("ICRS")
    >>> _ = axes[1].plot(icrs.ra.degree, icrs.pm_ra_cosdec.value, linestyle="none", marker=".")
    >>> _ = axes[1].set_ylabel(
    ...     rf"$\mu_\alpha \, \cos\delta$ [{icrs.pm_ra_cosdec.unit.to_string('latex_inline')}]"
    ... )
    >>> _ = axes[1].grid(visible=True)
    >>> _ = axes[2].set_title("ICRS")
    >>> _ = axes[2].plot(icrs.ra.degree, icrs.pm_dec.value, linestyle="none", marker=".")
    >>> _ = axes[2].set_xlabel("RA [deg]")
    >>> _ = axes[2].set_ylabel(rf"$\mu_\delta$ [{icrs.pm_dec.unit.to_string('latex_inline')}]")
    >>> _ = axes[2].grid(visible=True)
    >>> plt.draw()

..
  EXAMPLE END


.. _sphx_glr_generated_examples_coordinates_rv-to-gsr.py:

Convert a radial velocity to the Galactic Standard of Rest (GSR)
================================================================

..
  EXAMPLE START
  Convert a radial velocity to the Galactic Standard of Rest (GSR)

Radial or line-of-sight velocities of sources are often reported in a
Heliocentric or Solar-system barycentric reference frame. A common
transformation incorporates the projection of the Sun's motion along the
line-of-sight to the target, hence transforming it to a Galactic rest frame
instead (sometimes referred to as the Galactic Standard of Rest, GSR). This
transformation depends on the assumptions about the orientation of the Galactic
frame relative to the bary- or Heliocentric frame. It also depends on the
assumed solar velocity vector. Here we will demonstrate how to perform this
transformation using a sky position and barycentric radial-velocity.

Use the latest convention for the Galactocentric coordinates:

>>> import astropy.coordinates as coord
>>> _ = coord.galactocentric_frame_defaults.set("latest")

For this example, let's work with the coordinates and barycentric radial
velocity of the star HD 155967, as obtained from
`Simbad <https://simbad.unistra.fr/simbad/>`_:

>>> from astropy import units as u
>>> icrs = coord.SkyCoord(
...     ra=258.58356362 * u.deg,
...     dec=14.55255619 * u.deg,
...     radial_velocity=-16.1 * u.km / u.s,
...     frame="icrs",
... )

Next, we need to decide on the velocity of the Sun in the assumed GSR frame.
We will use the same velocity vector as used in the
`~astropy.coordinates.Galactocentric` frame, and convert it to a
`~astropy.coordinates.CartesianRepresentation` object using the
``.to_cartesian()`` method of the
`~astropy.coordinates.CartesianDifferential` object ``galcen_v_sun``:

>>> v_sun = coord.Galactocentric().galcen_v_sun.to_cartesian()

We now need to get a unit vector in the assumed Galactic frame from the sky
position in the ICRS frame above. We will use this unit vector to project the
solar velocity onto the line-of-sight:

>>> gal = icrs.transform_to(coord.Galactic)
>>> cart_data = gal.data.to_cartesian()
>>> unit_vector = cart_data / cart_data.norm()

Now we project the solar velocity using this unit vector:

>>> v_proj = v_sun.dot(unit_vector)

Finally, we add the projection of the solar velocity to the radial velocity
to get a GSR radial velocity:

>>> rv_gsr = icrs.radial_velocity + v_proj
>>> print(rv_gsr)  # doctest: +FLOAT_CMP
123.30460087379765 km / s

We could wrap this in a function so we can control the solar velocity and
reuse the above code:

>>> def rv_to_gsr(c, v_sun=None):
...     """Transform a barycentric radial velocity to the Galactic Standard of Rest
...     (GSR).
...
...     Parameters
...     ----------
...     c : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
...         The radial velocity, associated with a sky coordinates, to be
...         transformed.
...     v_sun : `~astropy.units.Quantity`, optional
...         The 3D velocity of the solar system barycenter in the GSR frame.
...         Defaults to the same solar motion as in the
...         `~astropy.coordinates.Galactocentric` frame.
...
...     Returns
...     -------
...     v_gsr : `~astropy.units.Quantity`
...         The input radial velocity transformed to a GSR frame.
...     """
...     if v_sun is None:
...         v_sun = coord.Galactocentric().galcen_v_sun.to_cartesian()
...
...     gal = c.transform_to(coord.Galactic)
...     cart_data = gal.data.to_cartesian()
...     unit_vector = cart_data / cart_data.norm()
...
...     v_proj = v_sun.dot(unit_vector)
...
...     return c.radial_velocity + v_proj

>>> rv_gsr = rv_to_gsr(icrs)
>>> print(rv_gsr)  # doctest: +FLOAT_CMP
123.30460087379765 km / s

..
  EXAMPLE END
