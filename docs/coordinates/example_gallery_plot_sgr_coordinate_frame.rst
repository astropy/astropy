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
    >>> axes[0].set_title("Sagittarius")  # doctest: +IGNORE_OUTPUT
    >>> axes[0].plot(
    ...     sgr.Lambda.wrap_at(180 * u.deg).radian,
    ...     sgr.Beta.radian,
    ...     linestyle="none",
    ...     marker=".",
    ... )  # doctest: +IGNORE_OUTPUT
    >>> axes[0].grid(visible=True)  # doctest: +IGNORE_OUTPUT
    >>> axes[1].set_title("ICRS")  # doctest: +IGNORE_OUTPUT
    >>> axes[1].plot(
    ...     icrs.ra.wrap_at(180 * u.deg).radian, icrs.dec.radian, linestyle="none", marker="."
    ... )  # doctest: +IGNORE_OUTPUT
    >>> axes[1].grid(visible=True)  # doctest: +IGNORE_OUTPUT

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
    >>> axes[0].set_title("Sagittarius")  # doctest: +IGNORE_OUTPUT
    >>> axes[0].plot(
    ...     sgr.Lambda.degree, sgr.pm_Lambda_cosBeta.value, linestyle="none", marker="."
    ... )  # doctest: +IGNORE_OUTPUT
    >>> axes[0].set_xlabel(r"$\Lambda$ [deg]")  # doctest: +IGNORE_OUTPUT
    >>> axes[0].set_ylabel(
    ...     rf"$\mu_\Lambda \, \cos B$ [{sgr.pm_Lambda_cosBeta.unit.to_string('latex_inline')}]"
    ... )  # doctest: +IGNORE_OUTPUT
    >>> axes[0].grid(visible=True)  # doctest: +IGNORE_OUTPUT
    >>> axes[1].set_title("ICRS")  # doctest: +IGNORE_OUTPUT
    >>> axes[1].plot(icrs.ra.degree, icrs.pm_ra_cosdec.value, linestyle="none", marker=".")  # doctest: +IGNORE_OUTPUT
    >>> axes[1].set_ylabel(
    ...     rf"$\mu_\alpha \, \cos\delta$ [{icrs.pm_ra_cosdec.unit.to_string('latex_inline')}]"
    ... )  # doctest: +IGNORE_OUTPUT
    >>> axes[1].grid(visible=True)  # doctest: +IGNORE_OUTPUT
    >>> axes[2].set_title("ICRS")  # doctest: +IGNORE_OUTPUT
    >>> axes[2].plot(icrs.ra.degree, icrs.pm_dec.value, linestyle="none", marker=".")  # doctest: +IGNORE_OUTPUT
    >>> axes[2].set_xlabel("RA [deg]")  # doctest: +IGNORE_OUTPUT
    >>> axes[2].set_ylabel(rf"$\mu_\delta$ [{icrs.pm_dec.unit.to_string('latex_inline')}]")  # doctest: +IGNORE_OUTPUT
    >>> axes[2].grid(visible=True)  # doctest: +IGNORE_OUTPUT
    >>> plt.draw()

..
  EXAMPLE END
