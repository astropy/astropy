# -*- coding: utf-8 -*-
"""
==========================================================
Create a new coordinate class (for the Sagittarius stream)
==========================================================

This document describes in detail how to subclass and define a custom spherical
coordinate frame, as discussed in :ref:`astropy-coordinates-design` and the
docstring for `~astropy.coordinates.BaseCoordinateFrame`. In this example, we
will define a coordinate system defined by the plane of orbit of the Sagittarius
Dwarf Galaxy (hereafter Sgr; as defined in Majewski et al. 2003).  The Sgr
coordinate system is often referred to in terms of two angular coordinates,
:math:`\Lambda,B`.

To do this, wee need to define a subclass of
`~astropy.coordinates.BaseCoordinateFrame` that knows the names and units of the
coordinate system angles in each of the supported representations.  In this case
we support `~astropy.coordinates.SphericalRepresentation` with "Lambda" and
"Beta". Then we have to define the transformation from this coordinate system to
some other built-in system. Here we will use Galactic coordinates, represented
by the `~astropy.coordinates.Galactic` class.

See Also
--------

* The `gala package <http://gala.adrian.pw/>`_, which defines a number of
  Astropy coordinate frames for stellar stream coordinate systems.
* Majewski et al. 2003, "A Two Micron All Sky Survey View of the Sagittarius
  Dwarf Galaxy. I. Morphology of the Sagittarius Core and Tidal Arms",
  https://arxiv.org/abs/astro-ph/0304198
* Law & Majewski 2010, "The Sagittarius Dwarf Galaxy: A Model for Evolution in a
  Triaxial Milky Way Halo", https://arxiv.org/abs/1003.1132
* David Law's Sgr info page http://www.stsci.edu/~dlaw/Sgr/

-------------------

*By: Adrian Price-Whelan, Erik Tollerud*

*License: BSD*

-------------------

"""

##############################################################################
# Make `print` work the same in all versions of Python, set up numpy,
# matplotlib, and use a nicer set of plot parameters:

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)


##############################################################################
# Import the packages necessary for coordinates

from astropy.coordinates import frame_transform_graph
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
import astropy.coordinates as coord
import astropy.units as u

##############################################################################
# The first step is to create a new class, which we'll call
# ``Sagittarius`` and make it a subclass of
# `~astropy.coordinates.BaseCoordinateFrame`:

class Sagittarius(coord.BaseCoordinateFrame):
    """
    A Heliocentric spherical coordinate system defined by the orbit
    of the Sagittarius dwarf galaxy, as described in
        http://adsabs.harvard.edu/abs/2003ApJ...599.1082M
    and further explained in
        http://www.stsci.edu/~dlaw/Sgr/.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    Lambda : `Angle`, optional, must be keyword
        The longitude-like angle corresponding to Sagittarius' orbit.
    Beta : `Angle`, optional, must be keyword
        The latitude-like angle corresponding to Sagittarius' orbit.
    distance : `Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
    pm_Lambda_cosBeta : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion along the stream in ``Lambda`` (including the
        ``cos(Beta)`` factor) for this object (``pm_Beta`` must also be given).
    pm_Beta : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Declination for this object (``pm_ra_cosdec`` must
        also be given).
    radial_velocity : :class:`~astropy.units.Quantity`, optional, must be keyword
        The radial velocity of this object.

    """

    default_representation = coord.SphericalRepresentation
    default_differential = coord.SphericalCosLatDifferential

    frame_specific_representation_info = {
        coord.SphericalRepresentation: [
            coord.RepresentationMapping('lon', 'Lambda'),
            coord.RepresentationMapping('lat', 'Beta'),
            coord.RepresentationMapping('distance', 'distance')]
    }

##############################################################################
# Breaking this down line-by-line, we define the class as a subclass of
# `~astropy.coordinates.BaseCoordinateFrame`. Then we include a descriptive
# docstring.  The final lines are class-level attributes that specify the
# default representation for the data, default differential for the velocity
# information, and mappings from the attribute names used by representation
# objects to the names that are to be used by the ``Sagittarius`` frame. In this
# case we override the names in the spherical representations but don't do
# anything with other representations like cartesian or cylindrical.
#
# Next we have to define the transformation from this coordinate system to some
# other built-in coordinate system; we will use Galactic coordinates. We can do
# this by defining functions that return transformation matrices, or by simply
# defining a function that accepts a coordinate and returns a new coordinate in
# the new system. Because the transformation to the Sagittarius coordinate
# system is just a spherical rotation from Galactic coordinates, we'll just
# define a function that returns this matrix. We'll start by constructing the
# transformation matrix using pre-deteremined Euler angles and the
# ``rotation_matrix`` helper function:

SGR_PHI = (180 + 3.75) * u.degree # Euler angles (from Law & Majewski 2010)
SGR_THETA = (90 - 13.46) * u.degree
SGR_PSI = (180 + 14.111534) * u.degree

# Generate the rotation matrix using the x-convention (see Goldstein)
D = rotation_matrix(SGR_PHI, "z")
C = rotation_matrix(SGR_THETA, "x")
B = rotation_matrix(SGR_PSI, "z")
A = np.diag([1.,1.,-1.])
SGR_MATRIX = matrix_product(A, B, C, D)

##############################################################################
# Since we already constructed the transformation (rotation) matrix above, and
# the inverse of a rotation matrix is just its transpose, the required
# transformation functions are very simple:

@frame_transform_graph.transform(coord.StaticMatrixTransform, coord.Galactic, Sagittarius)
def galactic_to_sgr():
    """ Compute the transformation matrix from Galactic spherical to
        heliocentric Sgr coordinates.
    """
    return SGR_MATRIX

##############################################################################
# The decorator ``@frame_transform_graph.transform(coord.StaticMatrixTransform,
# coord.Galactic, Sagittarius)``  registers this function on the
# ``frame_transform_graph`` as a coordinate transformation. Inside the function,
# we simply return the previously defined rotation matrix.
#
# We then register the inverse transformation by using the transpose of the
# rotation matrix (which is faster to compute than the inverse):

@frame_transform_graph.transform(coord.StaticMatrixTransform, Sagittarius, coord.Galactic)
def sgr_to_galactic():
    """ Compute the transformation matrix from heliocentric Sgr coordinates to
        spherical Galactic.
    """
    return matrix_transpose(SGR_MATRIX)

##############################################################################
# Now that we've registered these transformations between ``Sagittarius`` and
# `~astropy.coordinates.Galactic`, we can transform between *any* coordinate
# system and ``Sagittarius`` (as long as the other system has a path to
# transform to `~astropy.coordinates.Galactic`). For example, to transform from
# ICRS coordinates to ``Sagittarius``, we would do:

icrs = coord.ICRS(280.161732*u.degree, 11.91934*u.degree)
sgr = icrs.transform_to(Sagittarius)
print(sgr)

##############################################################################
# Or, to transform from the ``Sagittarius`` frame to ICRS coordinates (in this
# case, a line along the ``Sagittarius`` x-y plane):

sgr = Sagittarius(Lambda=np.linspace(0, 2*np.pi, 128)*u.radian,
                  Beta=np.zeros(128)*u.radian)
icrs = sgr.transform_to(coord.ICRS)
print(icrs)

##############################################################################
# As an example, we'll now plot the points in both coordinate systems:

fig, axes = plt.subplots(2, 1, figsize=(8, 10),
                         subplot_kw={'projection': 'aitoff'})

axes[0].set_title("Sagittarius")
axes[0].plot(sgr.Lambda.wrap_at(180*u.deg).radian, sgr.Beta.radian,
             linestyle='none', marker='.')

axes[1].set_title("ICRS")
axes[1].plot(icrs.ra.wrap_at(180*u.deg).radian, icrs.dec.radian,
             linestyle='none', marker='.')

plt.show()

##############################################################################
# This particular transformation is just a spherical rotation, which is a
# special case of an Affine transformation with no vector offset. The
# transformation of velocity components is therefore natively supported as
# well:

sgr = Sagittarius(Lambda=np.linspace(0, 2*np.pi, 128)*u.radian,
                  Beta=np.zeros(128)*u.radian,
                  pm_Lambda_cosBeta=np.random.uniform(-5, 5, 128)*u.mas/u.yr,
                  pm_Beta=np.zeros(128)*u.mas/u.yr)
icrs = sgr.transform_to(coord.ICRS)
print(icrs)

fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

axes[0].set_title("Sagittarius")
axes[0].plot(sgr.Lambda.degree,
             sgr.pm_Lambda_cosBeta.value,
             linestyle='none', marker='.')
axes[0].set_xlabel(r"$\Lambda$ [deg]")
axes[0].set_ylabel(r"$\mu_\Lambda \, \cos B$ [{0}]"
                   .format(sgr.pm_Lambda_cosBeta.unit.to_string('latex_inline')))

axes[1].set_title("ICRS")
axes[1].plot(icrs.ra.degree, icrs.pm_ra_cosdec.value,
             linestyle='none', marker='.')
axes[1].set_ylabel(r"$\mu_\alpha \, \cos\delta$ [{0}]"
                   .format(icrs.pm_ra_cosdec.unit.to_string('latex_inline')))

axes[2].set_title("ICRS")
axes[2].plot(icrs.ra.degree, icrs.pm_dec.value,
             linestyle='none', marker='.')
axes[2].set_xlabel("RA [deg]")
axes[2].set_ylabel(r"$\mu_\delta$ [{0}]"
                   .format(icrs.pm_dec.unit.to_string('latex_inline')))

plt.show()
