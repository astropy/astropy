# -*- coding: utf-8 -*-
"""
================================================================
Convert a radial velocity to the Galactic Standard of Rest (GSR)
================================================================

Radial or line-of-sight velocities of sources are often reported in a
Heliocentric or Solar-system barycentric reference frame. A common
transformation incorporates the projection of the Sun's motion along the
line-of-sight to the target, hence transforming it to a Galactic rest frame
instead (sometimes referred to as the Galactic Standard of Rest, GSR). This
transformation depends on the assumptions about the orientation of the Galactic
frame relative to the bary- or Heliocentric frame. It also depends on the
assumed solar velocity vector. Here we'll demonstrate how to perform this
transformation using a sky position and barycentric radial-velocity.


*By: Adrian Price-Whelan*

*License: BSD*


"""

################################################################################
# Make print work the same in all versions of Python and import the required
# Astropy packages:
import astropy.units as u
import astropy.coordinates as coord

################################################################################
# Use the latest convention for the Galactocentric coordinates
coord.galactocentric_frame_defaults.set('latest')

################################################################################
# For this example, let's work with the coordinates and barycentric radial
# velocity of the star HD 155967, as obtained from
# `Simbad <https://simbad.u-strasbg.fr/simbad/>`_:
icrs = coord.SkyCoord(ra=258.58356362*u.deg, dec=14.55255619*u.deg,
                      radial_velocity=-16.1*u.km/u.s, frame='icrs')

################################################################################
# We next need to decide on the velocity of the Sun in the assumed GSR frame.
# We'll use the same velocity vector as used in the
# `~astropy.coordinates.Galactocentric` frame, and convert it to a
# `~astropy.coordinates.CartesianRepresentation` object using the
# ``.to_cartesian()`` method of the
# `~astropy.coordinates.CartesianDifferential` object ``galcen_v_sun``:
v_sun = coord.Galactocentric().galcen_v_sun.to_cartesian()

################################################################################
# We now need to get a unit vector in the assumed Galactic frame from the sky
# position in the ICRS frame above. We'll use this unit vector to project the
# solar velocity onto the line-of-sight:
gal = icrs.transform_to(coord.Galactic)
cart_data = gal.data.to_cartesian()
unit_vector = cart_data / cart_data.norm()

################################################################################
# Now we project the solar velocity using this unit vector:
v_proj = v_sun.dot(unit_vector)

################################################################################
# Finally, we add the projection of the solar velocity to the radial velocity
# to get a GSR radial velocity:
rv_gsr = icrs.radial_velocity + v_proj
print(rv_gsr)


################################################################################
# We could wrap this in a function so we can control the solar velocity and
# re-use the above code:
def rv_to_gsr(c, v_sun=None):
    """Transform a barycentric radial velocity to the Galactic Standard of Rest
    (GSR).

    The input radial velocity must be passed in as a

    Parameters
    ----------
    c : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
        The radial velocity, associated with a sky coordinates, to be
        transformed.
    v_sun : `~astropy.units.Quantity`, optional
        The 3D velocity of the solar system barycenter in the GSR frame.
        Defaults to the same solar motion as in the
        `~astropy.coordinates.Galactocentric` frame.

    Returns
    -------
    v_gsr : `~astropy.units.Quantity`
        The input radial velocity transformed to a GSR frame.

    """
    if v_sun is None:
        v_sun = coord.Galactocentric().galcen_v_sun.to_cartesian()

    gal = c.transform_to(coord.Galactic)
    cart_data = gal.data.to_cartesian()
    unit_vector = cart_data / cart_data.norm()

    v_proj = v_sun.dot(unit_vector)

    return c.radial_velocity + v_proj


rv_gsr = rv_to_gsr(icrs)
print(rv_gsr)
