# coding: utf-8

""" Astropy coordinate class for the Sagittarius coordinate system """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import numpy as np
from numpy import cos, sin

from astropy.coordinates import frame_transform_graph
from astropy.coordinates.angles import rotation_matrix
import astropy.coordinates as coord
import astropy.units as u


__all__ = ["Sagittarius"]


class Sagittarius(coord.BaseCoordinateFrame):
    """
    A Heliocentric spherical coordinate system defined by the orbit
    of the Sagittarius dwarf galaxy, as described in
        http://adsabs.harvard.edu/abs/2003ApJ...599.1082M
    and further explained in
        http://www.astro.virginia.edu/~srm4n/Sgr/.

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

    """
    default_representation = coord.SphericalRepresentation

    frame_specific_representation_info = {
        'spherical': [coord.RepresentationMapping('lon', 'Lambda'),
                      coord.RepresentationMapping('lat', 'Beta'),
                      coord.RepresentationMapping('distance', 'distance')],
        'unitspherical': [coord.RepresentationMapping('lon', 'Lambda'),
                          coord.RepresentationMapping('lat', 'Beta')]
    }

# Define the Euler angles (from Law & Majewski 2010)
phi = np.radians(180+3.75)
theta = np.radians(90-13.46)
psi = np.radians(180+14.111534)

# Generate the rotation matrix using the x-convention (see Goldstein)
D = rotation_matrix(phi, "z", unit=u.radian)
C = rotation_matrix(theta, "x", unit=u.radian)
B = rotation_matrix(psi, "z", unit=u.radian)
sgr_matrix = np.array(B.dot(C).dot(D))


# Galactic to Sgr coordinates
@frame_transform_graph.transform(coord.FunctionTransform, coord.Galactic, Sagittarius)
def galactic_to_sgr(gal_coord, sgr_frame):
    """ Compute the transformation from Galactic spherical to
        heliocentric Sgr coordinates.
    """

    l = np.atleast_1d(gal_coord.l.radian)
    b = np.atleast_1d(gal_coord.b.radian)

    X = cos(b)*cos(l)
    Y = cos(b)*sin(l)
    Z = sin(b)

    # Calculate X,Y,Z,distance in the Sgr system
    Xs, Ys, Zs = sgr_matrix.dot(np.array([X, Y, Z]))
    Zs = -Zs

    # Calculate the angular coordinates lambda,beta
    Lambda = np.arctan2(Ys, Xs)*u.radian
    Lambda[Lambda < 0] = Lambda[Lambda < 0] + 2.*np.pi*u.radian
    Beta = np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs))*u.radian

    return Sagittarius(Lambda=Lambda, Beta=Beta,
                       distance=gal_coord.distance)


# Sgr to Galactic coordinates
@frame_transform_graph.transform(coord.FunctionTransform, Sagittarius, coord.Galactic)
def sgr_to_galactic(sgr_coord, gal_frame):
    """ Compute the transformation from heliocentric Sgr coordinates to
        spherical Galactic.
    """
    L = np.atleast_1d(sgr_coord.Lambda.radian)
    B = np.atleast_1d(sgr_coord.Beta.radian)

    Xs = cos(B)*cos(L)
    Ys = cos(B)*sin(L)
    Zs = sin(B)
    Zs = -Zs

    X, Y, Z = sgr_matrix.T.dot(np.array([Xs, Ys, Zs]))

    l = np.arctan2(Y, X)*u.radian
    b = np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z))*u.radian

    l[l<0] += 2*np.pi*u.radian

    return coord.Galactic(l=l, b=b, distance=sgr_coord.distance)

if __name__ == "__main__":
    # Example use case for our newly defined coordinate class
    icrs = coord.ICRS(152.88572*u.degree, 11.57281*u.degree)
    sgr = icrs.transform_to(Sagittarius)
    print(sgr)
