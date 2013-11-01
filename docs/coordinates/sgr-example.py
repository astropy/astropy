# coding: utf-8

""" Astropy coordinate class for the Sagittarius coordinate system """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os, sys

# Third-party
import numpy as np
from numpy import cos, sin

import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import transformations
from astropy.coordinates.angles import rotation_matrix

__all__ = ["SgrCoordinates"]

class SgrCoordinates(coord.SphericalCoordinatesBase):
    """ A Heliocentric spherical coordinate system defined by the orbit
        of the Sagittarius dwarf galaxy, as described in
            http://adsabs.harvard.edu/abs/2003ApJ...599.1082M
        and further explained in
            http://www.astro.virginia.edu/~srm4n/Sgr/.

    """
    __doc__ = __doc__.format(params=coord.SphericalCoordinatesBase. \
            _init_docstring_param_templ.format(lonnm='Lambda', latnm='Beta'))

    def __init__(self, *args, **kwargs):
        super(SgrCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and \
            isinstance(args[0], coord.SphericalCoordinatesBase):

            newcoord = args[0].transform_to(self.__class__)
            self._lonangle = newcoord._lonangle
            self._latangle = newcoord._latangle
            self._distance = newcoord._distance
        else:
            super(SgrCoordinates, self).\
                _initialize_latlon('Lambda', 'Beta', args, kwargs)

    # strings used for making __repr__ work
    _repr_lon_name = 'Lambda'
    _repr_lat_name = 'Beta'

    # Default format for to_string
    _default_string_style = 'dmsdms'

    @property
    def Lambda(self):
        return self._lonangle

    @property
    def Beta(self):
        return self._latangle

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
@transformations.transform_function(coord.Galactic, SgrCoordinates)
def galactic_to_sgr(galactic_coord):
    """ Compute the transformation from Galactic spherical to
        heliocentric Sgr coordinates.
    """

    l = np.atleast_1d(galactic_coord.l.radian)
    b = np.atleast_1d(galactic_coord.b.radian)

    X = cos(b)*cos(l)
    Y = cos(b)*sin(l)
    Z = sin(b)

    # Calculate X,Y,Z,distance in the Sgr system
    Xs, Ys, Zs = sgr_matrix.dot(np.array([X, Y, Z]))
    Zs = -Zs

    # Calculate the angular coordinates lambda,beta
    Lambda = np.degrees(np.arctan2(Ys,Xs))
    Lambda[Lambda < 0] = Lambda[Lambda < 0] + 360
    Beta = np.degrees(np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs)))

    return SgrCoordinates(Lambda, Beta, distance=galactic_coord.distance,
                          unit=(u.degree, u.degree))

@transformations.transform_function(SgrCoordinates, coord.Galactic)
def sgr_to_galactic(sgr_coord):
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

    l = np.degrees(np.arctan2(Y,X))
    b = np.degrees(np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z)))

    if l<0:
        l += 360

    return coord.Galactic(l, b, distance=sgr_coord.distance,
                          unit=(u.degree, u.degree))

if __name__ == "__main__":
    # Example use case for our newly defined coordinate class
    ra_dec = coord.ICRS(152.88572, 11.57281, unit=(u.degree,
                                                              u.degree))
    sgr = ra_dec.transform_to(SgrCoordinates)
    print(sgr)
