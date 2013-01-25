# coding: utf-8

""" Astropy coordinate class for the Sagittarius coordinate system """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os, sys

# Third-party
import numpy as np
from numpy import radians, degrees, cos, sin

import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import transformations
from astropy.coordinates.angles import rotation_matrix

__all__ = ["SgrCoordinates"]

class SgrCoordinates(coord.SphericalCoordinatesBase):
    """ A spherical coordinate system defined by the orbit of the Sagittarius 
        dwarf galaxy, as described in 
            http://adsabs.harvard.edu/abs/2003ApJ...599.1082M 
        and further explained in
            http://www.astro.virginia.edu/~srm4n/Sgr/.

    """
    __doc__ = __doc__.format(params=coord.SphericalCoordinatesBase. \
                                          _init_docstring_param_templ. \
                                          format(lonnm='Lambda', latnm='Beta'))

    def __init__(self, *args, **kwargs):
        super(SgrCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and \
            isinstance(args[0], coord.SphericalCoordinatesBase):
            
            newcoord = args[0].transform_to(self.__class__)
            self.Lambda = newcoord.Lambda
            self.Beta = newcoord.Beta
            self._distance = newcoord._distance
        else:
            super(SgrCoordinates, self).\
                _initialize_latlon('Lambda', 'Beta', False, args, kwargs, 
                                   anglebounds=((0, 360), (-90,90)))

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, 
                                                        self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} Lambda={1:.5f} deg, Beta={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.Lambda.degrees,
                          self.Beta.degrees, diststr)

    @property
    def lonangle(self):
        return self.Lambda

    @property
    def latangle(self):
        return self.Beta

# Define the Euler angles (from Law & Majewski 2010)
phi = radians(180+3.75)
theta = radians(90-13.46)
psi = radians(180+14.111534)

# Generate the rotation matrix using the x-convention (see Goldstein)
D = rotation_matrix(phi, "z", degrees=False)
C = rotation_matrix(theta, "x", degrees=False)
B = rotation_matrix(psi, "z", degrees=False)
sgr_matrix = np.array(B.dot(C).dot(D))

# Galactic to Sgr coordinates
@transformations.transform_function(coord.GalacticCoordinates, SgrCoordinates)
def galactic_to_sgr(galactic_coord):
    """ Compute the transformation from Galactic spherical to Sgr coordinates. 
    """

    l = galactic_coord.l.radians
    b = galactic_coord.b.radians

    X = cos(b)*cos(l)
    Y = cos(b)*sin(l)
    Z = sin(b)

    # Calculate X,Y,Z,distance in the Sgr system
    Xs, Ys, Zs = sgr_matrix.dot(np.array([X, Y, Z]))

    Zs = -Zs

    # Calculate the angular coordinates lambda,beta
    Lambda = degrees(np.arctan2(Ys,Xs))
    if Lambda<0:
        Lambda += 360

    Beta = degrees(np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs)))

    return SgrCoordinates(Lambda, Beta, distance=galactic_coord.distance, 
                          unit=(u.degree, u.degree))

@transformations.transform_function(SgrCoordinates, coord.GalacticCoordinates)
def sgr_to_galactic(sgr_coord):
    """ Compute the transformation from Sgr coordinates to spherical Galactic. 
    """
    L = sgr_coord.Lambda.radians
    B = sgr_coord.Beta.radians

    Xs = cos(B)*cos(L)
    Ys = cos(B)*sin(L)
    Zs = sin(B)
    Zs = -Zs

    X, Y, Z = sgr_matrix.T.dot(np.array([Xs, Ys, Zs]))

    l = degrees(np.arctan2(Y,X))
    b = degrees(np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z)))

    if l<0:
        l += 360

    return coord.GalacticCoordinates(l, b, distance=sgr_coord.distance, 
                                     unit=(u.degree, u.degree))

if __name__ == "__main__":
    # Example use case for our newly defined coordinate class
    ra_dec = coord.ICRSCoordinates(152.88572, 11.57281, unit=(u.degree,
                                                              u.degree))
    sgr = ra_dec.transform_to(SgrCoordinates)
    print(sgr)