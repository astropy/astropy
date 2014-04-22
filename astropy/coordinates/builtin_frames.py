# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains the coordinate frames actually implemented by astropy.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library

# Dependencies
import numpy as np

# Project
from ..extern import six
from .. import units as u
from ..time import Time
from .angles import Angle
from .representation import SphericalRepresentation
from .baseframe import BaseCoordinateFrame, frame_transform_graph
from .transformations import StaticMatrixTransform, FunctionTransform, \
                             DynamicMatrixTransform


__all__ = ['ICRS', 'FK5', 'FK4', 'FK4NoETerms', 'Galactic', 'AltAz']

# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
_EQUINOX_J2000 = Time('J2000', scale='utc')
_EQUINOX_B1950 = Time('B1950', scale='tai')

class ICRS(BaseCoordinateFrame):
    """
    A coordinate in the ICRS.

    If you're looking for "J2000" coordinates, and aren't sure if you
    want to use this or `FK5`, you probably want to use ICRS.
    It's more well-defined as a catalog coordinate and is an inertial
    system.
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'ra': 'lon', 'dec': 'lat', 'distance': 'distance'}
    frame_attr_names = {}  # not necessary if empty, but this makes it clearer

    @property
    def equinox(self):
        """
        ICRS is by design very close to equatorial J2000, so we call this the
        equinox for ICRS.
        """
        return _EQUINOX_J2000

    @staticmethod
    def _icrs_to_fk5_matrix():
        """
        B-matrix from USNO circular 179.  Used by the ICRS->FK5 transformation
        functions.
        """
        from .angles import rotation_matrix

        eta0 = -19.9 / 3600000.
        xi0 = 9.1 / 3600000.
        da0 = -22.9 / 3600000.

        m1 = rotation_matrix(-eta0, 'x')
        m2 = rotation_matrix(xi0, 'y')
        m3 = rotation_matrix(da0, 'z')

        return m1 * m2 * m3
# define this because it only needs to be computed once
ICRS._ICRS_TO_FK5_J2000_MAT = ICRS._icrs_to_fk5_matrix()


class FK5(BaseCoordinateFrame):
    """
    docstr
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'ra': 'lon', 'dec': 'lat', 'distance': 'distance'}
    frame_attr_names = {'equinox': _EQUINOX_J2000}

    @staticmethod
    def _precession_matrix(oldequinox, newequinox):
        """
        Compute and return the precession matrix for FK5 based on Capitaine et
        al. 2003/IAU2006.  Used inside some of the transformation functions.

        Parameters
        ----------
        oldequinox : `~astropy.time.Time`
            The equinox to precess from.
        newequinox : `~astropy.time.Time`
            The equinox to precess to.

        Returns
        -------
        newcoord : array
            The precession matrix to transform to the new equinox
        """
        from .earth_orientation import precession_matrix_Capitaine

        return precession_matrix_Capitaine(oldequinox, newequinox)

# Has to be defined at module level, because `transform` needs an FK5 reference
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK5)
def fk5_to_fk5(fk5coord1, fk5frame2):
    return fk5coord1._precession_matrix(fk5coord1.equinox, fk5frame2.equinox)


class FK4(BaseCoordinateFrame):
    """
    docstr
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'ra': 'lon', 'dec': 'lat', 'distance': 'distance'}
    frame_attr_names = {'equinox': _EQUINOX_B1950, 'obstime': None}

    @staticmethod
    def _precession_matrix(oldequinox, newequinox):
        """
        Compute and return the precession matrix for FK4 using Newcomb's method.
        Used inside some of the transformation functions.

        Parameters
        ----------
        oldequinox : `~astropy.time.Time`
            The equinox to precess from.
        newequinox : `~astropy.time.Time`
            The equinox to precess to.

        Returns
        -------
        newcoord : array
            The precession matrix to transform to the new equinox
        """
        from .earth_orientation import _precession_matrix_besselian

        return _precession_matrix_besselian(oldequinox.byear, newequinox.byear)

@frame_transform_graph.transform(DynamicMatrixTransform, FK4, FK4)
def fk4_to_fk4(fk4coord1, fk4frame2):
    return fk4coord1._precession_matrix(fk4coord1.equinox, fk4frame2.equinox)


class FK4NoETerms(BaseCoordinateFrame):
    """
    docstr
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'ra': 'lon', 'dec': 'lat', 'distance': 'distance'}
    frame_attr_names = {'equinox': _EQUINOX_B1950, 'obstime': None}


class Galactic(BaseCoordinateFrame):
    """
    docstr
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'l': 'lon', 'b': 'lat', 'distance': 'distance'}
    frame_attr_names = {}

    # North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
    # transformations to/from FK4/5
    # These are from Reid & Brunthaler 2004
    _ngp_J2000 = FK5(ra=192.859508*u.degree, dec=27.128336*u.degree)
    _lon0_J2000 = Angle(122.932, u.degree)
    # These are from the IAU's definition of galactic coordinates
    _ngp_B1950 = FK4(ra=192.25*u.degree, dec=27.4*u.degree)
    _lon0_B1950 = Angle(123, u.degree)


class AltAz(BaseCoordinateFrame):
    """
    docstr
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'az': 'lon', 'alt': 'lat', 'distance': 'distance'}
    frame_attr_names = {}


#<--------------------------------transformations------------------------------>
# Transformations are defined here instead of in the classes themselves, because
# we need references to the various objects to give to the decorators.

# ICRS to/from FK5
@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, FK5)
def icrs_to_fk5(icrscoord, fk5frame):
    # ICRS equinox should always be J2000, but just in case, use attribute
    pmat = fk5frame._precession_matrix(icrscoord.equinox, fk5frame.equinox)
    return pmat * icrscoord._ICRS_TO_FK5_J2000_MAT


# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, ICRS)
def fk5_to_icrs(fk5coord, icrsframe):
    # ICRS equinox should always be J2000, but just in case, use attribute
    pmat = fk5coord._precession_matrix(fk5coord.equinox, icrsframe.equinox)
    return icrsframe._ICRS_TO_FK5_J2000_MAT.T * pmat



# Galactic to/from FK4/FK5
# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, Galactic)
def fk5_to_gal(fk5coord, galframe):
    from .angles import rotation_matrix

    #need precess to J2000 first
    pmat = fk5coord._precession_matrix(fk5coord.equinox, _EQUINOX_J2000)
    mat1 = rotation_matrix(180 - Galactic._lon0_J2000.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_J2000.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_J2000.ra.degree, 'z')

    return mat1 * mat2 * mat3 * pmat


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK5)
def _gal_to_fk5(galcoord, fk5frame):
    return fk5_to_gal(fk5frame, galcoord).T


@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, Galactic)
def fk4_to_gal(fk4coords, galframe):
    from .angles import rotation_matrix

    mat1 = rotation_matrix(180 - Galactic._lon0_B1950.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_B1950.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_B1950.ra.degree, 'z')
    matprec = fk4coords._precession_matrix_besselian(fk4coords.equinox, _EQUINOX_B1950)

    return mat1 * mat2 * mat3 * matprec


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK4NoETerms)
def gal_to_fk4(galcoords, fk4frame):
    return fk4_to_gal(fk4frame, galcoords).T
