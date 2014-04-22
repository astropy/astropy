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


class FK5(BaseCoordinateFrame):
    """
    docstr
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'ra': 'lon', 'dec': 'lat', 'distance': 'distance'}
    frame_attr_names = {'equinox': _EQUINOX_J2000}


class FK4(BaseCoordinateFrame):
    """
    docstr
    """

    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'ra': 'lon', 'dec': 'lat', 'distance': 'distance'}
    frame_attr_names = {'equinox': _EQUINOX_B1950, 'obstime': None}


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
@frame_transform_graph.transform(StaticMatrixTransform, ICRS, FK5)
def icrs_to_fk5():
    """
    B-matrix from USNO circular 179
    """
    from .angles import rotation_matrix

    eta0 = -19.9 / 3600000.
    xi0 = 9.1 / 3600000.
    da0 = -22.9 / 3600000.

    m1 = rotation_matrix(-eta0, 'x')
    m2 = rotation_matrix(xi0, 'y')
    m3 = rotation_matrix(da0, 'z')

    return m1 * m2 * m3


# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, ICRS)
def fk5_to_icrs(fk5c):
    from .earth_orientation import _precess_from_J2000_Capitaine

    pmat = _precess_from_J2000_Capitaine(fk5c.equinox.jyear).T

    # transpose gets equinox -> J2000
    fk5toicrsmat = icrs_to_fk5().T

    return fk5toicrsmat * pmat
