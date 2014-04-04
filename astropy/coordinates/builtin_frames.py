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
from .baseframe import BaseCoordinateFrame
from . import representation import SphericalRepresentation


__all__ = ['ICRS', 'FK5', 'FK4', 'FK4NoETerms', 'Galactic', 'AltAz']

class ICRS(BaseCoordinateFrame):
    """
    A coordinate in the ICRS.

    If you're looking for "J2000" coordinates, and aren't sure if you
    want to use this or `FK5`, you probably want to use ICRS.
    It's more well-defined as a catalog coordinate and is an inertial
    system.
    """
    
    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'lon': 'ra', 'lat': 'dec'}


class FK5(BaseCoordinateFrame):
    """
    docstr
    """
    
    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'lon': 'ra', 'lat': 'dec'}

class FK4(BaseCoordinateFrame):
    """
    docstr
    """
    
    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'lon': 'ra', 'lat': 'dec'}

class FK4NoETerms(BaseCoordinateFrame):
    """
    docstr
    """
    
    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'lon': 'ra', 'lat': 'dec'}

class Galactic(BaseCoordinateFrame):
    """
    docstr
    """
    
    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'lon': 'l', 'lat': 'b'}

class AltAz(BaseCoordinateFrame):
    """
    docstr
    """
    
    preferred_representation = SphericalRepresentation
    preferred_attr_names = {'lon': 'az', 'lat': 'alt'}
