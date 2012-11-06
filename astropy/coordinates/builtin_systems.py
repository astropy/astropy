# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the implementations of specific coordinate systems
and the conversions between them.
"""

from .angles import RA, Dec, Angle
from .coordsystems import SphericalCoordinatesBase
from .. import units as u

__all__ = ['ICRSCoordinates', 'FK5Coordinates', 'FK4Coordinates',
           'GalacticCoordinates', 'HorizontalCoordinates'
          ]

#<-----------------Coordinate definitions; transforms are below---------------->


class ICRSCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in the ICRS system.

    If you're looking for "J2000" coordinates, this is probably what you want;
    ICRS is better defined and is within a few microarcsec of J2000.


    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to ICRSCoordinates and used as this
    coordinate.

    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].convert_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(ICRSCoordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec


class FK5Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK5 system.

    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK5Coordinates` and used as this
    coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        super(FK5Coordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].convert_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK5Coordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec


class FK4Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK4 system.

    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK4Coordinates` and used as this
    coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        super(FK4Coordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].convert_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK4Coordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec


class GalacticCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in Galactic Coordinates.

    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `GalacticCoordinates` and
    used as this coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='l', latnm='b'))
    def __init__(self, *args, **kwargs):
        super(GalacticCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].convert_to(self.__class__)
            self.l = newcoord.l
            self.b = newcoord.b
            self._distance = newcoord._distance
        else:
            super(GalacticCoordinates, self)._initialize_latlong('l', 'b', False, args, kwargs)

    @property
    def longangle(self):
        return self.l

    @property
    def latangle(self):
        return self.b


class HorizontalCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in the Horizontal or "alt/az" system.

    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `HorizontalCoordinates` and used
    as this coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='az', latnm='alt'))
    def __init__(self, *args, **kwargs):
        super(HorizontalCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].convert_to(self.__class__)
            self.az = newcoord.az
            self.alt = newcoord.alt
            self._distance = newcoord._distance
        else:
            super(HorizontalCoordinates, self)._initialize_latlong('az', 'alt', False, args, kwargs)

    @property
    def longangle(self):
        return self.az

    @property
    def latangle(self):
        return self.alt

#<-----------------------------------Transforms-------------------------------->
