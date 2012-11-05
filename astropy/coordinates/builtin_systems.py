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

    .. note::


    {params}
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()
        super(ICRSCoordinates, self)._initialize_latlong('ra', 'dec', args,
                                                         kwargs)

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
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()
        super(ICRSCoordinates, self)._initialize_latlong('ra', 'dec', args,
                                                         kwargs)

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
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()
        super(ICRSCoordinates, self)._initialize_latlong('ra', 'dec', args,
                                                         kwargs)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec


class GalacticCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in Galactic Coordinates

    {params}
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='l', latnm='b'))
    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()
        super(ICRSCoordinates, self)._initialize_latlong('l', 'b', args,
                                                         kwargs)

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
    """.format(params=SphericalCoordinatesBase._init_docstring_templ.format(longnm='az', latnm='alt'))
    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()
        super(ICRSCoordinates, self)._initialize_latlong('az', 'alt', args,
                                                         kwargs)

    @property
    def longangle(self):
        return self.az

    @property
    def latangle(self):
        return self.alt

#<-----------------------------------Transforms-------------------------------->
