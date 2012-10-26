# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the base classes and framework for coordinate objects.
"""

from abc import ABCMeta, abstractproperty, abstractmethod

from .angles import RA, Dec, Angle
from .. import units as u

__all__ = ['SphericalCoordinatesBase', 'Coordinates', 'CartesianPoint'
          ]

class SphericalCoordinatesBase(object):
    """
    Abstract superclass for all coordinate classes representing points in three
    dimensions.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, **kwargs):
        """
        Subclasses must override this, but they should also call this to set up
        internal state.
        """
        self._distancequant = None
        self._cartpoint = None

    @abstractproperty
    def latangle(self):
        """
        The latitudinal/elevation angle for these coordinates as an
        `~astropy.coorinates.angles.Angle` object.

        .. note ::
            This should be overridden in subclasses as a read-only property that
            just returns an attribute a way to abstract the exact choice of
            names for the coordiantes. E.g., `ICRSCoordinates` implements this
            by doing ``return self.ra``.
        """
        pass

    @abstractproperty
    def longangle(self):
        """
        The longitudinal/azimuthal angle for these coordinates as an
        `~astropy.coorinates.angles.Angle` object.

        .. note ::
            This should be overridden in subclasses as a read-only property that
            just returns an attribute a way to abstract the exact choice of
            names for the coordinates. E.g., `ICRSCoordinates` implements this
            by doing ``return self.dec``.
        """

    @property
    def distance(self):
        """
        The radial distance for this coordinate object as an
        `~astropy.units.Quantity` object. It must have units of distance.

        Alternatively, this may be `None`, indicating an unknown/not given
        distance. Where necessary, this will be interpreted as a (dimensionless)
        unit sphere.
        """
        return self._distancequant

    @distance.setter
    def distance(self, val):
        if val is None:
            self._distance = None
        elif not hasattr(val, 'value') or not hasattr(val, 'unit'):
            raise TypeError('Spherical coordinate distances is not a Quantity-like object')
        elif not val.unit.is_equivalent(u.m):
            raise u.IncompatibleUnitError('')

        else:
            self._distance = u.Quantity(val)

    @property
    def x(self):
        self._make_cart()
        return self._cartpoint.x

    @property
    def y(self):
        self._make_cart()
        return self._cartpoint.y

    @property
    def z(self):
        self._make_cart()
        return self._cartpoint.z

    @property
    def cartesian(self):
        self._make_cart()
        return self._cartpoint

    def _make_cart(self, override=False):
        if override or self._cartpoint is None:
            from .transformations import spherical_to_cartesian

            if self._distancequant is None:
                r = 1
                runit = None
            else:
                r = self._distancequant.value
                runit = self._distancequant.unit
            x, y, z = spherical_to_cartesian(r, self.latangle, self.longangle)
            self._cartpoint = CartesianPoint(x, y, z, runit)

class CartesianPoint(object):
    """
    A cartesian representation of a point in three-dimensional space.

    Attributes
    ----------
    x : number or array
        The first cartesian coordinate.
    y : number or array
        The second cartesian coordinate.
    z : number or array
        The third cartesian coordinate.
    unit : `~astropy.units.UnitBase` object or None
        The physical unit of the coordinate values.
    """

    def __init__(self, x, y, z, unit=None):
        self.x = x
        self.y = y
        self.z = z
        self.unit = unit

    def to_spherical(self):
        """
        Converts returns the spherical representation of this point.

        Returns
        -------
        r : float or array
            The radial coordinate (in the same units as the inputs).
        lat : float or array
            The latitude in radians
        lng : float or array
            The longitude in radians

        """
        from .transformations import cartesian_to_spherical

        return cartesian_to_spherical(self.x, self.y, self.z)

class Coordinates(object):
    """
    A convenience factory class to create coordinate objects.

    This class can be used to create coordinate objects. The coordinate system is chosen
    based on the keywords used. For example, using the 'l' and 'b' keywords will return
    a `~astropy.coordinates.GalacticCoordinates` object. A "Coordinates" object cannot be
    created on its own.

    Parameters
    ----------
    (ra, dec) : `~astropy.coordinates.Angle`, str, float, int
        Right ascension and declination values. Returns an ICRSCoordinate object.
    (l, b) : `~astropy.coordinates.Angle`, str, float, int
        Galactic latitude and longitude. Returns a GalacticCoordinates object.
    (az, el) : `~astropy.coordinates.Angle`, str, float, int
        Azimuth and elevation values. Returns a HorizontaolCoordinates object.

    unit : `~astropy.units.Unit`, str, tuple
        Units must be provided for each of the angles provided. If the unit value can be
        determined from the Angle objects directly, `unit` does not need to be specified.
        If `unit` is a single value, it is applied to both of the given angles. If one angle
        requires a unit and the other does not, use `None` as a placeholder.
    """
    __meta__ = ABCMeta

    def __new__(self, *args, **kwargs):
        # coordinates, units=None, ra=None, dec=None, az=None, el=None, l=None, b=None):
        """
        Document me.
        """
        #units = kwargs["unit"] if "unit" in kwargs.keys() else list()
        try:
            units = kwargs["unit"]
            if isinstance(units, u.Unit) or isinstance(units, str):
                units = (units, units)
        except KeyError:
            units = list()

        # first see if the keywords suggest what kind of coordinate is being requested.
        if "ra" in kwargs.keys() or "dec" in kwargs.keys():
            try:
                ra = kwargs["ra"]
                dec = kwargs["dec"]
            except KeyError:
                raise ValueError("When an 'ra' or 'dec' value is provided, the "
                                 "other coordinate must also be given.")
            for kw in ["l", "b", "az", "el"]: # no others should be provided
                if kw in kwargs.keys():
                    raise ValueError("Conflicting coordinates were given.")
            ra = RA(ra, unit=units[0]) if len(units) > 0 else RA(ra)
            dec = Dec(dec, unit=units[1]) if len(units) > 1 else Dec(dec)
            return ICRSCoordinates(ra=ra, dec=dec)

        if "az" in kwargs.keys() or "el" in kwargs.keys():
            try:
                az = kwargs["az"]
                el = kwargs["el"]
            except KeyError:
                raise ValueError("When an 'az' or 'el' horizontal coordinates value "
                                 "is provided, the other coordinate must also be given.")
            for kw in ["ra", "dec", "l", "b"]: # no others should be provided
                if kw in kwargs.keys():
                    raise ValueError("Conflicting coordinates were given.")
            az = Angle(az, unit=units[0]) if len(units) > 0 else Angle(az)
            el = Angle(el, unit=units[1]) if len(units) > 1 else Angle(el)
            return HorizontalCoordinates(az=az, el=el)

        if "l" in kwargs.keys() or "b" in kwargs.keys():
            try:
                l = kwargs["l"]
                b = kwargs["b"]
            except KeyError:
                raise ValueError("When an 'l' or 'b' galactic coordinates value is "
                                 "provided, the other coordinate must also be given.")
            for kw in ["ra", "dec", "az", "el"]: # no others should be provided
                if kw in kwargs.keys():
                    raise ValueError("Conflicting coordinates were given.")
            l = Angle(l, unit=units[0]) if len(units) > 0 else Angle(l)
            b = Angle(b, unit=units[1]) if len(units) > 1 else Angle(b)
            return GalacticCoordinates(l=l, b=b)

        if len(args) == 1:
            x = args[0]

            if isinstance(x, str):
                raise ValueError("The coordinate system could not be determines from the value "
                                 "provided. Specify the system via keywords or use the "
                                 "corresponding class (e.g. GalacticCoordinate).")
            elif isinstance(x, list):
                return ValueError("Lists of coordinates are not yet supported")
            else:
                return ValueError("Could not create a Coordinate object from an object "
                                  "of type '{0}'.".format(type(x).__name__))
        if len(args) == 2:
            #a1, a2 = args[0:2]
            if isinstance(args[0], RA) and isinstance(args[1], Dec):
                return ICRSCoordinates(ra=args[0], dec=args[1])
            raise ValueError("Two angles were provided ('{0[0]}', '{0[1]}'), but the "
                             "coordinate system "
                             "was not provided. Specify the system via keywords or use the "
                             "corresponding class (e.g. GalacticCoordinate).".format(args))

        else:
            raise ValueError("Could not construct coordinates.")

        if False: # old code - still useful?
            # determine units
            if units is not None:
                if isinstance(units, u.Unit):
                    pass # great!
                elif isinstance(units, tuple):
                    if len(units) == 0 or len(units) > 2:
                        raise ValueError("The units parameter only accepts "
                                         "tuples with one or two values.")
                    else:
                        # validate them
                        for a in units:
                            if not isinstance(a, u.Unit):
                                raise ValueError("Units must be specified as u.degree, u.hour, etc.")
                    # units are valid
            else:
                #units were None - try to determine units from coordinate object
                if isinstance(coordinates, tuple):
                    if len(coordinates) != 2:
                        raise ValueError("Two coordinate values must be provided - '{0}' found.".format(len(coordinates)))
                    else:
                        # we have two values - the goal is to end up with two Angle
                        # objects.
                        pass
