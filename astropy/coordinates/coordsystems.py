# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the base classes and framework for coordinate objects.
"""

from abc import ABCMeta, abstractproperty, abstractmethod

from .angles import RA, Dec, Angle, AngularSeparation
from .. import units as u

__all__ = ['SphericalCoordinatesBase', 'Coordinates', 'Distance',
           'CartesianPoint'
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
        self._distance = None
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
        `~astropy.coordinates.coordsystems.Distance` object.

        If set as a tuple, the tuple will be passed into the
        `~astropy.coordinates.coordsystems.Distance` constructor.

        Alternatively, this may be `None`, indicating an unknown/not given
        distance. Where necessary, this object will be interpreted as angles on
        the unit sphere.
        """
        return self._distance

    @distance.setter
    def distance(self, val):
        if val is None:
            self._distance = None
        elif isinstance(val, tuple):
            self._distance = Distance(*val)
        elif isinstance(val, Distance):
            self._distance = val
        else:
            raise TypeError('Spherical coordinate distance must be a ')

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

            if self._distanceq is None:
                r = 1
                runit = None
            else:
                r = self._distance.value
                runit = self._distance.unit
            x, y, z = spherical_to_cartesian(r, self.latangle, self.longangle)
            self._cartpoint = CartesianPoint(x, y, z, runit)

    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        Parameters
        ----------
        other : `~astropy.coordinates.coordsystems.SphericalCoordinatesBase`
            The coordinate system to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.Angles.AngularSeparation`
            The on-sky separation between this and the `other` coordinate.
        """
        lat1 = self.latangle.radians
        long1 = other.latangle.radians
        lat2 = self.longangle.radians
        long2 = other.longangle.radians
        return AngularSeparation(lat1, long1, lat2, long2, u.radian)

#FIXME: make this subclass Quantity once Quantity is in master
class Distance(object):
    """
    A one-dimensional distance.

    Parameters
    ----------
    value : scalar
        The value of this distance
    unit : `~astropy.units.core.UnitBase`
        The units for this distance.  Must have dimensions of distance.


    Raises
    ------
    ValueError
        If the `unit` is not a distance.
    """
    def __init__(self, value, unit):
        if not unit.is_equivalent(u.m):
            raise ValueError('provided unit for Distance is not a length')
        self._value = value
        self._unit = unit

    @property
    def lightyear(self):
        """
        The value of this distance in light years
        """
        return self._unit.to(u.lyr, self._value)

    @property
    def pc(self):
        """
        The value of this distance in parsecs
        """
        return self._unit.to(u.parsec, self._value)

    @property
    def kpc(self):
        """
        The value of this distance in kiloparsecs
        """
        return self._unit.to(u.kpc, self._value)

    @property
    def Mpc(self):
        """
        The value of this distance in megaparsecs
        """
        return self._unit.to(u.Mpc, self._value)

    @property
    def au(self):
        """
        The value of this distance in astronomical units
        """
        return self._unit.to(u.au, self._value)

    @property
    def m(self):
        """
        The value of this distance in meters
        """
        return self._unit.to(u.m, self._value)

    @property
    def km(self):
        """
        The value of this distance in kilometers
        """
        return self._unit.to(u.km, self._value)

    @property
    def z(self):
        """
        The redshift for this distance assuming its physical distance is
        a luminosity distance.

        .. note::
            This uses the "current" cosmology to determine the appropriate
            distance to redshift conversions.  See `astropy.cosmology`
            for details on how to change this.

        """
        return self.compute_z()

    def compute_z(self, cosmology=None):
        """
        The redshift for this distance assuming its physical distance is
        a luminosity distance.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.cosmology` or None
            The cosmology to assume for this calculation, or None to use the
            current cosmology.

        """
        from ..cosmology import luminosity_distance
        from scipy import optimize

        #FIXME: array: need to make this calculation more vector-friendly

        f = lambda z, d, cos: (luminosity_distance(z, cos) - d)**2
        return optimize.brent(f, (self.Mpc, cosmology))


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
        Converts to the spherical representation of this point.

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
