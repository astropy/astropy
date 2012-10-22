# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the implementations of specific coordinate systems
and the conversions between them.
"""

from abc import ABCMeta

from .core import CoordinatesBase, RA, Dec, Angle
from .. import units as u

__all__ = ['Coordinates', 'ICRSCoordinates', 'GalacticCoordinates', 
           'HorizontalCoordinates']


class ICRSCoordinates(CoordinatesBase):
    """
    Object representing a coordinate in the ICRS system.
    
    Parameters
    ----------
    ra : `~astropy.coordinates.angle`, float, int, str
    dec : `~astropy.coordinates.angle`, float, int, str
    unit : tuple
        If the units cannot be determined from the angle values provided, they must
        be specified as a tuple in the 'unit' parameter. The first value in the tuple
        is paired with the first angle provided, and the second with the second angle.
        (If the unit is specified in the first angle but not the second, the first
        value in the tuple may be 'None'.)
    """
    def __init__(self, *args, **kwargs):
        
        # Initialize values.
        # _ra, _dec are what we parse as potential values that still need validation
        _ra = None
        _dec = None
        self.ra = None
        self.dec = None
        
        if "unit" in kwargs:
            units = kwargs["unit"]
            del kwargs["unit"]
        else:
             units = list()

        if isinstance(units, tuple) or isinstance(units, list):
            pass # good
        elif isinstance(units, u.Unit) or isinstance(units, str):
            # Only a single unit given, which is fine (assigned to 'ra').
            # The value, even if given as a tuple, is unpacked. Just make it
            # a tuple for consistency
            units = [units]
        else:
            raise ValueError("The value for units must be given as a tuple, e.g. "
                             "unit=(u.hour, u.degree). An object of type '{0}' "
                             "was given.".format(type(units).__name__))

        
        if len(args) == 0 and len(kwargs) == 0:
            raise ValueError("A coordinate object cannot be created without ra,dec values.")
        elif len(args) > 0 and len(kwargs) > 0:
            raise ValueError("The angle values can only be specified as keyword arguments "
                             "(e.g. ra=x, dec=y) or as a single value (e.g. a string) "
                             "not a combination.")
        elif len(args) == 0 and len(kwargs) > 0:
            # only "ra" and "dec" accepted as keyword arguments
            try:
                _ra = kwargs["ra"]
                _dec = kwargs["dec"]
            except KeyError:
                raise ValueError("When values are supplied as keyword arguments, both "
                                 "'ra' and 'dec' must be specified.")
            if isinstance(_ra, RA):
                self.ra = _ra
            if isinstance(_dec, Dec):
                self.dec = _dec

        elif len(args) == 1 and len(kwargs) == 0:
            # need to try to parge the coordinate from a single argument
            x = args[0]
            if isinstance(args[0], str):
                parsed = False
                if "," in x:
                    _ra, _dec = split(",")
                    parsed = True
                elif "\t" in x:
                    _ra, _dec = split("\t")
                    parsed = True
                elif len(x.split()) == 6:
                    _ra = " ".join(x.split()[0:3])
                    _dec = " ".join(x.split()[3:])
                    parsed = True
                elif len(x.split()) == 2:
                    _ra, _dec = x.split()
                    parsed = True
                
                if not parsed:
                    values = x.split()
                    i = 1
                    while i < len(values) and not parsed:
                        try:
                            self.ra = RA(" ".join(values[0:i]))
                            parsed = True
                        except:
                            i += 1
                    
                    if parsed == True:
                        self.dec = Dec(" ".join(values[i:]))
                
                if not parsed:
                    raise ValueError("Could not parse ra,dec values from the string provided: '{0}'.".format(x))
            else:
                raise ValueError("A coordinate cannot be created with a value of type "
                                 "'{0}'.".format(type(arg[0]).__name___))

        elif len(args) == 2 and len(kwargs) == 0:
            _ra = args[0]
            _dec = args[1]

        elif len(args) > 2 and len(kwargs) == 0:
            raise ValueError("More than two values were found where only ra and dec "
                             "were expected.")
        else:
            raise ValueError("Unable to create a coordinate using the values provided.")


#             # First try to see if RA, Dec objects were provided in the args.            
#             for arg in args:
#                 if isinstance(arg, RA):
#                     _ra = arg
#                 elif isinstance(arg, Dec):
#                     _dec = arg
#             
#             if None not in [_ra, _dec]:
#                 self.ra = _ra
#                 self.dec = _dec
#                 return
#             elif (_ra and not _dec) or (not _ra and _dec):
#                 raise ValueError("When an RA or Dec value is provided, the other "
#                                  "coordinate must also be given.")
# 
#             # see if the whole coordinate might be parseable from arg[0]
#         
#         try:
#             if isinstance(args[0], RA) and isinstance(args[1], Dec):
#                 _ra = args[0]
#                 _dec = args[1]
#             elif isinstance(args[1], RA) and isinstance(args[0], Dec):
#                 _ra = args[1]
#                 _dec = args[0]
#         except IndexError:
#             raise ValueError("Not enough parameters were provided.")
        
        if self.ra is None:
            self.ra = RA(_ra, unit=units[0]) if len(units) > 0 else RA(_ra)
        if self.dec is None:
            self.dec = Dec(_dec, unit=units[1]) if len(units) > 1 else Dec(_dec)

    @property
    def angle1(self):
        return self.ra
    
    @property
    def angle2(self):
        return self.dec

    @property
    def icrs(self):
        return self

    @property
    def galactic(self):
        raise NotImplementedError()
        
    @property
    def horizontal(self):
        raise NotImplementedError()

class GalacticCoordinates(CoordinatesBase):
    """ 
    Galactic coordinate (l,b) class.
    
    Parameters
    ----------
    l : `~astropy.coordinates.angle`, float, int, str
        galactic latitude
    b : `~astropy.coordinates.angle`, float, int, str
        galactic longitude
    unit : tuple
        Units associated with the l and b values provided. Only needed if the unit cannot be
        unambiguously determined from the given values.

    """
    def __init__(self, *args, **kwargs):
        
        # initialize values
        # _ra, _dec are what we parse as potential values that still need validation
        _l = None
        _b = None
        self.l = None
        self.b = None
        
        if "unit" in kwargs:
            units = kwargs["unit"]
            del kwargs["unit"]
        else:
            units = list()
            
        if isinstance(units, tuple) or isinstance(units, list):
            pass # good
        elif isinstance(units, u.Unit) or isinstance(units, str):
            # Only a single unit given, which is fine (assigned to 'ra').
            # The value, even if given as a tuple, is unpacked. Just make it
            # a tuple for consistency
            units = [units]
        else:
            raise ValueError("The value for units must be given as a tuple, e.g. "
                             "unit=(u.hour, u.degree). An object of type '{0}' "
                             "was given.".format(type(units).__name__))

        if len(args) == 0 and len(kwargs) == 0:
            raise ValueError("A coordinate object cannot be created without l,b values.")
        elif len(args) > 0 and len(kwargs) > 0:
            raise ValueError("The angle values can only be specified as keyword arguments "
                             "(e.g. l=x, b=y) or as a single value (e.g. a string) "
                             "not a combination.")
                             
        if len(args) == 0 and len(kwargs) > 0:
            # only "l" and "b" accepted as keyword arguments
            try:
                _l = kwargs["l"]
                _b = kwargs["b"]
            except KeyError:
                raise ValueError("When values are supplied as keyword arguments, both "
                                 "'l' and 'b' must be specified.")
            if isinstance(_l, Angle):
                self.l = _l
            if isinstance(_b, Angle):
                self.b = _b

        elif len(args) > 0:
            # make sure someone isn't using RA/Dec objects
            for arg in args:
                if isinstance(arg, RA) or isinstance(arg, Dec):
                    raise TypeError("The class {0} doesn't accept RA or Dec values; "
                                     "use Angle objects instead.".format(type(self).__name__))
    
            if len(args) == 1 and len(kwargs) == 0:
                # need to try to parse the coordinate from a single argument
                x = args[0]
                if isinstance(args[0], str):
                    parsed = False
                    if "," in x:
                        _l, _b = split(",")
                        parsed = True
                    elif "\t" in x:
                        _l, _b = split("\t")
                        parsed = True
                    elif len(x.split()) == 6:
                        _l = " ".join(x.split()[0:3])
                        _b = " ".join(x.split()[3:])
                        parsed = True
                    elif len(x.split()) == 2:
                        _l, _b = x.split()
                        parsed = True
    
                    if not parsed:
                        values = x.split()
                        i = 1
                        while i < len(values) and not parsed:
                            try:
                                self.l = Angle(" ".join(x.values[0:i]))
                                parsed = True
                            except:
                                i += 1
                        
                        if parsed == True:
                            self.b = Angle(" ".join(x.values[i:]))
                    
                    if not parsed:
                        raise ValueError("Could not parse l,b values from the string provided: '{0}'.".format(x))
                else:
                    raise ValueError("A coordinate cannot be created with a value of type "
                                     "'{0}'.".format(type(arg[0]).__name___))
    
            elif len(args) == 2 and len(kwargs) == 0:
                _l = args[0]
                _b = args[1]
    
            elif len(args) > 2 and len(kwargs) == 0:
                raise ValueError("More than two values were found where only ra and dec "
                                 "were expected.")
        else:
            raise ValueError("Unable to create a coordinate using the values provided.")

        if self.l is None:
            self.l = Angle(_l, unit=units[0]) if len(units) > 0 else Angle(_l)
        if self.b is None:
            self.b = Angle(_b, unit=units[1]) if len(units) > 1 else Angle(_b)

    @property
    def angle1(self):
        return self.l
    
    @property
    def angle2(self):
        return self.b

    @property
    def icrs(self):
        raise NotImplementedError()

    @property
    def galactic(self):
        return self

    @property
    def horizontal(self):
        raise NotImplementedError()

class HorizontalCoordinates(CoordinatesBase):
    """ 
    Horizontal coordinate (az,el) class.
    """
    @property
    def angle1(self):
        return self.az
    
    @property
    def angle2(self):
        return self.el

    @property
    def icrs(self):
        raise NotImplementedError()

    @property
    def galactic(self):
        raise NotImplementedError()
        
    @property
    def horizontal(self):
        return self

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
