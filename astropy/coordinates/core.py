# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the fundamental classes used for representing
coordinates in astropy.
"""

import math
from types import *
from abc import ABCMeta
#from abc import abstractmethod
from abc import abstractproperty

import numpy as np

import conversions as convert
from errors import *
from .. import units as u

__all__ = ['Angle', 'RA', 'Dec', 'Coordinates', 'ICRSCoordinates', 'GalacticCoordinates', 'HorizontalCoordinates']

twopi = math.pi * 2.0 # no need to calculate this all the time

class Angle(object):
    """ This class represents an Angle. 
        
        Units must be specified by the units parameter.
        Degrees and hours both accept either a string like '15:23:14.231,' or a
        decimal representation of the value, e.g. 15.387.
        
        Parameters
        ----------
        angle : float, int, str
            The angle value
        units : {'degrees', 'radians', 'hours'}

    """
    
    def __init__(self, angle, unit=None, bounds=[-360,360]):
        
        self._bounds = bounds
        
        if isinstance(angle, type(self)):
            angle = angle.radians
            unit = u.radian
        
        # short circuit arrays for now
        if isinstance(angle, list):
            raise TypeError("Angles as lists are not yet supported.")
        
        self.isArray = type(angle) in [list, np.ndarray]
        
        if angle == None:
            raise ValueError("The Angle class requires a unit")

        #angle_type = type(angle[0]) if self.isArray else type(angle)
        
        # -------------------------------
        # unit validation and angle value
        # -------------------------------
        if isinstance(unit, u.Unit):
            pass
        elif isinstance(unit, str):
            unit = u.Unit(unit)
        elif unit == None:
            # try to determine unit from the "angle" value
            if self.isArray:
                try:
                    angle = [x.lower() for x in angle]
                except AttributeError:
                        # If units are not specified as a parameter, the only chance
                        # to determine them is in a string value - if it's not a string,
                        # then there's not enough information to create an Angle.
                        raise ValueError("Could not parse an angle value in the array provided"
                                         "- units could not be determined.".format(angle[idx]))
                for idx, a in enumerate(angle):
                    a_unit = None
                    # order is important here - longest name first
                    for unitStr in ["degrees", "degree", "deg", "°"]:
                        if unitStr in a:
                            a_unit = u.radian
                            a = angle.replace(unitStr, "")
                            angle[idx] = math.radians(convert.parseDegrees(a))
                            break
                    if unit == None:
                        for unitStr in ["hours", "hour", "hr"]:
                            if unitStr in a:
                                a_unit = u.radian
                                a = angle.replace(unitStr, "")
                                angle[idx] = math.radians(convert.parseHours(a)*15.)
                                break
                    if unit == None:
                        for unitStr in ["radians", "radian", "rad"]:
                            if unitStr in angle:
                                a_unit = u.radian
                                a = angle.replace(unitStr, "")
                                angle[idx] = convert.parseRadians(a)
                                break
                    if a_unit == None:
                        raise ValueError("Could not parse the angle value '{0}' "
                                         "- units could not be determined.".format(angle[idx]))
                unit = u.radian
                
            else: # single value
                if isinstance(angle, str):
                    angle = angle.lower()
                else:
                    raise ValueError("Could not parse the angle value '{0}' "
                                     "- units could not be determined.".format(angle))
                angle = angle.strip()
                for unitStr in ["degrees", "degree", "deg"]:
                    if unitStr in angle:
                        unit = u.degree
                        angle = angle.replace(unitStr, "")
                if unit == None:
                    for unitStr in ["hours", "hour", "hr"]:
                        if unitStr in angle:
                            unit = u.hour
                            angle = angle.replace(unitStr, "")
                if unit == None:
                    for unitStr in ["radians", "radian", "rad"]:
                        if unitStr in angle:
                            unit = u.radian
                            angle = angle.replace(unitStr, "")
                if unit == None:
                    if "h" in angle:
                        unit = u.hour
                    elif "d" in angle or "°" in angle:
                        unit = u.degree

        if unit == None:
            raise ValueError("The unit parameter should be an object from the "
                             "astropy.unit module (e.g. 'from astropy import units as u',"
                             "then use 'u.degree').")
        
        if self.isArray:
            pass # already performed conversions to radians above
        else:
            if unit == u.degree:
                self._radians = math.radians(convert.parseDegrees(angle))
            elif unit == u.radian:
                self._radians = float(angle)
            elif unit == u.hour:
                self._radians = convert.hoursToRadians(convert.parseHours(angle))
            else:
                raise IllegalUnitsError("The unit value provided was not one of u.degree, u.hour, u.radian'.")

        # ---------------
        # bounds checking
        # ---------------
        # TODO: handle arrays
        # handle bounds units, convert to radians
        if bounds == None:
           pass # no range checking performed
        else:
            try:
                if unit == u.radian:
                    lower_bound = bounds[0]
                    upper_bound = bounds[1]
                elif unit == u.degree:
                    lower_bound = math.radians(bounds[0])
                    upper_bound = math.radians(bounds[1])
                elif unit == u.hour:
                    lower_bound = math.radians(bounds[0]*15.)
                    upper_bound = math.radians(bounds[1]*15.)
                # invalid units handled above
            except TypeError:
                raise TypeError("Bounds specified for Angle must be a two element list, "
                                "e.g. [0,360] (was given '{0}').".format(type(bounds).__name__))
            
            # bounds check
            if lower_bound < self._radians < upper_bound:
                pass # looks good
            else:
                if self._radians > upper_bound:
                    while True:
                        self._radians -= twopi
                        if self._radians < lower_bound:
                            raise RangeError("The angle given falls outside of the specified bounds.")
                        elif lower_bound < self._radians < upper_bound:
                            break
                
                if self._radians < lower_bound:
                    while True:
                        self._radians += twopi
                        if self._radians > upper_bound:
                            raise RangeError("The angle given falls outside of the specified bounds.")
                        elif lower_bound < self._radians < upper_bound:
                            break
    
    @property
    def bounds(self):
        """" Returns the angle's bounds, an immutable property. """
        return self._bounds
    
    @property
    def degrees(self):
        """ Returns the angle's value in degrees (read-only property). """
        return math.degrees(self.radians) # converts radians to degrees

    @property
    def radians(self):
        """ Returns the angle's value in radians (read-only property). """
        return self._radians

    @property
    def hours(self):
        """ Returns the angle's value in hours (read-only property). """
        return convert.radiansToHours(self.radians)

    @property
    def hms(self):
        """ Returns the angle's value in hours, and print as an (h,m,s) tuple (read-only property). """
        return convert.radiansToHMS(self.radians)

    @property
    def dms(self):
        """ Returns the angle's value in degrees, and print as an (d,m,s) tuple (read-only property). """
        return convert.radiansToDMS(self.radians)

    def string(self, unit=u.degree, decimal=False, sep=" ", precision=5, pad=False): 
        """ Returns a string representation of the angle.
        
            Parameters
            ----------
            units : str
                Specifies the units, value should be one of the allowed units values (see: `Angle`)
            decimal : bool
                Specifies whether to return a sexagesimal representation (e.g. a tuple 
                (hours, minutes, seconds)), or decimal
            sep : str
                The separator between numbers in a sexagesimal representation, e.g. 12:41:11.1241
                where the separator is ":". Also accepts 2 or 3 separators, e.g. 12h41m11.1241s would be sep="hms",
                or 11-21:17.124 would be sep="-:"
        """
        
        if isinstance(unit, u.Unit):
            pass # great!
        elif isinstance(unit, str):
            unit = unit.lower()
            if unit == "degrees":
                unit = u.degree
            elif unit == "hours":
                unit = u.hour
            elif unit == "radians":
                unit = u.radian
            else:
                raise IllegalUnitsError("The unit value provided was not one "
                                        "of u.degree, u.hour, u.radian'.")
        else:
                raise IllegalUnitsError("The unit value provided was not one "
                                        "of u.degree, u.hour, u.radian'.")
        
        if unit == u.degree:
            if decimal:
                return ("{0:0." + str(precision) + "f}").format(self.degrees)
            else:
                return convert.degreesToString(self.degrees, precision=precision, sep=sep, pad=pad)
                
        elif unit == u.radian:
            return str(self.radians)
                
        elif unit == u.hour:
            if decimal:
                return ("{0:0." + str(precision) + "f}").format(self.hours)
            else:
                return convert.hoursToString(self.hours, precision=precision, sep=sep, pad=pad)
        else:
            raise IllegalUnitsError(unit)

    # ----------------------------------------------------------------------------
    # Emulating numeric types
    # -----------------------
    # Ref: http://docs.python.org/reference/datamodel.html#emulating-numeric-types
    # ----------------------------------------------------------------------------

    # Addition
    def __add__(self, other):
        if isinstance(other, type(self)):
            if self.bounds != other.bounds:
                raise ValueError("An {0} object can only be subtracted from another "
                                 "{0} object.".format(type(self).__name__))
            else:
                return Angle(self.radians + other.radians, unit=u.radian, bounds=self.bounds)
        else:
            raise NotImplementedError("An {0} object can only be added to another "
                                      "{0} object.".format(type(self).__name__))
    
    # Subtraction
    def __sub__(self, other):
        if isinstance(other, type(self)):
            if self.bounds != other.bounds:
                raise ValueError("An {0} object can only be subtracted from another "
                                 "{0} object.".format(type(self).__name__))
            else:
                return Angle(self.radians - other.radians, unit=u.radian, bounds=self.bounds)
        else:
            raise NotImplementedError("An {0} object can only be subtracted from another "
                                      "{0} object.".format(type(self).__name__))
    
    # Multiplication
    def __mul__(self, other):
        if isinstance(other, type(self)):
            raise NotImplementedError("Multiplication is not supported between two {0} "
                                      "objects ".format(type(self).__name__))
        elif type(other) in [float, int]:
            return Angle(self.radians*other, unit=u.radian) 
        else:
            raise NotImplementedError("An {0} object can only be multiplied by a float or integer.".format(type(self).__name__))
    
    # Division
    def __div__(self, other):
        if isinstance(other, type(self)):
            raise NotImplementedError("Division is not supported between two {0} "
                                      "objects.".format(type(self).__name__))
        elif type(other) in [float, int]:
            return Angle(self.radians/other, unit=u.radian)
        else:
            raise NotImplementedError("An {0} object can only be divided by a float or integer.".format(type(self).__name__))
    
    def __truediv__(self, other):
        raise NotImplementedError("Division is not supported between two {0} "
                                  "objects.".format(type(self).__name__))

    def __neg__(self):
        return Angle(-self.radians, unit=u.radian)

    def __eq__(self, other):
        # Ref: http://stackoverflow.com/questions/3049101/floating-point-equality-in-python-and-in-general
        #return self.radians == other.radians
        #return abs(self.radians - other.radians) < 1e-11
        #
        if isinstance(other, type(self)):
            # abs(x - y) <= nulps * spacing(max(abs(x), abs(y)))
            #print(abs(self.radians - other.radians), np.spacing(max(abs(self.radians), abs(other.radians))))
            return abs(self.radians - other.radians) <= 4 * np.spacing(max(abs(self.radians), abs(other.radians)))
        elif (other == None):
            return False
        else:
             raise NotImplementedError("An {0} object can only be compared to another {0} "
                                       "object.".format(type(self).__name__))
        
    def __ne__(self, other):
        #return not self.radians == other.radians
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, type(self)):
            return self.radians < other.radians
        else:
             raise NotImplementedError("An {0} object can only be compared to another {0} "
                                       "object.".format(type(self).__name__))
    
    def __gt__(self, other):
        if isinstance(other, type(self)):
            return self.radians > other.radians
        else:
             raise NotImplementedError("An {0} object can only be compared to another {0} "
                                       "object.".format(type(self).__name__))
    
    def __ge__(self, other):
        if isinstance(other, type(self)):
            return self.radians >= other.radians
        else:
             raise NotImplementedError("An {0} object can only be compared to another {0} "
                                       "object.".format(type(self).__name__))
        
    def __le__(self, other):
        if isinstance(other, type(self)):
            return self.radians <= other.radians
        else:
             raise NotImplementedError("An {0} object can only be compared to another {0} "
                                       "object.".format(type(self).__name__))
    
    def __abs__(self):
        return Angle(abs(self.radians), unit=u.radian)
    
class RA(Angle):
    """ Represents a J2000 Right Ascension 
    
        Accepts a right ascension angle value. The angle parameter accepts degrees, 
        hours, or radians.
        Degrees and hours both accept either a string like '15:23:14.231,' or a
        decimal representation of the value, e.g. 15.387.
        
        Parameters
        ----------
        angle : float, int
            The angle value
        unit : an instance of astropy.Unit or an equivalent string
            The unit of the angle value
    
    """
    
    def __init__(self, angle, unit=None):
   
        # This block attempts to determine the validity of the unit,
        # particularly with regard to the specifics of RA.
        # After this block, the normal Angle initializer handles most of the
        # validation/creation.
        
        if isinstance(angle, type(self)):
            return super(RA, self).__init__(angle.radians, unit=u.radian, bounds=(0,360))
                    
        if unit == u.hour:
            pass # to Angle initializer
            #self._radians = math.radians(decimal_hours * 15.)
        elif unit == u.degree:
            pass # to Angle initializer
            #decimal_degrees = convert.parseDegrees(angle)
            #self._radians = math.radians(decimal_degrees)
        elif unit == u.radian:
            pass # to Angle initializer
            #self._radians = convert.parseRadians(angle)
        elif unit == None:
            # Try to figure out the unit if we can.
            if isinstance(angle, str):
                # Try to deduce the units from hints in the string.
                # Further, enforce absolute bounds here, i.e. don't let
                # Angle +-2π to see if the angle falls in the bounds.
                if "d" in angle or "°" in angle:
                    # If in the form "12d32m53s", look for the "d" and assume degrees.
                    angle = math.radians(convert.parseDegrees(angle))
                    if 0 < angle < twopi:
                        unit = u.radian
                    else:
                        raise RangeError("The provided angle was assumed to be in degrees, but was out of the range (0,360) degrees.")
                elif "h" in angle:
                    # Same for "12h32m53s" for hours.
                    #self._radians = math.radians(convert.parseHours(angle)*15.0)
                    unit = u.hour
                else:
                    # could be in a form: "54:43:26" -
                    # if so AND the resulting decimal value is > 24 or < -24, assume degrees
                    decimal_value = convert.parseDegrees(angle)
                    if decimal_value > 24:
                        unit = u.degree
                    elif 0 <= decimal_value <= 24.0:
                        raise ValueError("No units were specified, and the angle value was ambiguous between hours and degrees.")
                    elif decimal_value < 0:
                        raise RangeError("No units were specified; could not assume any since the value was less than zero.")
            elif isinstance(angle, tuple):
                if len(angle) == 3 and 0 <= angle[0] < 24.0:
                    raise ValueError("No units were specified, and the angle value was ambiguous between hours and degrees.")
                else:
                    unit = u.degree
            else:
                raise ValueError("Angle values of type {0} not supported.".format(type(angle).__name__))
        if unit == None:
            raise ValueError("Units must be specified for RA, one of u.degree, u.hour, or u.radian.")
        
        # By here, the unit should be defined.
        super(RA, self).__init__(angle, unit=unit, bounds=(0,360))
        
    def hourAngle(self, lst, unit=None):
        """ Given a Local Sidereal Time (LST), calculate the hour angle for this RA
        
            Parameters
            ----------
            lst : float, str, `Angle`
                A Local Sidereal Time (LST)
            unit : str
                The units of the LST, if not an `Angle` object or datetime.datetime object
                .. note::
                    * if lst is **not** an `Angle`-like object, you can specify the units by passing a `unit` parameter into the call
                    * this function currently returns an `Angle` object
        """
        # TODO : this should return an HA() object, and accept an Angle or LST object
        if not isinstance(lst, Angle):
            lst = Angle(lst, units)
        
        return Angle(lst.radians - self.radians, unit=u.radian)
    
    def lst(self, hourAngle, unit=u.hour):
        """ Given an Hour Angle, calculate the Local Sidereal Time (LST) for this RA
    
            Parameters
            ----------
            ha :  float, str, `Angle`
                An Hour Angle
            units : str
                The units of the ha, if not an `Angle` object
            
            .. note:: if ha is *not* an Angle object, you can specify the units by passing a 'units' parameter into the call
                'units' can be radians, degrees, or hours
                this function always returns an Angle object
        """
        # TODO : I guess this should return an HA() object, and accept an Angle or LST object
        if not isinstance(ha, Angle):
            ha = Angle(ha, unit)
            
        return Angle(ha.radians + self.radians, units=u.radian)

class Dec(Angle):
    """ Represents a J2000 Declination """
    
    def __init__(self, angle, unit=u.degree):
        """ Accepts a Declination angle value. The angle parameter accepts degrees, 
            hours, or radians and the default units are hours.
            Degrees and hours both accept either a string like '15:23:14.231,' or a
            decimal representation of the value, e.g. 15.387.
            Bounds are fixed at [-90,90]
        
            Parameters
            ----------
            angle : float, int
                The angular value
            units : str 
                The units of the specified declination
            
        """
        if isinstance(angle, type(self)):
            return super(Dec, self).__init__(angle.radians, unit=u.radian, bounds=(-90,90))


        super(Dec, self).__init__(angle, unit=unit, bounds=(-90,90))

class CoordinatesBase(object):
    """
    Abstract superclass for all coordinate classes (except the factory class 'Coordinates').
    """
    
    __metaclass__ = ABCMeta
    
    @abstractproperty
    def angle1(self):
        pass
    
    @abstractproperty
    def angle2(self):
        pass

    def __eq__(self, other):
        if isinstance(other, type(self)):
            angle1_eq = abs(self.angle1.radians - other.angle1.radians) <= 4 * np.spacing(max(abs(self.angle1.radians), abs(other.angle1.radians)))
            angle2_eq = abs(self.angle2.radians - other.angle2.radians) <= 4 * np.spacing(max(abs(self.angle2.radians), abs(other.angle2.radians)))
            return angle1_eq and angle2_eq
        #elif isinstance(other, type(GalacticCoordinates)):
            ## TODO: support comparisons to other coordinate systems.
        #elif isinstance(other, type(HorizontalCoordinates)):
            ## TODO
        elif (other == None):
            return False
        else:
            raise NotImplementedError("An {0} object can only be compared to another {0} "
                                      "object.".format(type(self).__name__))
    def __ne__(self, other):
        return not self.__eq__(other)

class ICRSCoordinates(CoordinatesBase):
    """
    RA/Dec coordinate class.
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
        
class GalacticCoordinates(CoordinatesBase):
    """ 
    Galactic coordinate (l,b) class.
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
        if len(args) > 0:
            # make sure someone isn't using RA/Dec objects
            for arg in args:
                if isinstance(arg, RA) or isinstance(arg, Dec):
                    raise TypeError("The class {0} doesn't accept RA or Dec values; "
                                     "use Angle objects instead.".format(type(self).__name__))
            if len(args) == 0 and len(kwargs) > 0:
                # only "l" and "b" accepted as keyword arguments
                try:
                    _l = kwargs["l"]
                    _b = kwargs["b"]
                except KeyError:
                    raise ValueError("When values are supplied as keyword arguments, both "
                                     "'l' and 'b' must be specified.")
                if isinstance(_ra, Angle):
                    self.l = _l
                if isinstance(_dec, Angle):
                    self.b = _b
    
            elif len(args) == 1 and len(kwargs) == 0:
                # need to try to parge the coordinate from a single argument
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
                                print " ".join(x.values[0:i])
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

class Coordinates(object):
    """
    Document me.
    
    This is a factory class that will instantiate the appropriate CoordinateBase subclass.
    A "Coordinates" object cannot be created on its own.
    """
    __meta__ = ABCMeta
    
    def __new__(self, *args, **kwargs):
        # coordinates, units=None, ra=None, dec=None, az=None, el=None, l=None, b=None):
        """
        Document me.
        """
        units = kwargs["units"] if "units" in kwargs.keys() else list()
        
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
