# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the fundamental classes used for representing
coordinates in astropy.
"""

import math
from types import *

import numpy as np

import angle_utilities as util
from errors import *
from .. import units as u

__all__ = ['Angle', 'RA', 'Dec', 'AngularSeparation']

twopi = math.pi * 2.0  # no need to calculate this all the time


class Angle(object):
    """ This class represents an angle.

        Units must be specified by the units parameter.
        Degrees and hours both accept either a string like '15:23:14.231,' or a
        decimal representation of the value, e.g. 15.387.

        Parameters
        ----------
        angle : float, int, str
            The angle value.
        unit : `~astropy.units` (preferred), str
            The unit of the value specified for the angle. It is preferred that
            the unit be an object from the `~astropy.units` package, e.g.
            "from astropy import units as u; u.degree". Also accepts any string that the Unit class
            maps to "degrees", "radians", "hours".
        bounds : tuple
            A tuple indicating the upper and lower value that the new angle object may
            have.

    """

    def __init__(self, angle, unit=None, bounds=(-360, 360)):

        self._bounds = bounds

        if isinstance(angle, type(self)):
            angle = angle.radians
            unit = u.radian

        # short circuit arrays for now
        if isinstance(angle, list):
            raise TypeError("Angles as lists are not yet supported.")

        self.is_array = type(angle) in [list, np.ndarray]

        if angle == None:
            raise ValueError("The Angle class requires a unit")

        #angle_type = type(angle[0]) if self.is_array else type(angle)

        # -------------------------------
        # unit validation and angle value
        # -------------------------------
        if isinstance(unit, u.Unit):
            pass
        elif isinstance(unit, str):
            unit = u.Unit(unit)
        elif unit == None:
            # try to determine unit from the "angle" value
            if self.is_array:
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
                            angle[idx] = math.radians(util.parse_degrees(a))
                            break
                    if unit == None:
                        for unitStr in ["hours", "hour", "hr"]:
                            if unitStr in a:
                                a_unit = u.radian
                                a = angle.replace(unitStr, "")
                                angle[idx] = math.radians(util.parse_hours(a) * 15.)
                                break
                    if unit == None:
                        for unitStr in ["radians", "radian", "rad"]:
                            if unitStr in angle:
                                a_unit = u.radian
                                a = angle.replace(unitStr, "")
                                angle[idx] = util.parse_radians(a)
                                break
                    if a_unit == None:
                        raise ValueError("Could not parse the angle value '{0}' "
                                         "- units could not be determined.".format(angle[idx]))
                unit = u.radian

            else:  # single value
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

        if self.is_array:
            pass  # already performed conversions to radians above
        else:
            if unit == u.degree:
                self._radians = math.radians(util.parse_degrees(angle))
            elif unit == u.radian:
                self._radians = float(angle)
            elif unit == u.hour:
                self._radians = util.hours_to_radians(util.parse_hours(angle))
            else:
                raise ValueError("The unit value provided was not one of u.degree, u.hour, u.radian'.")

        # ---------------
        # bounds checking
        # ---------------
        # TODO: handle arrays
        # handle bounds units, convert to radians
        if bounds == None:
            Ppass  # no range checking performed
        else:
            try:
                if unit == u.radian:
                    lower_bound = bounds[0]
                    upper_bound = bounds[1]
                elif unit == u.degree:
                    lower_bound = math.radians(bounds[0])
                    upper_bound = math.radians(bounds[1])
                elif unit == u.hour:
                    lower_bound = math.radians(bounds[0] * 15.)
                    upper_bound = math.radians(bounds[1] * 15.)
                # invalid units handled above
            except TypeError:
                raise TypeError("Bounds specified for Angle must be a two element list, "
                                "e.g. [0,360] (was given '{0}').".format(type(bounds).__name__))

            # bounds check
            if lower_bound < self._radians < upper_bound:
                pass  # looks good
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
        """" The angle's bounds, an immutable property. """
        return self._bounds

    @property
    def degrees(self):
        """ The angle's value in degrees (read-only property). """
        return math.degrees(self.radians)  # converts radians to degrees

    @property
    def radians(self):
        """ The angle's value in radians (read-only property). """
        return self._radians

    @property
    def hours(self):
        """ The angle's value in hours (read-only property). """
        return util.radians_to_hours(self.radians)

    @property
    def hms(self):
        """ The angle's value in hours, and print as an (h,m,s) tuple (read-only property). """
        return util.radians_to_hms(self.radians)

    @property
    def dms(self):
        """ The angle's value in degrees, and print as an (d,m,s) tuple (read-only property). """
        return util.radians_to_dms(self.radians)

    #TODO: Check with @demitri and @adrn re: this vs __str__ vs fromat vs to_string and so on
    def string(self, unit=u.degree, decimal=False, sep=" ", precision=5, pad=False):
        """ A string representation of the angle.

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

            Returns
            -------
            strrepr : str
                A string representation of the angle.

        """

        if isinstance(unit, u.Unit):
            pass  # great!
        elif isinstance(unit, str):
            unit = unit.lower()
            if unit == "degrees":
                unit = u.degree
            elif unit == "hours":
                unit = u.hour
            elif unit == "radians":
                unit = u.radian
            else:
                raise ValueError("The unit value provided was not one of u.degree, u.hour, u.radian'.")
        else:
                raise ValueError("The unit value provided was not one of u.degree, u.hour, u.radian'.")

        if unit == u.degree:
            if decimal:
                return ("{0:0." + str(precision) + "f}").format(self.degrees)
            else:
                return util.degrees_to_string(self.degrees, precision=precision, sep=sep, pad=pad)

        elif unit == u.radian:
            return str(self.radians)

        elif unit == u.hour:
            if decimal:
                return ("{0:0." + str(precision) + "f}").format(self.hours)
            else:
                return util.hours_to_string(self.hours, precision=precision, sep=sep, pad=pad)
        else:
            raise ValueError(unit)

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
            return Angle(self.radians * other, unit=u.radian)
        else:
            raise NotImplementedError("An {0} object can only be multiplied by a float or integer.".format(type(self).__name__))

    # Division
    def __div__(self, other):
        if isinstance(other, type(self)):
            raise NotImplementedError("Division is not supported between two {0} "
                                      "objects.".format(type(self).__name__))
        elif type(other) in [float, int]:
            return Angle(self.radians / other, unit=u.radian)
        else:
            raise NotImplementedError("An {0} object can only be divided by a float or integer.".format(type(self).__name__))

    def __truediv__(self, other):
        raise NotImplementedError("Division is not supported between two {0} "
                                  "objects.".format(type(self).__name__))

    def __neg__(self):
        return Angle(-self.radians, unit=u.radian)

    def __eq__(self, other):
        if isinstance(other, Angle):
            return self.radians == other.radians
        if other == None:
            return False
        else:
            raise NotImplementedError("To compare {0} objects, compare their "
                                      "float values directly.".format(type(self).__name__))

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

    def __repr__(self):
        return "<{0}.{1} {2:.5f} deg>".format(__name__, type(self).__name__, self.degrees)


class RA(Angle):
    """ An object that represents a J2000 right ascension angle.

        This object can be created from a numeric value along with a unit. If the
        value specified is greater than "24", then a unit of degrees is assumed. Bounds
        are fixed to [0,360] degrees.

        Parameters
        ----------
        angle : float, int, str
            The angle value
        unit : `~astropy.units` (preferred), str
            The unit of the value specified for the angle. It is preferred that
            the unit be an object from the `~astropy.units` package, e.g.
            "from astropy import units as u; u.degree". Also accepts any string that the Unit class
            maps to "degrees", "radians", "hours". If not unit value is provided and
            the value of the angle is greater than "24.0", then the units are assumed
            to be degrees, otherwise, an exception is raised.
    """

    def __init__(self, angle, unit=None):

        # This block attempts to determine the validity of the unit,
        # particularly with regard to the specifics of RA.
        # After this block, the normal Angle initializer handles most of the
        # validation/creation.

        if isinstance(angle, type(self)):
            return super(RA, self).__init__(angle.radians, unit=u.radian, bounds=(0, 360))

        if unit == u.hour:
            pass  # to Angle initializer
            #self._radians = math.radians(decimal_hours * 15.)
        elif unit == u.degree:
            pass  # to Angle initializer
            #decimal_degrees = util.parse_degrees(angle)
            #self._radians = math.radians(decimal_degrees)
        elif unit == u.radian:
            pass  # to Angle initializer
            #self._radians = util.parse_radians(angle)
        elif unit == None:
            # Try to figure out the unit if we can.
            if isinstance(angle, float) or isinstance(angle, int):
                if angle > 24:
                    unit = u.degree
                else:
                    raise ValueError("No units were specified, and the angle value was ambiguous between hours and degrees.")
            elif isinstance(angle, str):
                # Try to deduce the units from hints in the string.
                # Further, enforce absolute bounds here, i.e. don't let
                # Angle +-2π to see if the angle falls in the bounds.
                if "d" in angle or "°" in angle:
                    # If in the form "12d32m53s", look for the "d" and assume degrees.
                    angle = math.radians(util.parse_degrees(angle))
                    if 0 < angle < twopi:
                        unit = u.radian
                    else:
                        raise RangeError("The provided angle was assumed to be in degrees, but was out of the range (0,360) degrees.")
                elif "h" in angle:
                    # Same for "12h32m53s" for hours.
                    #self._radians = math.radians(util.parse_hours(angle)*15.0)
                    unit = u.hour
                else:
                    # could be in a form: "54:43:26" -
                    # if so AND the resulting decimal value is > 24 or < -24, assume degrees
                    decimal_value = util.parse_degrees(angle)
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
        super(RA, self).__init__(angle, unit=unit, bounds=(0, 360))

    def hour_angle(self, lst, unit=None):
        """ Given a local sidereal time (LST), returns the hour angle for this RA.

            Parameters
            ----------
            lst : float, str, `~astropy.coordinates.angle`
                A local sidereal time (LST)
            unit : str
                The units of the LST, if not an `~astropy.coordinates.angle` object or datetime.datetime object
                .. note::
                    * if lst is **not** an `~astropy.coordinates.angle`-like object, you can specify the units by passing a `unit` parameter into the call
                    * this function currently returns an `Angle` object
        """
        # TODO : this should return an HA() object, and accept an Angle or LST object
        if not isinstance(lst, Angle):
            lst = Angle(lst, unit=unit)

        return Angle(lst.radians - self.radians, unit=u.radian)

    def lst(self, hour_angle, unit=u.hour):
        """
        Given an hour angle, calculate the local sidereal time (LST), returning an `~astropy.coordinates.Angle` object.

        Parameters
        ----------
        ha :  float, str, `~astropy.coordinates.angle`
            An hour angle
        unit : `~astropy.units` (preferred), str
            The unit of the value specified for the hour angle if it cannot be determined
            from the value provided. It is preferred that the unit be an object from the
            `~astropy.units` package, e.g. "from astropy import units as u; u.degree".
            Also accepts any string that the Unit class maps to "degrees", "radians", "hours".
        """
        # TODO : I guess this should return an HA() object, and accept an Angle or LST object
        if not isinstance(ha, Angle):
            ha = Angle(ha, unit)

        return Angle(ha.radians + self.radians, units=u.radian)


class Dec(Angle):
    """
    Represents a J2000 declination value.

    This object can be created from a numeric value along with a unit, or else a
    string in any commonly represented format, e.g. "12 43 23.53", "-32d52m29s".
    Unless otherwise specified via the 'unit' parameter, degrees are assumed.
    Bounds are fixed to [-90,90] degrees.

    Parameters
    ----------
    angle : float, int, str
        The angle value
    units : `~astropy.units` (preferred), str
        The units of the specified angle
    """

    def __init__(self, angle, unit=u.degree):
        if isinstance(angle, type(self)):
            return super(Dec, self).__init__(angle.radians, unit=u.radian, bounds=(-90, 90))

        super(Dec, self).__init__(angle, unit=unit, bounds=(-90, 90))


class AngularSeparation(Angle):
    """
    An on-sky separation between two directions.

    Parameters
    ----------
    lat1 : float
        The value of the first latitudinal/elevation angle.
    long1 : float
        The value of the first longitudinal/azimuthal angle.
    lat2 : float
        The value of the second latitudinal/elevation angle.
    long2 : float
        The value of the second longitudinal/azimuthal angle.
    units : `~astropy.units`
        The units of the given angles.

    .. note::


    """
    def __init__(self, lat1, long1, lat2, long2, units):

        units = u.Unit(units)
        lat1 = units.to(u.radian, lat1)
        if 0 == long1 == lat2 == long2:
            sepval = lat1
        else:
            long1 = units.to(u.radian, long1)
            lat2 = units.to(u.radian, lat2)
            long2 = units.to(u.radian, long2)

            sepval = self._haversine_dist_atan(lat1, long1, lat2, long2)

        super(AngularSeparation, self).__init__(sepval, u.radian)

    @staticmethod
    def _small_angle_dist(lat1, long1, lat2, long2):
        """
        Euclidean distance - only valid on sphere in the small-angle
        approximation.
        """

        dlat = lat2 - lat1
        dlong = long2 - long1

        return (dlat ** 2 + dlong ** 2) ** 0.5

    @staticmethod
    def _sphere_dist(lat1, long1, lat2, long2):
        """
        Simple formula for distance on a sphere: numerically unstable
        for small distances

        inputs must be in radians
        """
        #FIXME: array: use numpy functions
        from math import acos, sin, cos

        dlong = long2 - long1
        return acos(sin(lat1) * sin(-lat2) + cos(lat1) * cos(-lat2) * dlong)

    @staticmethod
    def _haversine_dist(lat1, long1, lat2, long2):
        """
        Haversine formula for distance on a sphere: more stable at poles

        inputs must be in radians
        """
        #FIXME: array: use numpy functions
        from math import asin, sin, cos

        sdlat = sin((lat2 - lat1) / 2)
        sdlong = sin((long2 - long1) / 2)
        coslats = cos(lat1) * cos(lat2)

        return 2 * asin((sdlat ** 2 + coslats * sdlong ** 2) ** 0.5)

    @staticmethod
    def _haversine_dist_atan(lat1, long1, lat2, long2):
        """
        Haversine formula for distance on a sphere: more stable at poles.
        This version uses arctan instead of arcsin and thus does better
        with sign convnentions.

        inputs must be in radians
        """
        #FIXME: array: use numpy functions
        from math import atan2, sin, cos

        sdlat = sin((lat2 - lat1) / 2)
        sdlong = sin((long2 - long1) / 2)
        coslats = cos(lat1) * cos(lat2)

        numerator = sdlat ** 2 + coslats * sdlong ** 2

        return atan2(numerator, 1 - numerator)

    @staticmethod
    def _vicenty_dist(lat1, long1, lat2, long2):
        """
        Vincenty formula for distance on a sphere: stable at poles and
        antipodes but more complex/computationally expensive

        inputs must be in radians
        """
        #FIXME: array: use numpy functions
        from math import atan2, sin, cos

        dlong = long2 - long1

        num1 = cos(lat2) * sin(dlong)
        num2 = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlong)
        denominator = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlong)

        return atan2((num1 ** 2 + num2 ** 2) ** 0.5, denominator)

    def __add__(self, other):
        raise TypeError('+ is ambiguous for AngularSeparation objects; not supported')

    def __radd__(self, other):
        raise TypeError('+ is ambiguous for AngularSeparation objects; not supported')

    def __sub__(self, other):
        raise TypeError('- is ambiguous for AngularSeparation objects; not supported')

    def __rsub__(self, other):
        raise TypeError('- is ambiguous for AngularSeparation objects; not supported')

    @property
    def arcmins(self):
        """
        The value of this separation in arcminutes.
        """
        return self.degrees * 60.

    @property
    def arcsecs(self):
        """
        The value of this separation in arcseconds.
        """
        return self.degrees * 3600.
