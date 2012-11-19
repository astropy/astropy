# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the fundamental classes used for representing
coordinates in astropy.
"""

import math

import numpy as np

from . import angle_utilities as util
from .errors import *
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
        angle : float, int, str, tuple
            The angle value. If a tuple, will be interpreted as (h, m s) or
            (d, m, s) depending on `unit`.
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
        from ..utils import isiterable

        self._bounds = bounds

        if isinstance(angle, type(self)):
            angle = angle.radians
            unit = u.radian

        self.is_array = (isiterable(angle) and
                         not isinstance(angle, basestring) and
                         not isinstance(angle, tuple))

        # short circuit arrays for now
        if self.is_array:
            raise NotImplementedError("Angles as arrays are not yet supported.")

        if angle == None:
            raise UnitsError("The Angle class requires a unit")

        # -------------------------------
        # unit validation and angle value
        # -------------------------------
        if isinstance(unit, u.UnitBase):
            pass
        elif isinstance(unit, basestring):
            unit = u.Unit(unit)
        elif unit is None:
            # try to determine unit from the "angle" value
            if self.is_array:
                try:
                    angle = [x.lower() for x in angle]
                except AttributeError:
                        # If units are not specified as a parameter, the only chance
                        # to determine them is in a string value - if it's not a string,
                        # then there's not enough information to create an Angle.
                        raise UnitsError("Could not parse an angle value in the array provided"
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
                    if unit is None:
                        for unitStr in ["hours", "hour", "hr"]:
                            if unitStr in a:
                                a_unit = u.radian
                                a = angle.replace(unitStr, "")
                                angle[idx] = math.radians(util.parse_hours(a) * 15.)
                                break
                    if unit is None:
                        for unitStr in ["radians", "radian", "rad"]:
                            if unitStr in angle:
                                a_unit = u.radian
                                a = angle.replace(unitStr, "")
                                angle[idx] = util.parse_radians(a)
                                break
                    if a_unit is None:
                        raise UnitsError("Could not parse the angle value '{0}' "
                                         "- units could not be determined.".format(angle[idx]))
                unit = u.radian

            else:  # single value
                if isinstance(angle, basestring):
                    angle = angle.lower()
                else:
                    raise UnitsError("Could not parse the angle value '{0}' "
                                     "- units could not be determined.".format(angle))
                angle = angle.strip()
                for unitStr in ["degrees", "degree", "deg"]:
                    if unitStr in angle:
                        unit = u.degree
                        angle = angle.replace(unitStr, "")
                if unit is None:
                    for unitStr in ["hours", "hour", "hr"]:
                        if unitStr in angle:
                            unit = u.hour
                            angle = angle.replace(unitStr, "")
                if unit is None:
                    for unitStr in ["radians", "radian", "rad"]:
                        if unitStr in angle:
                            unit = u.radian
                            angle = angle.replace(unitStr, "")
                if unit is None:
                    if "h" in angle:
                        unit = u.hour
                    elif "d" in angle or "°" in angle:
                        unit = u.degree

        if unit is None:
            raise UnitsError("The unit parameter should be an object from the "
                             "astropy.unit module (e.g. 'from astropy import units as u',"
                             "then use 'u.degree').")

        if self.is_array:
            pass  # already performed conversions to radians above
        else:
            if unit is u.degree:
                self._radians = math.radians(util.parse_degrees(angle))
            elif unit is u.radian:
                self._radians = float(angle)
            elif unit is u.hour:
                self._radians = util.hours_to_radians(util.parse_hours(angle))
            else:
                raise UnitsError("The unit value provided was not one of u.degree, u.hour, u.radian'.")

        # ---------------
        # bounds checking
        # ---------------
        # TODO: handle arrays
        # handle bounds units, convert to radians
        if bounds == None:
            pass  # no range checking performed
        else:
            try:
                if unit is u.radian:
                    lower_bound = bounds[0]
                    upper_bound = bounds[1]
                elif unit is u.degree:
                    lower_bound = math.radians(bounds[0])
                    upper_bound = math.radians(bounds[1])
                elif unit is u.hour:
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
                            raise BoundsError("The angle given falls outside of the specified bounds.")
                        elif lower_bound < self._radians < upper_bound:
                            break

                if self._radians < lower_bound:
                    while True:
                        self._radians += twopi
                        if self._radians > upper_bound:
                            raise BoundsError("The angle given falls outside of the specified bounds.")
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

        if isinstance(unit, u.UnitBase):
            pass  # great!
        elif isinstance(unit, basestring):
            unit = unit.lower()
            if unit is "degrees":
                unit = u.degree
            elif unit is "hours":
                unit = u.hour
            elif unit is "radians":
                unit = u.radian
            else:
                raise UnitsError("The unit value provided was not one of u.degree, u.hour, u.radian'.")
        else:
                raise UnitsError("The unit value provided was not one of u.degree, u.hour, u.radian'.")

        if unit is u.degree:
            if decimal:
                return ("{0:0." + str(precision) + "f}").format(self.degrees)
            else:
                return util.degrees_to_string(self.degrees, precision=precision, sep=sep, pad=pad)

        elif unit is u.radian:
            return str(self.radians)

        elif unit is u.hour:
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
                msg = "Can't add angles because bounds don't match: {0} and {1}"
                raise ValueError(msg.format(self.bounds, other.bounds))
            else:
                return Angle(self.radians + other.radians, unit=u.radian, bounds=self.bounds)
        else:
            raise NotImplementedError("An {0} object can only be added to another "
                                      "{0} object.".format(type(self).__name__))

    def __radd__(self, other):
        return self.__add__(other)

    # Subtraction
    def __sub__(self, other):
        if isinstance(other, type(self)):
            if self.bounds != other.bounds:
                msg = "Can't add angles because bounds don't match: {0} and {1}"
                raise ValueError(msg.format(self.bounds, other.bounds))
            else:
                return Angle(self.radians - other.radians, unit=u.radian, bounds=self.bounds)
        else:
            raise NotImplementedError("An {0} object can only be subtracted from another "
                                      "{0} object.".format(type(self).__name__))

    def __rsub__(self, other):
        if isinstance(other, type(self)):
            if self.bounds != other.bounds:
                msg = "Can't add angles because bounds don't match: {0} and {1}"
                raise ValueError(msg.format(self.bounds, other.bounds))
            else:
                return Angle(other.radians - self.radians, unit=u.radian, bounds=self.bounds)
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

    def __rmul__(self, other):
        return self.__mul__(other)

    # Division
    def __div__(self, other):
        if isinstance(other, type(self)):
            raise NotImplementedError("Division is not supported between two {0} "
                                      "objects.".format(type(self).__name__))
        elif type(other) in [float, int]:
            return Angle(self.radians / other, unit=u.radian)
        else:
            raise NotImplementedError("An {0} object can only be divided by a float or integer.".format(type(self).__name__))

    def __rdiv__(self, other):
        if isinstance(other, type(self)):
            raise NotImplementedError("Division is not supported between two {0} "
                                      "objects.".format(type(self).__name__))
        elif type(other) in [float, int]:
            return Angle(other / self.radians, unit=u.radian)
        else:
            raise NotImplementedError("An {0} object can only be divided by a float or integer.".format(type(self).__name__))

    def __truediv__(self, other):
        return self.__div__(other)

    def __rtruediv__(self, other):
        return self.__rdiv__(other)

    # other operations

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
        return "<{0} {1:.5f} deg>".format(type(self).__name__, self.degrees)


class RA(Angle):
    """ An object that represents a right ascension angle.

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

        if isinstance(angle, Angle):
            return super(RA, self).__init__(angle.radians, unit=u.radian, bounds=(0, 360))

        if unit is u.hour:
            pass  # to Angle initializer
            #self._radians = math.radians(decimal_hours * 15.)
        elif unit is u.degree:
            pass  # to Angle initializer
            #decimal_degrees = util.parse_degrees(angle)
            #self._radians = math.radians(decimal_degrees)
        elif unit is u.radian:
            pass  # to Angle initializer
            #self._radians = util.parse_radians(angle)
        elif unit is None:
            # Try to figure out the unit if we can.
            if isinstance(angle, float) or isinstance(angle, int):
                if angle > 24:
                    unit = u.degree
                else:
                    raise UnitsError("No units were specified, and the angle value was ambiguous between hours and degrees.")
            elif isinstance(angle, basestring):
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
                        raise UnitsError("No units were specified, and the angle value was ambiguous between hours and degrees.")
                    elif decimal_value < 0:
                        raise RangeError("No units were specified; could not assume any since the value was less than zero.")
            elif isinstance(angle, tuple):
                if len(angle) == 3 and 0 <= angle[0] < 24.0:
                    raise UnitsError("No units were specified, and the angle value was ambiguous between hours and degrees.")
                else:
                    unit = u.degree
            else:
                raise ValueError("Angle values of type {0} not supported.".format(type(angle).__name__))
        if unit is None:
            raise UnitsError("Units must be specified for RA, one of u.degree, u.hour, or u.radian.")

        # By here, the unit should be defined.
        super(RA, self).__init__(angle, unit=unit, bounds=(0, 360))

    def hour_angle(self, lst):
        """ Computes the hour angle for this RA given a local sidereal
        time (LST).

        Parameters
        ----------
        lst : `~astropy.coordinates.angle.Angle`, `~astropy.time.Time`
            A local sidereal time (LST).

        Returns
        -------
        hour_angle : `~astropy.coordinates.angle.Angle`
            The hour angle for this RA at the LST `lst`.
        """
        if hasattr(lst, 'mjd'):
            lst = Angle(np.remainder(lst.mjd, 1), unit=u.hour)

        return Angle(lst.radians - self.radians, unit=u.radian, bounds=(0, 360.))

    def lst(self, hour_angle):
        """
        Calculates the local sidereal time (LST) if this RA is at a
        particular hour angle.

        Parameters
        ----------
        hour_angle :  `~astropy.coordinates.angle.Angle`
            An hour angle.

        Returns
        -------
        lst : `~astropy.coordinates.angle.Angle`
            The local siderial time as an angle.

        """
        return Angle(hour_angle.radians + self.radians, unit=u.radian, bounds=(0, 360.))


class Dec(Angle):
    """
    Represents a declination value.

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
        if isinstance(angle, Angle):
            return super(Dec, self).__init__(angle.radians, unit=u.radian, bounds=(-90, 90))

        super(Dec, self).__init__(angle, unit=unit, bounds=(-90, 90))


class AngularSeparation(Angle):
    """
    An on-sky separation between two directions.

    Parameters
    ----------
    lat1 : float
        The value of the first latitudinal/elevation angle.
    lon1 : float
        The value of the first longitudinal/azimuthal angle.
    lat2 : float
        The value of the second latitudinal/elevation angle.
    lon2 : float
        The value of the second longitudinal/azimuthal angle.
    units : `~astropy.units`
        The units of the given angles.


    """
    def __init__(self, lat1, lon1, lat2, lon2, units):

        units = u.Unit(units)
        lat1 = units.to(u.radian, lat1)
        if 0 == lon1 == lat2 == lon2:
            sepval = lat1
        else:
            lon1 = units.to(u.radian, lon1)
            lat2 = units.to(u.radian, lat2)
            lon2 = units.to(u.radian, lon2)

            sepval = util.vicenty_sphere_dist(lat1, lon1, lat2, lon2)

        super(AngularSeparation, self).__init__(sepval, u.radian)

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

#<----------------------------------Rotations---------------------------------->


def rotation_matrix(angle, axis='z', degrees=True):
    """
    Generate a 3x3 cartesian rotation matrix in for rotation about
    a particular axis.

    Parameters
    ----------
    angle : scalar
        The amount of rotation this matrix should represent. In degrees
        if `degrees` is True, otherwise radians.
    axis : str or 3-sequence
        Either 'x','y', 'z', or a (x,y,z) specifying an axis to rotate
        about. If 'x','y', or 'z', the rotation sense is
        counterclockwise looking down the + axis (e.g. positive
        rotations obey left-hand-rule).
    degrees : bool
        If True the input angle is degrees, otherwise radians.

    Returns
    -------
    rmat: `numpy.matrix`
        A unitary rotation matrix.
    """
    from math import sin, cos, radians, sqrt

    if degrees:
        angle = radians(angle)

    if axis == 'z':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c, s, 0),
                          (-s, c, 0),
                          ( 0, 0, 1)))
    elif axis == 'y':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c, 0,-s),
                          ( 0, 1, 0),
                          ( s, 0, c)))
    elif axis == 'x':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( 1, 0, 0),
                          ( 0, c, s),
                          ( 0,-s, c)))
    else:
        x,y,z = axis
        w = cos(angle/2)

        #normalize
        if w == 1:
            x=y=z=0
        else:
            l = sqrt((x*x + y*y + z*z)/(1 - w*w))
            x /= l
            y /= l
            z /= l

        wsq = w*w
        xsq = x*x
        ysq = y*y
        zsq = z*z
        return np.matrix((( wsq+xsq-ysq-zsq, 2*x*y-2*w*z, 2*x*z+2*w*y),
                          ( 2*x*y+2*w*z, wsq-xsq+ysq-zsq,2*y*z-2*w*x),
                          ( 2*x*z-2*w*y, 2*y*z+2*w*x, wsq-xsq-ysq+zsq)))


def angle_axis(matrix, degrees=True):
    """
    Computes the angle of rotation and the rotation axis for a given rotation
    matrix.

    Parameters
    ----------
    matrix : array-like
        A 3 x 3 unitary rotation matrix.
    degrees : bool
        If True, output is in degrees.

    Returns
    -------
    angle : scalar
        The angle of rotation for this matrix. In degrees if `degrees is
        True, otherwise radians.
    axis : array (length 3)
        The axis of rotation for this matrix.

    """
    from math import sin, cos, acos, degrees, sqrt

    m = np.asmatrix(matrix)
    if m.shape != (3,3):
        raise ValueError('matrix is not 3x3')



    angle = acos((m[0,0] + m[1,1] + m[2,2] - 1)/2)
    denom = sqrt(2*((m[2,1]-m[1,2])+(m[0,2]-m[2,0])+(m[1,0]-m[0,1])))
    axis = np.array((m[2,1]-m[1,2],m[0,2]-m[2,0],m[1,0]-m[0,1]))/denom
    axis /= sqrt(np.sum(axis**2))

    if degrees:
        return degrees(angle),axis
    else:
        return angle,axis
