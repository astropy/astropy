# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the fundamental classes used for representing
coordinates in astropy.
"""
from __future__ import unicode_literals

import math

import numpy as np

from . import angle_utilities as util
from .errors import *
from .. import units as u
from ..utils import deprecated


__all__ = ['Angle', 'RA', 'Dec', 'AngularSeparation']


TWOPI = math.pi * 2.0  # no need to calculate this all the time


class Angle(u.Quantity):
    """
    An angle.

    An angle can be specified either as an array, scalar, tuple (see
    below), string, `~astropy.units.Quantity` or another
    `~astropy.coordinates.Angle`.

    If a string, it must be in one of the following formats:

        * ``'1:2:30.43 degrees'``
        * ``'1 2 0 hours'``
        * ``'1°2′3″'``
        * ``'1d2m3.4s'``
        * ``'-1h2m3s'``

    Parameters
    ----------
    angle : array, scalar, Quantity, Angle
        The angle value. If a tuple, will be interpreted as ``(h, m
        s)`` or ``(d, m, s)`` depending on `unit`. If a string, it
        will be interpreted following the rules described above.

        If `angle` is a sequence or array of strings, the resulting
        values will be in the given `unit`, or if None is provided,
        the unit will be taken from the first given value.

    unit : `~astropy.units.UnitBase`, str, optional
        The unit of the value specified for the angle.  This may be
        any string that `~astropy.units.Unit` understands, but it is
        better to give an actual unit object.  Must be an angular
        unit.

    bounds : 2-sequence, optional
        A length-2 sequence indicating the upper and lower value that
        the new angle object may have.  Each value may be in any of
        the forms accepted by the `angle` argument.  If no units are
        specified, (that is the value is a float, tuple, or string
        without an explicit unit), the unit is taken from the `unit`
        argument.  Pass `None` to perform no bounds checking.  By
        default the bounds are `(-360, 360)` degrees.

    dtype : ~numpy.dtype, optional
        See `~astropy.units.Quantity`.

    equivalencies : list of equivalence pairs, optional
        See `~astropy.units.Quantity`.

    Raises
    ------
    `~astropy.units.core.UnitsException`
        If a unit is not provided or it is not an angular unit.
    """
    def __new__(cls, angle, unit=None, bounds=(), dtype=None,
                equivalencies=[]):
        unit = cls._convert_unit_to_angle_unit(unit)
        if (unit is not None and
            not unit.is_equivalent(u.radian, equivalencies)):
            raise u.UnitsException(
                "Given unit {0} is not convertible to an angle".format(
                    unit))

        if isinstance(angle, u.Quantity):
            # This includes Angle subclasses as well
            if unit is not None:
                angle = angle.to(unit).value
            else:
                unit = angle.unit
                angle = angle.value

                unit = cls._convert_unit_to_angle_unit(unit)
        elif isinstance(angle, basestring):
            angle, unit = util.parse_angle(angle, unit)

        angle = cls._tuple_to_float(angle, unit)

        angle = np.asarray(angle)

        if angle.dtype.type in (np.string_, np.unicode_):
            # We need to modify this value from within
            # convert_string_to_angle, and the only way to do that
            # across Python 2.6 - 3.3 is to use this "store it in a
            # list" trick.
            determined_unit = [unit]

            def convert_string_to_angle(x):
                ang, new_unit = util.parse_angle(str(x), unit)
                if determined_unit[0] is None:
                    determined_unit[0] = new_unit
                    return cls._tuple_to_float(ang, unit)
                else:
                    return new_unit.to(
                        determined_unit[0], cls._tuple_to_float(ang, unit))

            convert_string_to_angle_ufunc = np.vectorize(
                convert_string_to_angle,
                otypes=[np.float_])
            angle = convert_string_to_angle_ufunc(angle)
            unit = determined_unit[0]

        elif angle.dtype.kind not in 'iuf':
            raise TypeError("Unsupported dtype '{0}'".format(angle.dtype))

        if unit is None:
            raise u.UnitsException("No unit was specified")

        bounded_angle, bounds = cls._bounds_check(angle, bounds, unit)

        self = super(Angle, cls).__new__(
            cls, bounded_angle, unit, dtype=dtype,
            equivalencies=equivalencies)

        self._bounds = bounds

        return self

    @staticmethod
    def _tuple_to_float(angle, unit):
        """
        Converts an angle represented as a 3-tuple into a floating
        point number in the given unit.
        """
        if isinstance(angle, tuple):
            # TODO: Numpy array of tuples?
            if unit is u.hourangle:
                util.check_hms_ranges(*angle)
                angle = util.hms_to_hours(*angle)
            elif unit is u.degree:
                angle = util.dms_to_degrees(*angle)
            else:
                raise u.UnitsException(
                    "Can not parse '{0}' as unit '{1}'".format(
                        angle, unit))
        return angle

    @staticmethod
    def _get_default_bounds():
        return (Angle(-360, u.degree, bounds=None),
                Angle(360, u.degree, bounds=None))

    @staticmethod
    def _bounds_check(angle, bounds, unit):
        def raise_error(original_angle, lower_angle, upper_angle, unit):
            raise BoundsError(
                "The angle(s) {0} falls outside of the specified "
                "bounds ({1}, {2}) (in {3})".format(
                    u.radian.to(unit, original_angle),
                    u.radian.to(unit, lower_angle),
                    u.radian.to(unit, upper_angle),
                    unit))

        if bounds is None:
            return angle, None

        if bounds == ():
            bounds = Angle._get_default_bounds()
        else:
            if len(bounds) != 2:
                raise ValueError(
                    "Bounds specified for Angle must be a two-element "
                    "sequence, e.g. [0, 360] (was given '{0}').".format(bounds))

            new_bounds = []
            for bound in bounds:
                if isinstance(bound, u.Quantity):
                    bound = Angle(bound, bounds=None)
                else:
                    bound = Angle(bound, unit, bounds=None)
                new_bounds.append(bound)
            bounds = tuple(new_bounds)

        lower_bound, upper_bound = bounds
        # Convert everything to radians, and keep a copy of the
        # original (unmoved) values for any resulting error messages.
        original_angle = unit.to(u.radian, angle)
        angle = np.array(original_angle)
        lower_angle = lower_bound.radian
        upper_angle = upper_bound.radian

        # TODO: This is perhaps a candidate for something to do in C.
        TWOPI = np.pi * 2.0

        if np.all((lower_angle < angle) & (angle < upper_angle)):
            pass
        else:
            too_big = angle > upper_angle
            while np.any(too_big):
                angle = np.where(too_big, angle - TWOPI, angle)
                if np.any(angle < lower_angle):
                    raise_error(original_angle, lower_angle, upper_angle, unit)

                too_big = angle > upper_angle

            too_small = angle < lower_angle
            while np.any(too_small):
                angle = np.where(too_small, angle + TWOPI, angle)
                if np.any(angle > upper_angle):
                    raise_error(original_angle, lower_angle, upper_angle, unit)

                too_small = angle < lower_angle

        return u.radian.to(unit, angle), bounds

    @staticmethod
    def _convert_unit_to_angle_unit(unit):
        if unit is not None:
            unit = u.Unit(unit)

            if unit is u.hour:
                unit = u.hourangle
        return unit

    def __quantity_view__(self, obj, unit):
        unit = self._convert_unit_to_angle_unit(unit)
        if unit is not None and unit.is_equivalent(u.radian):
            result = obj.view(Angle)
            result._bounds = self.bounds
            return result
        return super(Angle, self).__quantity_view__(
            obj, unit)

    def __quantity_instance__(self, val, unit, dtype=None, equivalencies=[]):
        unit = self._convert_unit_to_angle_unit(unit)
        if unit is not None and unit.is_equivalent(u.radian):
            return Angle(val, unit, bounds=self.bounds, dtype=dtype,
                         equivalencies=equivalencies)
        return super(Angle, self).__quantity_instance__(
            val, unit, dtype=dtype, equivalencies=equivalencies)

    def __array_wrap__(self, obj, context=None):
        obj = super(Angle, self).__array_wrap__(obj, context=context)

        if isinstance(obj, Angle):
            return Angle(obj.value, obj.unit, bounds=obj.bounds)

        return obj

    def _bounds_match_check(self, other, operation):
        if (isinstance(other, Angle) and
            self.bounds != other.bounds):
            default_bounds = self._get_default_bounds()
            if (self.bounds != default_bounds and
                other.bounds != default_bounds):
                msg = "Can't {0} angles because bounds don't match: {1} and {2}"
                raise BoundsError(msg.format(
                        operation, self.bounds, other.bounds))

    def __add__(self, other):
        self._bounds_match_check(other, 'add')
        return super(Angle, self).__add__(other)

    def __sub__(self, other):
        self._bounds_match_check(other, 'subtract')
        return super(Angle, self).__sub__(other)

    def __mul__(self, other):
        if isinstance(other, type(self)):
            raise TypeError(
                "multiplication is not supported between two {0} "
                "objects".format(
                    type(self).__name__))
        return super(Angle, self).__mul__(other)

    def __div__(self, other):
        if isinstance(other, type(self)):
            raise TypeError(
                "division is not supported between two {0} objects".format(
                    type(self).__name__))
        return super(Angle, self).__div__(other)

    __truediv__ = __div__

    @property
    def bounds(self):
        """
        The angle's bounds, an immutable property.  Returns a 2-tuple
        of `Angle` objects.
        """
        return self._bounds

    @property
    def hour(self):
        """
        The angle's value in hours (read-only property).
        """
        return self.hourangle

    @property
    def hms(self):
        """
        The angle's value in hours, as a ``(h, m, s)`` tuple
        (read-only property).
        """
        return util.hours_to_hms(self.hourangle)

    @property
    def dms(self):
        """
        The angle's value in degrees, as a ``(d, m, s)`` tuple
        (read-only property).
        """
        return util.degrees_to_dms(self.degree)

    def to_string(self, unit=None, decimal=False, sep='fromunit',
                  precision=5, alwayssign=False, pad=False):
        """ A string representation of the angle.

        Parameters
        ----------
        units : `~astropy.units.UnitBase`, optional
            Specifies the units.  Must be an angular unit.  If not
            provided, the unit used to initialize the angle will be
            used.

        decimal : bool, optional
            If True, a decimal respresentation will be used, otherwise
            the returned string will be in sexagesimal form.

        sep : str, optional
            The separator between numbers in a sexagesimal
            representation.  E.g., if it is ':', the result is
            "12:41:11.1241". Also accepts 2 or 3 separators. E.g.,
            ``sep='hms'`` would give the result "12h41m11.1241s", or
            sep='-:' would yield "11-21:17.124".  Alternatively, the
            special string 'fromunit' means 'dms' if the unit is
            degrees, or 'hms' if the unit is hours.

        precision : int, optional
            The level of decimal precision.  If `decimal` is True,
            this is the raw precision, otherwise it gives the
            precision of the last place of the sexagesimal
            representation (seconds).

        alwayssign : bool, optional
            If `True`, include the sign no matter what.  If `False`,
            only include the sign if it is negative.

        pad : bool, optional
            If `True`, include leading zeros when needed to ensure a
            fixed number of characters for sexagesimal representation.

        Returns
        -------
        strrepr : str
            A string representation of the angle.

        """
        if unit is None:
            unit = self.unit
        unit = self._convert_unit_to_angle_unit(unit)

        # Create an iterator so we can format each element of what
        # might be an array.
        if unit is u.degree:
            if decimal:
                values = self.degree
                func = ("{0:0." + str(precision) + "}").format
            else:
                if sep == 'fromunit':
                    sep = 'dms'
                values = self.degree
                func = lambda x: util.degrees_to_string(
                    x, precision=precision, sep=sep, pad=pad)

        elif unit is u.hourangle:
            if decimal:
                values = self.hour
                func = ("{0:0." + str(precision) + "}").format
            else:
                if sep == 'fromunit':
                    sep = 'hms'
                values = self.hour
                func = lambda x: util.hours_to_string(
                    x, precision=precision, sep=sep, pad=pad)

        elif unit.is_equivalent(u.radian):
            if decimal:
                values = self.to(unit).value
                func = ("{0:0." + str(precision) + "}").format
            elif sep == 'fromunit':
                values = self.to(unit).value
                func = ("{0:0." + str(precision) + "}" + unit.name).format
            else:
                raise ValueError(
                    "'{0}' can not be represented in sexagesimal "
                    "notation".format(
                        unit.name))

        else:
            raise u.UnitsException(
                "The unit value provided is not an angular unit.")

        def format(val):
            s = func(float(val))
            if alwayssign and not s.startswith('-'):
                s = '+' + s
            return s

        format_ufunc = np.vectorize(format, otypes=[np.object])
        return format_ufunc(values)

    @deprecated("0.3", name="format", alternative="to_string")
    def format(self, unit=u.degree, decimal=False, sep='fromunit', precision=5,
               alwayssign=False, pad=False):
        return self.to_string(
            unit=unit, decimal=decimal, sep=sep, precision=precision,
            alwayssign=alwayssign, pad=pad)

    def __str__(self):
        return str(self.to_string())


class RA(Angle):
    """
    An object that represents a right ascension angle.

    This object can be created from a numeric value along with a
    unit.

    Parameters
    ----------
    angle : array, scalar, str, tuple, Quantity, Angle
        The angle value. If a tuple, will be interpreted as `(h, m s)`
        or `(d, m, s)` depending on `unit`. If a string, it will be
        interpreted following the rules described in
        `~astropy.coordinates.Angle`.

    unit : `~astropy.units.UnitBase`, str, optional
        The unit of the value specified for the angle.  This may be
        any string that `~astropy.units.Unit` understands, but it is
        better to give an actual unit object.  Must be one an angular
        unit.

    Raises
    ------
    `~astropy.coordinates.errors.UnitsException`
        If a unit is not provided or it is not an angular unit.
    """

    def __new__(cls, angle, unit=None):
        return super(RA, cls).__new__(
            cls, angle, unit=unit, bounds=(
                Angle(0, u.degree, bounds=None),
                Angle(360, u.degree, bounds=None)))

    # The initializer as originally conceived allowed the unit to be unspecified
    # if it's bigger  than 24, because hours typically aren't past 24.
    # It is also then implicit that `RA` is usually in hours.
    # This is commented out for now because in discussion for v0.2 we decided to
    # keep things simple and just use the same behavior as Angle.
    # In v0.3 it may be either uncommented,
    # moved somewhere else, or eliminated completely
    # TODO: decide if this should stay in or permanently be removed
    #
    # def __init__(self, angle, unit=None):
    #
    #     # This block attempts to determine the validity of the unit,
    #     # particularly with regard to the specifics of RA.
    #     # After this block, the normal Angle initializer handles most of the
    #     # validation/creation.
    #
    #     if isinstance(angle, Angle):
    #         return super(RA, self).__init__(angle.radians, unit=u.radian, bounds=(0, 360))
    #
    #     if unit is u.hour:
    #         pass  # to Angle initializer
    #         # self._radians = math.radians(decimal_hours * 15.)
    #     elif unit is u.degree:
    #         pass  # to Angle initializer
    #         # decimal_degrees = util.parse_degrees(angle)
    #         # self._radians = math.radians(decimal_degrees)
    #     elif unit is u.radian:
    #         pass  # to Angle initializer
    #         # self._radians = util.parse_radians(angle)
    #
    #
    #     elif unit is None:
    #         # Try to figure out the unit if we can.
    #         if isinstance(angle, float) or isinstance(angle, int):
    #             if angle > 24:
    #                 unit = u.degree
    #             else:
    #                 raise u.UnitsException("No units were specified, and the angle value was ambiguous between hours and degrees.")
    #         elif isinstance(angle, basestring):
    #             # Try to deduce the units from hints in the string.
    #             # Further, enforce absolute bounds here, i.e. don't let
    #             # Angle +-2π to see if the angle falls in the bounds.
    #             if "d" in angle or "°" in angle:
    #                 # If in the form "12d32m53s", look for the "d" and assume degrees.
    #                 angle = math.radians(util.parse_degrees(angle))
    #                 if 0 < angle < TWOPI:
    #                     unit = u.radian
    #                 else:
    #                     raise RangeError("The provided angle was assumed to be in degrees, but was out of the range (0,360) degrees.")
    #             elif "h" in angle:
    #                 # Same for "12h32m53s" for hours.
    #                 # self._radians = math.radians(util.parse_hours(angle)*15.0)
    #                 unit = u.hour
    #             else:
    #                 # could be in a form: "54:43:26" -
    #                 # if so AND the resulting decimal value is > 24 or < -24, assume degrees
    #                 decimal_value = util.parse_degrees(angle)
    #                 if decimal_value > 24:
    #                     unit = u.degree
    #                 elif 0 <= decimal_value <= 24.0:
    #                     raise u.UnitsException("No units were specified, and the angle value was ambiguous between hours and degrees.")
    #                 elif decimal_value < 0:
    #                     raise RangeError("No units were specified; could not assume any since the value was less than zero.")
    #         elif isinstance(angle, tuple):
    #             if len(angle) == 3 and 0 <= angle[0] < 24.0:
    #                 raise u.UnitsException("No units were specified, and the angle value was ambiguous between hours and degrees.")
    #             else:
    #                 unit = u.degree
    #         else:
    #             raise ValueError("Angle values of type {0} not supported.".format(type(angle).__name__))
    #
    #     if unit is None:
    #         raise u.UnitsException("Units must be specified for RA, one of u.degree, u.hour, or u.radian.")
    #
    #     # By here, the unit should be defined.
    #     super(RA, self).__init__(angle, unit=unit, bounds=(0, 360))

    def hour_angle(self, lst):
        """
        Computes the hour angle for this RA given a local sidereal
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

        return Angle(
            lst.radian - self.radian, unit=u.radian, bounds=(0, TWOPI))

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
        return Angle(
            hour_angle.radian + self.radian, unit=u.radian, bounds=(0, TWOPI))


class Dec(Angle):
    """
    Represents a declination value.

    This object can be created from a numeric value along with a unit,
    or else a string in any commonly represented format, e.g. "12 43
    23.53", "-32d52m29s".  Unless otherwise specified via the 'unit'
    parameter, degrees are assumed.  Bounds are fixed to [-90,90]
    degrees.

    Parameters
    ----------
    angle : array, scalar, str, tuple, Quantity, Angle
        The angle value. If a tuple, will be interpreted as `(h, m s)`
        or `(d, m, s)` depending on `unit`. If a string, it will be
        interpreted following the rules described in
        `~astropy.coordinates.Angle`.

    unit : `~astropy.units.UnitBase`, str, optional
        The unit of the value specified for the angle.  This may be
        any string that `~astropy.units.Unit` understands, but it is
        better to give an actual unit object.  Must be one an angular
        unit.

    Raises
    ------
    `~astropy.units.core.UnitsException`
        If a unit is not provided or it is not an angular unit.
    """
    def __new__(cls, angle, unit=None):
        return super(Dec, cls).__new__(
            cls, angle, unit=unit,
            bounds=(Angle(-90, u.degree, bounds=None),
                    Angle(90, u.degree, bounds=None)))

    # TODO: do here whatever is decided for the "smart" RA initializer above
    #
    # def __init__(self, angle, unit=u.degree):
    #     super(RA, self).__init__(angle, unit=unit, bounds=(0, 360))
    #
    #     if isinstance(angle, Angle):
    #         return super(Dec, self).__init__(angle.radians, unit=u.radian, bounds=(-90, 90))
    #
    #     super(Dec, self).__init__(angle, unit=unit, bounds=(-90, 90))


class AngularSeparation(Angle):
    """
    An on-sky separation between two directions.

    .. note::
        This is computed using the Vincenty great circle distance
        formula, and hence should be numerically stable even for
        near antipodal points.

    Parameters
    ----------
    lon1 : float
        The value of the first longitudinal/azimuthal angle.

    lat1 : float
        The value of the first latitudinal/elevation angle.

    lon2 : float
        The value of the second longitudinal/azimuthal angle.

    lat2 : float
        The value of the second latitudinal/elevation angle.

    unit : `~astropy.units`
        The unit of the given angles.
    """
    def __new__(cls, lon1, lat1, lon2, lat2, unit,
        _supresslatlonswap_warning=False):  # TODO: remove this parameter in v0.4
        # TODO: remove this warning in v0.4
        if not _supresslatlonswap_warning:
            from warnings import warn
            from ..utils.exceptions import AstropyBackwardsIncompatibleChangeWarning
            warn(AstropyBackwardsIncompatibleChangeWarning('The ordering of '
                ' the AngularSeparation initializer angles was changed '
                'from lat1, lon1, lat2, lon2 in v0.2 to "lon1, lat1, lon2, '
                'lat2" in v0.3.  You MUST update your code to swap lat/lon '
                'if you are not using keywords, or you will get the wrong '
                'result.'))

        unit = cls._convert_unit_to_angle_unit(unit)
        lon1 = unit.to(u.radian, lon1)
        if lat1 or lat2 or lon2:
            lat1 = unit.to(u.radian, lat1)
            lon2 = unit.to(u.radian, lon2)
            lat2 = unit.to(u.radian, lat2)

            sepval = util.vincenty_sphere_dist(lon1, lat1, lon2, lat2)
        else:  # this is the case where lat1, lat2, and lon2 are all 0 or None or False
            sepval = lon1

        self = super(AngularSeparation, cls).__new__(cls, sepval, u.radian)

        return self

    def __quantity__(self):
        unit = self._convert_unit_to_angle_unit(unit)
        if unit is not None and unit.is_equivalent(u.radian):
            return AngularSeparation
        return Quantity

    def __add__(self, other):
        raise TypeError('+ is ambiguous for AngularSeparation objects; not supported')

    def __radd__(self, other):
        raise TypeError('+ is ambiguous for AngularSeparation objects; not supported')

    def __sub__(self, other):
        raise TypeError('- is ambiguous for AngularSeparation objects; not supported')

    def __rsub__(self, other):
        raise TypeError('- is ambiguous for AngularSeparation objects; not supported')


#<----------------------------------Rotations---------------------------------->


def rotation_matrix(angle, axis='z', unit=None):
    """
    Generate a 3x3 cartesian rotation matrix in for rotation about
    a particular axis.

    Parameters
    ----------
    angle : convertible to Angle
        The amount of rotation this matrix should represent.

    axis : str or 3-sequence
        Either 'x','y', 'z', or a (x,y,z) specifying an axis to rotate
        about. If 'x','y', or 'z', the rotation sense is
        counterclockwise looking down the + axis (e.g. positive
        rotations obey left-hand-rule).

    unit : UnitBase, optional
        If `angle` does not have associated units, they are in this
        unit.  If neither are provided, it is assumed to be degrees.

    Returns
    -------
    rmat: `numpy.matrix`
        A unitary rotation matrix.
    """
    # TODO: This doesn't handle arrays of angles

    from numpy import sin, cos, radians, sqrt

    if unit is None:
        unit = u.degree

    angle = Angle(angle, unit=unit)

    if axis == 'z':
        s = sin(angle)
        c = cos(angle)
        return np.matrix(((c, s, 0),
                          (-s, c, 0),
                          (0, 0, 1)))
    elif axis == 'y':
        s = sin(angle)
        c = cos(angle)
        return np.matrix(((c, 0, -s),
                          (0, 1, 0),
                          (s, 0, c)))
    elif axis == 'x':
        s = sin(angle)
        c = cos(angle)
        return np.matrix(((1, 0, 0),
                          (0, c, s),
                          (0, -s, c)))
    else:
        x, y, z = axis
        w = cos(angle / 2)

        # normalize
        if w == 1:
            x = y = z = 0
        else:
            l = sqrt((x * x + y * y + z * z) / (1 - w * w))
            x /= l
            y /= l
            z /= l

        wsq = w * w
        xsq = x * x
        ysq = y * y
        zsq = z * z
        return np.matrix(((wsq + xsq - ysq - zsq, 2 * x * y - 2 * w * z, 2 * x * z + 2 * w * y),
                          (2 * x * y + 2 * w * z, wsq - xsq + ysq - zsq, 2 * y * z - 2 * w * x),
                          (2 * x * z - 2 * w * y, 2 * y * z + 2 * w * x, wsq - xsq - ysq + zsq)))


def angle_axis(matrix, unit=None):
    """
    Computes the angle of rotation and the rotation axis for a given rotation
    matrix.

    Parameters
    ----------
    matrix : array-like
        A 3 x 3 unitary rotation matrix.

    unit : UnitBase
        The output unit.  If `None`, the output unit is degrees.

    Returns
    -------
    angle : Angle
        The angle of rotation for this matrix.

    axis : array (length 3)
        The axis of rotation for this matrix.
    """
    # TODO: This doesn't handle arrays of angles

    from numpy import sin, cos, acos, degrees, sqrt

    m = np.asmatrix(matrix)
    if m.shape != (3, 3):
        raise ValueError('matrix is not 3x3')

    angle = acos((m[0, 0] + m[1, 1] + m[2, 2] - 1) / 2)
    denom = sqrt(2 * ((m[2, 1] - m[1, 2]) + (m[0, 2] - m[2, 0]) + (m[1, 0] - m[0, 1])))
    axis = np.array((m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1])) / denom
    axis /= sqrt(np.sum(axis ** 2))

    angle = Angle(angle, u.radian)
    if unit is None:
        unit = u.degree
    return angle.to(unit), axis
