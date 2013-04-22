# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

''' This module defines custom errors and exceptions used in astropy.coordinates.

    Note: Recently switched from a custom defined UnitsError to
    astropy.units.core.UnitsException. This will be renamed to
    astropy.units.core.UnitsError, and we'll have to rename
    everything in astropy.coordinates again...
'''

__all__ = ['RangeError', 'BoundsError', 'IllegalHourError',
           'IllegalMinuteError', 'IllegalSecondError', 'ConvertError',
           'IllegalHourWarning', 'IllegalMinuteWarning', 'IllegalSecondWarning']


class RangeError(ValueError):
    """
    Raised when some part of an angle is out of its valid range.
    """


class BoundsError(RangeError):
    """
    Raised when an angle is outside of its user-specified bounds.
    """


class IllegalHourError(RangeError):
    """
    Raised when an hour value is not in the range [0,24).

    Usage::

        if not 0 <= hr < 24:
        raise IllegalHourError(hour)

    Parameters
    ----------
    hour : int, float
    """
    def __init__(self, hour):
        self.hour = hour

    def __str__(self):
        return "An invalid value for 'hours' was found ('{0}'); must be in the range [0,24).".format(self.hour)


class IllegalHourWarning(Warning):
    """
    Raised when an hour value is 24.

    Parameters
    ----------
    hour : int, float
    """
    def __init__(self, hour):
        self.hour = hour

    def __str__(self):
        return "'hour' was found  to be '{0}'; it is best to in the range [0,24).".format(self.minute)


class IllegalMinuteError(RangeError):
    """
    Raised when an minute value is not in the range [0,60].

    Usage:
        if not 0 <= min < 60:
            raise IllegalMinuteError(minute)

    Parameters
    ----------
    minute : int, float
    """
    def __init__(self, minute):
        self.minute = minute

    def __str__(self):
        return "An invalid value for 'minute' was found ('{0}'); should be in the range [0,60).".format(self.minute)


class IllegalMinuteWarning(Warning):
    """
    Raised when a minute value is 60.

    Parameters
    ----------
    minute : int, float
    """
    def __init__(self, minute):
        self.minute = minute

    def __str__(self):
        return "'minute' was found  to be '{0}'; it is best to in the range [0,60).".format(self.minute)


class IllegalSecondError(RangeError):
    """
    Raised when an second value (time) is not in the range [0,60].

    Usage:
        if not 0 <= sec < 60:
            raise IllegalSecondError(second)

    Parameters
    ----------
    second : int, float
    """
    def __init__(self, second):
        self.second = second

    def __str__(self):
        return "An invalid value for 'second' was found ('{0}'); should be in the range [0,60).".format(self.second)


class IllegalSecondWarning(Warning):
    """
    Raised when a second value is 60.

    Parameters
    ----------
    second : int, float
    """
    def __init__(self, second):
        self.second = second

    def __str__(self):
        return "'second' was found  to be '{0}'; it is best to in the range [0,60).".format(self.second)


# TODO: consider if this should be used to `units`?
class UnitsError(ValueError):
    """
    Raised if units are missing or invalid.
    """


class ConvertError(Exception):
    """
    Raised if a coordinate system cannot be converted to another
    """
