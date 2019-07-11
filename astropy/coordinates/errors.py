# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

''' This module defines custom errors and exceptions used in astropy.coordinates.
'''

from astropy.utils.exceptions import AstropyWarning

__all__ = ['RangeError', 'BoundsError', 'IllegalHourError',
           'IllegalMinuteError', 'IllegalSecondError', 'ConvertError',
           'IllegalHourWarning', 'IllegalMinuteWarning', 'IllegalSecondWarning',
           'UnknownSiteException']


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

    Parameters
    ----------
    hour : int, float

    Examples
    --------

    .. code-block:: python

        if not 0 <= hr < 24:
           raise IllegalHourError(hour)
    """
    def __init__(self, hour):
        self.hour = hour

    def __str__(self):
        return f"An invalid value for 'hours' was found ('{self.hour}'); must be in the range [0,24)."


class IllegalHourWarning(AstropyWarning):
    """
    Raised when an hour value is 24.

    Parameters
    ----------
    hour : int, float
    """
    def __init__(self, hour, alternativeactionstr=None):
        self.hour = hour
        self.alternativeactionstr = alternativeactionstr

    def __str__(self):
        message = f"'hour' was found  to be '{self.hour}', which is not in range (-24, 24)."
        if self.alternativeactionstr is not None:
            message += ' ' + self.alternativeactionstr
        return message


class IllegalMinuteError(RangeError):
    """
    Raised when an minute value is not in the range [0,60].

    Parameters
    ----------
    minute : int, float

    Examples
    --------

    .. code-block:: python

        if not 0 <= min < 60:
            raise IllegalMinuteError(minute)

    """
    def __init__(self, minute):
        self.minute = minute

    def __str__(self):
        return f"An invalid value for 'minute' was found ('{self.minute}'); should be in the range [0,60)."


class IllegalMinuteWarning(AstropyWarning):
    """
    Raised when a minute value is 60.

    Parameters
    ----------
    minute : int, float
    """
    def __init__(self, minute, alternativeactionstr=None):
        self.minute = minute
        self.alternativeactionstr = alternativeactionstr

    def __str__(self):
        message = f"'minute' was found  to be '{self.minute}', which is not in range [0,60)."
        if self.alternativeactionstr is not None:
            message += ' ' + self.alternativeactionstr
        return message


class IllegalSecondError(RangeError):
    """
    Raised when an second value (time) is not in the range [0,60].

    Parameters
    ----------
    second : int, float

    Examples
    --------

    .. code-block:: python

        if not 0 <= sec < 60:
            raise IllegalSecondError(second)
    """
    def __init__(self, second):
        self.second = second

    def __str__(self):
        return f"An invalid value for 'second' was found ('{self.second}'); should be in the range [0,60)."


class IllegalSecondWarning(AstropyWarning):
    """
    Raised when a second value is 60.

    Parameters
    ----------
    second : int, float
    """
    def __init__(self, second, alternativeactionstr=None):
        self.second = second
        self.alternativeactionstr = alternativeactionstr

    def __str__(self):
        message = f"'second' was found  to be '{self.second}', which is not in range [0,60)."
        if self.alternativeactionstr is not None:
            message += ' ' + self.alternativeactionstr
        return message


# TODO: consider if this should be used to `units`?
class UnitsError(ValueError):
    """
    Raised if units are missing or invalid.
    """


class ConvertError(Exception):
    """
    Raised if a coordinate system cannot be converted to another
    """


class UnknownSiteException(KeyError):
    def __init__(self, site, attribute, close_names=None):
        message = f"Site '{site}' not in database. Use {attribute} to see available sites."
        if close_names:
            message += " Did you mean one of: '{}'?'".format("', '".join(close_names))
        self.site = site
        self.attribute = attribute
        self.close_names = close_names
        return super().__init__(message)
