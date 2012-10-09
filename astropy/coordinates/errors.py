# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

''' This module defines custom errors and exceptions used in astropy.coordinates. '''

__all__ = ['RangeError', 'IllegalUnitsError', 'IllegalHourError', 'IllegalMinuteError', 'IllegalSecondError']

class RangeError(Exception):
	pass

class IllegalUnitsError(Exception):
    """
    Usage:
        if units not in VALIDUNITS:
            raise IllegalUnitsError("")
    """
    def __init__(self, units):
        self.units = units
    def __str__(self):
        return "The units specified must be one of the following: {0}. You specified: \'{1}\'".format(",".join(globals.VALIDUNITS), self.units)

	
class IllegalHourError(Exception):
    """
    Usage:
        if not 0 <= hr < 24:
            raise IllegalHourError(hour)
    """
    def __init__(self, hour):
        self.hour = hour
    def __str__(self):
        return "An invalid value for 'hours' was found ('{0}'); must be in the range [0,24).".format(self.second)

class IllegalMinuteError(Exception):
    """
    Usage:
        if not 0 <= min < 60:
            raise IllegalMinuteError(minute)
    """
    def __init__(self, minute):
        self.minute = minute
    def __str__(self):
        return "An invalid value for 'minute' was found ('{0}'); must be in the range [0,60).".format(self.minute)

class IllegalSecondError(Exception):
    """
    Usage:
        if not 0 <= sec < 60:
            raise IllegalSecondError(second)
    """
    def __init__(self, second):
        self.second = second
    def __str__(self):
        return "An invalid value for 'seconds' was found ('{0}'); must be in the range [0,60).".format(self.second)
        
class CoordinatesConversionError(Exception):
	"""
	"""
	pass
