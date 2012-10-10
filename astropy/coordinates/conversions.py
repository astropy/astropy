#!/usr/bin/python
# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains utility functions that are for internal use of
astropy.coordinates.core. Mainly they are conversions from one format
of data to another.
"""
import re
import math
import datetime as py_datetime

import core
from errors import *

__all__ = ['parseDegrees']

def _checkHourRange(hrs):
    ''' Checks that the given value is in the range (-24,24). '''
    if not -24 < hrs < 24:
        raise IllegalHourError(hrs)

def _checkMinuteRange(min):
    ''' Checks that the given value is in the range [0,60). '''
    if not 0 <= min < 60:
        raise IllegalMinuteError(min) #"Error: minutes not in range [0,60) ({0}).".format(min))

def _checkSecondRange(sec):
    ''' Checks that the given value is in the range [0,60). '''
    if not 0 <= sec < 60:
        raise IllegalSecondError(sec)#"Error: seconds not in range [0,60) ({0}).".format(sec))

def checkHMSRanges(h, m, s):
    _checkHourRange(h)
    _checkMinuteRange(m)
    _checkSecondRange(s)
    return None

def parseDegrees(degrees, outputDMS=False):
    """ Parses an input "degrees" value into decimal degrees or a 
        degree,arcminute,arcsecond tuple.
        
        Convert degrees given in any parseable format (float, string, or Angle) into 
        degrees, arcminutes, and arcseconds components or decimal degrees.
        
        Parameters
        ----------
        degrees : float, int, str
            If a string, accepts values in these formats:
                * [+|-]DD:MM:SS.sss (string), e.g. +12:52:32.423 or -12:52:32.423
                * DD.dddd (float, string), e.g. 12.542326
                * DD MM SS.sss (string, array), e.g. +12 52 32.423
            Whitespace may be spaces and/or tabs.
        outputDMS : bool
            If True, returns a tuple of (degree, arcminute, arcsecond)
    
        Returns degrees in decimal form unless the keyword "outputDMS" is True, in which
        case returns a tuple: (d, m, s).
        
    """
    
    # either a string or a float
    x = degrees
    
    if isinstance(x, float) or isinstance(x, int):
        parsedDegrees = float(x)
        #parsedDMS = degreesToDMS(parsedDegrees)

    elif isinstance(x, str):
        x = x.strip()
        
        string_parsed = False
        
        # See if the string is just a float or int value.
        try:
            parsedDegrees = float(x)
            string_parsed = True
        except ValueError:
            pass

        if not string_parsed:
            # Look for a pattern where d,m,s is specified
            div = '[:|/|\t|\-|\sDdMmSs]{1,2}' # accept these as (one or more repeated) delimiters: :, whitespace, /
            pattr = '^([+-]{0,1}\d{1,3})' + div + '(\d{1,2})' + div + '(\d{1,2}[\.0-9]+)' + '[Ss]{0,1}' + '$'
    
            try:
                elems = re.search(pattr, x).groups()
                parsedDegrees = dmsToDegrees(int(elems[0]), int(elems[1]), float(elems[2]))
                string_parsed = True
            except AttributeError:
                # regular expression did not match - try again below
                # make sure to let things like IllegalMinuteError, etc. through
                pass
                
        if not string_parsed:
            # look for a pattern where only d,m is specified
            pattr = '^([+-]{0,1}\d{1,3})' + div + '(\d{1,2})' + '[Mm]{0,1}' + '$'
            try:
                elems = re.search(pattr, x).groups()
                parsedDegrees = dmsToDegrees(int(elems[0]), int(elems[1]), 0.0)
                string_parsed = True
            except AttributeError:
                # regular expression did not match - try again below
                # make sure to let things like IllegalMinuteError, etc. through
                pass

        if not string_parsed:
            # look for a '°' symbol
            for unitStr in ["degrees", "degree", "deg", "d", "°"]:
                x = x.replace(unitStr, '')
                try:
                    parsedDegrees = float(x)
                    string_parsed = True
                except ValueError:
                    pass
                    
        if not string_parsed:
            raise ValueError("convert.parseDegrees: Invalid input string! ('{0}')".format(x))

    elif isinstance(x, core.Angle):
        parsedDegrees = x.degrees
        #parsedDMS = degreesToDMS(parsedDegrees)
        
    elif isinstance(x, tuple):
        parsedDegrees = dmsToDegrees(*x)
        #parsedDMS = x
        
    else:
        raise ValueError("convert.parseDegrees: could not parse value of {0}.".format(type(x)))
    
    return degreesToDMS(parsedDegrees) if outputDMS else parsedDegrees

def parseHours(hours, outputHMS=False):
    """ Parses an input "hour" value to decimal hours or an hour, minute, second tuple.
        
        Convert hours given in any parseable format (float, string, tuple, list, or Angle) into 
        hour, minute, and seconds components or decimal hours.
        
        Parameters
        ----------
        hours : float, str
            If a string, accepts values in these formats:
                * HH:MM:SS.sss (string), e.g. 12:52:32.423
                * HH.dddd (float, string), e.g. 12.542326
                * HH MM SS.sss (string, array), e.g. 12 52 32.423
            Whitespace may be spaces and/or tabs.
        outputHMS : bool
            If True, returns a tuple of (hour, minute, second)

    """
    
    # either a string or a float
    x = hours

    if isinstance(x, float) or isinstance(x, int):
        parsedHours = x
        parsedHMS = hoursToHMS(parsedHours)
    
    elif isinstance(x, str):
        x = x.strip()
        
        try:
            parsedHours = float(x)
            parsedHMS = hoursToHMS(parsedHours)
        except ValueError:

            string_parsed = False
            div = '[:|/|\t|\-|\sHhMmSs]{1,2}' # accept these as (one or more repeated) delimiters: :, whitespace, /
            
            # First look for a pattern where h,m,s is specified
            pattr = '^([+-]{0,1}\d{1,2})' + div + '(\d{1,2})' + div + '(\d{1,2}[\.0-9]+)' + '[Ss]{0,1}' + '$'

            try:
                elems = re.search(pattr, x).groups()
                string_parsed = True
            except:
                pass # try again below
                #raise ValueError("convert.parseHours: Invalid input string, can't parse to HMS. ({0})".format(x))
            
            if string_parsed:
                h, m, s = float(elems[0]), int(elems[1]), float(elems[2])
                parsedHours = hmsToHours(h, m, s)
                parsedHMS = (h, m, s)
            
            else:
                
                # look for a pattern where only d,m is specified
                pattr = '^([+-]{0,1}\d{1,2})' + div + '(\d{1,2})' + '[Mm]{0,1}' + '$'
                
                try:
                    elems = re.search(pattr, x).groups()
                    string_parsed = True
                except:
                    raise ValueError("convert.parseHours: Invalid input string, can't parse to HMS. ({0})".format(x))
                h, m, s = float(elems[0]), int(elems[1]), 0.0
                parsedHours = hmsToHours(h, m, s)
                parsedHMS = (h, m, s)

    elif isinstance(x, core.Angle):
        parsedHours = x.hours
        parsedHMS = hoursToHMS(parsedHours)
    
    elif isinstance(x, py_datetime.datetime):
        parsedHours = datetimeToDecimalTime(x)
        parsedHMS = hoursToHMS(parsedHours)
        
    elif isinstance(x, tuple):
        if len(x) == 3:
            parsedHours = hmsToHours(*x)
            parsedHMS = x
        else:
            raise ValueError("{0}.{1}: Incorrect number of values given, expected (h,m,s), got: {2}".format(os.path.basename(__file__), stack()[0][3], x))

    elif isinstance(x, list):
        if len(x) == 3:
            try:
                h = float(x[0])
                m = float(x[1])
                s = float(x[2])
            except ValueError:
                raise ValueError("{0}.{1}: Array values ([h,m,s] expected) could not be coerced into floats. {2}".format(os.path.basename(__file__), stack()[0][3], x))

            parsedHours = hmsToHours(h, m, s)
            parsedHMS = (h, m, s)
            if outputHMS:
                return (h, m, s)
            else:
                return hmsToHours(h, m, s)

        else:
            raise ValueError("{0}.{1}: Array given must contain exactly three elements ([h,m,s]), provided: {2}".format(os.path.basename(__file__), stack()[0][3], x)) # current filename/method should be made into a convenience method
    
    else:
        raise ValueError("parseHours: could not parse value of type {0}.".format(type(x).__name__))
    
    if outputHMS:
        return parsedHMS
    else:
        return parsedHours

def parseRadians(radians):
    """ Parses an input "radians" value into a float number.
        
        Convert radians given in any parseable format (float or Angle) into float radians.
        
        ..Note::
            This function is mostly for consistency with the other "parse" functions, like
            parseHours and parseDegrees. 
        
        Parameters
        ----------
        radians : float, int, Angle
            The input angle.
    """
    x = radians
    
    if type(x) in [float, int]:
        return float(x)
    elif isinstance(x, core.Angle):
        return x.radians
    else:
        raise ValueError("parseRadians: could not parse value of type {0}.".format(type(x).__name__))

def degreesToDMS(d):
    """ Convert any parseable degree value (see: parseDegrees) into a 
        degree,arcminute,arcsecond tuple 
    """    
    sign = math.copysign(1.0, d)
        
    (df, d) = math.modf(abs(d)) # (degree fraction, degree)
    (mf, m) = math.modf(df * 60.) # (minute fraction, minute)
    s = mf * 60.
    
    _checkMinuteRange(m)
    _checkSecondRange(s)

    return (float(sign * d), int(m), s)

def dmsToDegrees(d, m, s):
    """ Convert degrees, arcminute, arcsecond to a float degrees value. """

    _checkMinuteRange(m)
    _checkSecondRange(s)

   # determine sign
    sign = math.copysign(1.0, d)

    try:
        d = int(abs(d))
        m = int(abs(m))
        s = float(abs(s))
    except ValueError:
        raise ValueError("convert.dmsToDegrees: dms values ({0[0]},{0[1]},{0[2]}) could not be converted to numbers.".format(d,m,s))

    return sign * (d + m/60. + s/3600.)

def hmsToHours(h, m, s):
    """ Convert hour, minute, second to a float hour value. """

    checkHMSRanges(h, m, s);

    try:
        h = int(h)
        m = int(m)
        s = float(s)
    except ValueError:
        raise ValueError("convert.HMStoHours: HMS values ({0[0]},{0[1]},{0[2]}) could not be converted to numbers.".format(h,m,s))
    return h + m/60. + s/3600.

def hmsToDegrees(h, m, s):
    """ Convert hour, minute, second to a float degrees value. """
    return hmsToHours(h, m, s)*15.

def hmsToRadians(h, m, s):
    """ Convert hour, minute, second to a float radians value. """
    return math.radians(hmsToDegrees(h, m, s))

def hmsToDMS(h, m, s):
    """ Convert degrees, arcminutes, arcseconds to an hour, minute, second tuple """
    return degreesToDMS(hmsToDegrees(h, m, s))

def hoursToDecimal(h):
    """ Convert any parseable hour value (see: parseHours) into a float value """
    return parseHours(h, outputHMS=False)

def hoursToRadians(h):
    """ Convert an angle in Hours to Radians """
    return math.radians(h*15.)

def hoursToHMS(h):
    """ Convert any parseable hour value (see: parseHours) into an hour,minute,second tuple """
    sign = math.copysign(1.0, h)
        
    (hf, h) = math.modf(abs(h)) # (degree fraction, degree)
    (mf, m) = math.modf(hf * 60.) # (minute fraction, minute)
    s = mf * 60.
    
    checkHMSRanges(h,m,s) # throws exception if out of range
    
    return (float(sign*h), int(m), s)

def radiansToDegrees(r):
    """ Convert an angle in Radians to Degrees """
    try:
        r = float(r)
    except ValueError:
        raise ValueError("convert.radiansToHours: degree value ({0[0]}) could not be converted to a float.".format(r))
    return math.degrees(r)

def radiansToHours(r):
    """ Convert an angle in Radians to Hours """
    return radiansToDegrees(r) / 15.

def radiansToHMS(r):
    """ Convert an angle in Radians to an hour,minute,second tuple """
    hours = radiansToHours(r)
    return hoursToHMS(hours)
 
def radiansToDMS(r):
    """ Convert an angle in Radians to an degree,arcminute,arcsecond tuple """
    degrees = math.degrees(r)
    return degreesToDMS(degrees)  

def hoursToString(h, precision=5, pad=False, sep=("h", "m", "s")):
    """ Takes a decimal hour value and returns a string formatted as hms with separator
        specified by the 'sep' parameter. 
        
        More detailed description here!
    """
    if pad:
        hPad = ":02"
    else:
        hPad = ""
        
    if len(sep) == 1:
        literal = "{0" + hPad + "}"+ str(sep) + "{1:02d}" + str(sep) + "{2:0" + str(precision+3) + "." + str(precision) + "f}"
    elif len(sep) == 2:
        literal = "{0" + hPad + "}"+ str(sep[0]) + "{1:02d}" + str(sep[1]) + "{2:0" + str(precision+3) + "." + str(precision) + "f}"
    elif len(sep) == 3:
        literal = "{0" + hPad + "}"+ str(sep[0]) + "{1:02d}" + str(sep[1]) + "{2:0" + str(precision+3) + "." + str(precision) + "f}" + str(sep[2])
    else:
        raise ValueError("Invalid separator specification for converting angle to string.")
    
    (h,m,s) = hoursToHMS(h)
    h = "-{0}".format(int(h)) if math.copysign(1,h) == -1 else int(h)
    return literal.format(h,m,s)

def degreesToString(d, precision=5, pad=False, sep=":"):
    """ Takes a decimal hour value and returns a string formatted as dms with separator
        specified by the 'sep' parameter. 
    """
    if pad:
        dPad = ":02"
    else:
        dPad = ""
        
    if len(sep) == 1:
        literal = "{0" + dPad + "}" + str(sep) + "{1:02d}" + str(sep) + "{2:0" + str(precision+3) + "." + str(precision) + "f}"
    elif len(sep) == 2:
        literal = "{0" + dPad + "}"+ str(sep[0]) + "{1:02d}" + str(sep[1]) + "{2:0" + str(precision+3) + "." + str(precision) + "f}"
    elif len(sep) == 3:
        literal = "{0" + dPad + "}"+ str(sep[0]) + "{1:02d}" + str(sep[1]) + "{2:0" + str(precision+3) + "." + str(precision) + "f}" + str(sep[2])
    else:
        raise ValueError("Invalid separator specification for converting angle to string.")

    d,m,s = degreesToDMS(d)
    d = "-{0}".format(int(d)) if math.copysign(1,d) == -1 else int(d)
    return literal.format(d,m,s)


