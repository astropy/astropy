# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

''' This module defines custom errors and exceptions used in astropy.coordinates. '''

__all__ = ['RangeError', 'IllegalUnitsError']

class RangeError(Exception):
	pass

class IllegalUnitsError(Exception):
	pass