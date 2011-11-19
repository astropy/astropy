# Licensed under a 3-clause BSD style license - see LICENSE.rst  
"""This module contains functions to standardize access to configuration files
for astropy and affiliated modules.

.. note::
    The configuration system makes use of the 'configobj' pacakge, which stores
    configuration in a text format like that used in the standard library
    `ConfigParser`. More information and documentation for confobj can be found
    at `http://www.voidspace.org.uk/python/configobj.html`_.
"""

from __future__ import division
from ..extern.configobj import configobj

__all__ = []