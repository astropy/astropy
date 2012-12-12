# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains...
"""

from __future__ import division, print_function

# Standard library
import re
import sys
import urllib
import urllib2

# Third party
import numpy as np

# Astropy
from .builtin_systems import ICRSCoordinates
from .. import units as u

__all__ = ["get_icrs_coordinates"]

sesame_url = "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame?{name}"

def get_icrs_coordinates(name):
    """ """

    try:
        # Retrieve ascii name resolve data from CDS
        resp = urllib2.urlopen(sesame_url.format(name=urllib.quote(name)))
    except urllib2.URLError:
        raise OSError("Unable to connect to name resolve web service. Check your internet connection, and try again.")

    resp_data = resp.read()

    pattr = re.compile(r"%J\s*([0-9\.]+)\s*([\+\-\.0-9]+)")
    matched = pattr.search(resp_data)

    if matched == None:
        raise ValueError("Unable to find coordinates for name '{0}'".format(name))

    ra,dec = matched.groups()

    return ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))