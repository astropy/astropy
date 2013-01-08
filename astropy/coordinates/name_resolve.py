# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains...
"""

from __future__ import division, print_function

# Standard library
import re
import sys
import httplib
import urllib
import urllib2

# Third party
import numpy as np

# Astropy
from .builtin_systems import ICRSCoordinates
from .. import units as u

__all__ = ["get_icrs_coordinates"]

SESAME_URL = "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/{db}?{name}"
MIRROR_URL = "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/{db}?{name}"
ALLOWED_DATABASES = ['ned', 'simbad', 'vizier', 'all']
NAME_RESOLVE_TIMEOUT = 30 # seconds

class NameResolveException(Exception):
    pass

def get_icrs_coordinates(name, database='all'):
    """ Retrieve an ICRSCoordinates object by using an online name resolving
        service to retrieve coordinates for the specified name.

        Parameters
        ----------
        name : str
            The name of the object to get coordinates for, e.g. m42.
        database : str (optional)
            Specify which database to search. Can be 'ned', 'simbad', 'vizier', or 'all.'

        Returns
        -------
        coord : SphericalCoordinatesBase
            An `ICRSCoordinates` instance for the object name specified.
    """

    if database.lower() not in ALLOWED_DATABASES:
        raise ValueError("Invalid database name: {0}. Allowed databases: {1}".format(database, ",".join(ALLOWED_DATABASES)))

    # The web API just takes the first letter of the database name
    db = database.upper()[0]
    url = SESAME_URL.format(name=urllib.quote(name), db=db)
    try:
        # Retrieve ascii name resolve data from CDS
        resp = urllib2.urlopen(url, timeout=NAME_RESOLVE_TIMEOUT)
    except:
        raise NameResolveException("Unable to retrieve coordinates for name '{0}'".format(name))

    resp_data = resp.read()

    pattr = re.compile(r"%J\s*([0-9\.]+)\s*([\+\-\.0-9]+)")
    matched = pattr.search(resp_data.decode('utf-8'))

    if matched == None:
        if db == "A":
            err = "Unable to find coordinates for name '{0}'".format(name)
        else:
            err = "Unable to find coordinates for name '{0}' in database {1}".format(name, database)

        raise NameResolveException(err)

    ra,dec = matched.groups()

    return ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))