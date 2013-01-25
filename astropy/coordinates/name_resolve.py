# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains convenience functions for getting a coordinate object
for a named object by querying SESAME and getting the first returned result.
Note that this is intended to be a convenience, and is very simple. If you
need precise coordinates for an object you should find the appropriate
reference for that measurement and input the coordinates manually.
"""

from __future__ import division, print_function

# Standard library
import os
import re
import sys
import httplib
import urllib
import urllib2
from urlparse import urlparse
import socket

# Third party
import numpy as np

# Astropy
from ..config import ConfigurationItem
from .builtin_systems import ICRSCoordinates
from .. import units as u

__all__ = ["get_icrs_coordinates"]

SESAME_URL = ConfigurationItem("sesame_url", 
                        ["http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/",
                         "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/"],
                        "The URL to Sesame's web-queryable database.",
                        cfgtype='string_list')

SESAME_DATABASE = ConfigurationItem("sesame_database", ['all', 'simbad', 'ned', 
                                    'vizier'],
                                    "This specifies the default database that "
                                    "SESAME will query when using the name "
                                    "resolve mechanism in the coordinates "
                                    "subpackage. Default is to search all "
                                    "databases, but this can be 'all', "
                                    "'simbad', 'ned', or 'vizier'.")
                                    
NAME_RESOLVE_TIMEOUT = ConfigurationItem('name_resolve_timeout', 5,
                                         "This is the maximum time to wait "
                                         "for a response from a name resolve "
                                         "query to SESAME in seconds.")

class NameResolveError(Exception):
    pass

def get_icrs_coordinates(name):
    """ 
    Retrieve an ICRSCoordinates object by using an online name resolving
    service to retrieve coordinates for the specified name. By default,
    this will search all available databases until a match is found. If
    you would like to specify the database, use the configuration item
    `name_resolve.SESAME_DATABASE` . You can also specify a list of servers 
    to use for querying Sesame using the configuration item 
    `name_resolve.SESAME_URL`. This will try each one in order until a valid
    response is returned. By default, this list includes the main Sesame 
    host and a mirror at vizier. A final configuration item, 
    `name_resolve.NAME_RESOLVE_TIMEOUT`, is the number of seconds to wait
    for a response from the server before giving up. By default this is 
    5 seconds.

    Parameters
    ----------
    name : str
        The name of the object to get coordinates for, e.g. M42.

    Returns
    -------
    coord : SphericalCoordinatesBase
        An `ICRSCoordinates` instance for the object name specified.
        
    """
    
    database = SESAME_DATABASE()
    # The web API just takes the first letter of the database name
    db = database.upper()[0]
    
    # Make sure we don't have duplicates in the url list
    urls = []
    domains = []
    for url in SESAME_URL():
        domain = urlparse(url).netloc
        
        # Check for duplicates
        if domain not in domains:
            domains.append(domain)
            
            # Add the query to the end of the url, add to url list
            fmt_url = os.path.join(url, "{db}?{name}")
            fmt_url = fmt_url.format(name=urllib.quote(name), db=db)
            urls.append(fmt_url)
    
    for url in urls:
        try:
            # Retrieve ascii name resolve data from CDS
            resp = urllib2.urlopen(url, timeout=NAME_RESOLVE_TIMEOUT())
            break
        except urllib2.URLError, e:
            # This catches a timeout error, see:
            #   http://stackoverflow.com/questions/2712524/handling-urllib2s-timeout-python
            if isinstance(e.reason, socket.timeout):
                # If it was a timeout, try with the next URL
                continue
            else:
                raise NameResolveError("Unable to retrieve coordinates for name "
                                       "'{0}'".format(name))
    
    # All Sesame URL's timed out...
    else:
        raise NameResolveError("All Sesame queries timed out. Unable to "
                               "retrieve coordinates.")

    resp_data = resp.read()

    pattr = re.compile(r"%J\s*([0-9\.]+)\s*([\+\-\.0-9]+)")
    matched = pattr.search(resp_data.decode('utf-8'))

    if matched == None:
        if db == "A":
            err = "Unable to find coordinates for name '{0}'".format(name)
        else:
            err = "Unable to find coordinates for name '{0}' in database {1}"\
                  .format(name, database)

        raise NameResolveError(err)

    ra,dec = matched.groups()

    return ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))