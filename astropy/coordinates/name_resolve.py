# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains convenience functions for getting a coordinate object
for a named object by querying SESAME and getting the first returned result.
Note that this is intended to be a convenience, and is very simple. If you
need precise coordinates for an object you should find the appropriate
reference for that measurement and input the coordinates manually.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
import os
import re
import socket

# Astropy
from ..extern.six.moves import urllib
from .. import units as u
from .sky_coordinate import SkyCoord
from ..utils import data
from ..utils.state import ScienceState

__all__ = ["get_icrs_coordinates"]


class sesame_url(ScienceState):
    """
    The URL(s) to Sesame's web-queryable database.
    """
    _value = ["http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/",
              "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/"]

    @classmethod
    def validate(cls, value):
        # TODO: Implement me
        return value


class sesame_database(ScienceState):
    """
    This specifies the default database that SESAME will query when
    using the name resolve mechanism in the coordinates
    subpackage. Default is to search all databases, but this can be
    'all', 'simbad', 'ned', or 'vizier'.
    """
    _value = 'all'

    @classmethod
    def validate(cls, value):
        if value not in ['all', 'simbad', 'ned', 'vizier']:
            raise ValueError("Unknown database '{0}'".format(value))
        return value


class NameResolveError(Exception):
    pass


def _parse_response(resp_data):
    """
    Given a string response from SESAME, parse out the coordinates by looking
    for a line starting with a J, meaning ICRS J2000 coordinates.

    Parameters
    ----------
    resp_data : str
        The string HTTP response from SESAME.

    Returns
    -------
    ra : str
        The string Right Ascension parsed from the HTTP response.
    dec : str
        The string Declination parsed from the HTTP response.
    """

    pattr = re.compile(r"%J\s*([0-9\.]+)\s*([\+\-\.0-9]+)")
    matched = pattr.search(resp_data.decode('utf-8'))

    if matched is None:
        return None, None
    else:
        ra, dec = matched.groups()
        return ra, dec


def get_icrs_coordinates(name):
    """
    Retrieve an ICRS object by using an online name resolving service to
    retrieve coordinates for the specified name. By default, this will
    search all available databases until a match is found. If you would like
    to specify the database, use the science state
    ``astropy.coordinates.name_resolve.sesame_database``. You can also
    specify a list of servers to use for querying Sesame using the science
    state ``astropy.coordinates.name_resolve.sesame_url``. This will try
    each one in order until a valid response is returned. By default, this
    list includes the main Sesame host and a mirror at vizier.  The
    configuration item `astropy.utils.data.Conf.remote_timeout` controls the
    number of seconds to wait for a response from the server before giving
    up.

    Parameters
    ----------
    name : str
        The name of the object to get coordinates for, e.g. ``'M42'``.

    Returns
    -------
    coord : `astropy.coordinates.ICRS` object
        The object's coordinates in the ICRS frame.

    """

    database = sesame_database.get()
    # The web API just takes the first letter of the database name
    db = database.upper()[0]

    # Make sure we don't have duplicates in the url list
    urls = []
    domains = []
    for url in sesame_url.get():
        domain = urllib.parse.urlparse(url).netloc

        # Check for duplicates
        if domain not in domains:
            domains.append(domain)

            # Add the query to the end of the url, add to url list
            fmt_url = os.path.join(url, "{db}?{name}")
            fmt_url = fmt_url.format(name=urllib.parse.quote(name), db=db)
            urls.append(fmt_url)

    exceptions = []
    for url in urls:
        try:
            # Retrieve ascii name resolve data from CDS
            resp = urllib.request.urlopen(url, timeout=data.conf.remote_timeout)
            resp_data = resp.read()
            break
        except urllib.error.URLError as e:
            exceptions.append(e)
            continue
        except socket.timeout as e:
            # There are some cases where urllib2 does not catch socket.timeout
            # especially while receiving response data on an already previously
            # working request
            e.reason = "Request took longer than the allowed {:.1f} " \
                       "seconds".format(data.conf.remote_timeout)
            exceptions.append(e)
            continue

    # All Sesame URL's failed...
    else:
        messages = ["{url}: {e.reason}".format(url=url, e=e)
                    for url, e in zip(urls, exceptions)]
        raise NameResolveError("All Sesame queries failed. Unable to "
                               "retrieve coordinates. See errors per URL "
                               "below: \n {}".format("\n".join(messages)))

    ra, dec = _parse_response(resp_data)

    if ra is None and dec is None:
        if db == "A":
            err = "Unable to find coordinates for name '{0}'".format(name)
        else:
            err = "Unable to find coordinates for name '{0}' in database {1}"\
                  .format(name, database)

        raise NameResolveError(err)

    # Return SkyCoord object
    sc = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    return sc
