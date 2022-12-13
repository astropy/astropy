# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains convenience functions for getting a coordinate object
for a named object by querying SESAME and getting the first returned result.
Note that this is intended to be a convenience, and is very simple. If you
need precise coordinates for an object you should find the appropriate
reference for that measurement and input the coordinates manually.
"""

# Standard library
import os
import re
import socket
import urllib.error
import urllib.parse
import urllib.request

# Astropy
from astropy import units as u
from astropy.utils import data
from astropy.utils.data import download_file, get_file_contents
from astropy.utils.state import ScienceState

from .sky_coordinate import SkyCoord

__all__ = ["get_icrs_coordinates"]


class sesame_url(ScienceState):
    """
    The URL(s) to Sesame's web-queryable database.
    """

    _value = [
        "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/",
        "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/",
    ]

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

    _value = "all"

    @classmethod
    def validate(cls, value):
        if value not in ["all", "simbad", "ned", "vizier"]:
            raise ValueError(f"Unknown database '{value}'")
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
    matched = pattr.search(resp_data)

    if matched is None:
        return None, None
    else:
        ra, dec = matched.groups()
        return ra, dec


def get_icrs_coordinates(name, parse=False, cache=False):
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
    parse : bool
        Whether to attempt extracting the coordinates from the name by
        parsing with a regex. For objects catalog names that have
        J-coordinates embedded in their names eg:
        'CRTS SSS100805 J194428-420209', this may be much faster than a
        sesame query for the same object name. The coordinates extracted
        in this way may differ from the database coordinates by a few
        deci-arcseconds, so only use this option if you do not need
        sub-arcsecond accuracy for coordinates.
    cache : bool, str, optional
        Determines whether to cache the results or not. Passed through to
        `~astropy.utils.data.download_file`, so pass "update" to update the
        cached value.

    Returns
    -------
    coord : `astropy.coordinates.ICRS` object
        The object's coordinates in the ICRS frame.

    """

    # if requested, first try extract coordinates embedded in the object name.
    # Do this first since it may be much faster than doing the sesame query
    if parse:
        from . import jparser

        if jparser.search(name):
            return jparser.to_skycoord(name)
        else:
            # if the parser failed, fall back to sesame query.
            pass
            # maybe emit a warning instead of silently falling back to sesame?

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
            resp_data = get_file_contents(
                download_file(url, cache=cache, show_progress=False)
            )
            break
        except urllib.error.URLError as e:
            exceptions.append(e)
            continue
        except socket.timeout as e:
            # There are some cases where urllib2 does not catch socket.timeout
            # especially while receiving response data on an already previously
            # working request
            e.reason = (
                "Request took longer than the allowed "
                f"{data.conf.remote_timeout:.1f} seconds"
            )
            exceptions.append(e)
            continue

    # All Sesame URL's failed...
    else:
        messages = [f"{url}: {e.reason}" for url, e in zip(urls, exceptions)]
        raise NameResolveError(
            "All Sesame queries failed. Unable to retrieve coordinates. See errors per"
            f" URL below: \n {os.linesep.join(messages)}"
        )

    ra, dec = _parse_response(resp_data)

    if ra is None or dec is None:
        if db == "A":
            err = f"Unable to find coordinates for name '{name}' using {url}"
        else:
            err = (
                f"Unable to find coordinates for name '{name}' in database"
                f" {database} using {url}"
            )

        raise NameResolveError(err)

    # Return SkyCoord object
    sc = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs")
    return sc
