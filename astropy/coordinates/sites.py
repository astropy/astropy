# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Observatories accessible by the `sites` module originate from the IRAF
Observatory Database, and are stored in ``data/observatories.json``.
Longitudes are listed with positive to the West.

Additions or corrections to the observatory list can be submitted via Pull
Request to the [astropy-data GitHub repository](https://github.com/astropy/astropy),
updating the ``location.json`` file.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
from difflib import get_close_matches

from ..utils.data import get_pkg_data_contents, get_file_contents
from ..extern import six
from .earth import EarthLocation

__all__ = ['get_site', 'get_site_names']

# Observatory database and list of names:
_builtin_site_dict = {'_site_db_uninitialized': True}
_builtin_site_names = []


def _parse_sites_json(jsondb, names, sitedict):
    for site in jsondb:
        location = EarthLocation.from_geodetic(jsondb[site]['longitude'],
                                               jsondb[site]['latitude'],
                                               jsondb[site]['elevation_meters'])
        names.append(jsondb[site]['name'])
        for alias in jsondb[site]['aliases']:
            sitedict[alias.lower()] = location

def _get_builtin_sites():
    """
    Load observatory database from data/observatories.json and parse them into
    a dictionary in memory.
    """
    if _builtin_site_dict.get('_site_db_uninitialized', False):
        jsondb = json.loads(get_pkg_data_contents('data/observatories.json'))
        _parse_sites_json(jsondb, _builtin_site_names, _builtin_site_dict)
        del _builtin_site_dict['_site_db_uninitialized']

    return _builtin_site_dict, _builtin_site_names

def _get_downloaded_sites():
    """
    Load observatory database from data.astropy.org and parse into data files.
    """
    sitedict = {}
    names = []

    jsondb = json.loads(get_file_contents('http://data.astropy.org/locations.json'))
    _parse_sites_json(jsondb, names, sitedict)

    return sitedict, names

def get_site(site_name, online):
    """
    Return an `~astropy.coordinates.EarthLocation` object for known observatory.

    Use `~astropy.coordinates.get_site_names` to see a list of available
    observatories.

    Parameters
    ----------
    site_name : str
        Name of the observatory.
    online : bool
        Use the online registry of observatories instead of the version included
        with astropy.  Requires an active internet connection.

    Returns
    -------
    `~astropy.coordinates.EarthLocation`
        The location of the observatory.

    See Also
    --------
    get_site_names : the list of sites that this function can access
    """
    if online:
        site_db, site_names = _get_downloaded_sites()
    else:
        site_db, site_names = _get_builtin_sites()

    if site_name.lower() not in site_db:
        # If site name not found, find close matches and suggest them in error
        close_names = get_close_matches(site_name, site_db)
        close_names = sorted(close_names, key=lambda x: len(x))
        if close_names:
            errmsg = ('Site not in database. Use ``get_site_names()`` '
                      'to see available sites. Did you mean one of: "{0}"?')
            errmsg = errmsg.format("', '".join(close_names))
        else:
            errmsg = 'Site "{0}" not in database.'.format(site_name)
        raise KeyError(errmsg)

    return site_db[site_name.lower()]


def get_site_names(show_aliases, online):
    """
    Get list of names of observatories for use with
    `~astropy.coordinates.get_site`.

    Parameters
    ----------
    show_aliases : bool
        If True, show the full list observatory names and aliases, or just the
        list of names if False.
    online : bool
        Use the online registry of observatories instead of the version included
        with astropy.  Requires an active internet connection.

    Returns
    -------
    names : list of str
        List of valid observatory names

    See Also
    --------
    get_site : gets the `~astropy.coordinates.EarthLocation` for one of the
               sites this returns.
    """
    if online:
        site_db, site_names = _get_downloaded_sites()
    else:
        site_db, site_names = _get_builtin_sites()

    if show_aliases:
        return sorted(site_db.keys())
    else:
        return sorted(site_names)
