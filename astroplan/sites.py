# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Observatories accessible by the `sites` module originate from the IRAF
Observatory Database, and are stored in astroplan/data/observatories.json.
Longitudes are listed with positive to the West.

Additions and corrections to the observatory list can be submitted via Pull
Request to the [astroplan GitHub repository](https://github.com/astroplanners/astroplan).

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.utils.data import get_pkg_data_contents
from difflib import get_close_matches
import json
from astropy.coordinates import EarthLocation

__all__ = ['get_site', 'get_site_names', 'add_site']

# Observatory database and list of names:
_site_db = None
_site_names = []

def _load_sites():
    """
    Load observatory database from astroplan/data/observatories.json
    """
    global _site_db, _site_names
    _site_db = dict()
    db = json.loads(get_pkg_data_contents('data/observatories.json'))
    for site in db:
        location = EarthLocation.from_geodetic(db[site]['longitude'],
                                   db[site]['latitude'],
                                   db[site]['elevation_meters'])
        _site_names.append(db[site]['name'])
        for alias in db[site]['aliases']:
            _site_db[alias.lower()] = location

def get_site(site_name):
    """
    Return an `~astropy.coordinates.EarthLocation` object for known observatory.

    Use `~astroplan.core.get_site_names` to see a list of available
    observatories.

    Parameters
    ----------
    site_name : str
        Name of the observatory.

    Returns
    -------
    `~astropy.coordinates.EarthLocation`
        The location of the observatory.
    """
    if _site_db is None:
        _load_sites()

    if site_name.lower() not in _site_db.keys():
        # If site name not found, find close matches and suggest them in error
        close_names = get_close_matches(site_name, _site_db.keys())
        close_names = sorted(close_names, key=lambda x: len(x))
        if len(close_names) > 0:
            errmsg = ("Site not in database. Use `astroplan.get_site_names()` "
                      "to see available sites. Did you mean: '{}'?".format(
                      "', '".join(close_names)))
        else:
            errmsg = 'Site not in database.'
        raise KeyError(errmsg)

    return _site_db[site_name.lower()]

def get_site_names(full_list=True):
    """
    Get list of names of observatories for use with `~astroplan.core.get_site`.

    Parameters
    ----------
    full_list : bool
        Show full list observatory names and aliases (True), or just the list
        of names (False)? Default to True.

    Returns
    -------
    List of observatory names (strings)
    """
    if _site_db is None:
        _load_sites()

    if full_list:
        return sorted(_site_db.keys())
    else:
        return sorted(_site_names)

def add_site(site_name, location):
    """
    Add a site to the list of available observatories.

    Parameters
    ----------
    site_name : string
        Name of the observatory

    location : `~astropy.coordinates.EarthLocation`
        Location of the observatory
    """
    if _site_db is None:
        _load_sites()

    if not isinstance(location, EarthLocation):
        raise ValueError('Location must be an EarthLocation.')

    if site_name.lower() not in _site_db.keys():
        _site_db[site_name.lower()] = location
        _site_names.append(site_name)
    else:
        raise KeyError('The site "{}" already exists at (longitude,latitude,'
                       'elevation)={}'.format(site_name,
                                     _site_db[site_name.lower()].to_geodetic()))
