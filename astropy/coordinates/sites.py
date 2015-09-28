# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Observatories accessible by the `sites` module originate from the IRAF
Observatory Database, and are stored in ``data/observatories.json``.
Longitudes are listed with positive to the West.

Additions or corrections to the observatory list can be submitted via Pull
Request to the [astropy GitHub repository](https://github.com/astropy/astropy),
updating the ``observatories.json`` file.

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
from difflib import get_close_matches

from ..utils.data import get_pkg_data_contents
from ..extern import six
from .earth import EarthLocation

__all__ = ['get_site', 'get_site_names', 'add_site', 'remove_site']

# Observatory database and list of names:
_site_db = None
_site_names = []


def _load_sites():
    """
    Load observatory database from data/observatories.json and parse them into
    a dictionary in memory.
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

    Use `~astropy.coordinates.get_site_names` to see a list of available
    observatories.

    Parameters
    ----------
    site_name : str
        Name of the observatory.

    Returns
    -------
    `~astropy.coordinates.EarthLocation`
        The location of the observatory.

    See Also
    --------
    get_site_names : the list of sites that this function can access
    add_site : adds a new site to the list
    remove_site : remove a site by name from the list of sites
    """
    if _site_db is None:
        _load_sites()

    if site_name.lower() not in _site_db:
        # If site name not found, find close matches and suggest them in error
        close_names = get_close_matches(site_name, _site_db)
        close_names = sorted(close_names, key=lambda x: len(x))
        if close_names:
            errmsg = ('Site not in database. Use ``get_site_names()`` '
                      'to see available sites. Did you mean one of: "{0}"?')
            errmsg = errmsg.format("', '".join(close_names))
        else:
            errmsg = 'Site "{0}" not in database.'.format(site_name)
        raise KeyError(errmsg)

    return _site_db[site_name.lower()]


def get_site_names(show_aliases=True):
    """
    Get list of names of observatories for use with
    `~astropy.coordinates.get_site`.

    Parameters
    ----------
    show_aliases : bool
        If True, show the full list observatory names and aliases, or just the
        list of names if False.

    Returns
    -------
    names : list of str
        List of valid observatory names

    See Also
    --------
    get_site : gets the `~astropy.coordinates.EarthLocation` for one of the
               sites this returns.
    get_site_names : the list of sites names that this function can access
    remove_site : remove a site by name from the list of sites

    """
    if _site_db is None:
        _load_sites()

    if show_aliases:
        return sorted(_site_db.keys())
    else:
        return sorted(_site_names)


def add_site(site_names, location):
    """
    Add a site to the list of available observatories.

    Parameters
    ----------
    site_names : string or list of strings
        Name of the observatory.  If a list, all names will be aliases for the
        same ``location``, and the first will be taken as the canonical name.
    location : `~astropy.coordinates.EarthLocation`
        Location of the observatory

    See Also
    --------
    get_site : gets the `~astropy.coordinates.EarthLocation` for one of the
               sites this returns.
    get_site_names : the list of site names that this function will add to
    remove_site : remove a site by name from the list of sites

    """
    if _site_db is None:
        _load_sites()

    if not isinstance(location, EarthLocation):
        raise ValueError('Location must be an EarthLocation.')

    if isinstance(site_names, six.string_types):
        site_names = [site_names]

    firstnamedone = False
    for name in site_names:
        if name.lower() not in _site_db.keys():
            _site_db[name.lower()] = location
            if not firstnamedone:
                _site_names.append(name)
                firstnamedone = True
        else:
            raise KeyError('The site "{0}" already exists at (longitude,latitude,'
                           'elevation)={1}'.format(name, _site_db[name.lower()].to_geodetic()))


def remove_site(site_name):
    """
    Removes a site (and all its alias) from the list of available observatories.

    Parameters
    ----------
    site_name : string
        Name of the observatory (or one of its aliases) to remove.

    Raises
    ------
    KeyError
        If the name is not in the database

    See Also
    --------
    get_site : gets the `~astropy.coordinates.EarthLocation` for one of the
               sites this returns.
    get_site_names : the list of site names that this function will add to
    add_site : adds a new site to the list
    """
    lname = site_name.lower()
    if lname in _site_db:
        remloc = _site_db.pop(lname)

        # now take care of aliases
        namestorem = []
        for name, loc in six.iteritems(_site_db):
            if remloc is loc:
                namestorem.append(name)
        for nm in namestorem:
            del _site_db[nm]

        # now go through and make sure none of them are in _site_names
        for nm in namestorem:
            lnm = nm.lower()
            if lnm in _site_names:
                _site_names.remove(lnm)
    else:
        msg = 'Site name "{0}" not in the database of sites'
        raise KeyError(msg.format(lname))
