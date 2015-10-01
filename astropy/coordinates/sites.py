# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Observatories accessible without internet access originate from the IRAF
Observatory Database, and are stored in ``data/observatories.json``.  This is
inteded mainly as a fallback file, and the online file is where new changes
should go.

Additions or corrections to the observatory list can be submitted via Pull
Request to the [astropy-data GitHub repository](https://github.com/astropy/astropy-data),
updating the ``location.json`` file.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
from difflib import get_close_matches

from ..utils.data import get_pkg_data_contents, get_file_contents
from .earth import EarthLocation
from .. import units as u


class SiteRegistry(object):
    """
    A bare-bones registry of EarthLocation objects.
    """
    def __init__(self):
        # the keys to this are always lower-case
        self._site_dict = {}
        # these can be whatever case is appropriate
        self._site_names = []

    def get_site(self, site_name):
        """
        Returns an EarthLocation for a known site in this registry.

        Parameters
        ----------
        site_name : str
            Name of the observatory (case-insensitive).

        Returns
        -------
        site : `~astropy.coordinates.EarthLocation`
            The location of the observatory.
        """
        if site_name.lower() not in self._site_dict:
            # If site name not found, find close matches and suggest them in error
            close_names = get_close_matches(site_name, self._site_dict)
            close_names = sorted(close_names, key=lambda x: len(x))
            if close_names:
                errmsg = ('Site not in database. Use ``get_site_names()`` '
                          'to see available sites. Did you mean one of: "{0}"?')
                errmsg = errmsg.format("', '".join(close_names))
            else:
                errmsg = 'Site "{0}" not in database.'.format(site_name)
            raise KeyError(errmsg)

        return self._site_dict[site_name.lower()]

    def get_site_names(self):
        """
        Returns the names in this registry

        Returns
        -------
        site : list of str
            The names of the sites in this registry
        """
        return sorted(self._site_names)

    def add_site(self, names, locationobj):
        """
        Adds a location to the registry.

        Parameters
        ----------
        names : list of str
            All the names this site should go under
        locationobj: `~astropy.coordinates.EarthLocation`
            The actual site object
        """
        for name in names:
            self._site_dict[name.lower()] = locationobj
            self._site_names.append(name)

    @classmethod
    def from_json(cls, jsondb):
        reg = cls()
        for site in jsondb:
            location = EarthLocation.from_geodetic(jsondb[site]['longitude'] * u.Unit(jsondb[site]['longitude_unit']),
                                                   jsondb[site]['latitude'] * u.Unit(jsondb[site]['latitude_unit']),
                                                   jsondb[site]['elevation'] * u.Unit(jsondb[site]['elevation_unit']))
            location.info.name = jsondb[site]['name']

            reg.add_site([site] + jsondb[site]['aliases'], location)
        return reg


def get_builtin_sites():
    """
    Load observatory database from data/observatories.json and parse them into
    a SiteRegistry.
    """
    jsondb = json.loads(get_pkg_data_contents('data/observatories.json'))
    return SiteRegistry.from_json(jsondb)


def get_downloaded_sites(jsonurl='http://data.astropy.org/locations.json'):
    """
    Load observatory database from data.astropy.org and parse into a SiteRegistry
    """
    jsondb = json.loads(get_file_contents(jsonurl))
    return SiteRegistry.from_json(jsondb)
