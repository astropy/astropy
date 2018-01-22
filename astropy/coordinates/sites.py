# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Currently the only site accessible without internet access is the Royal
Greenwich Observatory, as an example (and for testing purposes).  In future
releases, a canonical set of sites may be bundled into astropy for when the
online registry is unavailable.

Additions or corrections to the observatory list can be submitted via Pull
Request to the [astropy-data GitHub repository](https://github.com/astropy/astropy-data),
updating the ``location.json`` file.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
from difflib import get_close_matches
from collections import Mapping

from ..utils.data import get_pkg_data_contents, get_file_contents
from .earth import EarthLocation
from .errors import UnknownSiteException
from .. import units as u


class SiteRegistry(Mapping):
    """
    A bare-bones registry of EarthLocation objects.

    This acts as a mapping (dict-like object) but with the important caveat that
    it's always transforms its inputs to lower-case.  So keys are always all
    lower-case, and even if you ask for something that's got mixed case, it will
    be interpreted as the all lower-case version.
    """
    def __init__(self):
        # the keys to this are always lower-case
        self._lowercase_names_to_locations = {}
        # these can be whatever case is appropriate
        self._names = []

    def __getitem__(self, site_name):
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
        if site_name.lower() not in self._lowercase_names_to_locations:
            # If site name not found, find close matches and suggest them in error
            close_names = get_close_matches(site_name, self._lowercase_names_to_locations)
            close_names = sorted(close_names, key=len)

            raise UnknownSiteException(site_name, "the 'names' attribute", close_names=close_names)

        return self._lowercase_names_to_locations[site_name.lower()]

    def __len__(self):
        return len(self._lowercase_names_to_locations)

    def __iter__(self):
        return iter(self._lowercase_names_to_locations)

    def __contains__(self, site_name):
        return site_name.lower() in self._lowercase_names_to_locations

    @property
    def names(self):
        """
        The names in this registry.  Note that these are *not* exactly the same
        as the keys: keys are always lower-case, while `names` is what you
        should use for the actual readable names (which may be case-sensitive)

        Returns
        -------
        site : list of str
            The names of the sites in this registry
        """
        return sorted(self._names)

    def add_site(self, names, locationobj):
        """
        Adds a location to the registry.

        Parameters
        ----------
        names : list of str
            All the names this site should go under
        locationobj : `~astropy.coordinates.EarthLocation`
            The actual site object
        """
        for name in names:
            self._lowercase_names_to_locations[name.lower()] = locationobj
            self._names.append(name)

    @classmethod
    def from_json(cls, jsondb):
        reg = cls()
        for site in jsondb:
            site_info = jsondb[site]
            location = EarthLocation.from_geodetic(site_info['longitude'] * u.Unit(site_info['longitude_unit']),
                                                   site_info['latitude'] * u.Unit(site_info['latitude_unit']),
                                                   site_info['elevation'] * u.Unit(site_info['elevation_unit']))
            location.info.name = site_info['name']

            reg.add_site([site] + site_info['aliases'], location)
        reg._loaded_jsondb = jsondb
        return reg


def get_builtin_sites():
    """
    Load observatory database from data/observatories.json and parse them into
    a SiteRegistry.
    """
    jsondb = json.loads(get_pkg_data_contents('data/sites.json'))
    return SiteRegistry.from_json(jsondb)


def get_downloaded_sites(jsonurl=None):
    """
    Load observatory database from data.astropy.org and parse into a SiteRegistry
    """

    # we explicitly set the encoding because the default is to leave it set by
    # the users' locale, which may fail if it's not matched to the sites.json
    if jsonurl is None:
        content = get_pkg_data_contents('coordinates/sites.json', encoding='UTF-8')
    else:
        content = get_file_contents(jsonurl, encoding='UTF-8')

    jsondb = json.loads(content)
    return SiteRegistry.from_json(jsondb)
