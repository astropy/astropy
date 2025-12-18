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

import json
from collections.abc import Iterator, Mapping
from difflib import get_close_matches
from typing import Any, Final, Self

from astropy import units as u
from astropy.utils.data import get_file_contents, get_pkg_data_contents

from .earth import EarthLocation
from .errors import UnknownSiteException

# An EarthLocation instance for the tests that need one.
_GREENWICH: Final = EarthLocation(
    lon=-0.001475 * u.deg, lat=51.477811 * u.deg, height=46 * u.m
)


class SiteRegistry(Mapping):
    """
    A bare-bones registry of EarthLocation objects.

    This acts as a mapping (dict-like object) but with the important caveat that
    it's always transforms its inputs to lower-case.  So keys are always all
    lower-case, and even if you ask for something that's got mixed case, it will
    be interpreted as the all lower-case version.
    """

    def __init__(self) -> None:
        # the keys to this are always lower-case
        self._lowercase_names_to_locations: dict[str, EarthLocation] = {}
        # these can be whatever case is appropriate
        self._names: list[str] = []

    def __getitem__(self, site_name: str) -> EarthLocation:
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
        try:
            return self._lowercase_names_to_locations[site_name.lower()]
        except KeyError:
            raise UnknownSiteException(
                site=site_name,
                attribute="the 'names' attribute",
                close_names=sorted(
                    get_close_matches(site_name, self._lowercase_names_to_locations),
                    key=len,
                ),
            ) from None
        except AttributeError:
            raise TypeError(f"site name {site_name!r} is not a 'str'") from None

    def __len__(self) -> int:
        return len(self._lowercase_names_to_locations)

    def __iter__(self) -> Iterator[str]:
        return iter(self._lowercase_names_to_locations)

    @property
    def names(self) -> list[str]:
        """
        The names in this registry.  Note that these are *not* exactly the same
        as the keys: keys are always lower-case, while `names` is what you
        should use for the actual readable names (which may be case-sensitive).

        Returns
        -------
        site : list of str
            The names of the sites in this registry
        """
        return sorted(self._names)

    def add_site(self, names: list[str], locationobj: EarthLocation) -> None:
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
    def from_json(cls, jsondb: Mapping[str, dict[str, Any]]) -> Self:
        reg = cls()
        for site in jsondb:
            site_info = jsondb[site].copy()
            location = EarthLocation.from_geodetic(
                site_info.pop("longitude") * u.Unit(site_info.pop("longitude_unit")),
                site_info.pop("latitude") * u.Unit(site_info.pop("latitude_unit")),
                site_info.pop("elevation") * u.Unit(site_info.pop("elevation_unit")),
            )
            name = site_info.pop("name")
            location.info.name = name
            aliases = [alias for alias in site_info.pop("aliases") if alias]
            if name not in aliases and name != site:
                aliases.append(name)
            location.info.meta = site_info  # whatever is left

            reg.add_site([site] + aliases, location)
        return reg


def get_downloaded_sites(jsonurl: str | None = None) -> SiteRegistry:
    """
    Load observatory database from data.astropy.org and parse into a SiteRegistry.
    """
    # we explicitly set the encoding because the default is to leave it set by
    # the users' locale, which may fail if it's not matched to the sites.json
    if jsonurl is None:
        content = get_pkg_data_contents("coordinates/sites.json", encoding="UTF-8")
    else:
        content = get_file_contents(jsonurl, encoding="UTF-8")

    jsondb = json.loads(content)
    return SiteRegistry.from_json(jsondb)
