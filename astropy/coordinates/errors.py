# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""This module defines custom errors and exceptions used in astropy.coordinates.
"""

__all__ = ["ConvertError", "UnknownSiteException"]


# TODO: consider if this should be used to `units`?
class UnitsError(ValueError):
    """
    Raised if units are missing or invalid.
    """


class ConvertError(Exception):
    """
    Raised if a coordinate system cannot be converted to another.
    """


class UnknownSiteException(KeyError):
    def __init__(self, site, attribute, close_names=None):
        message = (
            f"Site '{site}' not in database. Use {attribute} to see available sites."
            f" If '{site}' exists in the online astropy-data repository, use the"
            " 'refresh_cache=True' option to download the latest version."
        )
        if close_names:
            message += " Did you mean one of: '{}'?'".format("', '".join(close_names))
        self.site = site
        self.attribute = attribute
        self.close_names = close_names
        return super().__init__(message)
