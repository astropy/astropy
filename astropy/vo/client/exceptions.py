# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Exceptions related to Virtual Observatory (VO)."""
from __future__ import absolute_import, division, print_function, unicode_literals


__all__ = ['BaseVOError', 'VOSError', 'MissingCatalog', 'DuplicateCatalogName',
           'DuplicateCatalogURL', 'InvalidAccessURL', 'ConeSearchError']


class BaseVOError(Exception):  # pragma: no cover
    """Base class for VO exceptions."""
    pass


class VOSError(BaseVOError):  # pragma: no cover
    """General VO service exception."""
    pass


class MissingCatalog(VOSError):  # pragma: no cover
    """VO catalog is missing."""
    pass


class DuplicateCatalogName(VOSError):  # pragma: no cover
    """VO catalog of the same title already exists."""
    pass


class DuplicateCatalogURL(VOSError):  # pragma: no cover
    """VO catalog of the same access URL already exists."""
    pass


class InvalidAccessURL(VOSError):  # pragma: no cover
    """Invalid access URL."""
    pass


class ConeSearchError(BaseVOError):  # pragma: no cover
    """General Cone Search exception."""
    pass
