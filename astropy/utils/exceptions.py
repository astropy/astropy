# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains errors/exceptions and warnings of general use for
astropy. Exceptions that are specific to a given subpackage should *not* be
here, but rather in the particular subpackage. Exception is the _erfa module
as we rather have the users import those exceptions from here.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


class AstropyWarning(Warning):
    """
    The base warning class from which all Astropy warnings should inherit.

    Any warning inheriting from this class is handled by the Astropy logger.
    """


class AstropyUserWarning(UserWarning, AstropyWarning):
    """
    The primary warning class for Astropy.

    Use this if you do not need a specific sub-class.
    """


class AstropyDeprecationWarning(AstropyWarning):
    """
    A warning class to indicate a deprecated feature.
    """


class AstropyPendingDeprecationWarning(PendingDeprecationWarning, AstropyWarning):
    """
    A warning class to indicate a soon-to-be deprecated feature.
    """


class AstropyBackwardsIncompatibleChangeWarning(AstropyWarning):
    """
    A warning class indicating a change in astropy that is incompatible
    with previous versions.

    The suggested procedure is to issue this warning for the version in
    which the change occurs, and remove it for all following versions.
    """


class ErfaError(ValueError):
    """
    A class for errors triggered by ERFA functions (status codes < 0)

    Note: this class should *not* be referenced by fully-qualified name, because
    it may move to ERFA in a future version.  In a future such move it will
    still be imported here as an alias, but the true namespace of the class may
    change.
    """


class ErfaWarning(AstropyUserWarning):
    """
    A class for warnings triggered by ERFA functions (status codes > 0)

    Note: this class should *not* be referenced by fully-qualified name, because
    it may move to ERFA in a future version.  In a future such move it will
    still be imported here as an alias, but the true namespace of the class may
    change.
    """
