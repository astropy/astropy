# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains errors/exceptions and warnings of general use for
astropy. Exceptions that are specific to a given subpackage should *not*
be here, but rather in the particular subpackage.
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
