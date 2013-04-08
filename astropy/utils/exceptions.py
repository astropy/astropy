# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains errors/exceptions and warnings of general use for
astropy. Exceptions that are specific to a given subpackage should *not*
be here, but rather in the particular subpackage.
"""


class AstropyBackwardsIncompatibleChangeWarning(Warning):
    """
    A warning class indicating a change in astropy that is incompatible
    with previous versions.

    The suggested procedure is to issue this warning for the version in
    which the change occurs, and remove it for all following versions.
    """
