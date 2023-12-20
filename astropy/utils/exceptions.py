# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains errors/exceptions and warnings of general use for
astropy. Exceptions that are specific to a given subpackage should *not* be
here, but rather in the particular subpackage.
"""

__all__ = [
    "AstropyWarning",
    "AstropyUserWarning",
    "AstropyDeprecationWarning",
    "AstropyPendingDeprecationWarning",
    "AstropyBackwardsIncompatibleChangeWarning",
    "DuplicateRepresentationWarning",
    "NoValue",
]


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


class DuplicateRepresentationWarning(AstropyWarning):
    """
    A warning class indicating a representation name was already registered.
    """


class _NoValue:
    """Special keyword value.

    This class may be used as the default value assigned to a
    deprecated keyword in order to check if it has been given a user
    defined value.
    """

    def __repr__(self):
        return "astropy.utils.exceptions.NoValue"


NoValue = _NoValue()


def __getattr__(name: str):
    if name in ("ErfaError", "ErfaWarning"):
        import warnings

        warnings.warn(
            f"Importing {name} from astropy.utils.exceptions was deprecated "
            "in version 6.1 and will stop working in a future version. "
            f"Instead, please use\nfrom erfa import {name}\n\n",
            category=AstropyDeprecationWarning,
            stacklevel=1,
        )

        import erfa

        return getattr(erfa, name)

    raise AttributeError(f"Module {__name__!r} has no attribute {name!r}.")
