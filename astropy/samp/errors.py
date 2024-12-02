# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines custom errors and exceptions used in `astropy.samp`.
"""

import xmlrpc.client as xmlrpc

from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["SAMPClientError", "SAMPHubError", "SAMPProxyError", "SAMPWarning"]


class SAMPWarning(AstropyUserWarning):
    """
    SAMP-specific Astropy warning class.
    """


class SAMPHubError(Exception):
    """
    SAMP Hub exception.
    """


class SAMPClientError(Exception):
    """
    SAMP Client exceptions.
    """


class SAMPProxyError(xmlrpc.Fault):
    """
    SAMP Proxy Hub exception.
    """
