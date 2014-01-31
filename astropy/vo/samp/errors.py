# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines custom errors and exceptions used in `astropy.vo.samp`.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern.six.moves import xmlrpc_client as xmlrpc
from ...utils.exceptions import AstropyUserWarning


__all__ = ['SAMPWarning', 'SAMPHubError', 'SAMPClientError', 'SAMPProxyError']


class SAMPWarning(AstropyUserWarning):
    """
    SAMP-specific Astropy warning class
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
    SAMP Proxy Hub exception
    """
