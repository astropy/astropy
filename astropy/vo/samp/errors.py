# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines custom errors and exceptions used in `astropy.vo.samp`.
"""

from ...extern import six

if six.PY3:
    import xmlrpc.client as xmlrpc
else:
    import xmlrpclib as xmlrpc


__all__ = ['SAMPHubError', 'SAMPClientError', 'SAMPProxyError']


class SAMPHubError(Exception):
    """
    SAMP Hub exception.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class SAMPClientError(Exception):
    """
    SAMP Client exceptions.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class SAMPProxyError(xmlrpc.Fault):
    """
    SAMP Proxy Hub exception
    """
