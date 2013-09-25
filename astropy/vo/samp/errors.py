# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines custom errors and exceptions used in `astropy.vo.samp`."""

import platform
PYTHON_VERSION = float(platform.python_version()[:3])
if PYTHON_VERSION >= 3.0:
    import xmlrpc.client as xmlrpc
else:
    import xmlrpclib as xmlrpc


__all__ = ['SAMPHubError', 'SAMPClientError', 'SAMPProxyError']


class SAMPHubError(Exception):
    """SAMP Hub exception."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class SAMPClientError(Exception):
    """SAMP Client exceptions."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

#: SAMP Proxy Hub exceptions (overwites xmlrpc.Fault).
SAMPProxyError = xmlrpc.Fault
