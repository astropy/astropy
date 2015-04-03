# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import ssl
import socket

from ...extern import six

if six.PY3:
    from xmlrpc.server import SimpleXMLRPCRequestHandler
else:
    from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler

from ...extern.six.moves import xmlrpc_client as xmlrpc
from .standard_profile import ThreadingXMLRPCServer

__all__ = []


from ...extern.six.moves.http_client import HTTPConnection, HTTPS_PORT


class HTTPSConnection(HTTPConnection):
    """
    This class allows communication via SSL.
    """

    default_port = HTTPS_PORT

    def __init__(self, host, port=None, key_file=None, cert_file=None,
                 cert_reqs=ssl.CERT_NONE, ca_certs=None,
                 ssl_version=None):

        HTTPConnection.__init__(self, host, port)

        self.key_file = key_file
        self.cert_file = cert_file
        self.cert_reqs = cert_reqs
        self.ca_certs = ca_certs
        self.ssl_version = ssl_version

    def connect(self):
        "Connect to a host on a given (SSL) port."

        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((self.host, self.port))
        # We have to explicitly not pass the ssl_version to
        # `ssl.wrap_socket` if it's None.
        kwargs = {
            'server_side': False,
            'certfile': self.cert_file,
            'keyfile': self.key_file,
            'cert_reqs': self.cert_reqs,
            'ca_certs': self.ca_certs,
        }
        if self.ssl_version is not None:
            kwargs['ssl_version'] = self.ssl_version
        else:
            kwargs['ssl_version'] = ssl.PROTOCOL_TLSv1
        sslconn = ssl.wrap_socket(sock, **kwargs)
        self.sock = sslconn


if six.PY2:
    from ...extern.six.moves.http_client import HTTP

    class HTTPS(HTTP):
        """
        Facility class fo HTTP communication.
        """

        _connection_class = HTTPSConnection

        def __init__(self, host='', port=None, key_file=None, cert_file=None,
                     cert_reqs=ssl.CERT_NONE, ca_certs=None,
                     ssl_version=None):

            # provide a default host, pass the X509 cert info

            # urf. compensate for bad input.
            if port == 0:
                port = None

            self._setup(self._connection_class(host, port, key_file,
                                               cert_file, cert_reqs,
                                               ca_certs, ssl_version))

            # we never actually use these for anything, but we keep them
            # here for compatibility with post-1.5.2 CVS.
            self.key_file = key_file
            self.cert_file = cert_file

        def getresponse(self, buffering=False):
            "Get the response from the server."
            return self._conn.getresponse(buffering)


class SafeTransport(xmlrpc.Transport):
    """
    Handles an HTTPS transaction to an XML-RPC server.
    """

    def __init__(self, key_file=None, cert_file=None,
                 cert_reqs=ssl.CERT_NONE, ca_certs=None,
                 ssl_version=None, strict=None,
                 use_datetime=0):

        xmlrpc.Transport.__init__(self, use_datetime)
        self._connection = (None, None)
        self.key_file = key_file
        self.cert_file = cert_file
        self.cert_reqs = cert_reqs
        self.ca_certs = ca_certs
        self.ssl_version = ssl_version

    def make_connection(self, host):

        if self._connection and host == self._connection[0]:
            return self._connection[1]

        # create a HTTPS connection object from a host descriptor
        # host may be a string, or a (host, x509-dict) tuple
        host, extra_headers, x509 = self.get_host_info(host)
        if six.PY2:
            return HTTPS(host, None, self.key_file, self.cert_file,
                         self.cert_reqs, self.ca_certs, self.ssl_version)
        else:
            self._connection = host, HTTPSConnection(host, None, **(x509 or {}))
            return self._connection[1]


class SecureXMLRPCServer(ThreadingXMLRPCServer):
    """
    An XMLRPC server supporting secure sockets connections.
    """

    def __init__(self, addr, key_file, cert_file, cert_reqs, ca_certs, ssl_version,
                 log=None, requestHandler=SimpleXMLRPCRequestHandler,
                 logRequests=True, allow_none=True, encoding=None):
        """
        Secure XML-RPC server.

        It it very similar to SimpleXMLRPCServer but it uses HTTPS for
        transporting XML data.
        """
        self.key_file = key_file
        self.cert_file = cert_file
        self.cert_reqs = cert_reqs
        self.ca_certs = ca_certs
        self.ssl_version = ssl_version
        self.allow_reuse_address = True

        ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                       logRequests, allow_none, encoding)

    def get_request(self):
        # override this to wrap socket with SSL
        sock, addr = self.socket.accept()
        # We have to explicitly not pass the ssl_version to
        # `ssl.wrap_socket` if it's None.
        kwargs = {
            'server_side': True,
            'certfile': self.cert_file,
            'keyfile': self.key_file,
            'cert_reqs': self.cert_reqs,
            'ca_certs': self.ca_certs
        }
        if self.ssl_version is not None:
            kwargs['ssl_version'] = self.ssl_version
        else:
            kwargs['ssl_version'] = ssl.PROTOCOL_TLSv1
        sslconn = ssl.wrap_socket(sock, **kwargs)
        return sslconn, addr


def get_ssl_version_name(ssl_version):
    if ssl_version is None:
        # create_default_context added after the OpenSSL bugfix in
        # Python 2.7.9 etc.  It's the best way to get the default SSL
        # protocol from Python ... otherwise, we just assume it's the
        # old default of PROTOCOL_SSLv23.
        if hasattr(ssl, 'create_default_context'):
            context = ssl.create_default_context()
            ssl_version = context.protocol
        else:
            ssl_version = ssl.PROTOCOL_TLSv1

    # get_protocol_name is an undocumented method
    if hasattr(ssl, 'get_protocol_name'):
        return ssl.get_protocol_name(ssl_version)
    else:
        # Not all versions of Python support all protocols,
        # so we have to only accept those that are present.
        for protocol in ['SSLv2', 'SSLv23', 'SSLv3', 'TLSv1', 'TLSv1_1',
                         'TLSv1_2']:
            value = getattr(ssl, 'PROTOCOL_' + protocol, None)
            if ssl_version == value:
                return protocol
        return '<unknown>'
