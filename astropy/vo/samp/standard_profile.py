# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import traceback
import warnings

from ...extern.six.moves import xmlrpc_client as xmlrpc
from ...extern.six.moves import socketserver
from ...extern import six

if six.PY3:
    from xmlrpc.server import SimpleXMLRPCRequestHandler, SimpleXMLRPCServer
else:
    from SimpleXMLRPCServer import (SimpleXMLRPCRequestHandler,
                                    SimpleXMLRPCServer)

from .constants import SAMP_ICON
from .errors import SAMPWarning

__all__ = []


class SAMPSimpleXMLRPCRequestHandler(SimpleXMLRPCRequestHandler):
    """
    XMLRPC handler of Standard Profile requests.
    """

    def do_GET(self):

        if self.path == '/samp/icon':
            self.send_response(200, 'OK')
            self.send_header('Content-Type', 'image/png')
            self.end_headers()
            self.wfile.write(SAMP_ICON)

    if sys.version_info[:2] >= (2, 7):
        def do_POST(self):
            """
            Handles the HTTP POST request.

            Attempts to interpret all HTTP POST requests as XML-RPC calls,
            which are forwarded to the server's ``_dispatch`` method for
            handling.
            """

            # Check that the path is legal
            if not self.is_rpc_path_valid():
                self.report_404()
                return

            try:
                # Get arguments by reading body of request.
                # We read this in chunks to avoid straining
                # socket.read(); around the 10 or 15Mb mark, some platforms
                # begin to have problems (bug #792570).
                max_chunk_size = 10 * 1024 * 1024
                size_remaining = int(self.headers["content-length"])
                L = []
                while size_remaining:
                    chunk_size = min(size_remaining, max_chunk_size)
                    L.append(self.rfile.read(chunk_size))
                    size_remaining -= len(L[-1])
                data = b''.join(L)

                params, method = xmlrpc.loads(data)

                if method == "samp.webhub.register":
                    params = list(params)
                    params.append(self.client_address)
                    if 'Origin' in self.headers:
                        params.append(self.headers.get('Origin'))
                    else:
                        params.append('unknown')
                    params = tuple(params)
                    data = xmlrpc.dumps(params, methodname=method)

                elif method in ('samp.hub.notify', 'samp.hub.notifyAll',
                                'samp.hub.call', 'samp.hub.callAll',
                                'samp.hub.callAndWait'):

                    user = "unknown"

                    if method == 'samp.hub.callAndWait':
                        params[2]["host"] = self.address_string()
                        params[2]["user"] = user
                    else:
                        params[-1]["host"] = self.address_string()
                        params[-1]["user"] = user

                    data = xmlrpc.dumps(params, methodname=method)

                data = self.decode_request_content(data)
                if data is None:
                    return  # response has been sent

                # In previous versions of SimpleXMLRPCServer, _dispatch
                # could be overridden in this class, instead of in
                # SimpleXMLRPCDispatcher. To maintain backwards compatibility,
                # check to see if a subclass implements _dispatch and dispatch
                # using that method if present.
                response = self.server._marshaled_dispatch(
                    data, getattr(self, '_dispatch', None), self.path
                )
            except Exception as e:
                # This should only happen if the module is buggy
                # internal error, report as HTTP server error
                self.send_response(500)

                # Send information about the exception if requested
                if hasattr(self.server, '_send_traceback_header') and \
                   self.server._send_traceback_header:
                    self.send_header("X-exception", str(e))
                    trace = traceback.format_exc()
                    trace = str(trace.encode('ASCII', 'backslashreplace'), 'ASCII')
                    self.send_header("X-traceback", trace)

                self.send_header("Content-length", "0")
                self.end_headers()
            else:
                # got a valid XML RPC response
                self.send_response(200)
                self.send_header("Content-type", "text/xml")
                if self.encode_threshold is not None:
                    if len(response) > self.encode_threshold:
                        q = self.accept_encodings().get("gzip", 0)
                        if q:
                            try:
                                response = xmlrpc.gzip_encode(response)
                                self.send_header("Content-Encoding", "gzip")
                            except NotImplementedError:
                                pass
                self.send_header("Content-length", str(len(response)))
                self.end_headers()
                self.wfile.write(response)

    else:

        def do_POST(self):
            """
            Handles the HTTP POST request.

            Attempts to interpret all HTTP POST requests as XML-RPC calls,
            which are forwarded to the server's ``_dispatch`` method for
            handling.
            """

            # Check that the path is legal
            if not self.is_rpc_path_valid():
                self.report_404()
                return

            try:
                # Get arguments by reading body of request.
                # We read this in chunks to avoid straining
                # socket.read(); around the 10 or 15Mb mark, some platforms
                # begin to have problems (bug #792570).
                max_chunk_size = 10 * 1024 * 1024
                size_remaining = int(self.headers["content-length"])
                L = []
                while size_remaining:
                    chunk_size = min(size_remaining, max_chunk_size)
                    L.append(self.rfile.read(chunk_size))
                    size_remaining -= len(L[-1])
                data = b''.join(L)

                params, method = xmlrpc.loads(data)

                if method == "samp.webhub.register":
                    params = list(params)
                    params.append(self.client_address)
                    if 'Origin' in self.headers:
                        params.append(self.headers.get('Origin'))
                    else:
                        params.append('unknown')
                    params = tuple(params)
                    data = xmlrpc.dumps(params, methodname=method)

                elif method in ('samp.hub.notify', 'samp.hub.notifyAll',
                                'samp.hub.call', 'samp.hub.callAll',
                                'samp.hub.callAndWait'):

                    user = "unknown"

                    if method == 'samp.hub.callAndWait':
                        params[2]["host"] = self.address_string()
                        params[2]["user"] = user
                    else:
                        params[-1]["host"] = self.address_string()
                        params[-1]["user"] = user

                    data = xmlrpc.dumps(params, methodname=method)

                # In previous versions of SimpleXMLRPCServer, _dispatch
                # could be overridden in this class, instead of in
                # SimpleXMLRPCDispatcher. To maintain backwards compatibility,
                # check to see if a subclass implements _dispatch and dispatch
                # using that method if present.
                response = self.server._marshaled_dispatch(
                    data, getattr(self, '_dispatch', None)
                )
            except Exception as e:  # This should only happen if the module is buggy
                # internal error, report as HTTP server error
                self.send_response(500)

                # Send information about the exception if requested
                if hasattr(self.server, '_send_traceback_header') and \
                   self.server._send_traceback_header:
                    self.send_header("X-exception", str(e))
                    self.send_header("X-traceback", traceback.format_exc())

                self.end_headers()
            else:
                # got a valid XML RPC response
                self.send_response(200)
                self.send_header("Content-Type", "text/xml")
                self.send_header("Content-Length", str(len(response)))
                self.end_headers()
                self.wfile.write(response)

                # shut down the connection
                self.wfile.flush()
                self.connection.shutdown(1)


class ThreadingXMLRPCServer(socketserver.ThreadingMixIn, SimpleXMLRPCServer):
    """
    Asynchronous multithreaded XMLRPC server.
    """

    def __init__(self, addr, log=None,
                 requestHandler=SAMPSimpleXMLRPCRequestHandler,
                 logRequests=True, allow_none=True, encoding=None):
        self.log = log
        SimpleXMLRPCServer.__init__(self, addr, requestHandler,
                                    logRequests, allow_none, encoding)

    def handle_error(self, request, client_address):
        if self.log is None:
            socketserver.BaseServer.handle_error(self, request, client_address)
        else:
            warnings.warn("Exception happened during processing of request "
                          "from %s: %s" % (client_address, sys.exc_info()[1]),
                          SAMPWarning)
