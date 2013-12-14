# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utility functions and classes"""
from __future__ import print_function, division
import base64
import datetime
import hashlib
import inspect
import platform
import socket
import sys
import traceback
from ...extern.six.moves import queue
from ...extern.six.moves import socketserver
from ...extern.six import StringIO

try:
    import bsddb
except ImportError:
    BDB_SUPPORT = False
else:
    BDB_SUPPORT = True

try:
    import ssl
except ImportError:
    SSL_SUPPORT = False
else:
    SSL_SUPPORT = True

try:
    from ..extern.six.moves import tkinter as tk
    HAS_TKINTER = True
except:
    HAS_TKINTER = False


PYTHON_VERSION = float(platform.python_version()[:3])

if PYTHON_VERSION >= 3.0:
    import http.client
    from http.client import HTTPConnection, HTTPS_PORT
    # from http.server import *
else:
    import httplib
    from httplib import HTTPConnection, HTTPS_PORT, HTTP
    # from SimpleHTTPServer import *

if PYTHON_VERSION >= 3.0:
    import urllib.parse as urlparse
    import urllib.error
    import urllib.request
    try:
        from urllib.parse import parse_qs
    except ImportError:
        from cgi import parse_qs
    from urllib.error import URLError
    from urllib.request import urlopen
else:
    import urllib2
    import urlparse
    try:
        from urlparse import parse_qs
    except ImportError:
        from cgi import parse_qs
    from urllib2 import urlopen, URLError

if PYTHON_VERSION >= 3.0:
    import xmlrpc.client as xmlrpc
    from xmlrpc.server import SimpleXMLRPCRequestHandler, SimpleXMLRPCServer
else:
    import xmlrpclib as xmlrpc
    from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler, SimpleXMLRPCServer

from .constants import SAMP_STATUS_ERROR, SAMPY_ICON


def internet_on():
    try:
        response = urlopen('http://google.com',timeout=1)
        return True
    except URLError as err:
        pass
    return False

__all__ = ["SAMPLog", "SAMPMsgReplierWrapper"]

__doctest_skip__ = ['.']


class _ServerProxyPoolMethod:
    # some magic to bind an XML-RPC method to an RPC server.
    # supports "nested" methods (e.g. examples.getStateName)

    def __init__(self, proxies, name):
        self.__proxies = proxies
        self.__name = name

    def __getattr__(self, name):
        return _ServerProxyPoolMethod(self.__proxies, "%s.%s" % (self.__name, name))

    def __call__(self, *args, **kwrds):
        proxy = self.__proxies.get()
        try:
            response = eval("proxy.%s(*args, **kwrds)" % self.__name)
        except:
            self.__proxies.put(proxy)
            raise
        self.__proxies.put(proxy)
        return response


class ServerProxyPool(object):

    """A thread-safe pool of `xmlrpc.ServerProxy` objects."""

    def __init__(self, size, proxy_class, *args, **keywords):

        self._proxies = queue.Queue(size)
        for i in range(size):
            self._proxies.put(proxy_class(*args, **keywords))

    def __getattr__(self, name):
        # magic method dispatcher
        return _ServerProxyPoolMethod(self._proxies, name)


if HAS_TKINTER:
    class WebProfilePopupDialogue(tk.Tk):

        """Pop-up dialog (Tkinter backend) which asks the user to consent
        or reject a connection through the WebProfile.
        """

        def __init__(self, queue, screenName=None, baseName=None, className='Tk',
                     useTk=1, sync=0, use=None):

            tk.Tk.__init__(self, screenName, baseName, className,
                           useTk, sync, use)

            self._queue = queue

            self.title("SAMP Hub")
            self._text = tk.Label(self, font=("Helvetica", 14),
                                  fg="red", justify=tk.CENTER)
            self._text.pack(padx=5, pady=5, expand=1, fill=tk.X,)

            a = tk.Button(self, text="CONSENT", command=self._consent)
            a.pack(padx=5, pady=5, expand=1, fill=tk.X, side=tk.LEFT)

            r = tk.Button(self, text="REJECT", command=self._reject)
            r.pack(padx=5, pady=5, expand=1, fill=tk.X, side=tk.RIGHT)

            self.protocol("WM_DELETE_WINDOW", self._reject)

            self.withdraw()

        def showPopup(self, request):

            samp_name = "unknown"

            if isinstance(request[0], str):
        # To support the old protocol version
                samp_name = request[0]
            else:
                samp_name = request[0]["samp.name"]

            text = \
                """A Web application which declares to be

  Name: %s
  Origin: %s

  is requesting to be registered with the SAMP Hub.
  Pay attention that if you permit its registration, such
  application will acquire all current user privileges, like
  file read/write.

  Do you give your consent?""" % (samp_name, request[2])

            self._text.configure(text=text)
            self.deiconify()

        def _consent(self):
            self._queue.put(True)
            self.withdraw()

        def _reject(self):
            self._queue.put(False)
            self.withdraw()

else:

    class WebProfilePopupDialogue(object):

        """Terminal dialog (no backend) which asks the user to consent
        or reject a connection through the WebProfile
        """

        def __init__(self, queue):

            self._queue = queue

        def showPopup(self, request):

            samp_name = "unknown"

            if isinstance(request[0], str):
                # To support the old protocol version
                samp_name = request[0]
            else:
                samp_name = request[0]["samp.name"]

            self.cls()

            text = \
                """A Web application which declares to be

  Name: %s
  Origin: %s

  is requesting to be registered with the SAMP Hub.
  Pay attention that if you permit its registration, such
  application will acquire all current user privileges, like
  file read/write.

  Do you give your consent? [yes|no]
  """ % (samp_name, request[2])

            print (text)
            self.beep()

            answer = raw_input("[no] >>> ")
            if answer.lower() in ["yes", "y"]:
                print("OK!")
                self._consent()
            else:
                print("REJECTED!")
                self._reject()

        def _consent(self):
            self._queue.put(True)

        def _reject(self):
            self._queue.put(False)

        def update(self):
            pass


class SAMPMsgReplierWrapper(object):

    """Decorator class/function that allows to automatically grab
    errors and returned maps (if any) from a function bound
    to a SAMP call (or notify).


      Parameters
      ----------
      cli : `SAMPIntegratedClient` or `SAMPClient`
          SAMP client instance.
          Decorator initialization, accepting the instance of the
          client that receives the call or notification.
      """

    def __init__(self, cli):
        self.cli = cli

    def __call__(self, f):

        def wrapped_f(*args):

            if ((inspect.ismethod(f) and f.__func__.__code__.co_argcount == 6)
                     or (inspect.isfunction(f) and f.__code__.co_argcount == 5)
                     or args[2] is None):

                # It is a notification
                f(*args)

            else:
                # It's a call
                try:
                    result = f(*args)
                    if result:
                        self.cli.hub.reply(self.cli.getPrivateKey(), args[2],
                                           {"samp.status": SAMP_STATUS_ERROR,
                                            "samp.result": result})
                except:
                    err = StringIO()
                    traceback.print_exc(file=err)
                    txt = err.getvalue()
                    self.cli.hub.reply(self.cli.getPrivateKey(), args[2],
                                       {"samp.status": SAMP_STATUS_ERROR,
                                        "samp.result": {"txt": txt}})

        return wrapped_f


class SAMPLog(object):

    """SAMP Hub logging class.

    Provides methods for gracefully print SAMPy logging messages.
    """

    #: Disable logging at all
    OFF = 0
    #: Log error messages only
    ERROR = 1
    #: Log errors and warnings
    WARNING = 2
    #: Log info, warning and error messages
    INFO = 3
    #: Log everything for debugging
    DEBUG = 4

    def __init__(self, level=INFO, stdout=sys.stdout, stderr=sys.stderr, prefix="SAMP"):
        """
        Log class constructor.

        @param level: logging level (one among L{OFF}, L{ERROR}, L{WARNING},
        L{INFO} and L{DEBUG}). By default it is set to L{INFO}.
        @type level: int

        @param stdout: file-like output device. By default it is set to sys.stdout.
        @type stdout: file

        @param stderr: file-like error device. By default it is set to sys.stderr.
        @type stderr: file

        @param prefix: prefix string to logging messages ("SAMP" by default)
        @type prefix: string
        """
        self._level = level
        self._stdout = stdout
        self._stderr = stderr
        self._prefix = prefix

    def setLevel(self, level):
        """
        Set the logging level.

        @param level: one among L{OFF}, L{ERROR}, L{WARNING}, L{INFO} and L{DEBUG}.
        @type level: int
        """
        self._level = level

    def getLevel(self):
        """
        Return the current logging level.

        @return: the current logging level. See L{setLevel}.
        @rtype: int
        """
        return self._level

    def translateLevel(self, level):
        """
        Translate a logging level from the numeric form to a string form and vice versa.
        For example: L{ERROR} is translated in C{"ERROR"} string and vice versa.

        @param level: the logging level to be translated
        @type level: int/str

        @return: the logging level traslated from the numeric form to the string form and vice versa.
        @rtype: int/str
        """

        if isinstance(level, int):
            if level == SAMPLog.OFF:
                return "OFF"
            elif level == SAMPLog.INFO:
                return "INFO"
            elif level == SAMPLog.ERROR:
                return "ERROR"
            elif level == SAMPLog.WARNING:
                return "WARNING"
            elif level == SAMPLog.DEBUG:
                return "DEBUG"
        elif isinstance(level, str):
            if level.upper() == "OFF":
                return SAMPLog.OFF
            elif level.upper() == "INFO":
                return SAMPLog.INFO
            elif level.upper() == "ERROR":
                return SAMPLog.ERROR
            elif level.upper() == "WARNING":
                return SAMPLog.WARNING
            elif level.upper() == "DEBUG":
                return SAMPLog.DEBUG
        else:
            return None

    def info(self, message):
        """
        Log an INFO message.

        @param message: the message to be logged
        @type message: string
        """
        if self._level >= self.INFO:
            self._stdout.write('[%s] Info    (%s): %s\n' %
                               (self._prefix, datetime.datetime.now().isoformat(), message))
            self._stdout.flush()

    def error(self, message):
        """
        Log an ERROR message.

        @param message: the message to be logged
        @type message: string
        """
        if self._level >= self.ERROR:
            self._stderr.write('[%s] Error   (%s): %s\n' %
                               (self._prefix, datetime.datetime.now().isoformat(), message))
            self._stderr.flush()

    def warning(self, message):
        """
        Log a WARNING message.

        @param message: the message to be logged
        @type message: string
        """
        if self._level >= self.WARNING:
            self._stderr.write('[%s] Warning (%s): %s\n' %
                               (self._prefix, datetime.datetime.now().isoformat(), message))
            self._stderr.flush()

    def debug(self, message):
        """
        Log a DEBUG message.

        @param message: the message to be logged
        @type message: string
        """
        if self._level >= self.DEBUG:
            self._stdout.write('[%s] Debug   (%s): %s\n' %
                               (self._prefix, datetime.datetime.now().isoformat(), message))
            self._stdout.flush()


class SAMPSimpleXMLRPCRequestHandler(SimpleXMLRPCRequestHandler):

    """XMLRPC handler of Standar Profile requests (internal use only)"""

    def do_GET(self):

        if self.path == '/sampy/icon':
            self.send_response(200, 'OK')
            self.send_header('Content-Type', 'image/png')
            self.end_headers()
            self.wfile.write(base64.decodestring(SAMPY_ICON))

    if PYTHON_VERSION >= 2.7:

        def do_POST(self):
            """Handles the HTTP POST request.

            Attempts to interpret all HTTP POST requests as XML-RPC calls,
            which are forwarded to the server's `_dispatch` method for handling.
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
                max_chunk_size = 10*1024*1024
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

                    if 'Authorization' in self.headers:
                        # handle Basic authentication
                        (enctype, encstr) = self.headers.get('Authorization').split()
                        user, password = base64.standard_b64decode(encstr).split(':')

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
            except Exception as e:  # This should only happen if the module is buggy
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

    elif PYTHON_VERSION >= 2.6 and PYTHON_VERSION < 2.7:

        def do_POST(self):
            """Handles the HTTP POST request.

            Attempts to interpret all HTTP POST requests as XML-RPC calls,
            which are forwarded to the server's `_dispatch` method for handling.
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
                max_chunk_size = 10*1024*1024
                size_remaining = int(self.headers["content-length"])
                L = []
                while size_remaining:
                    chunk_size = min(size_remaining, max_chunk_size)
                    L.append(self.rfile.read(chunk_size))
                    size_remaining -= len(L[-1])
                data = ''.join(L)

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

                    if 'Authorization' in self.headers:
                        # handle Basic authentication
                        (enctype, encstr) = self.headers.get('Authorization').split()
                        user, password = base64.standard_b64decode(encstr).split(':')

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

    else:

        def do_POST(self):
            """Handles the HTTP POST request.

            Attempts to interpret all HTTP POST requests as XML-RPC calls,
            which are forwarded to the server's `_dispatch` method for handling.
            """

            # Check that the path is legal
            if not self.is_rpc_path_valid():
                self.report_404()
                return

            try:
                # Get arguments by reading body of request        self.connection.close().
                # We read this in chunks to avoid straining
                # socket.read(); around the 10 or 15Mb mark, some platforms
                # begin to have problems (bug #792570).
                max_chunk_size = 10*1024*1024
                size_remaining = int(self.headers["content-length"])
                L = []
                while size_remaining:
                    chunk_size = min(size_remaining, max_chunk_size)
                    L.append(self.rfile.read(chunk_size))
                    size_remaining -= len(L[-1])
                data = ''.join(L)

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

                    if 'Authorization' in self.headers:
                        # handle Basic authentication
                        (enctype, encstr) = self.headers.get('Authorization').split()
                        user, password = base64.standard_b64decode(encstr).split(':')

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
            except:  # This should only happen if the module is buggy
                # internal error, report as HTTP server error
                self.send_response(500)
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

    """Asynchronous multithreaded XMLRPC server (internal use only)"""

    def __init__(self, addr, log=None, requestHandler=SAMPSimpleXMLRPCRequestHandler,
                 logRequests=True, allow_none=True, encoding=None):
        self.log = log
        SimpleXMLRPCServer.__init__(self, addr, requestHandler,
                                    logRequests, allow_none, encoding)

    def handle_error(self, request, client_address):
        if self.log == None:
            socketserver.BaseServer.handle_error(self, request, client_address)
        else:
            self.log.warning("Exception happened during processing of request from %s: %s" % (client_address, sys.exc_info()[1]))


class WebProfileRequestHandler(SAMPSimpleXMLRPCRequestHandler):

    """Handler of XMLRPC requests performed through the WebProfile (internal use
    only)
    """

    def _send_CORS_header(self):

        if not self.headers.get('Origin') is None:

            method = self.headers.get('Access-Control-Request-Method')
            if method and self.command == "OPTIONS":
                # Preflight method
                self.send_header('Content-Length', '0')
                self.send_header('Access-Control-Allow-Origin', self.headers.get('Origin'))
                self.send_header('Access-Control-Allow-Methods', method)
                self.send_header('Access-Control-Allow-Headers', 'Content-Type')
                self.send_header('Access-Control-Allow-Credentials', 'true')
            else:
                # Simple method
                self.send_header('Access-Control-Allow-Origin', self.headers.get('Origin'))
                self.send_header('Access-Control-Allow-Headers', 'Content-Type')
                self.send_header('Access-Control-Allow-Credentials', 'true')

    def end_headers(self):
        self._send_CORS_header()
        SAMPSimpleXMLRPCRequestHandler.end_headers(self)

    def _serve_cross_domain_xml(self):

        cross_domain = False

        if self.path == "/crossdomain.xml":
            # Adobe standard
            response = """<?xml version='1.0'?>
<!DOCTYPE cross-domain-policy SYSTEM "http://www.adobe.com/xml/dtds/cross-domain-policy.dtd">
<cross-domain-policy>
  <site-control permitted-cross-domain-policies="all"/>
  <allow-access-from domain="*"/>
  <allow-http-request-headers-from domain="*" headers="*"/>
</cross-domain-policy>"""

            self.send_response(200, 'OK')
            self.send_header('Content-Type', 'text/x-cross-domain-policy')
            self.send_header("Content-Length", str(len(response)))
            self.end_headers()
            self.wfile.write(response)
            self.wfile.flush()
            cross_domain = True

        elif self.path == "/clientaccesspolicy.xml":
            # Microsoft standard
            response = """<?xml version='1.0'?>
<access-policy>
  <cross-domain-access>
    <policy>
      <allow-from>
        <domain uri="*"/>
      </allow-from>
      <grant-to>
        <resource path="/" include-subpaths="true"/>
      </grant-to>
    </policy>
  </cross-domain-access>
</access-policy>"""

            self.send_response(200, 'OK')
            self.send_header('Content-Type', 'text/xml')
            self.send_header("Content-Length", str(len(response)))
            self.end_headers()
            self.wfile.write(response)
            self.wfile.flush()
            cross_domain = True

        return cross_domain

    def do_POST(self):
        if self._serve_cross_domain_xml():
            return

        return SAMPSimpleXMLRPCRequestHandler.do_POST(self)

    def do_HEAD(self):

        if not self.is_http_path_valid():
            self.report_404()
            return

        if self._serve_cross_domain_xml():
            return

    def do_OPTIONS(self):

        self.send_response(200, 'OK')
        self.end_headers()

    def do_GET(self):

        if not self.is_http_path_valid():
            self.report_404()
            return

        split_path = self.path.split('?')

        if split_path[0] in ['/translator/%s' % clid for clid in self.server.clients]:
            # Request of a file proxying
            urlpath = parse_qs(split_path[1])
            try:
                proxyfile = urlopen(urlpath["ref"][0])
                self.send_response(200, 'OK')
                self.end_headers()
                self.wfile.write(proxyfile.read())
                proxyfile.close()
            except:
                self.report_404()
                return

        if self._serve_cross_domain_xml():
            return

    def is_http_path_valid(self):

        valid_paths = ["/clientaccesspolicy.xml", "/crossdomain.xml"] + ['/translator/%s' % clid for clid in self.server.clients]
        return self.path.split('?')[0] in valid_paths


class WebProfileXMLRPCServer(ThreadingXMLRPCServer):

    """XMLRPC server supporting the SAMP Web Profile"""

    def __init__(self, addr, log=None, requestHandler=WebProfileRequestHandler,
                 logRequests=True, allow_none=True, encoding=None):

        self.clients = []
        ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                       logRequests, allow_none, encoding)

    def add_client(self, client_id):
        self.clients.append(client_id)

    def remove_client(self, client_id):
        try:
            self.clients.remove(client_id)
        except:
            pass


if SSL_SUPPORT:

    class HTTPSConnection(HTTPConnection):

        """This class allows communication via SSL (client side - internal use
        only).
        """

        default_port = HTTPS_PORT

        def __init__(self, host, port=None, key_file=None, cert_file=None,
                     cert_reqs=ssl.CERT_NONE, ca_certs=None,
                     ssl_version=ssl.PROTOCOL_SSLv3, strict=None):

            HTTPConnection.__init__(self, host, port, strict)

            self.key_file = key_file
            self.cert_file = cert_file
            self.cert_reqs = cert_reqs
            self.ca_certs = ca_certs
            self.ssl_version = ssl_version

        def connect(self):
            "Connect to a host on a given (SSL) port."

            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect((self.host, self.port))
            sslconn = ssl.wrap_socket(sock, server_side=False,
                                      certfile=self.cert_file,
                                      keyfile=self.key_file,
                                      cert_reqs=self.cert_reqs,
                                      ca_certs=self.ca_certs,
                                      ssl_version=self.ssl_version)
            self.sock = sslconn

    if PYTHON_VERSION < 3.0:

        class HTTPS(HTTP):

            """Facility class fo HTTP communication (internal use only)"""

            _connection_class = HTTPSConnection

            def __init__(self, host='', port=None, key_file=None, cert_file=None,
                         cert_reqs=ssl.CERT_NONE, ca_certs=None,
                         ssl_version=ssl.PROTOCOL_SSLv3):

                # provide a default host, pass the X509 cert info

                # urf. compensate for bad input.
                if port == 0:
                    port = None

                self._setup(self._connection_class(host, port, key_file,
                                                   cert_file, cert_reqs,
                                                   ca_certs, ssl_version, None))

                # we never actually use these for anything, but we keep them
                # here for compatibility with post-1.5.2 CVS.
                self.key_file = key_file
                self.cert_file = cert_file

            def getresponse(self, buffering=False):
                "Get the response from the server."
                return self._conn.getresponse(buffering)

    class SafeTransport(xmlrpc.Transport):

        """Handles an HTTPS transaction to an XML-RPC server. (internal use only)"""

        def __init__(self, key_file=None, cert_file=None,
                     cert_reqs=ssl.CERT_NONE, ca_certs=None,
                     ssl_version=ssl.PROTOCOL_SSLv3, strict=None,
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
            if PYTHON_VERSION < 3.0:
                return HTTPS(host, None, self.key_file, self.cert_file,
                             self.cert_reqs, self.ca_certs, self.ssl_version)
            else:
                self._connection = host, http.client.HTTPSConnection(chost,
                                                                     None, **(x509 or {}))
                return self._connection[1]

    class SecureXMLRPCServer(ThreadingXMLRPCServer):

        """An XMLRPC server supporting secure sockets connections (internal use only)
        """

        def __init__(self, addr, keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                     log=None, requestHandler=SimpleXMLRPCRequestHandler,
                     logRequests=True, allow_none=True, encoding=None):
            """
            Secure XML-RPC server.

            It it very similar to SimpleXMLRPCServer but it uses HTTPS for transporting XML data.
            """
            self.keyfile = keyfile
            self.certfile = certfile
            self.cert_reqs = cert_reqs
            self.ca_certs = ca_certs
            self.ssl_version = ssl_version
            self.allow_reuse_address = True

            ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                           logRequests, allow_none, encoding)

        def get_request(self):
            # override this to wrap socket with SSL
            sock, addr = self.socket.accept()
            sslconn = ssl.wrap_socket(sock, server_side=True,
                                      certfile=self.certfile,
                                      keyfile=self.keyfile,
                                      cert_reqs=self.cert_reqs,
                                      ca_certs=self.ca_certs,
                                      ssl_version=self.ssl_version)
            return sslconn, addr


if BDB_SUPPORT:

    class BasicAuthSimpleXMLRPCRequestHandler(SAMPSimpleXMLRPCRequestHandler):

        """XML-RPC Request Handler for Basic Authentication support. (internal use only)

        Paramters
        ---------
        auth_file : str
            Authentication file path. It is a Berkeley DB file in Hash
            format containing a set of key=value pairs of the form:
            `<user name>=md5(<password>)<group 1>,<group 2>,<group 3>,...`.

        access_restrict : dict
            Dictionary containing the restriction rules for authentication.
            If the access must be restricted to a specific user then `access_restrict` is a dictionary
            containing `{"user"; <user name>}`. If the access must be restricted to the
            users belonging to a certain group, the `access_restrict` is a dictionary containing
            `{"group"; <group name>}`. An additional key can be present: `"admin": <administrator user>`.
            It defines the name of the administrator user with full access permission.
        """

        def __init__(self, request, client_address, server, auth_file, access_restrict=None):
            self.db = bsddb.hashopen(auth_file, "r")
            self.access_restrict = access_restrict
            SimpleXMLRPCRequestHandler.__init__(self, request, client_address, server)
            self.db.close()

        def checkId(self, id, pwd):

            if id in self.db.keys():

                pwdhash = self.db[id][0:16]
                groups = self.db[id][16:]
                pwd = hashlib.md5(pwd.encode('utf-8')).digest()

                if self.access_restrict != None:

                    # ADMIN TEST
                    if "admin" in self.access_restrict:
                        admin = self.access_restrict["admin"]
                        if admin in self.db:
                            adminpwdhash = self.db[admin][0:16]
                            if admin == id and adminpwdhash == pwd:
                                return True

                    # TEST USER RESTRICTION
                    if "user" in self.access_restrict:
                        if self.access_restrict["user"] == id and pwdhash == pwd:
                            return True
                        else:
                            return False

                    # TEST GROUP RESTRICTION
                    if "group" in self.access_restrict:
                        if self.access_restrict["group"] in groups.split(",") and pwdhash == pwd:
                            return True
                        else:
                            return False
                else:
                    if pwdhash == pwd:
                        return True
                    else:
                        return False
            else:
                return False

        def authenticate_client(self):
            validuser = False

            if 'Authorization' in self.headers:
                # handle Basic authentication
                (enctype, encstr) = self.headers.get('Authorization').split()
                (user, password) = base64.standard_b64decode(encstr).split(':')
                validuser = self.checkId(user, password)

            return validuser

        def do_POST(self):

            if self.authenticate_client():
                SAMPSimpleXMLRPCRequestHandler.do_POST(self)
            else:
                self.report_401()

        def report_401(self):
            # Report a 401 error
            self.send_response(401)
            self.send_header("WWW-Authenticate", "Basic realm=\"Protected access\"")
            self.end_headers()
            # shut down the connection
            self.connection.shutdown(1)
            self.connection.close()

    class BasicAuthXMLRPCServer(ThreadingXMLRPCServer):

        """XML-RPC server with Basic Authentication support. (internal use only).

        Parameters
        ----------
        auth_file : str
            Authentication file path. It is a Berkeley DB file in Hash
            format containing a set of key=value pairs of the form:
            `<user name>=md5(<password>)<group 1>,<group 2>,<group 3>,...`.

        access_restrict : dict
            Dictionary containing the restriction rules for authentication.
            If the access must be restricted to a specific user then access_restrict is a dictionary
            containing `{"user"; <user name>}`. If the access must be restricted to the
            users belonging to a certain group, the access_restrict is a dictionary containing
            `{"group"; <group name>}`. An additional key can be present: `"admin": <administrator user>`.
            It defines the name of the administrator user with full access permission.
        """

        def __init__(self, addr, auth_file, access_restrict=None, log=None,
                     requestHandler=BasicAuthSimpleXMLRPCRequestHandler,
                     logRequests=True, allow_none=True, encoding=None):

            self.auth_file = auth_file
            self.access_restrict = access_restrict

            ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                           logRequests, allow_none, encoding)

        def finish_request(self, request, client_address):
            if self.auth_file != None and self.RequestHandlerClass == BasicAuthSimpleXMLRPCRequestHandler:
                self.RequestHandlerClass(request, client_address, self,
                                         self.auth_file, self.access_restrict)
            else:
                ThreadingXMLRPCServer.finish_request(self, request, client_address)

if SSL_SUPPORT and BDB_SUPPORT:

    class BasicAuthSecureXMLRPCServer(ThreadingXMLRPCServer):

        """XML-RPC server with Basic Authentication support, secure socket
        connections and multithreaded. (internal use only)"""

        def __init__(self, addr, keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                     auth_file, access_restrict=None, log=None,
                     requestHandler=BasicAuthSimpleXMLRPCRequestHandler,
                     logRequests=True, allow_none=True, encoding=None):

            self.keyfile = keyfile
            self.certfile = certfile
            self.cert_reqs = cert_reqs
            self.ca_certs = ca_certs
            self.ssl_version = ssl_version
            self.allow_reuse_address = True
            self.auth_file = auth_file
            self.access_restrict = access_restrict

            ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                           logRequests, allow_none, encoding)

        def get_request(self):
            # override this to wrap socket with SSL
            sock, addr = self.socket.accept()
            sslconn = ssl.wrap_socket(sock, server_side=True,
                                      certfile=self.certfile,
                                      keyfile=self.keyfile,
                                      cert_reqs=self.cert_reqs,
                                      ca_certs=self.ca_certs,
                                      ssl_version=self.ssl_version)
            return sslconn, addr

        def finish_request(self, request, client_address):
            if self.auth_file != None and self.RequestHandlerClass == BasicAuthSimpleXMLRPCRequestHandler:
                self.RequestHandlerClass(request, client_address, self,
                                         self.auth_file, self.access_restrict)
            else:
                ThreadingXMLRPCServer.finish_request(self, request, client_address)


class _HubAsClient(object):

    def __init__(self, handler):
        self._handler = handler

    def __getattr__(self, name):
        # magic method dispatcher
        return _HubAsClientMethod(self._handler, name)


class _HubAsClientMethod(object):

    def __init__(self, send, name):
        self.__send = send
        self.__name = name

    def __getattr__(self, name):
        return _HubAsClientMethod(self.__send, "%s.%s" % (self.__name, name))

    def __call__(self, *args):
        return self.__send(self.__name, args)
