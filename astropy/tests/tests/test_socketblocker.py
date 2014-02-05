import urllib
from ..helper import pytest
import BaseHTTPServer
import SimpleHTTPServer
from threading import Thread
import time


def test_outgoing_fails():
    with pytest.raises(IOError):
        urllib.urlopen('http://www.astropy.org')


def run_while_true(server_class=BaseHTTPServer.HTTPServer,
                   handler_class=SimpleHTTPServer.SimpleHTTPRequestHandler):
    """
    This assumes that keep_running() is a function of no arguments which
    is tested initially and after each request.  If its return value
    is true, the server continues.
    """
    server_address = ('localhost', 8000)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()
 
server = Thread(target=run_while_true)
server.setDaemon(True)


def test_localconnect_succeeds():
    """
    Ensure that connections to localhost are allowed, since these are genuinely
    not remotedata.
    """

    server.start()
    time.sleep(0.1)

    urllib.urlopen('http://localhost:8000').close()
    urllib.urlopen('http://127.0.0.1:8000').close()
