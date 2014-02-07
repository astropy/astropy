from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..helper import pytest
from astropy.extern.six.moves import BaseHTTPServer, SimpleHTTPServer
from astropy.extern.six.moves.urllib.request import urlopen
from threading import Thread
import time


def test_outgoing_fails():
    with pytest.raises(IOError):
        urlopen('http://www.astropy.org')


class StoppableHTTPServer(BaseHTTPServer.HTTPServer):
    def __init__(self, *args):
        super(StoppableHTTPServer, self).__init__(*args)
        self.stop = False

    def handle_request(self):
        self.stop = True
        super(StoppableHTTPServer, self).handle_request()
     
    def serve_forever(self):
        """
        Serve until stop set, which will happen if any request is handled
        """
        while not self.stop:
            self.handle_request()

def run_while_true(server_class=BaseHTTPServer.HTTPServer,
                   handler_class=SimpleHTTPServer.SimpleHTTPRequestHandler):
    """
    This assumes that keep_running() is a function of no arguments which
    is tested initially and after each request.  If its return value
    is true, the server continues.
    """
    server_address = ('localhost', 8000)
    httpd = server_class(server_address, handler_class)

    # "serve forever" really means "serve until first request"
    httpd.serve_forever()


@pytest.mark.parametrize(('localhost'),('localhost','127.0.0.1'))
def test_localconnect_succeeds(localhost):
    """
    Ensure that connections to localhost are allowed, since these are genuinely
    not remotedata.
    """

    server = Thread(target=run_while_true)
    server.setDaemon(True)

    server.start()
    time.sleep(0.1)

    urlopen('http://{}:8000'.format(localhost)).close()
