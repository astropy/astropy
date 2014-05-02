from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..helper import pytest
from ..disable_internet import no_internet
from astropy.extern.six.moves import BaseHTTPServer, SimpleHTTPServer
from astropy.extern.six.moves.urllib.request import urlopen
from threading import Thread
import time


def test_outgoing_fails():
    with pytest.raises(IOError):
        with no_internet():
            urlopen('http://www.astropy.org')


class StoppableHTTPServer(BaseHTTPServer.HTTPServer,object):
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


@pytest.mark.parametrize(('localhost'), ('localhost', '127.0.0.1'))
def test_localconnect_succeeds(localhost):
    """
    Ensure that connections to localhost are allowed, since these are genuinely
    not remotedata.
    """

    # port "0" means find open port
    # see http://stackoverflow.com/questions/1365265/on-localhost-how-to-pick-a-free-port-number
    httpd = StoppableHTTPServer(('localhost', 0),
                                SimpleHTTPServer.SimpleHTTPRequestHandler)

    port = httpd.socket.getsockname()[1]

    server = Thread(target=httpd.serve_forever)
    server.setDaemon(True)

    server.start()
    time.sleep(0.1)

    urlopen('http://{localhost:s}:{port:d}'.format(localhost=localhost,port=port)).close()
