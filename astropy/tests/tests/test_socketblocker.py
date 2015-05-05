from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import time

from threading import Thread

from ..helper import pytest
from ..disable_internet import no_internet
from astropy.extern.six.moves import BaseHTTPServer, SimpleHTTPServer
from astropy.extern.six.moves.urllib.request import urlopen


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


PY3_4 = sys.version_info[:2] >= (3, 4)


# Used for the below test--inline functions aren't pickleable
# by multiprocessing?
def _square(x):
    return x ** 2


@pytest.mark.skipif('not PY3_4 or sys.platform == "win32"')
def test_multiprocessing_forkserver():
    """
    Test that using multiprocessing with forkserver works.  Perhaps
    a simpler more direct test would be to just open some local
    sockets and pass something through them.

    Regression test for https://github.com/astropy/astropy/pull/3713
    """

    import multiprocessing
    ctx = multiprocessing.get_context('forkserver')
    pool = ctx.Pool(1)
    result = pool.map(_square, [1, 2, 3, 4, 5])
    pool.close()
    pool.join()
    assert result == [1, 4, 9, 16, 25]
