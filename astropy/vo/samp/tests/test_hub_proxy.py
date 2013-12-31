import os
import tempfile

from ..hub_proxy import SAMPHubProxy
from ..hub import SAMPHubServer
from ..client import SAMPClient

# By default, tests should not use the internet.
from ..utils import ALLOW_INTERNET
ALLOW_INTERNET.set(False)



class TestHubProxy(object):

    def setup_method(self, method):

        fileobj, self.lockfile = tempfile.mkstemp()

        self.hub = SAMPHubServer(web_profile=False,
                                 lockfile=self.lockfile)
        self.hub.start()

        os.environ['SAMP_HUB'] = "std-lockurl:file://" + os.path.abspath(self.lockfile)

        self.proxy = SAMPHubProxy()
        self.proxy.connect()

        del os.environ['SAMP_HUB']  # hacky

    def teardown_method(self, method):

        if self.proxy.is_connected:
            self.proxy.disconnect()

        self.hub.stop()

        if os.path.exists(self.lockfile):
            os.remove(self.lockfile)

    def test_is_connected(self):
        assert self.proxy.is_connected

    def test_disconnect(self):
        self.proxy.disconnect()

    def test_get_running_hubs(self):
        SAMPHubProxy.get_running_hubs()

    def test_ping(self):
        self.proxy.ping()

    def test_registration(self):
        result = self.proxy.register(self.proxy.lockfile["samp.secret"])
        self.proxy.unregister(result['samp.private-key'])
