from ..hub_proxy import SAMPHubProxy
from ..hub import SAMPHubServer

from .. import conf

def setup_module(module):
    conf.use_internet = False


class TestHubProxy(object):

    def setup_method(self, method):

        self.hub = SAMPHubServer(web_profile=False, mode='multiple', pool_size=1)
        self.hub.start()

        self.proxy = SAMPHubProxy()
        self.proxy.connect(hub=self.hub, pool_size=1)

    def teardown_method(self, method):

        if self.proxy.is_connected:
            self.proxy.disconnect()

        self.hub.stop()

    def test_is_connected(self):
        assert self.proxy.is_connected

    def test_disconnect(self):
        self.proxy.disconnect()

    def test_ping(self):
        self.proxy.ping()

    def test_registration(self):
        result = self.proxy.register(self.proxy.lockfile["samp.secret"])
        self.proxy.unregister(result['samp.private-key'])


def test_custom_lockfile(tmpdir):

    lockfile = tmpdir.join('.samptest').realpath().strpath

    hub = SAMPHubServer(web_profile=False, lockfile=lockfile, pool_size=1)
    hub.start()

    proxy = SAMPHubProxy()
    proxy.connect(hub=hub, pool_size=1)

    hub.stop()
