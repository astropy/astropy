from ....tests.helper import pytest
from ..hub import SAMPHubServer
from ..integrated_client import SAMPIntegratedClient

class Receiver(object):
    def __init__(self):
        self.received = False
    def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
        self.params = params
        self.received = True

class TestIntegratedClient(object):

    def setup_method(self, method):

        self.hub = SAMPHubServer(web_profile=False)
        self.hub.start()

        self.client1 = SAMPIntegratedClient()
        self.client1.connect()
        self.client1_id = self.client1.get_public_id()

        self.client2 = SAMPIntegratedClient()
        self.client2.connect()
        self.client2_id = self.client2.get_public_id()

        self.metadata1 = {"samp.name": "Client 1",
                          "samp.description.text": "Client 1 Description",
                          "client.version": "1.1"}

        self.metadata2 = {"samp.name": "Client 2",
                          "samp.description.text": "Client 2 Description",
                          "client.version": "1.2"}

    def teardown_method(self, method):
        self.client1.disconnect()
        self.client2.disconnect()
        self.hub.stop()

    def test_ping(self):
        self.client1.ping()
        self.client2.ping()

    def test_registered(self):

        assert self.client1_id not in self.client1.get_registered_clients()  # ensure docstring makes this clear
        assert self.client2_id in self.client1.get_registered_clients()

        assert self.client1_id in self.client2.get_registered_clients()
        assert self.client2_id not in self.client2.get_registered_clients()  # ensure docstring makes this clear

    def test_metadata(self):

        assert self.client1.get_metadata(self.client1_id) == {}
        assert self.client1.get_metadata(self.client2_id) == {}
        assert self.client2.get_metadata(self.client1_id) == {}
        assert self.client2.get_metadata(self.client2_id) == {}

        self.client1.declare_metadata(self.metadata1)

        assert self.client1.get_metadata(self.client1_id) == self.metadata1
        assert self.client2.get_metadata(self.client1_id) == self.metadata1
        assert self.client1.get_metadata(self.client2_id) == {}
        assert self.client2.get_metadata(self.client2_id) == {}

        self.client2.declare_metadata(self.metadata2)

        assert self.client1.get_metadata(self.client1_id) == self.metadata1
        assert self.client2.get_metadata(self.client1_id) == self.metadata1
        assert self.client1.get_metadata(self.client2_id) == self.metadata2
        assert self.client2.get_metadata(self.client2_id) == self.metadata2

    def test_subscriptions(self):

        assert self.client1.get_subscribed_clients('table.load.votable') == {}

        rec = Receiver()
        self.client1.declare_subscriptions({'table.load.votable':rec})

        assert self.client2.get_subscribed_clients('table.load.votable') == {self.client1_id:{}}
        assert self.client1.get_subscribed_clients('table.load.votable') == {}  # a bit strange, make sure this is clear in docstring

        assert 'table.load.votable' in self.client1.get_subscriptions(self.client1_id)
        assert 'table.load.votable' in self.client2.get_subscriptions(self.client1_id)

    def test_is_connected(self):
        assert self.client1.is_connected
        assert self.client2.is_connected

    @pytest.mark.xfail  # need a better exception than current 'Fault'
    def test_notify(self):
        message = {}
        with pytest.raises(SAMPError):
            self.client1.notify(self.client2_id, message)

    def test_notify(self):

        self.client2.declare_subscriptions({'table.load.votable':Receiver()})

        # Standard notify
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.notify(self.client2_id, message)

        # Easy notify
        self.client1.enotify(self.client2_id, "table.load.votable")

    def test_notify_all(self):

        self.client2.declare_subscriptions({'table.load.votable':Receiver()})

        # Standard notify
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.notify_all(message)

        # Easy notify
        self.client1.enotify_all("table.load.votable")

    def test_call(self):

        self.client2.declare_subscriptions({'table.load.votable':Receiver()})

        # Standard call
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.call(self.client2_id, 'test_tag', message)

        # Easy call
        self.client1.ecall(self.client2_id, 'test_tag', "table.load.votable")

    def test_call_all(self):

        self.client2.declare_subscriptions({'table.load.votable':Receiver()})

        # Standard call
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.call_all('test_tag', message)

        # Easy call
        self.client1.ecall_all('test_tag', "table.load.votable")

    def test_call_and_wait(self):

        self.client2.declare_subscriptions({'table.load.votable':Receiver()})

        # Standard call_and_wait
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.call_and_wait(self.client2_id, 'test_tag', message)

        # Easy call_and_wait
        self.client1.ecall_and_wait(self.client2_id, 'test_tag', "table.load.votable")
