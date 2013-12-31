import os
import time
import tempfile

from ....tests.helper import pytest

from ..hub import SAMPHubServer
from ..integrated_client import SAMPIntegratedClient
from ..constants import SAMP_STATUS_OK

# By default, tests should not use the internet.
from ..utils import ALLOW_INTERNET
ALLOW_INTERNET.set(False)

# TODO:
# - reply/ereply
# - bind_receive_response
# - make proxy accept a lockfile name


class Receiver(object):

    def __init__(self):
        self.received = False

    def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
        self.params = params
        self.received = True


class TestIntegratedClient(object):

    def setup_method(self, method):

        fileobj, self.lockfile = tempfile.mkstemp()

        self.hub = SAMPHubServer(web_profile=False,
                                 lockfile=self.lockfile)
        self.hub.start()

        os.environ['SAMP_HUB'] = "std-lockurl:file://" + os.path.abspath(self.lockfile)

        self.client1 = SAMPIntegratedClient()
        self.client1.connect()
        self.client1_id = self.client1.get_public_id()

        self.client2 = SAMPIntegratedClient()
        self.client2.connect()
        self.client2_id = self.client2.get_public_id()

        del os.environ['SAMP_HUB']  # hacky

        self.metadata1 = {"samp.name": "Client 1",
                          "samp.description.text": "Client 1 Description",
                          "client.version": "1.1"}

        self.metadata2 = {"samp.name": "Client 2",
                          "samp.description.text": "Client 2 Description",
                          "client.version": "1.2"}

    def teardown_method(self, method):

        if self.client1.is_connected:
            self.client1.disconnect()
        if self.client2.is_connected:
            self.client2.disconnect()

        self.hub.stop()

        if os.path.exists(self.lockfile):
            os.remove(self.lockfile)

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
        self.client1.declare_subscriptions({'table.load.votable': rec})

        assert self.client2.get_subscribed_clients('table.load.votable') == {self.client1_id: {}}
        assert self.client1.get_subscribed_clients('table.load.votable') == {}  # a bit strange, make sure this is clear in docstring

        assert 'table.load.votable' in self.client1.get_subscriptions(self.client1_id)
        assert 'table.load.votable' in self.client2.get_subscriptions(self.client1_id)

    def test_is_connected(self):
        assert self.client1.is_connected
        assert self.client2.is_connected

    @pytest.mark.xfail  # need a better exception than current 'Fault'
    def test_no_mtype(self):
        message = {}
        with pytest.raises(SAMPError):
            self.client1.notify(self.client2_id, message)

    def test_notify(self):

        self.client2.declare_subscriptions({'table.load.votable': Receiver()})

        # Standard notify
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.notify(self.client2_id, message)

        # Easy notify
        self.client1.enotify(self.client2_id, "table.load.votable")

    def test_notify_extra(self):
        """
        Ensures that passing extra_kws works correctly when passed to _format_easy_msg
        """
        self.client2.declare_subscriptions({'table.load.votable': Receiver()})
        self.client1.enotify(self.client2_id, "table.load.votable",
                             extra_kws={'simple.example': 'test'})

    def test_notify_all(self):

        self.client2.declare_subscriptions({'table.load.votable': Receiver()})

        # Standard notify
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.notify_all(message)

        # Easy notify
        self.client1.enotify_all("table.load.votable")

    def test_call(self):

        self.client2.declare_subscriptions({'table.load.votable': Receiver()})

        # Standard call
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.call(self.client2_id, 'test_tag', message)

        # Easy call
        self.client1.ecall(self.client2_id, 'test_tag', "table.load.votable")

    def test_call_all(self):

        self.client2.declare_subscriptions({'table.load.votable': Receiver()})

        # Standard call
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.call_all('test_tag', message)

        # Easy call
        self.client1.ecall_all('test_tag', "table.load.votable")

    @pytest.mark.xfail
    def test_call_and_wait(self):

        self.client2.declare_subscriptions({'table.load.votable': Receiver()})

        # Standard call_and_wait
        message = {}
        message['samp.mtype'] = "table.load.votable"
        self.client1.call_and_wait(self.client2_id, 'test_tag', message)

        # Easy call_and_wait
        self.client1.ecall_and_wait(self.client2_id, 'test_tag', "table.load.votable")

    def test_bind_receive_notification(self):

        class TestReceiver(object):

            def test_receive_notification(self, private_key, sender_id, mtype,
                                          params, extra):
                self.private_key = private_key
                self.sender_id = sender_id
                self.mtype = mtype
                self.params = params
                self.extra = extra

        rec = TestReceiver()

        self.client2.bind_receive_notification('test.message', rec.test_receive_notification)

        self.client1.enotify(self.client2_id, "test.message", a=1, b='a')

        private_key = self.client2.get_private_key()

        # self.client2.unbind_receive_notification("test.message")  # TODO: fix

        self.client2.disconnect()  # TODO: should not be needed

        time.sleep(1.)

        assert rec.private_key == private_key
        assert rec.sender_id == self.client1_id
        assert rec.mtype == 'test.message'
        assert rec.params == {'a': 1, 'b': 'a'}

    def test_bind_receive_call(self):

        class TestReceiver(object):

            def test_receive_call(self, private_key, sender_id, msg_id, mtype,
                                  params, extra):
                self.private_key = private_key
                self.sender_id = sender_id
                self.msg_id = msg_id
                self.mtype = mtype
                self.params = params
                self.extra = extra

        rec = TestReceiver()

        self.client2.bind_receive_call('test.call', rec.test_receive_call)

        self.client1.ecall(self.client2_id, "message.id", "test.call", a=1, b='a')

        private_key = self.client2.get_private_key()

        # self.client2.unbind_receive_call("test.call")  # TODO: fix

        self.client2.disconnect()  # TODO: should not be needed

        time.sleep(1.)

        assert rec.private_key == private_key
        assert rec.sender_id == self.client1_id
        assert rec.msg_id.endswith('message.id')
        assert rec.mtype == 'test.call'
        assert rec.params == {'a': 1, 'b': 'a'}

    def test_ecall_and_wait(self):

        class TestReceiver(object):

            def test_receive_call(self, private_key, sender_id, msg_id, mtype,
                                  params, extra):
                self.private_key = private_key
                self.sender_id = sender_id
                self.msg_id = msg_id
                self.mtype = mtype
                self.params = params
                self.extra = extra
                self.client.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
                                  "samp.result": {"txt": "test"}})

        rec = TestReceiver()
        rec.client = self.client2

        self.client2.bind_receive_call('test.call', rec.test_receive_call)

        result = self.client1.ecall_and_wait(self.client2_id, "test.call", timeout=5., a=1, b='a')

        private_key = self.client2.get_private_key()

        # self.client2.unbind_receive_call("test.call")  # TODO: fix

        self.client2.disconnect()  # TODO: should not be needed

        time.sleep(1.)

        assert rec.private_key == private_key
        assert rec.sender_id == self.client1_id
        assert rec.mtype == 'test.call'
        assert rec.params == {'a': 1, 'b': 'a'}

        assert result['samp.status'] == SAMP_STATUS_OK
        assert result['samp.result'] == {'txt': 'test'}

    def test_del(self):
        self.client1.__del__()
        self.client2.__del__()
