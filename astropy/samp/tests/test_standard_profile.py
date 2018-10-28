import ssl
import tempfile

import pytest

from ...utils.data import get_pkg_data_filename

from ..hub import SAMPHubServer
from ..integrated_client import SAMPIntegratedClient
from ..errors import SAMPProxyError

# By default, tests should not use the internet.
from .. import conf

from .test_helpers import random_params, Receiver, assert_output, TEST_REPLY


def setup_module(module):
    conf.use_internet = False


class TestStandardProfile:

    @property
    def hub_init_kwargs(self):
        return {}

    @property
    def client_init_kwargs(self):
        return {}

    @property
    def client_connect_kwargs(self):
        return {}

    def setup_method(self, method):

        self.tmpdir = tempfile.mkdtemp()

        self.hub = SAMPHubServer(web_profile=False, mode='multiple', pool_size=1,
                                 **self.hub_init_kwargs)
        self.hub.start()

        self.client1 = SAMPIntegratedClient(**self.client_init_kwargs)
        self.client1.connect(hub=self.hub, pool_size=1, **self.client_connect_kwargs)

        self.client2 = SAMPIntegratedClient(**self.client_init_kwargs)
        self.client2.connect(hub=self.hub, pool_size=1, **self.client_connect_kwargs)

    def teardown_method(self, method):

        if self.client1.is_connected:
            self.client1.disconnect()
        if self.client2.is_connected:
            self.client2.disconnect()

        self.hub.stop()

    def test_main(self):

        self.client1_id = self.client1.get_public_id()
        self.client2_id = self.client2.get_public_id()

        self.metadata1 = {"samp.name": "Client 1",
                          "samp.description.text": "Client 1 Description",
                          "client.version": "1.1"}

        self.metadata2 = {"samp.name": "Client 2",
                          "samp.description.text": "Client 2 Description",
                          "client.version": "1.2"}

        # Check that the clients are connected

        assert self.client1.is_connected
        assert self.client2.is_connected

        # Check that ping works

        self.client1.ping()
        self.client2.ping()

        # Check that get_registered_clients works as expected.

        assert self.client1_id not in self.client1.get_registered_clients()
        assert self.client2_id in self.client1.get_registered_clients()
        assert self.client1_id in self.client2.get_registered_clients()
        assert self.client2_id not in self.client2.get_registered_clients()

        # Check that get_metadata works as expected

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

        # Check that, without subscriptions, sending a notification from one
        # client to another raises an error.

        message = {}
        message['samp.mtype'] = "table.load.votable"
        message['samp.params'] = {}

        with pytest.raises(SAMPProxyError):
            self.client1.notify(self.client2_id, message)

        # Check that there are no currently active subscriptions

        assert self.client1.get_subscribed_clients('table.load.votable') == {}
        assert self.client2.get_subscribed_clients('table.load.votable') == {}

        # We now test notifications and calls

        rec1 = Receiver(self.client1)
        rec2 = Receiver(self.client2)

        self.client2.bind_receive_notification('table.load.votable',
                                               rec2.receive_notification)

        self.client2.bind_receive_call('table.load.votable',
                                       rec2.receive_call)

        self.client1.bind_receive_response('test-tag', rec1.receive_response)

        # Check resulting subscriptions

        assert self.client1.get_subscribed_clients('table.load.votable') == {self.client2_id: {}}
        assert self.client2.get_subscribed_clients('table.load.votable') == {}

        assert 'table.load.votable' in self.client1.get_subscriptions(self.client2_id)
        assert 'table.load.votable' in self.client2.get_subscriptions(self.client2_id)

        # Once we have finished with the calls and notifications, we will
        # check the data got across correctly.

        # Test notify

        params = random_params(self.tmpdir)
        self.client1.notify(self.client2.get_public_id(),
                            {'samp.mtype': 'table.load.votable',
                             'samp.params': params})

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        params = random_params(self.tmpdir)
        self.client1.enotify(self.client2.get_public_id(),
                             "table.load.votable", **params)

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        # Test notify_all

        params = random_params(self.tmpdir)
        self.client1.notify_all({'samp.mtype': 'table.load.votable',
                                 'samp.params': params})

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        params = random_params(self.tmpdir)
        self.client1.enotify_all("table.load.votable", **params)

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        # Test call

        params = random_params(self.tmpdir)
        self.client1.call(self.client2.get_public_id(), 'test-tag',
                            {'samp.mtype': 'table.load.votable',
                             'samp.params': params})

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        params = random_params(self.tmpdir)
        self.client1.ecall(self.client2.get_public_id(), 'test-tag',
                           "table.load.votable", **params)

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        # Test call_all

        params = random_params(self.tmpdir)
        self.client1.call_all('tag1',
                              {'samp.mtype': 'table.load.votable',
                               'samp.params': params})

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        params = random_params(self.tmpdir)
        self.client1.ecall_all('tag2',
                               "table.load.votable", **params)

        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        # Test call_and_wait

        params = random_params(self.tmpdir)
        result = self.client1.call_and_wait(self.client2.get_public_id(),
                                            {'samp.mtype': 'table.load.votable',
                                             'samp.params': params}, timeout=5)

        assert result == TEST_REPLY
        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        params = random_params(self.tmpdir)
        result = self.client1.ecall_and_wait(self.client2.get_public_id(),
                                             "table.load.votable", timeout=5, **params)

        assert result == TEST_REPLY
        assert_output('table.load.votable', self.client2.get_private_key(),
                      self.client1_id, params, timeout=60)

        # TODO: check that receive_response received the right data
