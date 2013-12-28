# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Integration tests, i.e. tests involving several classes"""

import os
import time
import tempfile

from ....tests.helper import remote_data
from ... import samp


@remote_data
def test_SAMPMsgReplierWrapper():
    """Test that SAMPMsgReplierWrapper can be instantiated"""
    cli = samp.SAMPIntegratedClient()
    samp.SAMPMsgReplierWrapper(cli)


@remote_data
def test_SAMPClient_connect():
    """Test that SAMPClient can connect and register"""

    fileobj, lockfile = tempfile.mkstemp()

    hub = samp.SAMPHubServer(web_profile=False,
                             lockfile=lockfile)

    os.environ['SAMP_HUB'] = "std-lockurl:file://" + os.path.abspath(lockfile)

    proxy = samp.SAMPHubProxy()

    try:
        hub.start()
    except:
        print("Another hub is running.")

    proxy.connect()

    cli = samp.SAMPClient(proxy, name="Client", description="Test Client")
    cli.start()
    cli.register()

    metadata = {"cli.version":"0.01"}
    cli.declare_metadata(metadata)

    cli.unregister()
    cli.stop()
    proxy.disconnect()

    hub.stop()

    del os.environ['SAMP_HUB']  # hacky

    if os.path.exists(lockfile):
        os.remove(lockfile)


@remote_data
class TestSAMPCommunication(object):

    """Test SAMP client-server communications.
    """

    def setup_class(self):

        self.hub = samp.SAMPHubServer(web_profile=False)
        self.hub.start()

        self.myhub1 = samp.SAMPHubProxy()
        self.myhub1.connect()

        self.myhub2 = samp.SAMPHubProxy()
        self.myhub2.connect()

    def teardown_class(self):

        self.myhub1.disconnect()
        self.myhub2.disconnect()
        self.hub.stop()

    def test_2_clients(self):

        # Create a client that uses
        # the passed Hub Proxy
        cli1 = samp.SAMPClient(self.myhub1, name="Client 1",
                               description="Test Client 1")
        # Create another client
        cli2 = samp.SAMPClient(self.myhub2, name="Client 2",
                               description="Test Client 2")
        # Create metadata dictionaries
        metadata1 = {"cli1.version":"0.01"}
        metadata2 = {"cli2.version":"0.25"}

        # Start and register clients
        cli1.start()
        cli1.register()
        cli2.start()
        cli2.register()
        # Declare metadata
        cli1.declare_metadata(metadata1)
        cli2.declare_metadata(metadata2)

        print("\nCLI1", cli1.get_private_key(), cli1.get_public_id(), "\n")
        print("\nCLI2", cli2.get_private_key(), cli2.get_public_id(), "\n")

        # Function called when a notification is received
        def test_receive_notification(private_key, sender_id, mtype,
                                      params, extra):
            print("Notification:", private_key, sender_id, mtype, params,
                  extra, "\n\n")

        # Function called when a call is received
        def test_receive_call(private_key, sender_id, msg_id, mtype, params,
                              extra):
            print("Call:", private_key, sender_id, msg_id, mtype, params,
                  extra, "\n\n")
            self.myhub1.reply(cli1.get_private_key(), msg_id,
                              {"samp.status": samp.SAMP_STATUS_OK,
                               "samp.result": {"txt": "printed"}})

        # Function called when a response is received
        def test_receive_response(private_key, sender_id, msg_id, response):
            print("Response:", private_key, sender_id, msg_id, response, "\n\n")

        # Subscribe Client 1 to "samp.*" and "samp.app.*" MType and bind it to
        # the related functions
        cli1.bind_receive_notification("samp.app.*", test_receive_notification)
        cli1.bind_receive_call("samp.app.*", test_receive_call)

        # Bind Client 2 message-tags received to suitable functions
        cli2.bind_receive_response("my-dummy-print", test_receive_response)
        cli2.bind_receive_response("my-dummy-print-specific", test_receive_response)

        # Client 2 notifies to All "samp.app.echo" MType using myhub
        self.myhub2.notify_all(cli2.get_private_key(),
                              {"samp.mtype": "samp.app.echo",
                               "samp.params": {"txt": "Hello world!"}})

        time.sleep(1)

        # Client 2 calls to All "samp.app.echo" MType using "my-dummy-print"
        # as message-tag
        print(self.myhub2.call_all(cli2.get_private_key(), "my-dummy-print",
                                  {"samp.mtype": "samp.app.echo",
                                   "samp.params": {"txt": "Hello world!"}}),
              "\n\n")

        time.sleep(1)

        # Client 2 calls "samp.app.echo" MType on Client 1 tagging it as
        # "my-dummy-print-specific"
        try:
            print(cli2.hub.call(cli2.get_private_key(),
                                cli1.get_public_id(),
                                "my-dummy-print-specific",
                                {"samp.mtype": "samp.app.echo",
                                 "samp.params": {"txt": "Hello Cli 1!"}}),
                  "\n\n")

        except samp.SAMPProxyError as e:
            print("Error (%s): %s" % (e.faultCode, e.faultString))

        time.sleep(1)

        # Function called to test synchronous calls
        def test_receive_sync_call(private_key, sender_id, msg_id, mtype,
                                   params, extra):
            print("SYNC Call:", sender_id, msg_id, mtype, params, extra, "\n\n")
            time.sleep(1)
            self.myhub1.reply(cli1.get_private_key(), msg_id,
                              {"samp.status": samp.SAMP_STATUS_OK,
                               "samp.result": {"txt": "printed sync"}})
            return ""

        # Bind test MType for sync calls
        cli1.bind_receive_call("samp.test", test_receive_sync_call)

        now = time.time()

        print("SYNCRO --->\n\n")
        try:
            # Sync call
            print(self.myhub2.call_and_wait(cli2.get_private_key(),
                                          cli1.get_public_id(),
                                          {"samp.mtype": "samp.test",
                                           "samp.params": {"txt":
                                                           "Hello SYNCRO Cli 1!"}},
                                          "10"), "\n\n")
        except samp.SAMPProxyError as e:
            # If timeout expires than a SAMPProxyError is returned
            print("Error (%s): %s" % (e.faultCode, e.faultString))

        print("<------SYNCRO (%d) \n\n" % (time.time() - now))

        time.sleep(6)

        cli1.unregister()
        cli1.stop()
        cli2.unregister()
        cli2.stop()