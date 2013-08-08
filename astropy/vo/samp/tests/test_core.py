from __future__ import print_function, division

import time

from ....tests.helper import pytest
from ... import samp


def test_SAMPHubError():
    """Test that SAMPHubError can be instantiated"""
    samp.SAMPHubError("test")
    assert True
    
def test_SAMPClientError():
    """Test that SAMPClientError can be instantiated"""
    samp.SAMPClientError("test")
    
def test_SAMPProxyError():
    """Test that SAMPProxyError can be instantiated"""
    samp.SAMPProxyError("test", "any")
    assert True

def test_SAMPLog():
    """Test that SAMPLog can be instantiated"""
    samp.SAMPLog()
    assert True

def test_SAMPHubServer():
    """Test that SAMPHub can be instantiated"""
    samp.SAMPHubServer()
    assert True
    
def test_SAMPHubProxy():
    """Test that SAMPHubProxy can be instantiated"""
    samp.SAMPHubProxy()
    assert True
    
def test_SAMPClient():
    """Test that SAMPClient can be instantiated"""
    proxy = samp.SAMPHubProxy()
    samp.SAMPClient(proxy)
    assert True
  
def test_SAMPIntegratedClient():
    """Test that SAMPIntegratedClient can be instantiated"""
    samp.SAMPIntegratedClient()
    assert True
    
def test_SAMPMsgReplierWrapper():
    """Test that SAMPMsgReplierWrapper can be instantiated"""
    cli = samp.SAMPIntegratedClient()
    samp.SAMPMsgReplierWrapper(cli)
    assert True
    
    
def test_SAMPHubServer_run():
    """Test that SAMPHub can be run"""
    
    try:
        hub = samp.SAMPHubServer(web_profile=False)
        hub.start()
        time.sleep(1)
        hub.stop()
    except:
        assert False
        
    assert True
    
    
def test_SAMPClient_connect():
    """Test that SAMPClient can connect and register"""
    
    hub = samp.SAMPHubServer(web_profile=False)
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
    cli.declareMetadata(metadata)
    
    cli.unregister()
    cli.stop()
    proxy.disconnect()

    hub.stop()
    
    
    
class TestSAMPCommunication(object):
    """Test SAMP client-server communications.
    """
    def setup_class(self):
    
        self.hub = samp.SAMPHubServer()
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
        cli1.declareMetadata(metadata1)
        cli2.declareMetadata(metadata2)
        
        print("\nCLI1", cli1.getPrivateKey(), cli1.getPublicId(), "\n")
        print("\nCLI2", cli2.getPrivateKey(), cli2.getPublicId(), "\n")
        
        
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
            myhub1.reply(cli1.getPrivateKey(), msg_id,
                         {"samp.status": samp.SAMP_STATUS_OK,
                          "samp.result": {"txt": "printed"}})
        
        # Function called when a response is received
        def test_receive_response(private_key, sender_id, msg_id, response):
            print("Response:", private_key, sender_id, msg_id, response, "\n\n")
        
        # Subscribe Client 1 to "samp.*" and "samp.app.*" MType and bind it to
        # the related functions
        cli1.bindReceiveNotification("samp.app.*", test_receive_notification)
        cli1.bindReceiveCall("samp.app.*", test_receive_call)
        
        # Bind Client 2 message-tags received to suitable functions
        cli2.bindReceiveResponse("my-dummy-print", test_receive_response)
        cli2.bindReceiveResponse("my-dummy-print-specific", test_receive_response)
        
        # Client 2 notifies to All "samp.app.echo" MType using myhub
        self.myhub2.notifyAll(cli2.getPrivateKey(),
                              {"samp.mtype": "samp.app.echo",
                               "samp.params": {"txt": "Hello world!"}})
        
        time.sleep(1)
        
        # Client 2 calls to All "samp.app.echo" MType using "my-dummy-print"
        # as message-tag
        print(self.myhub2.callAll(cli2.getPrivateKey(), "my-dummy-print",
                                  {"samp.mtype": "samp.app.echo",
                                   "samp.params": {"txt": "Hello world!"}}),
              "\n\n")
        
        time.sleep(1)
        
        # Client 2 calls "samp.app.echo" MType on Client 1 tagging it as 
        # "my-dummy-print-specific"
        try:
            print(cli2.hub.call(cli2.getPrivateKey(),
                                cli1.getPublicId(),
                                "my-dummy-print-specific",
                                {"samp.mtype": "samp.app.echo",
                                 "samp.params": {"txt": "Hello Cli 1!"}}),
                  "\n\n")
            
        except SAMPProxyError as e:
            print("Error (%s): %s" % (e.faultCode, e.faultString))
        
          
        time.sleep(1)
          
        # Function called to test synchronous calls
        def test_receive_sync_call(private_key, sender_id, msg_id, mtype,
                                   params, extra):
            print("SYNC Call:", sender_id, msg_id, mtype, params, extra, "\n\n")
            time.sleep(1)
            self.myhub1.reply(cli1.getPrivateKey(), msg_id,
                              {"samp.status": samp.SAMP_STATUS_OK,
                               "samp.result": {"txt": "printed sync"}})
            return ""
        
        # Bind test MType for sync calls
        cli1.bindReceiveCall("samp.test", test_receive_sync_call)  
        
        now = time.time()
        
        print("SYNCRO --->\n\n")
        try:
            # Sync call
            print(self.myhub2.callAndWait(cli2.getPrivateKey(), 
                                          cli1.getPublicId(),
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
        
        assert True
        
    