# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ....tests.helper import remote_data
from ... import samp

@remote_data
def test_SAMPClient():
    """Test that SAMPClient can be instantiated"""
    proxy = samp.SAMPHubProxy()
    samp.SAMPClient(proxy)
    assert True

  
@remote_data
def test_SAMPIntegratedClient():
    """Test that SAMPIntegratedClient can be instantiated"""
    samp.SAMPIntegratedClient()
    assert True

    
@remote_data
def test_SAMPHubProxy():
    """Test that SAMPHubProxy can be instantiated"""
    samp.SAMPHubProxy()
    assert True
