# Licensed under a 3-clause BSD style license - see LICENSE.rst
import time
from ....tests.helper import remote_data
from ... import samp

@remote_data
def test_SAMPHubServer():
    """Test that SAMPHub can be instantiated"""
    samp.SAMPHubServer(web_profile=False)
    assert True

    
@remote_data
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
