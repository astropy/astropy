# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from ..hub_proxy import SAMPHubProxy
from ..client import SAMPClient
from ..integrated_client import SAMPIntegratedClient
from ..hub import SAMPHubServer

# By default, tests should not use the internet.
from .. import conf

def setup_module(module):
    conf.use_internet = False


def test_SAMPHubProxy():
    """Test that SAMPHubProxy can be instantiated"""
    SAMPHubProxy()


def test_SAMPClient():
    """Test that SAMPClient can be instantiated"""
    proxy = SAMPHubProxy()
    SAMPClient(proxy)


def test_SAMPIntegratedClient():
    """Test that SAMPIntegratedClient can be instantiated"""
    SAMPIntegratedClient()


# Not yet included because it seems that there is a problem with
# instantiating multiple hubs.
#def test_SAMPHubServer():
#    """Test that SAMPHubServer can be instantiated"""
#    my_hub = SAMPHubServer()
#    my_hub.start()
#    my_hub.stop()


@pytest.fixture
def samp_hub(request):
    my_hub = SAMPHubServer()
    my_hub.start()
    request.addfinalizer(my_hub.stop)


def test_reconnect(samp_hub):
    """Test that SAMPIntegratedClient can reconnect."""
    my_client = SAMPIntegratedClient()
    my_client.connect()
    my_client.disconnect()
    my_client.connect()

