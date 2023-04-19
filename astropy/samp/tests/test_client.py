# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

# By default, tests should not use the internet.
from astropy.samp import SAMPWarning, conf
from astropy.samp.client import SAMPClient
from astropy.samp.hub import SAMPHubServer
from astropy.samp.hub_proxy import SAMPHubProxy
from astropy.samp.integrated_client import SAMPIntegratedClient


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


@pytest.fixture
def samp_hub(request):
    """A fixture that can be used by client tests that require a HUB."""
    my_hub = SAMPHubServer()
    my_hub.start()
    request.addfinalizer(my_hub.stop)


def test_SAMPIntegratedClient_notify_all(samp_hub):
    """Test that SAMP returns a warning if no receiver got the message."""
    client = SAMPIntegratedClient()
    client.connect()
    message = {"samp.mtype": "coverage.load.moc.fits"}
    with pytest.warns(SAMPWarning):
        client.notify_all(message)
    client.disconnect()


def test_reconnect(samp_hub):
    """Test that SAMPIntegratedClient can reconnect.
    This is a regression test for bug [#2673]
    https://github.com/astropy/astropy/issues/2673
    """
    my_client = SAMPIntegratedClient()
    my_client.connect()
    my_client.disconnect()
    my_client.connect()
