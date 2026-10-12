# Licensed under a 3-clause BSD style license - see LICENSE.rst

import platform
import warnings

import pytest

# By default, tests should not use the internet.
from astropy.samp import SAMPProxyError, SAMPWarning, conf
from astropy.samp.client import SAMPClient
from astropy.samp.hub import SAMPHubServer
from astropy.samp.hub_proxy import SAMPHubProxy
from astropy.samp.integrated_client import SAMPIntegratedClient


def setup_module(module):
    conf.use_internet = False


@pytest.mark.skipif(platform.system() == "Darwin", reason="Takes too long on OSX")
def test_SAMPHubProxy():
    """Test that SAMPHubProxy can be instantiated"""
    SAMPHubProxy()


@pytest.mark.skipif(platform.system() == "Darwin", reason="Takes too long on OSX")
@pytest.mark.slow
def test_SAMPClient():
    """Test that SAMPClient can be instantiated"""
    proxy = SAMPHubProxy()
    SAMPClient(proxy)


@pytest.mark.skipif(platform.system() == "Darwin", reason="Takes too long on OSX")
def test_SAMPIntegratedClient():
    """Test that SAMPIntegratedClient can be instantiated"""
    SAMPIntegratedClient()


@pytest.fixture
def samp_hub():
    """A fixture that can be used by client tests that require a HUB."""
    my_hub = SAMPHubServer()
    with warnings.catch_warnings():
        # Turn SAMPWarning into an error so that it can be seen in a different context
        # The error will be caught internally twice and ultimately raised as a SAMPProxyError
        warnings.simplefilter("error", category=SAMPWarning)
        my_hub.start()
        yield
        my_hub.stop()


@pytest.mark.skipif(platform.system() == "Darwin", reason="Takes too long on OSX")
@pytest.mark.filterwarnings("ignore:unclosed <socket:ResourceWarning")
def test_SAMPIntegratedClient_notify_all(samp_hub):
    """Test that SAMP returns a warning if no receiver got the message."""
    client = SAMPIntegratedClient()
    client.connect()
    message = {"samp.mtype": "coverage.load.moc.fits"}
    # SAMPWarning will be internally converted to a SAMPProxyError
    with pytest.raises(
        SAMPProxyError,
        match=r".*SAMPWarning.*:No client was able to receive this message",
    ):
        client.notify_all(message)
    client.disconnect()


@pytest.mark.skipif(platform.system() == "Darwin", reason="Takes too long on OSX")
def test_reconnect(samp_hub):
    """Test that SAMPIntegratedClient can reconnect.
    This is a regression test for bug [#2673]
    https://github.com/astropy/astropy/issues/2673
    """
    my_client = SAMPIntegratedClient()
    my_client.connect()
    my_client.disconnect()
    my_client.connect()
