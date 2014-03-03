# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..hub_proxy import SAMPHubProxy
from ..client import SAMPClient
from ..integrated_client import SAMPIntegratedClient

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
