# Licensed under a 3-clause BSD style license - see LICENSE.rst

# By default, tests should not use the internet.
from .. import conf
from ..errors import SAMPHubError, SAMPProxyError, SAMPClientError


def setup_module(module):
    conf.use_internet = False


def test_SAMPHubError():
    """Test that SAMPHubError can be instantiated"""
    SAMPHubError("test")


def test_SAMPClientError():
    """Test that SAMPClientError can be instantiated"""
    SAMPClientError("test")


def test_SAMPProxyError():
    """Test that SAMPProxyError can be instantiated"""
    SAMPProxyError("test", "any")
