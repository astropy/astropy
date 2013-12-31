# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..errors import SAMPHubError, SAMPClientError, SAMPProxyError

# By default, tests should not use the internet.
from ..utils import ALLOW_INTERNET
ALLOW_INTERNET.set(False)

def test_SAMPHubError():
    """Test that SAMPHubError can be instantiated"""
    SAMPHubError("test")


def test_SAMPClientError():
    """Test that SAMPClientError can be instantiated"""
    SAMPClientError("test")


def test_SAMPProxyError():
    """Test that SAMPProxyError can be instantiated"""
    SAMPProxyError("test", "any")
