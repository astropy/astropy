# Licensed under a 3-clause BSD style license - see LICENSE.rst
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
