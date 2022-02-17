# Licensed under a 3-clause BSD style license - see LICENSE.rst

import time

import pytest

from astropy.samp import conf
from astropy.samp.hub import SAMPHubServer


def setup_module(module):
    conf.use_internet = False


def test_SAMPHubServer():
    """Test that SAMPHub can be instantiated"""
    SAMPHubServer(web_profile=False, mode='multiple', pool_size=1)


def test_SAMPHubServer_run():
    """Test that SAMPHub can be run"""
    hub = SAMPHubServer(web_profile=False, mode='multiple', pool_size=1)
    hub.start()
    time.sleep(1)
    hub.stop()


@pytest.mark.slow
def test_SAMPHubServer_run_repeated():
    """
    Test that SAMPHub can be restarted after it has been stopped, including
    when web profile support is enabled.
    """

    hub = SAMPHubServer(web_profile=True, mode='multiple', pool_size=1)
    hub.start()
    time.sleep(1)
    hub.stop()
    time.sleep(1)
    hub.start()
    time.sleep(1)
    hub.stop()
