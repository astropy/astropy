"""
Test the web profile using Python classes that have been adapted to act like a
web client. We can only put a single test here because only one hub can run
with the web profile active, and the user might want to run the tests in
parallel.
"""

import threading
from urllib.request import Request, urlopen

import pytest

from astropy.samp import SAMPHubServer, SAMPIntegratedClient, conf
from astropy.samp.web_profile import CLIENT_ACCESS_POLICY, CROSS_DOMAIN
from astropy.tests.helper import CI
from astropy.utils.data import get_readable_fileobj

from .test_standard_profile import TestStandardProfile as BaseTestStandardProfile
from .web_profile_test_helpers import (
    AlwaysApproveWebProfileDialog,
    SAMPIntegratedWebClient,
)


def setup_module(module):
    conf.use_internet = False


@pytest.mark.skipif(CI, reason="flaky in CI")
class TestWebProfile(BaseTestStandardProfile):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmp_path):
        self.dialog = AlwaysApproveWebProfileDialog()
        t = threading.Thread(target=self.dialog.poll)
        t.start()

        self.tmpdir = str(tmp_path)
        lockfile = str(tmp_path / ".samp")

        self.hub = SAMPHubServer(
            web_profile_dialog=self.dialog, lockfile=lockfile, web_port=0, pool_size=1
        )
        self.hub.start()

        self.client1 = SAMPIntegratedClient()
        self.client1.connect(hub=self.hub, pool_size=1)
        self.client1_id = self.client1.get_public_id()
        self.client1_key = self.client1.get_private_key()

        self.client2 = SAMPIntegratedWebClient()
        self.client2.connect(web_port=self.hub._web_port, pool_size=2)
        self.client2_id = self.client2.get_public_id()
        self.client2_key = self.client2.get_private_key()

    def teardown_method(self):
        if self.client1.is_connected:
            self.client1.disconnect()
        if self.client2.is_connected:
            self.client2.disconnect()

        self.hub.stop()
        self.dialog.stop()

    # The full communication tests are run since TestWebProfile inherits
    # test_main from TestStandardProfile

    def test_web_profile(self):
        # Check some additional queries to the server

        with get_readable_fileobj(
            f"http://localhost:{self.hub._web_port}/crossdomain.xml"
        ) as f:
            assert f.read() == CROSS_DOMAIN

        with get_readable_fileobj(
            f"http://localhost:{self.hub._web_port}/clientaccesspolicy.xml"
        ) as f:
            assert f.read() == CLIENT_ACCESS_POLICY

        # Check headers

        req = Request(f"http://localhost:{self.hub._web_port}/crossdomain.xml")
        req.add_header("Origin", "test_web_profile")
        resp = urlopen(req)

        assert resp.getheader("Access-Control-Allow-Origin") == "test_web_profile"
        assert resp.getheader("Access-Control-Allow-Headers") == "Content-Type"
        assert resp.getheader("Access-Control-Allow-Credentials") == "true"
