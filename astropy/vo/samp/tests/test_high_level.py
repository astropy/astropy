import tempfile

import numpy as np

from ..hub import SAMPHubServer
from ..integrated_client import SAMPIntegratedClient
from ..high_level import send, receive, list_clients

from ....table import Table
from ....nddata import NDData
from ....io import fits


class TestHighLevel(object):

    def setup_method(self, method):

        self.tmpdir = tempfile.mkdtemp()

        self.hub = SAMPHubServer(web_profile=False, mode='multiple')
        self.hub.start()

        self.client1 = SAMPIntegratedClient()
        self.client1.connect(hub=self.hub)

    def test_high_level(self):

        # Test sending Table
        t = Table([[1, 2, 3], [4, 5, 6]])
        send(t, 'test_table_1', hub=self.hub)

        # Test sending NDData
        # d = NDData(np.ones((128, 128)))
        # send(d, 'test_image_1', hub=self.hub)

        # Test sending numpy regular array
        a = np.ones((128, 128))
        send(a, 'test_image_2', hub=self.hub)

        # Test sending numpy structured array
        x = np.array([(1, 2), (5, 6)], dtype=[('a', float), ('b', float)])
        send(x, 'test_table_2', hub=self.hub)

        # Test HDU objects

        hdu1 = fits.ImageHDU(a)
        send(hdu1, 'test_image_3', hub=self.hub)

        hdu2 = fits.PrimaryHDU(a)
        send(hdu2, 'test_image_4', hub=self.hub)

        hdu3 = fits.BinTableHDU(x)
        send(hdu3, 'test_table_3', hub=self.hub)

        hdu4 = fits.BinTableHDU(x)
        send(hdu4, 'test_table_4', hub=self.hub)
