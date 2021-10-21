import pytest
import numpy as np
from numpy.testing import assert_equal

from astropy.table import Table

da = pytest.importorskip('dask.array')


class TestDaskHandler:

    def setup_method(self, method):
        self.t = Table()
        self.t['a'] = da.arange(10)

    def test_get_column(self):
        assert isinstance(self.t['a'], da.Array)
        assert_equal(self.t['a'].compute(), np.arange(10))

    def test_slicing(self):
        sub = self.t[5:]
        assert isinstance(sub['a'], da.Array)
        assert_equal(sub['a'].compute(), np.arange(5, 10))

    def test_pformat(self):
        assert self.t.pformat_all() == [' a ', '---', '  0', '  1', '  2',
                                        '  3', '  4', '  5', '  6', '  7',
                                        '  8', '  9']
