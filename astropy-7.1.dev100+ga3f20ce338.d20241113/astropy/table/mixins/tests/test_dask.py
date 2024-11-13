import numpy as np
import pytest
from numpy.testing import assert_equal

from astropy.table import Table

da = pytest.importorskip("dask.array")


class TestDaskHandler:
    def setup_method(self, method):
        self.t = Table()
        self.t["a"] = da.arange(10)

    def test_add_row(self):
        self.t.add_row(self.t[0])
        assert_equal(self.t["a"].compute(), np.hstack([np.arange(10), 0]))

    def test_get_column(self):
        assert isinstance(self.t["a"], da.Array)
        assert_equal(self.t["a"].compute(), np.arange(10))

    def test_slicing_row_single(self):
        sub = self.t[5]
        assert isinstance(sub["a"], da.Array)
        assert not hasattr(sub["a"], "info")  # should be a plain dask array
        assert sub["a"].compute() == 5

    def test_slicing_row_range(self):
        sub = self.t[5:]
        assert isinstance(sub["a"], da.Array)
        assert hasattr(sub["a"], "info")  # should be a mixin column
        assert_equal(sub["a"].compute(), np.arange(5, 10))

    def test_slicing_column_range(self):
        sub = self.t[("a",)]
        assert isinstance(sub["a"], da.Array)
        assert hasattr(sub["a"], "info")  # should be a mixin column
        assert_equal(sub["a"].compute(), np.arange(10))

    def test_pformat(self):
        assert self.t.pformat() == [
            " a ",
            "---",
            "  0",
            "  1",
            "  2",
            "  3",
            "  4",
            "  5",
            "  6",
            "  7",
            "  8",
            "  9",
        ]

    def test_info_preserved(self):
        self.t["a"].info.description = "A dask column"

        sub = self.t[1:3]
        assert sub["a"].info.name == "a"
        assert sub["a"].info.description == "A dask column"

        col = self.t["a"].copy()
        assert col.info.name == "a"
        assert col.info.description == "A dask column"

        self.t.add_row(self.t[0])
        assert self.t["a"].info.name == "a"
        assert self.t["a"].info.description == "A dask column"
