# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.table import QTable, hstack, join, vstack
from astropy.utils.compat.optional_deps import HAS_H5PY
from astropy.utils.masked import Masked

from .test_masked import assert_masked_equal

FILE_FORMATS = ["ecsv", "fits"]
if HAS_H5PY:
    FILE_FORMATS.append("h5")


class MaskedArrayTableSetup:
    @classmethod
    def setup_arrays(self):
        self.a = np.array([3.0, 5.0, 0.0])
        self.mask_a = np.array([True, False, False])

    @classmethod
    def setup_class(self):
        self.setup_arrays()
        self.ma = Masked(self.a, mask=self.mask_a)
        self.ma.info.format = ".1f"
        self.t = QTable([self.ma], names=["ma"])


class MaskedQuantityTableSetup(MaskedArrayTableSetup):
    @classmethod
    def setup_arrays(self):
        self.a = np.array([3.0, 5.0, 0.0]) << u.m
        self.mask_a = np.array([True, False, False])


class TestMaskedArrayTable(MaskedArrayTableSetup):
    def test_table_initialization(self):
        assert_array_equal(self.t["ma"].unmasked, self.a)
        assert_array_equal(self.t["ma"].mask, self.mask_a)
        assert repr(self.t).splitlines()[-3:] == [
            "    ———",
            "    5.0",
            "    0.0",
        ]

    def test_info_basics(self):
        assert self.t["ma"].info.name == "ma"
        assert "serialize_method" in self.t["ma"].info.attr_names
        t2 = self.t.copy()
        t2["ma"].info.format = ".2f"
        t2["ma"].info.serialize_method["fits"] = "nonsense"
        assert repr(t2).splitlines()[-3:] == [
            "    ———",
            "   5.00",
            "   0.00",
        ]
        # Check that if we slice, things get copied over correctly.
        t3 = t2[:2]
        assert t3["ma"].info.name == "ma"
        assert t3["ma"].info.format == ".2f"
        assert "serialize_method" in t3["ma"].info.attr_names
        assert t3["ma"].info.serialize_method["fits"] == "nonsense"

    @pytest.mark.parametrize("file_format", FILE_FORMATS)
    def test_table_write(self, file_format, tmp_path):
        name = tmp_path / f"a.{file_format}"
        kwargs = {}
        if file_format == "h5":
            kwargs["path"] = "trial"
            kwargs["serialize_meta"] = True

        self.t.write(name, **kwargs)
        t2 = QTable.read(name)
        assert isinstance(t2["ma"], self.ma.__class__)
        assert np.all(t2["ma"] == self.ma)
        assert np.all(t2["ma"].mask == self.mask_a)
        if file_format == "fits":
            # Imperfect roundtrip through FITS native format description.
            assert self.t["ma"].info.format in t2["ma"].info.format
        else:
            assert t2["ma"].info.format == self.t["ma"].info.format

    @pytest.mark.parametrize("serialize_method", ["data_mask", "null_value"])
    def test_table_write_serialization(self, serialize_method, tmp_path):
        name = tmp_path / "test.ecsv"
        self.t.write(name, serialize_method=serialize_method)
        with open(name) as fh:
            lines = fh.readlines()

        t2 = QTable.read(name)
        assert isinstance(t2["ma"], self.ma.__class__)

        if serialize_method == "data_mask":
            # Will data_mask, we have data and mask columns and should
            # exactly round-trip.
            assert len(lines[-1].split()) == 2
            assert_masked_equal(t2["ma"], self.ma)
        else:
            # With null_value we have just a data column with null values
            # marked, so we lost information on the data below the mask.
            assert len(lines[-1].split()) == 1
            assert np.all(t2["ma"] == self.ma)
            assert np.all(t2["ma"].mask == self.mask_a)

    def test_non_existing_serialize_method(self, tmp_path):
        name = tmp_path / "bad.ecsv"
        with pytest.raises(ValueError, match="serialize method must be"):
            self.t.write(name, serialize_method="bad_serialize_method")


class TestMaskedQuantityTable(TestMaskedArrayTable, MaskedQuantityTableSetup):
    # Runs tests from TestMaskedArrayTable as well as some extra ones.
    def test_table_operations_requiring_masking(self):
        t1 = self.t
        t2 = QTable({"ma2": Masked([1, 2] * u.m)})
        t12 = hstack([t1, t2], join_type="outer")
        assert np.all(t12["ma"].mask == [True, False, False])
        # 'ma2' is shorter by one so we expect one True from hstack so length matches
        assert np.all(t12["ma2"].mask == [False, False, True])

        t12 = hstack([t1, t2], join_type="inner")
        assert np.all(t12["ma"].mask == [True, False])
        assert np.all(t12["ma2"].mask == [False, False])

        # Vstack tables with different column names. In this case we get masked
        # values
        t12 = vstack([t1, t2], join_type="outer")
        #  ma ma2
        #  m   m
        # --- ---
        #  ——  ——
        # 5.0  ——
        # 0.0  ——
        #  —— 1.0
        #  —— 2.0
        assert np.all(t12["ma"].mask == [True, False, False, True, True])
        assert np.all(t12["ma2"].mask == [True, True, True, False, False])

    def test_table_operations_requiring_masking_auto_promote(self):
        MaskedQuantity = Masked(u.Quantity)
        t1 = QTable({"ma1": [1, 2] * u.m})
        t2 = QTable({"ma2": [3, 4, 5] * u.m})
        t12 = hstack([t1, t2], join_type="outer")
        assert isinstance(t12["ma1"], MaskedQuantity)
        assert np.all(t12["ma1"].mask == [False, False, True])
        assert np.all(t12["ma1"] == [1, 2, 0] * u.m)
        assert not isinstance(t12["ma2"], MaskedQuantity)
        assert isinstance(t12["ma2"], u.Quantity)
        assert np.all(t12["ma2"] == [3, 4, 5] * u.m)

        t12 = hstack([t1, t2], join_type="inner")
        assert isinstance(t12["ma1"], u.Quantity)
        assert not isinstance(t12["ma1"], MaskedQuantity)
        assert isinstance(t12["ma2"], u.Quantity)
        assert not isinstance(t12["ma2"], MaskedQuantity)

        # Vstack tables with different column names. In this case we get masked
        # values
        t12 = vstack([t1, t2], join_type="outer")
        assert np.all(t12["ma1"].mask == [False, False, True, True, True])
        assert np.all(t12["ma2"].mask == [True, True, False, False, False])

        t1["a"] = [1, 2]
        t2["a"] = [1, 3, 4]
        t12 = join(t1, t2, join_type="outer")
        assert np.all(t12["ma1"].mask == [False, False, True, True])
        assert np.all(t12["ma2"].mask == [False, True, False, False])
