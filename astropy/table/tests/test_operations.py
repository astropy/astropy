# Licensed under a 3-clause BSD style license - see LICENSE.rst

from collections import OrderedDict
from contextlib import nullcontext

import numpy as np
import pytest

from astropy import table
from astropy import units as u
from astropy.coordinates import (
    BaseRepresentationOrDifferential,
    CartesianRepresentation,
    SkyCoord,
    StokesCoord,
    SphericalRepresentation,
    UnitSphericalRepresentation,
    search_around_3d,
)
from astropy.coordinates.earth import EarthLocation
from astropy.coordinates.tests.helper import skycoord_equal
from astropy.coordinates.tests.test_representation import representation_equal
from astropy.table import Column, MaskedColumn, QTable, Table, TableMergeError
from astropy.table.operations import _get_out_class, join_distance, join_skycoord
from astropy.time import Time, TimeDelta
from astropy.units.quantity import Quantity
from astropy.utils import metadata
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.metadata import MergeConflictError


def sort_eq(list1, list2):
    return sorted(list1) == sorted(list2)


def check_mask(col, exp_mask):
    """Check that col.mask == exp_mask"""
    if hasattr(col, "mask"):
        # Coerce expected mask into dtype of col.mask. In particular this is
        # needed for types like EarthLocation where the mask is a structured
        # array.
        exp_mask = np.array(exp_mask).astype(col.mask.dtype)
        out = np.all(col.mask == exp_mask)
    else:
        # With no mask the check is OK if all the expected mask values
        # are False (i.e. no auto-conversion to MaskedQuantity if it was
        # not required by the join).
        out = np.all(exp_mask == False)  # noqa: E712
    return out


class TestJoin:
    def _setup(self, t_cls=Table):
        lines1 = [
            " a   b   c ",
            "  0 foo  L1",
            "  1 foo  L2",
            "  1 bar  L3",
            "  2 bar  L4",
        ]
        lines2 = [
            " a   b   d ",
            "  1 foo  R1",
            "  1 foo  R2",
            "  2 bar  R3",
            "  4 bar  R4",
        ]
        self.t1 = t_cls.read(lines1, format="ascii")
        self.t2 = t_cls.read(lines2, format="ascii")
        self.t3 = t_cls(self.t2, copy=True)

        self.t1.meta.update(OrderedDict([("b", [1, 2]), ("c", {"a": 1}), ("d", 1)]))
        self.t2.meta.update(OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)]))
        self.t3.meta.update(OrderedDict([("b", 3), ("c", [1, 2]), ("d", 2), ("a", 1)]))

        self.meta_merge = OrderedDict(
            [
                ("b", [1, 2, 3, 4]),
                ("c", {"a": 1, "b": 1}),
                ("d", 1),
                ("a", 1),
            ]
        )

    def test_table_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.join(self.t1, self.t2, join_type="inner")
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self, operation_table_type):
        self._setup(operation_table_type)

        with pytest.warns(metadata.MergeConflictWarning) as w:
            out = table.join(self.t1, self.t3, join_type="inner")
        assert len(w) == 3

        assert out.meta == self.t3.meta

        with pytest.warns(metadata.MergeConflictWarning) as w:
            out = table.join(
                self.t1, self.t3, join_type="inner", metadata_conflicts="warn"
            )
        assert len(w) == 3

        assert out.meta == self.t3.meta

        out = table.join(
            self.t1, self.t3, join_type="inner", metadata_conflicts="silent"
        )

        assert out.meta == self.t3.meta

        with pytest.raises(MergeConflictError):
            out = table.join(
                self.t1, self.t3, join_type="inner", metadata_conflicts="error"
            )

        with pytest.raises(ValueError):
            out = table.join(
                self.t1, self.t3, join_type="inner", metadata_conflicts="nonsense"
            )

    def test_both_unmasked_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Basic join with default parameters (inner join on common keys)
        t12 = table.join(t1, t2)
        assert type(t12) is operation_table_type
        assert type(t12["a"]) is type(t1["a"])
        assert type(t12["b"]) is type(t1["b"])
        assert type(t12["c"]) is type(t1["c"])
        assert type(t12["d"]) is type(t2["d"])
        assert t12.masked is False
        assert sort_eq(
            t12.pformat(),
            [
                " a   b   c   d ",
                "--- --- --- ---",
                "  1 foo  L2  R1",
                "  1 foo  L2  R2",
                "  2 bar  L4  R3",
            ],
        )
        # Table meta merged properly
        assert t12.meta == self.meta_merge

    def test_both_unmasked_left_right_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Left join
        t12 = table.join(t1, t2, join_type="left")
        assert t12.has_masked_columns is True
        assert t12.masked is False
        for name in ("a", "b", "c"):
            assert type(t12[name]) is Column
        assert type(t12["d"]) is MaskedColumn
        assert sort_eq(
            t12.pformat(),
            [
                " a   b   c   d ",
                "--- --- --- ---",
                "  0 foo  L1  --",
                "  1 bar  L3  --",
                "  1 foo  L2  R1",
                "  1 foo  L2  R2",
                "  2 bar  L4  R3",
            ],
        )

        # Right join
        t12 = table.join(t1, t2, join_type="right")
        assert t12.has_masked_columns is True
        assert t12.masked is False
        assert sort_eq(
            t12.pformat(),
            [
                " a   b   c   d ",
                "--- --- --- ---",
                "  1 foo  L2  R1",
                "  1 foo  L2  R2",
                "  2 bar  L4  R3",
                "  4 bar  --  R4",
            ],
        )

        # Outer join
        t12 = table.join(t1, t2, join_type="outer")
        assert t12.has_masked_columns is True
        assert t12.masked is False
        assert sort_eq(
            t12.pformat(),
            [
                " a   b   c   d ",
                "--- --- --- ---",
                "  0 foo  L1  --",
                "  1 bar  L3  --",
                "  1 foo  L2  R1",
                "  1 foo  L2  R2",
                "  2 bar  L4  R3",
                "  4 bar  --  R4",
            ],
        )

        # Check that the common keys are 'a', 'b'
        t12a = table.join(t1, t2, join_type="outer")
        t12b = table.join(t1, t2, join_type="outer", keys=["a", "b"])
        assert np.all(t12a.as_array() == t12b.as_array())

    def test_both_unmasked_single_key_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Inner join on 'a' column
        t12 = table.join(t1, t2, keys="a")
        assert type(t12) is operation_table_type
        assert type(t12["a"]) is type(t1["a"])
        assert type(t12["b_1"]) is type(t1["b"])
        assert type(t12["c"]) is type(t1["c"])
        assert type(t12["b_2"]) is type(t2["b"])
        assert type(t12["d"]) is type(t2["d"])
        assert t12.masked is False
        assert sort_eq(
            t12.pformat(),
            [
                " a  b_1  c  b_2  d ",
                "--- --- --- --- ---",
                "  1 foo  L2 foo  R1",
                "  1 foo  L2 foo  R2",
                "  1 bar  L3 foo  R1",
                "  1 bar  L3 foo  R2",
                "  2 bar  L4 bar  R3",
            ],
        )

    def test_both_unmasked_single_key_left_right_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Left join
        t12 = table.join(t1, t2, join_type="left", keys="a")
        assert t12.has_masked_columns is True
        assert sort_eq(
            t12.pformat(),
            [
                " a  b_1  c  b_2  d ",
                "--- --- --- --- ---",
                "  0 foo  L1  --  --",
                "  1 foo  L2 foo  R1",
                "  1 foo  L2 foo  R2",
                "  1 bar  L3 foo  R1",
                "  1 bar  L3 foo  R2",
                "  2 bar  L4 bar  R3",
            ],
        )

        # Right join
        t12 = table.join(t1, t2, join_type="right", keys="a")
        assert t12.has_masked_columns is True
        assert sort_eq(
            t12.pformat(),
            [
                " a  b_1  c  b_2  d ",
                "--- --- --- --- ---",
                "  1 foo  L2 foo  R1",
                "  1 foo  L2 foo  R2",
                "  1 bar  L3 foo  R1",
                "  1 bar  L3 foo  R2",
                "  2 bar  L4 bar  R3",
                "  4  --  -- bar  R4",
            ],
        )

        # Outer join
        t12 = table.join(t1, t2, join_type="outer", keys="a")
        assert t12.has_masked_columns is True
        assert sort_eq(
            t12.pformat(),
            [
                " a  b_1  c  b_2  d ",
                "--- --- --- --- ---",
                "  0 foo  L1  --  --",
                "  1 foo  L2 foo  R1",
                "  1 foo  L2 foo  R2",
                "  1 bar  L3 foo  R1",
                "  1 bar  L3 foo  R2",
                "  2 bar  L4 bar  R3",
                "  4  --  -- bar  R4",
            ],
        )

    def test_masked_unmasked(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        self._setup(operation_table_type)
        t1 = self.t1
        t1m = operation_table_type(self.t1, masked=True)
        t2 = self.t2

        # Result table is never masked
        t1m2 = table.join(t1m, t2, join_type="inner")
        assert t1m2.masked is False

        # Result should match non-masked result
        t12 = table.join(t1, t2)
        assert np.all(t12.as_array() == np.array(t1m2))

        # Mask out some values in left table and make sure they propagate
        t1m["b"].mask[1] = True
        t1m["c"].mask[2] = True
        t1m2 = table.join(t1m, t2, join_type="inner", keys="a")
        assert sort_eq(
            t1m2.pformat(),
            [
                " a  b_1  c  b_2  d ",
                "--- --- --- --- ---",
                "  1  --  L2 foo  R1",
                "  1  --  L2 foo  R2",
                "  1 bar  -- foo  R1",
                "  1 bar  -- foo  R2",
                "  2 bar  L4 bar  R3",
            ],
        )

        t21m = table.join(t2, t1m, join_type="inner", keys="a")
        assert sort_eq(
            t21m.pformat(),
            [
                " a  b_1  d  b_2  c ",
                "--- --- --- --- ---",
                "  1 foo  R2  --  L2",
                "  1 foo  R2 bar  --",
                "  1 foo  R1  --  L2",
                "  1 foo  R1 bar  --",
                "  2 bar  R3 bar  L4",
            ],
        )

    def test_masked_masked(self, operation_table_type):
        self._setup(operation_table_type)
        """Two masked tables"""
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        t1 = self.t1
        t1m = operation_table_type(self.t1, masked=True)
        t2 = self.t2
        t2m = operation_table_type(self.t2, masked=True)

        # Result table is never masked but original column types are preserved
        t1m2m = table.join(t1m, t2m, join_type="inner")
        assert t1m2m.masked is False
        for col in t1m2m.itercols():
            assert type(col) is MaskedColumn

        # Result should match non-masked result
        t12 = table.join(t1, t2)
        assert np.all(t12.as_array() == np.array(t1m2m))

        # Mask out some values in both tables and make sure they propagate
        t1m["b"].mask[1] = True
        t1m["c"].mask[2] = True
        t2m["d"].mask[2] = True
        t1m2m = table.join(t1m, t2m, join_type="inner", keys="a")
        assert sort_eq(
            t1m2m.pformat(),
            [
                " a  b_1  c  b_2  d ",
                "--- --- --- --- ---",
                "  1  --  L2 foo  R1",
                "  1  --  L2 foo  R2",
                "  1 bar  -- foo  R1",
                "  1 bar  -- foo  R2",
                "  2 bar  L4 bar  --",
            ],
        )

    def test_classes(self):
        """Ensure that classes and subclasses get through as expected"""

        class MyCol(Column):
            pass

        class MyMaskedCol(MaskedColumn):
            pass

        t1 = Table()
        t1["a"] = MyCol([1])
        t1["b"] = MyCol([2])
        t1["c"] = MyMaskedCol([3])

        t2 = Table()
        t2["a"] = Column([1, 2])
        t2["d"] = MyCol([3, 4])
        t2["e"] = MyMaskedCol([5, 6])

        t12 = table.join(t1, t2, join_type="inner")
        for name, exp_type in (
            ("a", MyCol),
            ("b", MyCol),
            ("c", MyMaskedCol),
            ("d", MyCol),
            ("e", MyMaskedCol),
        ):
            assert type(t12[name] is exp_type)

        t21 = table.join(t2, t1, join_type="left")
        # Note col 'b' gets upgraded from MyCol to MaskedColumn since it needs to be
        # masked, but col 'c' stays since MyMaskedCol supports masking.
        for name, exp_type in (
            ("a", MyCol),
            ("b", MaskedColumn),
            ("c", MyMaskedCol),
            ("d", MyCol),
            ("e", MyMaskedCol),
        ):
            assert type(t21[name] is exp_type)

    def test_col_rename(self, operation_table_type):
        self._setup(operation_table_type)
        """
        Test auto col renaming when there is a conflict.  Use
        non-default values of uniq_col_name and table_names.
        """
        t1 = self.t1
        t2 = self.t2
        t12 = table.join(
            t1,
            t2,
            uniq_col_name="x_{table_name}_{col_name}_y",
            table_names=["L", "R"],
            keys="a",
        )
        assert t12.colnames == ["a", "x_L_b_y", "c", "x_R_b_y", "d"]

    def test_rename_conflict(self, operation_table_type):
        self._setup(operation_table_type)
        """
        Test that auto-column rename fails because of a conflict
        with an existing column
        """
        t1 = self.t1
        t2 = self.t2
        t1["b_1"] = 1  # Add a new column b_1 that will conflict with auto-rename
        with pytest.raises(TableMergeError):
            table.join(t1, t2, keys="a")

    def test_missing_keys(self, operation_table_type):
        self._setup(operation_table_type)
        """Merge on a key column that doesn't exist"""
        t1 = self.t1
        t2 = self.t2
        with pytest.raises(TableMergeError):
            table.join(t1, t2, keys=["a", "not there"])

    def test_bad_join_type(self, operation_table_type):
        self._setup(operation_table_type)
        """Bad join_type input"""
        t1 = self.t1
        t2 = self.t2
        with pytest.raises(ValueError):
            table.join(t1, t2, join_type="illegal value")

    def test_no_common_keys(self, operation_table_type):
        self._setup(operation_table_type)
        """Merge tables with no common keys"""
        t1 = self.t1
        t2 = self.t2
        del t1["a"]
        del t1["b"]
        del t2["a"]
        del t2["b"]
        with pytest.raises(TableMergeError):
            table.join(t1, t2)

    def test_masked_key_column(self, operation_table_type):
        self._setup(operation_table_type)
        """Merge on a key column that has a masked element"""
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        t1 = self.t1
        t2 = operation_table_type(self.t2, masked=True)
        table.join(t1, t2)  # OK
        t2["a"].mask[0] = True
        with pytest.raises(TableMergeError):
            table.join(t1, t2)

    def test_col_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t2.rename_column("d", "c")  # force col conflict and renaming
        meta1 = OrderedDict([("b", [1, 2]), ("c", {"a": 1}), ("d", 1)])
        meta2 = OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)])

        # Key col 'a', should first value ('cm')
        t1["a"].unit = "cm"
        t2["a"].unit = "m"
        # Key col 'b', take first value 't1_b'
        t1["b"].info.description = "t1_b"
        # Key col 'b', take first non-empty value 't1_b'
        t2["b"].info.format = "%6s"
        # Key col 'a', should be merged meta
        t1["a"].info.meta = meta1
        t2["a"].info.meta = meta2
        # Key col 'b', should be meta2
        t2["b"].info.meta = meta2

        # All these should pass through
        t1["c"].info.format = "%3s"
        t1["c"].info.description = "t1_c"

        t2["c"].info.format = "%6s"
        t2["c"].info.description = "t2_c"

        if operation_table_type is Table:
            ctx = pytest.warns(
                metadata.MergeConflictWarning,
                match=(
                    r"In merged column 'a' the 'unit' attribute does not match \(cm"
                    r" != m\)"
                ),
            )
        else:
            ctx = nullcontext()

        with ctx:
            t12 = table.join(t1, t2, keys=["a", "b"])

        assert t12["a"].unit == "m"
        assert t12["b"].info.description == "t1_b"
        assert t12["b"].info.format == "%6s"
        assert t12["a"].info.meta == self.meta_merge
        assert t12["b"].info.meta == meta2
        assert t12["c_1"].info.format == "%3s"
        assert t12["c_1"].info.description == "t1_c"
        assert t12["c_2"].info.format == "%6s"
        assert t12["c_2"].info.description == "t2_c"

    def test_join_multidimensional(self, operation_table_type):
        self._setup(operation_table_type)

        # Regression test for #2984, which was an issue where join did not work
        # on multi-dimensional columns.

        t1 = operation_table_type()
        t1["a"] = [1, 2, 3]
        t1["b"] = np.ones((3, 4))

        t2 = operation_table_type()
        t2["a"] = [1, 2, 3]
        t2["c"] = [4, 5, 6]

        t3 = table.join(t1, t2)

        np.testing.assert_allclose(t3["a"], t1["a"])
        np.testing.assert_allclose(t3["b"], t1["b"])
        np.testing.assert_allclose(t3["c"], t2["c"])

    def test_join_multidimensional_masked(self, operation_table_type):
        self._setup(operation_table_type)
        """
        Test for outer join with multidimensional columns where masking is required.
        (Issue #4059).
        """
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")

        a = table.MaskedColumn([1, 2, 3], name="a")
        a2 = table.Column([1, 3, 4], name="a")
        b = table.MaskedColumn(
            [
                [1, 2],
                [3, 4],
                [5, 6],
            ],
            name="b",
            mask=[
                [1, 0],
                [0, 1],
                [0, 0],
            ],
        )
        c = table.Column(
            [
                [1, 1],
                [2, 2],
                [3, 3],
            ],
            name="c",
        )
        t1 = operation_table_type([a, b])
        t2 = operation_table_type([a2, c])
        t12 = table.join(t1, t2, join_type="inner")

        assert np.all(
            t12["b"].mask
            == [
                [True, False],
                [False, False],
            ]
        )
        assert not hasattr(t12["c"], "mask")

        t12 = table.join(t1, t2, join_type="outer")
        assert np.all(
            t12["b"].mask
            == [
                [True, False],
                [False, True],
                [False, False],
                [True, True],
            ]
        )
        assert np.all(
            t12["c"].mask
            == [
                [False, False],
                [True, True],
                [False, False],
                [False, False],
            ]
        )

    def test_mixin_functionality(self, mixin_cols):
        col = mixin_cols["m"]
        cls_name = type(col).__name__
        len_col = len(col)
        idx = np.arange(len_col)
        t1 = table.QTable([idx, col], names=["idx", "m1"])
        t2 = table.QTable([idx, col], names=["idx", "m2"])
        # Set up join mismatches for different join_type cases
        t1 = t1[[0, 1, 3]]
        t2 = t2[[0, 2, 3]]

        # Test inner join, which works for all mixin_cols
        out = table.join(t1, t2, join_type="inner")
        assert len(out) == 2
        assert out["m2"].__class__ is col.__class__
        assert np.all(out["idx"] == [0, 3])
        if cls_name == "SkyCoord":
            # SkyCoord doesn't support __eq__ so use our own
            assert skycoord_equal(out["m1"], col[[0, 3]])
            assert skycoord_equal(out["m2"], col[[0, 3]])
        elif "Repr" in cls_name or "Diff" in cls_name:
            assert np.all(representation_equal(out["m1"], col[[0, 3]]))
            assert np.all(representation_equal(out["m2"], col[[0, 3]]))
        else:
            assert np.all(out["m1"] == col[[0, 3]])
            assert np.all(out["m2"] == col[[0, 3]])

        # Check for left, right, outer join which requires masking. Works for
        # the listed mixins classes.
        if isinstance(col, (Quantity, Time, TimeDelta)):
            out = table.join(t1, t2, join_type="left")
            assert len(out) == 3
            assert np.all(out["idx"] == [0, 1, 3])
            assert np.all(out["m1"] == t1["m1"])
            assert np.all(out["m2"] == t2["m2"])
            check_mask(out["m1"], [False, False, False])
            check_mask(out["m2"], [False, True, False])

            out = table.join(t1, t2, join_type="right")
            assert len(out) == 3
            assert np.all(out["idx"] == [0, 2, 3])
            assert np.all(out["m1"] == t1["m1"])
            assert np.all(out["m2"] == t2["m2"])
            check_mask(out["m1"], [False, True, False])
            check_mask(out["m2"], [False, False, False])

            out = table.join(t1, t2, join_type="outer")
            assert len(out) == 4
            assert np.all(out["idx"] == [0, 1, 2, 3])
            assert np.all(out["m1"] == col)
            assert np.all(out["m2"] == col)
            assert check_mask(out["m1"], [False, False, True, False])
            assert check_mask(out["m2"], [False, True, False, False])
        else:
            # Otherwise make sure it fails with the right exception message
            for join_type in ("outer", "left", "right"):
                with pytest.raises(NotImplementedError) as err:
                    table.join(t1, t2, join_type=join_type)
                assert "join requires masking" in str(
                    err.value
                ) or "join unavailable" in str(err.value)

    def test_cartesian_join(self, operation_table_type):
        t1 = Table(rows=[(1, "a"), (2, "b")], names=["a", "b"])
        t2 = Table(rows=[(3, "c"), (4, "d")], names=["a", "c"])
        t12 = table.join(t1, t2, join_type="cartesian")

        assert t1.colnames == ["a", "b"]
        assert t2.colnames == ["a", "c"]
        assert len(t12) == len(t1) * len(t2)
        assert str(t12).splitlines() == [
            "a_1  b  a_2  c ",
            "--- --- --- ---",
            "  1   a   3   c",
            "  1   a   4   d",
            "  2   b   3   c",
            "  2   b   4   d",
        ]

        with pytest.raises(ValueError, match="cannot supply keys for a cartesian join"):
            t12 = table.join(t1, t2, join_type="cartesian", keys="a")

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    def test_join_with_join_skycoord_sky(self):
        sc1 = SkyCoord([0, 1, 1.1, 2], [0, 0, 0, 0], unit="deg")
        sc2 = SkyCoord([0.5, 1.05, 2.1], [0, 0, 0], unit="deg")
        t1 = Table([sc1], names=["sc"])
        t2 = Table([sc2], names=["sc"])
        t12 = table.join(t1, t2, join_funcs={"sc": join_skycoord(0.2 * u.deg)})
        exp = [
            "sc_id   sc_1    sc_2  ",
            "      deg,deg deg,deg ",
            "----- ------- --------",
            "    1 1.0,0.0 1.05,0.0",
            "    1 1.1,0.0 1.05,0.0",
            "    2 2.0,0.0  2.1,0.0",
        ]
        assert str(t12).splitlines() == exp

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    @pytest.mark.parametrize("distance_func", ["search_around_3d", search_around_3d])
    def test_join_with_join_skycoord_3d(self, distance_func):
        sc1 = SkyCoord([0, 1, 1.1, 2] * u.deg, [0, 0, 0, 0] * u.deg, [1, 1, 2, 1] * u.m)
        sc2 = SkyCoord([0.5, 1.05, 2.1] * u.deg, [0, 0, 0] * u.deg, [1, 1, 1] * u.m)
        t1 = Table([sc1], names=["sc"])
        t2 = Table([sc2], names=["sc"])
        join_func = join_skycoord(np.deg2rad(0.2) * u.m, distance_func=distance_func)
        t12 = table.join(t1, t2, join_funcs={"sc": join_func})
        exp = [
            "sc_id     sc_1        sc_2    ",
            "       deg,deg,m   deg,deg,m  ",
            "----- ----------- ------------",
            "    1 1.0,0.0,1.0 1.05,0.0,1.0",
            "    2 2.0,0.0,1.0  2.1,0.0,1.0",
        ]
        assert str(t12).splitlines() == exp

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    def test_join_with_join_distance_1d(self):
        c1 = [0, 1, 1.1, 2]
        c2 = [0.5, 1.05, 2.1]
        t1 = Table([c1], names=["col"])
        t2 = Table([c2], names=["col"])
        join_func = join_distance(
            0.2, kdtree_args={"leafsize": 32}, query_args={"p": 2}
        )
        t12 = table.join(t1, t2, join_type="outer", join_funcs={"col": join_func})
        exp = [
            "col_id col_1 col_2",
            "------ ----- -----",
            "     1   1.0  1.05",
            "     1   1.1  1.05",
            "     2   2.0   2.1",
            "     3   0.0    --",
            "     4    --   0.5",
        ]
        assert str(t12).splitlines() == exp

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    def test_join_with_join_distance_1d_multikey(self):
        from astropy.table.operations import _apply_join_funcs

        c1 = [0, 1, 1.1, 1.2, 2]
        id1 = [0, 1, 2, 2, 3]
        o1 = ["a", "b", "c", "d", "e"]
        c2 = [0.5, 1.05, 2.1]
        id2 = [0, 2, 4]
        o2 = ["z", "y", "x"]
        t1 = Table([c1, id1, o1], names=["col", "id", "o1"])
        t2 = Table([c2, id2, o2], names=["col", "id", "o2"])
        join_func = join_distance(0.2)
        join_funcs = {"col": join_func}
        t12 = table.join(t1, t2, join_type="outer", join_funcs=join_funcs)
        exp = [
            "col_id col_1  id  o1 col_2  o2",
            "------ ----- --- --- ----- ---",
            "     1   1.0   1   b    --  --",
            "     1   1.1   2   c  1.05   y",
            "     1   1.2   2   d  1.05   y",
            "     2   2.0   3   e    --  --",
            "     2    --   4  --   2.1   x",
            "     3   0.0   0   a    --  --",
            "     4    --   0  --   0.5   z",
        ]
        assert str(t12).splitlines() == exp

        left, right, keys = _apply_join_funcs(t1, t2, ("col", "id"), join_funcs)
        assert keys == ("col_id", "id")

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    def test_join_with_join_distance_1d_quantity(self):
        c1 = [0, 1, 1.1, 2] * u.m
        c2 = [500, 1050, 2100] * u.mm
        t1 = QTable([c1], names=["col"])
        t2 = QTable([c2], names=["col"])
        join_func = join_distance(20 * u.cm)
        t12 = table.join(t1, t2, join_funcs={"col": join_func})
        exp = [
            "col_id col_1 col_2 ",
            "         m     mm  ",
            "------ ----- ------",
            "     1   1.0 1050.0",
            "     1   1.1 1050.0",
            "     2   2.0 2100.0",
        ]
        assert str(t12).splitlines() == exp

        # Generate column name conflict
        t2["col_id"] = [0, 0, 0]
        t2["col__id"] = [0, 0, 0]
        t12 = table.join(t1, t2, join_funcs={"col": join_func})
        exp = [
            "col___id col_1 col_2  col_id col__id",
            "           m     mm                 ",
            "-------- ----- ------ ------ -------",
            "       1   1.0 1050.0      0       0",
            "       1   1.1 1050.0      0       0",
            "       2   2.0 2100.0      0       0",
        ]
        assert str(t12).splitlines() == exp

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    def test_join_with_join_distance_2d(self):
        c1 = np.array([[0, 1, 1.1, 2], [0, 0, 1, 0]]).transpose()
        c2 = np.array([[0.5, 1.05, 2.1], [0, 0, 0]]).transpose()
        t1 = Table([c1], names=["col"])
        t2 = Table([c2], names=["col"])
        join_func = join_distance(
            0.2, kdtree_args={"leafsize": 32}, query_args={"p": 2}
        )
        t12 = table.join(t1, t2, join_type="outer", join_funcs={"col": join_func})
        exp = [
            "col_id   col_1       col_2   ",
            f'{t12["col_id"].dtype.name}  float64[2]  float64[2]',  # int32 or int64
            "------ ---------- -----------",
            "     1 1.0 .. 0.0 1.05 .. 0.0",
            "     2 2.0 .. 0.0  2.1 .. 0.0",
            "     3 0.0 .. 0.0    -- .. --",
            "     4 1.1 .. 1.0    -- .. --",
            "     5   -- .. --  0.5 .. 0.0",
        ]
        assert t12.pformat(show_dtype=True) == exp

    def test_keys_left_right_basic(self):
        """Test using the keys_left and keys_right args to specify different
        join keys. This takes the standard test case but renames column 'a'
        to 'x' and 'y' respectively for tables 1 and 2. Then it compares the
        normal join on 'a' to the new join on 'x' and 'y'."""
        self._setup()

        for join_type in ("inner", "left", "right", "outer"):
            t1 = self.t1.copy()
            t2 = self.t2.copy()
            # Expected is same as joining on 'a' but with names 'x', 'y' instead
            t12_exp = table.join(t1, t2, keys="a", join_type=join_type)
            t12_exp.add_column(t12_exp["a"], name="x", index=1)
            t12_exp.add_column(t12_exp["a"], name="y", index=len(t1.colnames) + 1)
            del t12_exp["a"]

            # Different key names
            t1.rename_column("a", "x")
            t2.rename_column("a", "y")
            keys_left_list = ["x"]  # Test string key name
            keys_right_list = [["y"]]  # Test list of string key names
            if join_type == "outer":
                # Just do this for the outer join (others are the same)
                keys_left_list.append([t1["x"].tolist()])  # Test list key column
                keys_right_list.append([t2["y"]])  # Test Column key column

            for keys_left, keys_right in zip(keys_left_list, keys_right_list):
                t12 = table.join(
                    t1,
                    t2,
                    keys_left=keys_left,
                    keys_right=keys_right,
                    join_type=join_type,
                )

                assert t12.colnames == t12_exp.colnames
                for col in t12.values_equal(t12_exp).itercols():
                    assert np.all(col)
                assert t12_exp.meta == t12.meta

    def test_keys_left_right_exceptions(self):
        """Test exceptions using the keys_left and keys_right args to specify
        different join keys.
        """
        self._setup()
        t1 = self.t1
        t2 = self.t2

        msg = r"left table does not have key column 'z'"
        with pytest.raises(ValueError, match=msg):
            table.join(t1, t2, keys_left="z", keys_right=["a"])

        msg = r"left table has different length from key \[1, 2\]"
        with pytest.raises(ValueError, match=msg):
            table.join(t1, t2, keys_left=[[1, 2]], keys_right=["a"])

        msg = r"keys arg must be None if keys_left and keys_right are supplied"
        with pytest.raises(ValueError, match=msg):
            table.join(t1, t2, keys_left="z", keys_right=["a"], keys="a")

        msg = r"keys_left and keys_right args must have same length"
        with pytest.raises(ValueError, match=msg):
            table.join(t1, t2, keys_left=["a", "b"], keys_right=["a"])

        msg = r"keys_left and keys_right must both be provided"
        with pytest.raises(ValueError, match=msg):
            table.join(t1, t2, keys_left=["a", "b"])

        msg = r"cannot supply join_funcs arg and keys_left / keys_right"
        with pytest.raises(ValueError, match=msg):
            table.join(t1, t2, keys_left=["a"], keys_right=["a"], join_funcs={})

    def test_join_structured_column(self):
        """Regression tests for gh-13271."""
        # Two tables with matching names, including a structured column.
        t1 = Table(
            [
                np.array([(1.0, 1), (2.0, 2)], dtype=[("f", "f8"), ("i", "i8")]),
                ["one", "two"],
            ],
            names=["structured", "string"],
        )
        t2 = Table(
            [
                np.array([(2.0, 2), (4.0, 4)], dtype=[("f", "f8"), ("i", "i8")]),
                ["three", "four"],
            ],
            names=["structured", "string"],
        )
        t12 = table.join(t1, t2, ["structured"], join_type="outer")
        assert t12.pformat() == [
            "structured [f, i] string_1 string_2",
            "----------------- -------- --------",
            "          (1., 1)      one       --",
            "          (2., 2)      two    three",
            "          (4., 4)       --     four",
        ]


class TestSetdiff:
    def _setup(self, t_cls=Table):
        lines1 = [" a   b ", "  0 foo ", "  1 foo ", "  1 bar ", "  2 bar "]
        lines2 = [" a   b ", "  0 foo ", "  3 foo ", "  4 bar ", "  2 bar "]
        lines3 = [
            " a   b   d ",
            "  0 foo  R1",
            "  8 foo  R2",
            "  1 bar  R3",
            "  4 bar  R4",
        ]
        self.t1 = t_cls.read(lines1, format="ascii")
        self.t2 = t_cls.read(lines2, format="ascii")
        self.t3 = t_cls.read(lines3, format="ascii")

    def test_default_same_columns(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t1, self.t2)
        assert type(out["a"]) is type(self.t1["a"])
        assert type(out["b"]) is type(self.t1["b"])
        assert out.pformat() == [" a   b ", "--- ---", "  1 bar", "  1 foo"]

    def test_default_same_tables(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t1, self.t1)

        assert type(out["a"]) is type(self.t1["a"])
        assert type(out["b"]) is type(self.t1["b"])
        assert out.pformat() == [
            " a   b ",
            "--- ---",
        ]

    def test_extra_col_left_table(self, operation_table_type):
        self._setup(operation_table_type)

        with pytest.raises(ValueError):
            table.setdiff(self.t3, self.t1)

    def test_extra_col_right_table(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t1, self.t3)

        assert type(out["a"]) is type(self.t1["a"])
        assert type(out["b"]) is type(self.t1["b"])
        assert out.pformat() == [
            " a   b ",
            "--- ---",
            "  1 foo",
            "  2 bar",
        ]

    def test_keys(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t3, self.t1, keys=["a", "b"])

        assert type(out["a"]) is type(self.t1["a"])
        assert type(out["b"]) is type(self.t1["b"])
        assert out.pformat() == [
            " a   b   d ",
            "--- --- ---",
            "  4 bar  R4",
            "  8 foo  R2",
        ]

    def test_missing_key(self, operation_table_type):
        self._setup(operation_table_type)

        with pytest.raises(ValueError):
            table.setdiff(self.t3, self.t1, keys=["a", "d"])


class TestVStack:
    def _setup(self, t_cls=Table):
        self.t1 = t_cls.read(
            [
                " a   b",
                " 0. foo",
                " 1. bar",
            ],
            format="ascii",
        )

        self.t2 = t_cls.read(
            [
                " a    b   c",
                " 2.  pez  4",
                " 3.  sez  5",
            ],
            format="ascii",
        )

        self.t3 = t_cls.read(
            [
                " a    b",
                " 4.   7",
                " 5.   8",
                " 6.   9",
            ],
            format="ascii",
        )
        self.t4 = t_cls(self.t1, copy=True, masked=t_cls is Table)

        # The following table has meta-data that conflicts with t1
        self.t5 = t_cls(self.t1, copy=True)

        self.t1.meta.update(OrderedDict([("b", [1, 2]), ("c", {"a": 1}), ("d", 1)]))
        self.t2.meta.update(OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)]))
        self.t4.meta.update(OrderedDict([("b", [5, 6]), ("c", {"c": 1}), ("e", 1)]))
        self.t5.meta.update(OrderedDict([("b", 3), ("c", "k"), ("d", 1)]))
        self.meta_merge = OrderedDict(
            [
                ("b", [1, 2, 3, 4, 5, 6]),
                ("c", {"a": 1, "b": 1, "c": 1}),
                ("d", 1),
                ("a", 1),
                ("e", 1),
            ]
        )

    def test_validate_join_type(self):
        self._setup()
        with pytest.raises(TypeError, match="Did you accidentally call vstack"):
            table.vstack(self.t1, self.t2)

    def test_stack_rows(self, operation_table_type):
        self._setup(operation_table_type)
        t2 = self.t1.copy()
        t2.meta.clear()
        out = table.vstack([self.t1, t2[1]])
        assert type(out["a"]) is type(self.t1["a"])
        assert type(out["b"]) is type(self.t1["b"])
        assert out.pformat() == [
            " a   b ",
            "--- ---",
            "0.0 foo",
            "1.0 bar",
            "1.0 bar",
        ]

    def test_stack_table_column(self, operation_table_type):
        self._setup(operation_table_type)
        t2 = self.t1.copy()
        t2.meta.clear()
        out = table.vstack([self.t1, t2["a"]])
        assert out.masked is False
        assert out.pformat() == [
            " a   b ",
            "--- ---",
            "0.0 foo",
            "1.0 bar",
            "0.0  --",
            "1.0  --",
        ]

    def test_table_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.vstack([self.t1, self.t2, self.t4], join_type="inner")
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self, operation_table_type):
        self._setup(operation_table_type)

        with pytest.warns(metadata.MergeConflictWarning) as w:
            out = table.vstack([self.t1, self.t5], join_type="inner")
        assert len(w) == 2

        assert out.meta == self.t5.meta

        with pytest.warns(metadata.MergeConflictWarning) as w:
            out = table.vstack(
                [self.t1, self.t5], join_type="inner", metadata_conflicts="warn"
            )
        assert len(w) == 2

        assert out.meta == self.t5.meta

        out = table.vstack(
            [self.t1, self.t5], join_type="inner", metadata_conflicts="silent"
        )

        assert out.meta == self.t5.meta

        with pytest.raises(MergeConflictError):
            out = table.vstack(
                [self.t1, self.t5], join_type="inner", metadata_conflicts="error"
            )

        with pytest.raises(ValueError):
            out = table.vstack(
                [self.t1, self.t5], join_type="inner", metadata_conflicts="nonsense"
            )

    def test_bad_input_type(self, operation_table_type):
        self._setup(operation_table_type)
        with pytest.raises(ValueError):
            table.vstack([])
        with pytest.raises(TypeError):
            table.vstack(1)
        with pytest.raises(TypeError):
            table.vstack([self.t2, 1])
        with pytest.raises(ValueError):
            table.vstack([self.t1, self.t2], join_type="invalid join type")

    def test_stack_basic_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        t12 = table.vstack([t1, t2], join_type="inner")
        assert t12.masked is False
        assert type(t12) is operation_table_type
        assert type(t12["a"]) is type(t1["a"])
        assert type(t12["b"]) is type(t1["b"])
        assert t12.pformat() == [
            " a   b ",
            "--- ---",
            "0.0 foo",
            "1.0 bar",
            "2.0 pez",
            "3.0 sez",
        ]

        t124 = table.vstack([t1, t2, t4], join_type="inner")
        assert type(t124) is operation_table_type
        assert type(t12["a"]) is type(t1["a"])
        assert type(t12["b"]) is type(t1["b"])
        assert t124.pformat() == [
            " a   b ",
            "--- ---",
            "0.0 foo",
            "1.0 bar",
            "2.0 pez",
            "3.0 sez",
            "0.0 foo",
            "1.0 bar",
        ]

    def test_stack_basic_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4
        t12 = table.vstack([t1, t2], join_type="outer")
        assert t12.masked is False
        assert t12.pformat() == [
            " a   b   c ",
            "--- --- ---",
            "0.0 foo  --",
            "1.0 bar  --",
            "2.0 pez   4",
            "3.0 sez   5",
        ]

        t124 = table.vstack([t1, t2, t4], join_type="outer")
        assert t124.masked is False
        assert t124.pformat() == [
            " a   b   c ",
            "--- --- ---",
            "0.0 foo  --",
            "1.0 bar  --",
            "2.0 pez   4",
            "3.0 sez   5",
            "0.0 foo  --",
            "1.0 bar  --",
        ]

    def test_stack_incompatible(self, operation_table_type):
        self._setup(operation_table_type)
        with pytest.raises(TableMergeError) as excinfo:
            table.vstack([self.t1, self.t3], join_type="inner")
        assert "The 'b' columns have incompatible types: {}".format(
            [self.t1["b"].dtype.name, self.t3["b"].dtype.name]
        ) in str(excinfo.value)

        with pytest.raises(TableMergeError) as excinfo:
            table.vstack([self.t1, self.t3], join_type="outer")
        assert "The 'b' columns have incompatible types:" in str(excinfo.value)

        with pytest.raises(TableMergeError):
            table.vstack([self.t1, self.t2], join_type="exact")

        t1_reshape = self.t1.copy()
        t1_reshape["b"].shape = [2, 1]
        with pytest.raises(TableMergeError) as excinfo:
            table.vstack([self.t1, t1_reshape])
        assert "have different shape" in str(excinfo.value)

    def test_vstack_one_masked(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        self._setup(operation_table_type)
        t1 = self.t1
        t4 = self.t4
        t4["b"].mask[1] = True
        t14 = table.vstack([t1, t4])
        assert t14.masked is False
        assert t14.pformat() == [
            " a   b ",
            "--- ---",
            "0.0 foo",
            "1.0 bar",
            "0.0 foo",
            "1.0  --",
        ]

    def test_col_meta_merge_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        # Key col 'a', should last value ('km')
        t1["a"].info.unit = "cm"
        t2["a"].info.unit = "m"
        t4["a"].info.unit = "km"

        # Key col 'a' format should take last when all match
        t1["a"].info.format = "%f"
        t2["a"].info.format = "%f"
        t4["a"].info.format = "%f"

        # Key col 'b', take first value 't1_b'
        t1["b"].info.description = "t1_b"

        # Key col 'b', take first non-empty value '%6s'
        t4["b"].info.format = "%6s"

        # Key col 'a', should be merged meta
        t1["a"].info.meta.update(
            OrderedDict([("b", [1, 2]), ("c", {"a": 1}), ("d", 1)])
        )
        t2["a"].info.meta.update(
            OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)])
        )
        t4["a"].info.meta.update(
            OrderedDict([("b", [5, 6]), ("c", {"c": 1}), ("e", 1)])
        )

        # Key col 'b', should be meta2
        t2["b"].info.meta.update(
            OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)])
        )

        if operation_table_type is Table:
            ctx = pytest.warns(metadata.MergeConflictWarning)
        else:
            ctx = nullcontext()

        with ctx as warning_lines:
            out = table.vstack([t1, t2, t4], join_type="inner")

        if operation_table_type is Table:
            assert len(warning_lines) == 2
            assert (
                "In merged column 'a' the 'unit' attribute does not match (cm != m)"
                in str(warning_lines[0].message)
            )
            assert (
                "In merged column 'a' the 'unit' attribute does not match (m != km)"
                in str(warning_lines[1].message)
            )
            # Check units are suitably ignored for a regular Table
            assert out.pformat() == [
                "   a       b   ",
                "   km          ",
                "-------- ------",
                "0.000000    foo",
                "1.000000    bar",
                "2.000000    pez",
                "3.000000    sez",
                "0.000000    foo",
                "1.000000    bar",
            ]
        else:
            # Check QTable correctly dealt with units.
            assert out.pformat() == [
                "   a       b   ",
                "   km          ",
                "-------- ------",
                "0.000000    foo",
                "0.000010    bar",
                "0.002000    pez",
                "0.003000    sez",
                "0.000000    foo",
                "1.000000    bar",
            ]
        assert out["a"].info.unit == "km"
        assert out["a"].info.format == "%f"
        assert out["b"].info.description == "t1_b"
        assert out["b"].info.format == "%6s"
        assert out["a"].info.meta == self.meta_merge
        assert out["b"].info.meta == OrderedDict(
            [("b", [3, 4]), ("c", {"b": 1}), ("a", 1)]
        )

    def test_col_meta_merge_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        # Key col 'a', should last value ('km')
        t1["a"].unit = "cm"
        t2["a"].unit = "m"
        t4["a"].unit = "km"

        # Key col 'a' format should take last when all match
        t1["a"].info.format = "%0d"
        t2["a"].info.format = "%0d"
        t4["a"].info.format = "%0d"

        # Key col 'b', take first value 't1_b'
        t1["b"].info.description = "t1_b"

        # Key col 'b', take first non-empty value '%6s'
        t4["b"].info.format = "%6s"

        # Key col 'a', should be merged meta
        t1["a"].info.meta.update(
            OrderedDict([("b", [1, 2]), ("c", {"a": 1}), ("d", 1)])
        )
        t2["a"].info.meta.update(
            OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)])
        )
        t4["a"].info.meta.update(
            OrderedDict([("b", [5, 6]), ("c", {"c": 1}), ("e", 1)])
        )

        # Key col 'b', should be meta2
        t2["b"].info.meta.update(
            OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)])
        )

        # All these should pass through
        t2["c"].unit = "m"
        t2["c"].info.format = "%6s"
        t2["c"].info.description = "t2_c"

        with pytest.warns(metadata.MergeConflictWarning) as warning_lines:
            out = table.vstack([t1, t2, t4], join_type="outer")

        assert len(warning_lines) == 2
        assert (
            "In merged column 'a' the 'unit' attribute does not match (cm != m)"
            in str(warning_lines[0].message)
        )
        assert (
            "In merged column 'a' the 'unit' attribute does not match (m != km)"
            in str(warning_lines[1].message)
        )
        assert out["a"].unit == "km"
        assert out["a"].info.format == "%0d"
        assert out["b"].info.description == "t1_b"
        assert out["b"].info.format == "%6s"
        assert out["a"].info.meta == self.meta_merge
        assert out["b"].info.meta == OrderedDict(
            [("b", [3, 4]), ("c", {"b": 1}), ("a", 1)]
        )
        assert out["c"].info.unit == "m"
        assert out["c"].info.format == "%6s"
        assert out["c"].info.description == "t2_c"

    def test_vstack_one_table(self, operation_table_type):
        self._setup(operation_table_type)
        """Regression test for issue #3313"""
        assert (self.t1 == table.vstack(self.t1)).all()
        assert (self.t1 == table.vstack([self.t1])).all()

    def test_mixin_functionality(self, mixin_cols):
        col = mixin_cols["m"]
        len_col = len(col)
        t = table.QTable([col], names=["a"])
        cls_name = type(col).__name__

        # Vstack works for these classes:
        if isinstance(
            col,
            (
                u.Quantity,
                Time,
                TimeDelta,
                SkyCoord,
                EarthLocation,
                BaseRepresentationOrDifferential,
                StokesCoord,
            ),
        ):
            out = table.vstack([t, t])
            assert len(out) == len_col * 2
            if cls_name == "SkyCoord":
                # Argh, SkyCoord needs __eq__!!
                assert skycoord_equal(out["a"][len_col:], col)
                assert skycoord_equal(out["a"][:len_col], col)
            elif "Repr" in cls_name or "Diff" in cls_name:
                assert np.all(representation_equal(out["a"][:len_col], col))
                assert np.all(representation_equal(out["a"][len_col:], col))
            else:
                assert np.all(out["a"][:len_col] == col)
                assert np.all(out["a"][len_col:] == col)
        else:
            with pytest.raises(NotImplementedError) as err:
                table.vstack([t, t])
            assert "vstack unavailable for mixin column type(s): {}".format(
                cls_name
            ) in str(err.value)

        # Check for outer stack which requires masking.  Only Time supports
        # this currently.
        t2 = table.QTable([col], names=["b"])  # different from col name for t
        if isinstance(col, (Time, TimeDelta, Quantity)):
            out = table.vstack([t, t2], join_type="outer")
            assert len(out) == len_col * 2
            assert np.all(out["a"][:len_col] == col)
            assert np.all(out["b"][len_col:] == col)
            assert check_mask(out["a"], [False] * len_col + [True] * len_col)
            assert check_mask(out["b"], [True] * len_col + [False] * len_col)
            # check directly stacking mixin columns:
            out2 = table.vstack([t, t2["b"]])
            assert np.all(out["a"] == out2["a"])
            assert np.all(out["b"] == out2["b"])
        else:
            with pytest.raises(NotImplementedError) as err:
                table.vstack([t, t2], join_type="outer")
            assert "vstack requires masking" in str(
                err.value
            ) or "vstack unavailable" in str(err.value)

    def test_vstack_different_representation(self):
        """Test that representations can be mixed together."""
        rep1 = CartesianRepresentation([1, 2] * u.km, [3, 4] * u.km, 1 * u.km)
        rep2 = SphericalRepresentation([0] * u.deg, [0] * u.deg, 10 * u.km)
        t1 = Table([rep1])
        t2 = Table([rep2])
        t12 = table.vstack([t1, t2])
        expected = CartesianRepresentation(
            [1, 2, 10] * u.km, [3, 4, 0] * u.km, [1, 1, 0] * u.km
        )
        assert np.all(representation_equal(t12["col0"], expected))

        rep3 = UnitSphericalRepresentation([0] * u.deg, [0] * u.deg)
        t3 = Table([rep3])
        with pytest.raises(ValueError, match="representations are inconsistent"):
            table.vstack([t1, t3])

    def test_vstack_structured_column(self):
        """Regression tests for gh-13271."""
        # Two tables with matching names, including a structured column.
        t1 = Table(
            [
                np.array([(1.0, 1), (2.0, 2)], dtype=[("f", "f8"), ("i", "i8")]),
                ["one", "two"],
            ],
            names=["structured", "string"],
        )
        t2 = Table(
            [
                np.array([(3.0, 3), (4.0, 4)], dtype=[("f", "f8"), ("i", "i8")]),
                ["three", "four"],
            ],
            names=["structured", "string"],
        )
        t12 = table.vstack([t1, t2])
        assert t12.pformat() == [
            "structured [f, i] string",
            "----------------- ------",
            "          (1., 1)    one",
            "          (2., 2)    two",
            "          (3., 3)  three",
            "          (4., 4)   four",
        ]

        # One table without the structured column.
        t3 = t2[("string",)]
        t13 = table.vstack([t1, t3])
        assert t13.pformat() == [
            "structured [f, i] string",
            "----------------- ------",
            "         (1.0, 1)    one",
            "         (2.0, 2)    two",
            "               --  three",
            "               --   four",
        ]


class TestDStack:
    def _setup(self, t_cls=Table):
        self.t1 = t_cls.read(
            [
                " a   b",
                " 0. foo",
                " 1. bar",
            ],
            format="ascii",
        )

        self.t2 = t_cls.read(
            [
                " a    b   c",
                " 2.  pez  4",
                " 3.  sez  5",
            ],
            format="ascii",
        )
        self.t2["d"] = Time([1, 2], format="cxcsec")

        self.t3 = t_cls(
            {
                "a": [[5.0, 6.0], [4.0, 3.0]],
                "b": [["foo", "bar"], ["pez", "sez"]],
            },
            names=("a", "b"),
        )

        self.t4 = t_cls(self.t1, copy=True, masked=t_cls is Table)

        self.t5 = t_cls(
            {
                "a": [[4.0, 2.0], [1.0, 6.0]],
                "b": [["foo", "pez"], ["bar", "sez"]],
            },
            names=("a", "b"),
        )
        self.t6 = t_cls.read(
            [
                " a    b   c",
                " 7.  pez  2",
                " 4.  sez  6",
                " 6.  foo  3",
            ],
            format="ascii",
        )

    def test_validate_join_type(self):
        self._setup()
        with pytest.raises(TypeError, match="Did you accidentally call dstack"):
            table.dstack(self.t1, self.t2)

    @staticmethod
    def compare_dstack(tables, out):
        for ii, tbl in enumerate(tables):
            for name, out_col in out.columns.items():
                if name in tbl.colnames:
                    # Columns always compare equal
                    assert np.all(tbl[name] == out[name][:, ii])

                    # If input has a mask then output must have same mask
                    if hasattr(tbl[name], "mask"):
                        assert np.all(tbl[name].mask == out[name].mask[:, ii])

                    # If input has no mask then output might have a mask (if other table
                    # is missing that column). If so then all mask values should be False.
                    elif hasattr(out[name], "mask"):
                        assert not np.any(out[name].mask[:, ii])

                else:
                    # Column missing for this table, out must have a mask with all True.
                    assert np.all(out[name].mask[:, ii])

    def test_dstack_table_column(self, operation_table_type):
        """Stack a table with 3 cols and one column (gets auto-converted to Table)."""
        self._setup(operation_table_type)
        t2 = self.t1.copy()
        out = table.dstack([self.t1, t2["a"]])
        self.compare_dstack([self.t1, t2[("a",)]], out)

    def test_dstack_basic_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail("Quantity columns do not support masking.")
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4
        t4["a"].mask[0] = True
        # Test for non-masked table
        t12 = table.dstack([t1, t2], join_type="outer")
        assert type(t12) is operation_table_type
        assert type(t12["a"]) is type(t1["a"])
        assert type(t12["b"]) is type(t1["b"])
        self.compare_dstack([t1, t2], t12)

        # Test for masked table
        t124 = table.dstack([t1, t2, t4], join_type="outer")
        assert type(t124) is operation_table_type
        assert type(t124["a"]) is type(t4["a"])
        assert type(t124["b"]) is type(t4["b"])
        self.compare_dstack([t1, t2, t4], t124)

    def test_dstack_basic_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        # Test for masked table
        t124 = table.dstack([t1, t2, t4], join_type="inner")
        assert type(t124) is operation_table_type
        assert type(t124["a"]) is type(t4["a"])
        assert type(t124["b"]) is type(t4["b"])
        self.compare_dstack([t1, t2, t4], t124)

    def test_dstack_multi_dimension_column(self, operation_table_type):
        self._setup(operation_table_type)
        t3 = self.t3
        t5 = self.t5
        t2 = self.t2
        t35 = table.dstack([t3, t5])
        assert type(t35) is operation_table_type
        assert type(t35["a"]) is type(t3["a"])
        assert type(t35["b"]) is type(t3["b"])
        self.compare_dstack([t3, t5], t35)

        with pytest.raises(TableMergeError):
            table.dstack([t2, t3])

    def test_dstack_different_length_table(self, operation_table_type):
        self._setup(operation_table_type)
        t2 = self.t2
        t6 = self.t6
        with pytest.raises(ValueError):
            table.dstack([t2, t6])

    def test_dstack_single_table(self):
        self._setup(Table)
        out = table.dstack(self.t1)
        assert np.all(out == self.t1)

    def test_dstack_representation(self):
        rep1 = SphericalRepresentation([1, 2] * u.deg, [3, 4] * u.deg, 1 * u.kpc)
        rep2 = SphericalRepresentation([10, 20] * u.deg, [30, 40] * u.deg, 10 * u.kpc)
        t1 = Table([rep1])
        t2 = Table([rep2])
        t12 = table.dstack([t1, t2])
        assert np.all(representation_equal(t12["col0"][:, 0], rep1))
        assert np.all(representation_equal(t12["col0"][:, 1], rep2))

    def test_dstack_skycoord(self):
        sc1 = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        sc2 = SkyCoord([10, 20] * u.deg, [30, 40] * u.deg)
        t1 = Table([sc1])
        t2 = Table([sc2])
        t12 = table.dstack([t1, t2])
        assert skycoord_equal(sc1, t12["col0"][:, 0])
        assert skycoord_equal(sc2, t12["col0"][:, 1])

    def test_dstack_structured_column(self):
        """Regression tests for gh-13271."""
        # Two tables with matching names, including a structured column.
        t1 = Table(
            [
                np.array([(1.0, 1), (2.0, 2)], dtype=[("f", "f8"), ("i", "i8")]),
                ["one", "two"],
            ],
            names=["structured", "string"],
        )
        t2 = Table(
            [
                np.array([(3.0, 3), (4.0, 4)], dtype=[("f", "f8"), ("i", "i8")]),
                ["three", "four"],
            ],
            names=["structured", "string"],
        )
        t12 = table.dstack([t1, t2])
        assert t12.pformat() == [
            "structured [f, i]     string   ",
            "------------------ ------------",
            "(1., 1) .. (3., 3) one .. three",
            "(2., 2) .. (4., 4)  two .. four",
        ]

        # One table without the structured column.
        t3 = t2[("string",)]
        t13 = table.dstack([t1, t3])
        assert t13.pformat() == [
            "structured [f, i]    string   ",
            "----------------- ------------",
            "   (1.0, 1) .. -- one .. three",
            "   (2.0, 2) .. --  two .. four",
        ]


class TestHStack:
    def _setup(self, t_cls=Table):
        self.t1 = t_cls.read(
            [
                " a    b",
                " 0. foo",
                " 1. bar",
            ],
            format="ascii",
        )

        self.t2 = t_cls.read(
            [
                " a    b   c",
                " 2.  pez  4",
                " 3.  sez  5",
            ],
            format="ascii",
        )

        self.t3 = t_cls.read(
            [
                " d    e",
                " 4.   7",
                " 5.   8",
                " 6.   9",
            ],
            format="ascii",
        )
        self.t4 = t_cls(self.t1, copy=True, masked=True)
        self.t4["a"].name = "f"
        self.t4["b"].name = "g"

        # The following table has meta-data that conflicts with t1
        self.t5 = t_cls(self.t1, copy=True)

        self.t1.meta.update(OrderedDict([("b", [1, 2]), ("c", {"a": 1}), ("d", 1)]))
        self.t2.meta.update(OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)]))
        self.t4.meta.update(OrderedDict([("b", [5, 6]), ("c", {"c": 1}), ("e", 1)]))
        self.t5.meta.update(OrderedDict([("b", 3), ("c", "k"), ("d", 1)]))
        self.meta_merge = OrderedDict(
            [
                ("b", [1, 2, 3, 4, 5, 6]),
                ("c", {"a": 1, "b": 1, "c": 1}),
                ("d", 1),
                ("a", 1),
                ("e", 1),
            ]
        )

    def test_validate_join_type(self):
        self._setup()
        with pytest.raises(TypeError, match="Did you accidentally call hstack"):
            table.hstack(self.t1, self.t2)

    def test_stack_same_table(self, operation_table_type):
        """
        From #2995, test that hstack'ing references to the same table has the
        expected output.
        """
        self._setup(operation_table_type)
        out = table.hstack([self.t1, self.t1])
        assert out.masked is False
        assert out.pformat() == [
            "a_1 b_1 a_2 b_2",
            "--- --- --- ---",
            "0.0 foo 0.0 foo",
            "1.0 bar 1.0 bar",
        ]

    def test_stack_rows(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack([self.t1[0], self.t2[1]])
        assert out.masked is False
        assert out.pformat() == [
            "a_1 b_1 a_2 b_2  c ",
            "--- --- --- --- ---",
            "0.0 foo 3.0 sez   5",
        ]

    def test_stack_columns(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack([self.t1, self.t2["c"]])
        assert type(out["a"]) is type(self.t1["a"])
        assert type(out["b"]) is type(self.t1["b"])
        assert type(out["c"]) is type(self.t2["c"])
        assert out.pformat() == [
            " a   b   c ",
            "--- --- ---",
            "0.0 foo   4",
            "1.0 bar   5",
        ]

    def test_table_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack([self.t1, self.t2, self.t4], join_type="inner")
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self, operation_table_type):
        self._setup(operation_table_type)

        with pytest.warns(metadata.MergeConflictWarning) as w:
            out = table.hstack([self.t1, self.t5], join_type="inner")
        assert len(w) == 2

        assert out.meta == self.t5.meta

        with pytest.warns(metadata.MergeConflictWarning) as w:
            out = table.hstack(
                [self.t1, self.t5], join_type="inner", metadata_conflicts="warn"
            )
        assert len(w) == 2

        assert out.meta == self.t5.meta

        out = table.hstack(
            [self.t1, self.t5], join_type="inner", metadata_conflicts="silent"
        )

        assert out.meta == self.t5.meta

        with pytest.raises(MergeConflictError):
            out = table.hstack(
                [self.t1, self.t5], join_type="inner", metadata_conflicts="error"
            )

        with pytest.raises(ValueError):
            out = table.hstack(
                [self.t1, self.t5], join_type="inner", metadata_conflicts="nonsense"
            )

    def test_bad_input_type(self, operation_table_type):
        self._setup(operation_table_type)
        with pytest.raises(ValueError):
            table.hstack([])
        with pytest.raises(TypeError):
            table.hstack(1)
        with pytest.raises(TypeError):
            table.hstack([self.t2, 1])
        with pytest.raises(ValueError):
            table.hstack([self.t1, self.t2], join_type="invalid join type")

    def test_stack_basic(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t3 = self.t3
        t4 = self.t4

        out = table.hstack([t1, t2], join_type="inner")
        assert out.masked is False
        assert type(out) is operation_table_type
        assert type(out["a_1"]) is type(t1["a"])
        assert type(out["b_1"]) is type(t1["b"])
        assert type(out["a_2"]) is type(t2["a"])
        assert type(out["b_2"]) is type(t2["b"])
        assert out.pformat() == [
            "a_1 b_1 a_2 b_2  c ",
            "--- --- --- --- ---",
            "0.0 foo 2.0 pez   4",
            "1.0 bar 3.0 sez   5",
        ]

        # stacking as a list gives same result
        out_list = table.hstack([t1, t2], join_type="inner")
        assert out.pformat() == out_list.pformat()

        out = table.hstack([t1, t2], join_type="outer")
        assert out.pformat() == out_list.pformat()

        out = table.hstack([t1, t2, t3, t4], join_type="outer")
        assert out.masked is False
        assert out.pformat() == [
            "a_1 b_1 a_2 b_2  c   d   e   f   g ",
            "--- --- --- --- --- --- --- --- ---",
            "0.0 foo 2.0 pez   4 4.0   7 0.0 foo",
            "1.0 bar 3.0 sez   5 5.0   8 1.0 bar",
            " --  --  --  --  -- 6.0   9  --  --",
        ]

        out = table.hstack([t1, t2, t3, t4], join_type="inner")
        assert out.masked is False
        assert out.pformat() == [
            "a_1 b_1 a_2 b_2  c   d   e   f   g ",
            "--- --- --- --- --- --- --- --- ---",
            "0.0 foo 2.0 pez   4 4.0   7 0.0 foo",
            "1.0 bar 3.0 sez   5 5.0   8 1.0 bar",
        ]

    def test_stack_incompatible(self, operation_table_type):
        self._setup(operation_table_type)
        # For join_type exact, which will fail here because n_rows
        # does not match
        with pytest.raises(TableMergeError):
            table.hstack([self.t1, self.t3], join_type="exact")

    def test_hstack_one_masked(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail()
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = operation_table_type(t1, copy=True, masked=True)
        t2.meta.clear()
        t2["b"].mask[1] = True
        out = table.hstack([t1, t2])
        assert out.pformat() == [
            "a_1 b_1 a_2 b_2",
            "--- --- --- ---",
            "0.0 foo 0.0 foo",
            "1.0 bar 1.0  --",
        ]

    def test_table_col_rename(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack(
            [self.t1, self.t2],
            join_type="inner",
            uniq_col_name="{table_name}_{col_name}",
            table_names=("left", "right"),
        )
        assert out.masked is False
        assert out.pformat() == [
            "left_a left_b right_a right_b  c ",
            "------ ------ ------- ------- ---",
            "   0.0    foo     2.0     pez   4",
            "   1.0    bar     3.0     sez   5",
        ]

    def test_col_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t3 = self.t3[:2]
        t4 = self.t4

        # Just set a bunch of meta and make sure it is the same in output
        meta1 = OrderedDict([("b", [1, 2]), ("c", {"a": 1}), ("d", 1)])
        t1["a"].unit = "cm"
        t1["b"].info.description = "t1_b"
        t4["f"].info.format = "%6s"
        t1["b"].info.meta.update(meta1)
        t3["d"].info.meta.update(
            OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)])
        )
        t4["g"].info.meta.update(
            OrderedDict([("b", [5, 6]), ("c", {"c": 1}), ("e", 1)])
        )
        t3["e"].info.meta.update(
            OrderedDict([("b", [3, 4]), ("c", {"b": 1}), ("a", 1)])
        )
        t3["d"].unit = "m"
        t3["d"].info.format = "%6s"
        t3["d"].info.description = "t3_c"

        out = table.hstack([t1, t3, t4], join_type="exact")

        for t in [t1, t3, t4]:
            for name in t.colnames:
                for attr in ("meta", "unit", "format", "description"):
                    assert getattr(out[name].info, attr) == getattr(t[name].info, attr)

        # Make sure we got a copy of meta, not ref
        t1["b"].info.meta["b"] = None
        assert out["b"].info.meta["b"] == [1, 2]

    def test_hstack_one_table(self, operation_table_type):
        self._setup(operation_table_type)
        """Regression test for issue #3313"""
        assert (self.t1 == table.hstack(self.t1)).all()
        assert (self.t1 == table.hstack([self.t1])).all()

    def test_mixin_functionality(self, mixin_cols):
        col1 = mixin_cols["m"]
        col2 = col1[2:4]  # Shorter version of col1
        t1 = table.QTable([col1])
        t2 = table.QTable([col2])

        cls_name = type(col1).__name__

        out = table.hstack([t1, t2], join_type="inner")
        assert type(out["col0_1"]) is type(out["col0_2"])
        assert len(out) == len(col2)

        # Check that columns are as expected.
        if cls_name == "SkyCoord":
            assert skycoord_equal(out["col0_1"], col1[: len(col2)])
            assert skycoord_equal(out["col0_2"], col2)
        elif "Repr" in cls_name or "Diff" in cls_name:
            assert np.all(representation_equal(out["col0_1"], col1[: len(col2)]))
            assert np.all(representation_equal(out["col0_2"], col2))
        else:
            assert np.all(out["col0_1"] == col1[: len(col2)])
            assert np.all(out["col0_2"] == col2)

        # Time class supports masking, all other mixins do not
        if isinstance(col1, (Time, TimeDelta, Quantity)):
            out = table.hstack([t1, t2], join_type="outer")
            assert len(out) == len(t1)
            assert np.all(out["col0_1"] == col1)
            assert np.all(out["col0_2"][: len(col2)] == col2)
            assert check_mask(out["col0_2"], [False, False, True, True])

            # check directly stacking mixin columns:
            out2 = table.hstack([t1, t2["col0"]], join_type="outer")
            assert np.all(out["col0_1"] == out2["col0_1"])
            assert np.all(out["col0_2"] == out2["col0_2"])
        else:
            with pytest.raises(NotImplementedError) as err:
                table.hstack([t1, t2], join_type="outer")
            assert "hstack requires masking" in str(err.value)


def test_unique(operation_table_type):
    t = operation_table_type.read(
        [
            " a b  c  d",
            " 2 b 7.0 0",
            " 1 c 3.0 5",
            " 2 b 6.0 2",
            " 2 a 4.0 3",
            " 1 a 1.0 7",
            " 2 b 5.0 1",
            " 0 a 0.0 4",
            " 1 a 2.0 6",
            " 1 c 3.0 5",
        ],
        format="ascii",
    )

    tu = operation_table_type(np.sort(t[:-1]))

    t_all = table.unique(t)
    assert sort_eq(t_all.pformat(), tu.pformat())
    t_s = t.copy()
    del t_s["b", "c", "d"]
    t_all = table.unique(t_s)
    assert sort_eq(
        t_all.pformat(),
        [
            " a ",
            "---",
            "  0",
            "  1",
            "  2",
        ],
    )

    key1 = "a"
    t1a = table.unique(t, key1)
    assert sort_eq(
        t1a.pformat(),
        [
            " a   b   c   d ",
            "--- --- --- ---",
            "  0   a 0.0   4",
            "  1   c 3.0   5",
            "  2   b 7.0   0",
        ],
    )
    t1b = table.unique(t, key1, keep="last")
    assert sort_eq(
        t1b.pformat(),
        [
            " a   b   c   d ",
            "--- --- --- ---",
            "  0   a 0.0   4",
            "  1   c 3.0   5",
            "  2   b 5.0   1",
        ],
    )
    t1c = table.unique(t, key1, keep="none")
    assert sort_eq(
        t1c.pformat(),
        [
            " a   b   c   d ",
            "--- --- --- ---",
            "  0   a 0.0   4",
        ],
    )

    key2 = ["a", "b"]
    t2a = table.unique(t, key2)
    assert sort_eq(
        t2a.pformat(),
        [
            " a   b   c   d ",
            "--- --- --- ---",
            "  0   a 0.0   4",
            "  1   a 1.0   7",
            "  1   c 3.0   5",
            "  2   a 4.0   3",
            "  2   b 7.0   0",
        ],
    )

    t2b = table.unique(t, key2, keep="last")
    assert sort_eq(
        t2b.pformat(),
        [
            " a   b   c   d ",
            "--- --- --- ---",
            "  0   a 0.0   4",
            "  1   a 2.0   6",
            "  1   c 3.0   5",
            "  2   a 4.0   3",
            "  2   b 5.0   1",
        ],
    )
    t2c = table.unique(t, key2, keep="none")
    assert sort_eq(
        t2c.pformat(),
        [
            " a   b   c   d ",
            "--- --- --- ---",
            "  0   a 0.0   4",
            "  2   a 4.0   3",
        ],
    )

    key2 = ["a", "a"]
    with pytest.raises(ValueError) as exc:
        t2a = table.unique(t, key2)
    assert exc.value.args[0] == "duplicate key names"

    with pytest.raises(ValueError) as exc:
        table.unique(t, key2, keep=True)
    assert exc.value.args[0] == "'keep' should be one of 'first', 'last', 'none'"

    t1_m = operation_table_type(t1a, masked=True)
    t1_m["a"].mask[1] = True

    with pytest.raises(ValueError) as exc:
        t1_mu = table.unique(t1_m)
    assert (
        exc.value.args[0] == "cannot use columns with masked values as keys; "
        "remove column 'a' from keys and rerun unique()"
    )

    t1_mu = table.unique(t1_m, silent=True)
    assert t1_mu.masked is False
    assert t1_mu.pformat() == [
        " a   b   c   d ",
        "--- --- --- ---",
        "  0   a 0.0   4",
        "  2   b 7.0   0",
        " --   c 3.0   5",
    ]

    with pytest.raises(ValueError):
        t1_mu = table.unique(t1_m, silent=True, keys="a")

    t1_m = operation_table_type(t, masked=True)
    t1_m["a"].mask[1] = True
    t1_m["d"].mask[3] = True

    # Test that multiple masked key columns get removed in the correct
    # order
    t1_mu = table.unique(t1_m, keys=["d", "a", "b"], silent=True)
    assert t1_mu.masked is False
    assert t1_mu.pformat() == [
        " a   b   c   d ",
        "--- --- --- ---",
        "  2   a 4.0  --",
        "  2   b 7.0   0",
        " --   c 3.0   5",
    ]


def test_vstack_bytes(operation_table_type):
    """
    Test for issue #5617 when vstack'ing bytes columns in Py3.
    This is really an upstream numpy issue numpy/numpy/#8403.
    """
    t = operation_table_type([[b"a"]], names=["a"])
    assert t["a"].itemsize == 1

    t2 = table.vstack([t, t])
    assert len(t2) == 2
    assert t2["a"].itemsize == 1


def test_vstack_unicode():
    """
    Test for problem related to issue #5617 when vstack'ing *unicode*
    columns.  In this case the character size gets multiplied by 4.
    """
    t = table.Table([["a"]], names=["a"])
    assert t["a"].itemsize == 4  # 4-byte / char for U dtype

    t2 = table.vstack([t, t])
    assert len(t2) == 2
    assert t2["a"].itemsize == 4


def test_join_mixins_time_quantity():
    """
    Test for table join using non-ndarray key columns.
    """
    tm1 = Time([2, 1, 2], format="cxcsec")
    q1 = [2, 1, 1] * u.m
    idx1 = [1, 2, 3]
    tm2 = Time([2, 3], format="cxcsec")
    q2 = [2, 3] * u.m
    idx2 = [10, 20]
    t1 = Table([tm1, q1, idx1], names=["tm", "q", "idx"])
    t2 = Table([tm2, q2, idx2], names=["tm", "q", "idx"])
    # Output:
    #
    # <Table length=4>
    #         tm            q    idx_1 idx_2
    #                       m
    #       object       float64 int64 int64
    # ------------------ ------- ----- -----
    # 0.9999999999969589     1.0     2    --
    #   2.00000000000351     1.0     3    --
    #   2.00000000000351     2.0     1    10
    #  3.000000000000469     3.0    --    20

    t12 = table.join(t1, t2, join_type="outer", keys=["tm", "q"])
    # Key cols are lexically sorted
    assert np.all(t12["tm"] == Time([1, 2, 2, 3], format="cxcsec"))
    assert np.all(t12["q"] == [1, 1, 2, 3] * u.m)
    assert np.all(t12["idx_1"] == np.ma.array([2, 3, 1, 0], mask=[0, 0, 0, 1]))
    assert np.all(t12["idx_2"] == np.ma.array([0, 0, 10, 20], mask=[1, 1, 0, 0]))


def test_join_mixins_not_sortable():
    """
    Test for table join using non-ndarray key columns that are not sortable.
    """
    sc = SkyCoord([1, 2], [3, 4], unit="deg,deg")
    t1 = Table([sc, [1, 2]], names=["sc", "idx1"])
    t2 = Table([sc, [10, 20]], names=["sc", "idx2"])

    with pytest.raises(TypeError, match="one or more key columns are not sortable"):
        table.join(t1, t2, keys="sc")


def test_join_non_1d_key_column():
    c1 = [[1, 2], [3, 4]]
    c2 = [1, 2]
    t1 = Table([c1, c2], names=["a", "b"])
    t2 = t1.copy()
    with pytest.raises(ValueError, match="key column 'a' must be 1-d"):
        table.join(t1, t2, keys="a")


def test_argsort_time_column():
    """Regression test for #10823."""
    times = Time(["2016-01-01", "2018-01-01", "2017-01-01"])
    t = Table([times], names=["time"])
    i = t.argsort("time")
    assert np.all(i == times.argsort())


def test_sort_indexed_table():
    """Test fix for #9473 and #6545 - and another regression test for #10823."""
    t = Table([[1, 3, 2], [6, 4, 5]], names=("a", "b"))
    t.add_index("a")
    t.sort("a")
    assert np.all(t["a"] == [1, 2, 3])
    assert np.all(t["b"] == [6, 5, 4])
    t.sort("b")
    assert np.all(t["b"] == [4, 5, 6])
    assert np.all(t["a"] == [3, 2, 1])

    times = ["2016-01-01", "2018-01-01", "2017-01-01"]
    tm = Time(times)
    t2 = Table([tm, [3, 2, 1]], names=["time", "flux"])
    t2.sort("flux")
    assert np.all(t2["flux"] == [1, 2, 3])
    t2.sort("time")
    assert np.all(t2["flux"] == [3, 1, 2])
    assert np.all(t2["time"] == tm[[0, 2, 1]])

    # Using the table as a TimeSeries implicitly sets the index, so
    # this test is a bit different from the above.
    from astropy.timeseries import TimeSeries

    ts = TimeSeries(time=times)
    ts["flux"] = [3, 2, 1]
    ts.sort("flux")
    assert np.all(ts["flux"] == [1, 2, 3])
    ts.sort("time")
    assert np.all(ts["flux"] == [3, 1, 2])
    assert np.all(ts["time"] == tm[[0, 2, 1]])


def test_get_out_class():
    c = table.Column([1, 2])
    mc = table.MaskedColumn([1, 2])
    q = [1, 2] * u.m

    assert _get_out_class([c, mc]) is mc.__class__
    assert _get_out_class([mc, c]) is mc.__class__
    assert _get_out_class([c, c]) is c.__class__
    assert _get_out_class([c]) is c.__class__

    with pytest.raises(ValueError):
        _get_out_class([c, q])

    with pytest.raises(ValueError):
        _get_out_class([q, c])


def test_masking_required_exception():
    """
    Test that outer join, hstack and vstack fail for a mixin column which
    does not support masking.
    """
    col = table.NdarrayMixin([0, 1, 2, 3])
    t1 = table.QTable([[1, 2, 3, 4], col], names=["a", "b"])
    t2 = table.QTable([[1, 2], col[:2]], names=["a", "c"])

    with pytest.raises(NotImplementedError) as err:
        table.vstack([t1, t2], join_type="outer")
    assert "vstack unavailable" in str(err.value)

    with pytest.raises(NotImplementedError) as err:
        table.hstack([t1, t2], join_type="outer")
    assert "hstack requires masking" in str(err.value)

    with pytest.raises(NotImplementedError) as err:
        table.join(t1, t2, join_type="outer")
    assert "join requires masking" in str(err.value)


def test_stack_columns():
    c = table.Column([1, 2])
    mc = table.MaskedColumn([1, 2])
    q = [1, 2] * u.m
    time = Time(["2001-01-02T12:34:56", "2001-02-03T00:01:02"])
    sc = SkyCoord([1, 2], [3, 4], unit="deg")
    cq = table.Column([11, 22], unit=u.m)

    t = table.hstack([c, q])
    assert t.__class__ is table.QTable
    assert t.masked is False
    t = table.hstack([q, c])
    assert t.__class__ is table.QTable
    assert t.masked is False

    t = table.hstack([mc, q])
    assert t.__class__ is table.QTable
    assert t.masked is False

    t = table.hstack([c, mc])
    assert t.__class__ is table.Table
    assert t.masked is False

    t = table.vstack([q, q])
    assert t.__class__ is table.QTable

    t = table.vstack([c, c])
    assert t.__class__ is table.Table

    t = table.hstack([c, time])
    assert t.__class__ is table.Table
    t = table.hstack([c, sc])
    assert t.__class__ is table.Table
    t = table.hstack([q, time, sc])
    assert t.__class__ is table.QTable

    with pytest.raises(ValueError):
        table.vstack([c, q])

    with pytest.raises(ValueError):
        t = table.vstack([q, cq])


def test_mixin_join_regression():
    # This used to trigger a ValueError:
    # ValueError: NumPy boolean array indexing assignment cannot assign
    # 6 input values to the 4 output values where the mask is true

    t1 = QTable()
    t1["index"] = [1, 2, 3, 4, 5]
    t1["flux1"] = [2, 3, 2, 1, 1] * u.Jy
    t1["flux2"] = [2, 3, 2, 1, 1] * u.Jy

    t2 = QTable()
    t2["index"] = [3, 4, 5, 6]
    t2["flux1"] = [2, 1, 1, 3] * u.Jy
    t2["flux2"] = [2, 1, 1, 3] * u.Jy

    t12 = table.join(t1, t2, keys=("index", "flux1", "flux2"), join_type="outer")

    assert len(t12) == 6
