# Licensed under a 3-clause BSD style license - see LICENSE.rst

from contextlib import nullcontext

import numpy as np
import pytest

from astropy import coordinates, time
from astropy import units as u
from astropy.table import Column, NdarrayMixin, QTable, Table, table_helpers, unique
from astropy.tests.helper import PYTEST_LT_8_0
from astropy.time import Time
from astropy.utils.compat import NUMPY_LT_1_22_1
from astropy.utils.exceptions import AstropyUserWarning


def sort_eq(list1, list2):
    return sorted(list1) == sorted(list2)


def test_column_group_by(T1q):
    """Test grouping a Column by various key types."""
    # T1q["a"] could be Column or Quantity, so force the object we want to group to be
    # Column. Then later we are using the "a" column as a grouping key.
    t1a = Column(T1q["a"])
    unit = T1q["a"].unit or 1

    # Group by a Column (i.e. numpy array)
    t1ag = t1a.group_by(T1q["a"])
    keys = t1ag.groups.keys
    assert np.all(t1ag.groups.indices == np.array([0, 1, 4, 8]))
    assert np.all(keys == np.array([0, 1, 2]) * unit)

    # Group by a Table and numpy structured array
    for t1ag, key_unit in (
        (t1a.group_by(T1q["a", "b"]), unit),
        (t1a.group_by(T1q["a", "b"].as_array()), 1),
    ):
        assert np.all(t1ag.groups.indices == np.array([0, 1, 3, 4, 5, 7, 8]))
        keys = t1ag.groups.keys
        assert keys.dtype.names == ("a", "b")
        assert np.all(keys["a"] == np.array([0, 1, 1, 2, 2, 2]) * key_unit)
        assert np.all(keys["b"] == np.array(["a", "a", "b", "a", "b", "c"]))


def test_column_group_by_no_argsort(T1b):
    t1a = T1b["a"]
    with pytest.raises(
        TypeError, match=r"keys input \(list\) must have an `argsort` method"
    ):
        # Pass a Python list with no argsort method
        t1a.group_by(list(range(len(t1a))))


def test_table_group_by(T1):
    """
    Test basic table group_by functionality for possible key types and for
    masked/unmasked tables.
    """
    for masked in (False, True):
        t1 = QTable(T1, masked=masked)
        # Group by a single column key specified by name
        tg = t1.group_by("a")
        assert np.all(tg.groups.indices == np.array([0, 1, 4, 8]))
        assert str(tg.groups) == "<TableGroups indices=[0 1 4 8]>"
        assert str(tg["a"].groups) == "<ColumnGroups indices=[0 1 4 8]>"

        # Sorted by 'a' and in original order for rest
        assert tg.pformat() == [
            " a   b   c   d   q ",
            "                 m ",
            "--- --- --- --- ---",
            "  0   a 0.0   4 4.0",
            "  1   b 3.0   5 5.0",
            "  1   a 2.0   6 6.0",
            "  1   a 1.0   7 7.0",
            "  2   c 7.0   0 0.0",
            "  2   b 5.0   1 1.0",
            "  2   b 6.0   2 2.0",
            "  2   a 4.0   3 3.0",
        ]
        assert tg.meta["ta"] == 1
        assert tg["c"].meta["a"] == 1
        assert tg["c"].description == "column c"

        # Group by a table column
        tg2 = t1.group_by(t1["a"])
        assert tg.pformat() == tg2.pformat()

        # Group by two columns spec'd by name
        for keys in (["a", "b"], ("a", "b")):
            tg = t1.group_by(keys)
            assert np.all(tg.groups.indices == np.array([0, 1, 3, 4, 5, 7, 8]))
            # Sorted by 'a', 'b' and in original order for rest
            assert tg.pformat() == [
                " a   b   c   d   q ",
                "                 m ",
                "--- --- --- --- ---",
                "  0   a 0.0   4 4.0",
                "  1   a 2.0   6 6.0",
                "  1   a 1.0   7 7.0",
                "  1   b 3.0   5 5.0",
                "  2   a 4.0   3 3.0",
                "  2   b 5.0   1 1.0",
                "  2   b 6.0   2 2.0",
                "  2   c 7.0   0 0.0",
            ]

        # Group by a Table
        tg2 = t1.group_by(t1["a", "b"])
        assert tg.pformat() == tg2.pformat()

        # Group by a structured array
        tg2 = t1.group_by(t1["a", "b"].as_array())
        assert tg.pformat() == tg2.pformat()

        # Group by a simple ndarray
        tg = t1.group_by(np.array([0, 1, 0, 1, 2, 1, 0, 0]))
        assert np.all(tg.groups.indices == np.array([0, 4, 7, 8]))
        assert tg.pformat() == [
            " a   b   c   d   q ",
            "                 m ",
            "--- --- --- --- ---",
            "  2   c 7.0   0 0.0",
            "  2   b 6.0   2 2.0",
            "  1   a 2.0   6 6.0",
            "  1   a 1.0   7 7.0",
            "  2   b 5.0   1 1.0",
            "  2   a 4.0   3 3.0",
            "  1   b 3.0   5 5.0",
            "  0   a 0.0   4 4.0",
        ]


def test_groups_keys(T1m: QTable):
    tg = T1m.group_by("a")
    unit = T1m["a"].unit or 1
    keys = tg.groups.keys
    assert keys.dtype.names == ("a",)
    assert np.all(keys["a"] == np.array([0, 1, 2]) * unit)

    tg = T1m.group_by(["a", "b"])
    keys = tg.groups.keys
    assert keys.dtype.names == ("a", "b")
    assert np.all(keys["a"] == np.array([0, 1, 1, 2, 2, 2]) * unit)
    assert np.all(keys["b"] == np.array(["a", "a", "b", "a", "b", "c"]))

    # Grouping by Column ignores column name
    tg = T1m.group_by(T1m["b"])
    keys = tg.groups.keys
    assert keys.dtype.names is None


def test_groups_keys_time(T1b: QTable):
    """Group a table with a time column using that column as a key."""
    T1b = T1b.copy()
    T1b["a"] = Time(T1b["a"], format="cxcsec")

    tg = T1b.group_by("a")
    keys = tg.groups.keys
    assert keys.dtype.names == ("a",)
    assert np.all(keys["a"] == Time(np.array([0, 1, 2]), format="cxcsec"))

    tg = T1b.group_by(["a", "b"])
    keys = tg.groups.keys
    assert keys.dtype.names == ("a", "b")
    assert np.all(keys["a"] == Time(np.array([0, 1, 1, 2, 2, 2]), format="cxcsec"))
    assert np.all(keys["b"] == np.array(["a", "a", "b", "a", "b", "c"]))


def test_groups_iterator(T1):
    tg = T1.group_by("a")
    for ii, group in enumerate(tg.groups):
        assert group.pformat() == tg.groups[ii].pformat()
        assert group["a"][0] == tg["a"][tg.groups.indices[ii]]


def test_grouped_copy(T1):
    """
    Test that copying a table or column copies the groups properly
    """
    for masked in (False, True):
        t1 = QTable(T1, masked=masked)
        tg = t1.group_by("a")
        tgc = tg.copy()
        assert np.all(tgc.groups.indices == tg.groups.indices)
        assert np.all(tgc.groups.keys == tg.groups.keys)

        tac = tg["a"].copy()
        assert np.all(tac.groups.indices == tg["a"].groups.indices)

        c1 = t1["a"].copy()
        gc1 = c1.group_by(t1["a"])
        gc1c = gc1.copy()
        assert np.all(gc1c.groups.indices == np.array([0, 1, 4, 8]))


def test_grouped_slicing(T1):
    """
    Test that slicing a table removes previous grouping
    """

    for masked in (False, True):
        t1 = QTable(T1, masked=masked)

        # Regular slice of a table
        tg = t1.group_by("a")
        tg2 = tg[3:5]
        assert np.all(tg2.groups.indices == np.array([0, len(tg2)]))
        assert tg2.groups.keys is None


def test_group_column_from_table(T1):
    """
    Group a column that is part of a table
    """
    cg = T1["c"].group_by(np.array(T1["a"]))
    assert np.all(cg.groups.keys == np.array([0, 1, 2]))
    assert np.all(cg.groups.indices == np.array([0, 1, 4, 8]))


def test_table_groups_mask_index(T1):
    """
    Use boolean mask as item in __getitem__ for groups
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked).group_by("a")

        t2 = t1.groups[np.array([True, False, True])]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[2].pformat()
        assert np.all(t2.groups.keys["a"] == np.array([0, 2]))


def test_table_groups_array_index(T1):
    """
    Use numpy array as item in __getitem__ for groups
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked).group_by("a")

        t2 = t1.groups[np.array([0, 2])]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[2].pformat()
        assert np.all(t2.groups.keys["a"] == np.array([0, 2]))


def test_table_groups_slicing(T1):
    """
    Test that slicing table groups works
    """

    for masked in (False, True):
        t1 = Table(T1, masked=masked).group_by("a")

        # slice(0, 2)
        t2 = t1.groups[0:2]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[1].pformat()
        assert np.all(t2.groups.keys["a"] == np.array([0, 1]))

        # slice(1, 2)
        t2 = t1.groups[1:2]
        assert len(t2.groups) == 1
        assert t2.groups[0].pformat() == t1.groups[1].pformat()
        assert np.all(t2.groups.keys["a"] == np.array([1]))

        # slice(0, 3, 2)
        t2 = t1.groups[0:3:2]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[2].pformat()
        assert np.all(t2.groups.keys["a"] == np.array([0, 2]))


def test_grouped_item_access(T1):
    """
    Test that column slicing preserves grouping
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked)

        # Regular slice of a table
        tg = t1.group_by("a")
        tgs = tg["a", "c", "d"]
        assert np.all(tgs.groups.keys == tg.groups.keys)
        assert np.all(tgs.groups.indices == tg.groups.indices)
        tgsa = tgs.groups.aggregate(np.sum)
        assert tgsa.pformat() == [
            " a   c    d ",
            "--- ---- ---",
            "  0  0.0   4",
            "  1  6.0  18",
            "  2 22.0   6",
        ]

        tgs = tg["c", "d"]
        assert np.all(tgs.groups.keys == tg.groups.keys)
        assert np.all(tgs.groups.indices == tg.groups.indices)
        tgsa = tgs.groups.aggregate(np.sum)
        assert tgsa.pformat() == [
            " c    d ",
            "---- ---",
            " 0.0   4",
            " 6.0  18",
            "22.0   6",
        ]


def test_mutable_operations(T1):
    """
    Operations like adding or deleting a row should removing grouping,
    but adding or removing or renaming a column should retain grouping.
    """
    for masked in (False, True):
        t1 = QTable(T1, masked=masked)

        # add row
        tg = t1.group_by("a")
        tg.add_row((0, "a", 3.0, 4, 4 * u.m))
        assert np.all(tg.groups.indices == np.array([0, len(tg)]))
        assert tg.groups.keys is None

        # remove row
        tg = t1.group_by("a")
        tg.remove_row(4)
        assert np.all(tg.groups.indices == np.array([0, len(tg)]))
        assert tg.groups.keys is None

        # add column
        tg = t1.group_by("a")
        indices = tg.groups.indices.copy()
        tg.add_column(Column(name="e", data=np.arange(len(tg))))
        assert np.all(tg.groups.indices == indices)
        assert np.all(tg["e"].groups.indices == indices)
        assert np.all(tg["e"].groups.keys == tg.groups.keys)

        # remove column (not key column)
        tg = t1.group_by("a")
        tg.remove_column("b")
        assert np.all(tg.groups.indices == indices)
        # Still has original key col names
        assert tg.groups.keys.dtype.names == ("a",)
        assert np.all(tg["a"].groups.indices == indices)

        # remove key column
        tg = t1.group_by("a")
        tg.remove_column("a")
        assert np.all(tg.groups.indices == indices)
        assert tg.groups.keys.dtype.names == ("a",)
        assert np.all(tg["b"].groups.indices == indices)

        # rename key column
        tg = t1.group_by("a")
        tg.rename_column("a", "aa")
        assert np.all(tg.groups.indices == indices)
        assert tg.groups.keys.dtype.names == ("a",)
        assert np.all(tg["aa"].groups.indices == indices)


def test_group_by_masked(T1):
    t1m = QTable(T1, masked=True)
    t1m["c"].mask[4] = True
    t1m["d"].mask[5] = True
    assert t1m.group_by("a").pformat() == [
        " a   b   c   d   q ",
        "                 m ",
        "--- --- --- --- ---",
        "  0   a  --   4 4.0",
        "  1   b 3.0  -- 5.0",
        "  1   a 2.0   6 6.0",
        "  1   a 1.0   7 7.0",
        "  2   c 7.0   0 0.0",
        "  2   b 5.0   1 1.0",
        "  2   b 6.0   2 2.0",
        "  2   a 4.0   3 3.0",
    ]


def test_group_by_errors(T1):
    """
    Appropriate errors get raised.
    """
    # Bad column name as string
    with pytest.raises(ValueError):
        T1.group_by("f")

    # Bad column names in list
    with pytest.raises(ValueError):
        T1.group_by(["f", "g"])

    # Wrong length array
    with pytest.raises(ValueError):
        T1.group_by(np.array([1, 2]))

    # Wrong type
    with pytest.raises(TypeError):
        T1.group_by(None)

    # Masked key column
    t1 = QTable(T1, masked=True)
    t1["a"].mask[4] = True
    with pytest.raises(ValueError):
        t1.group_by("a")


def test_groups_keys_meta(T1):
    """
    Make sure the keys meta['grouped_by_table_cols'] is working.
    """
    # Group by column in this table
    tg = T1.group_by("a")
    assert tg.groups.keys.meta["grouped_by_table_cols"] is True
    assert tg["c"].groups.keys.meta["grouped_by_table_cols"] is True
    assert tg.groups[1].groups.keys.meta["grouped_by_table_cols"] is True
    assert (
        tg["d"]
        .groups[np.array([False, True, True])]
        .groups.keys.meta["grouped_by_table_cols"]
        is True
    )

    # Group by external Table
    tg = T1.group_by(T1["a", "b"])
    assert tg.groups.keys.meta["grouped_by_table_cols"] is False
    assert tg["c"].groups.keys.meta["grouped_by_table_cols"] is False
    assert tg.groups[1].groups.keys.meta["grouped_by_table_cols"] is False

    # Group by external numpy array
    tg = T1.group_by(T1["a", "b"].as_array())
    assert not hasattr(tg.groups.keys, "meta")
    assert not hasattr(tg["c"].groups.keys, "meta")

    # Group by Column
    tg = T1.group_by(T1["a"])
    assert "grouped_by_table_cols" not in tg.groups.keys.meta
    assert "grouped_by_table_cols" not in tg["c"].groups.keys.meta


def test_table_aggregate(T1):
    """
    Aggregate a table
    """
    # Table with only summable cols
    t1 = T1["a", "c", "d"]
    tg = t1.group_by("a")
    tga = tg.groups.aggregate(np.sum)
    assert tga.pformat() == [
        " a   c    d ",
        "--- ---- ---",
        "  0  0.0   4",
        "  1  6.0  18",
        "  2 22.0   6",
    ]
    # Reverts to default groups
    assert np.all(tga.groups.indices == np.array([0, 3]))
    assert tga.groups.keys is None

    # metadata survives
    assert tga.meta["ta"] == 1
    assert tga["c"].meta["a"] == 1
    assert tga["c"].description == "column c"

    # Aggregate with np.sum with masked elements.  This results
    # in one group with no elements, hence a nan result and conversion
    # to float for the 'd' column.
    t1m = QTable(T1, masked=True)
    t1m["c"].mask[4:6] = True
    t1m["d"].mask[4:6] = True
    tg = t1m.group_by("a")

    if PYTEST_LT_8_0:
        ctx = nullcontext()
    else:
        ctx = pytest.warns(AstropyUserWarning, match="Cannot aggregate column")

    with pytest.warns(UserWarning, match="converting a masked element to nan"), ctx:
        tga = tg.groups.aggregate(np.sum)

    assert tga.pformat() == [
        " a   c    d    q  ",
        "               m  ",
        "--- ---- ---- ----",
        "  0  nan  nan  4.0",
        "  1  3.0 13.0 18.0",
        "  2 22.0  6.0  6.0",
    ]

    # Aggregate with np.sum with masked elements, but where every
    # group has at least one remaining (unmasked) element.  Then
    # the int column stays as an int.
    t1m = QTable(t1, masked=True)
    t1m["c"].mask[5] = True
    t1m["d"].mask[5] = True
    tg = t1m.group_by("a")
    tga = tg.groups.aggregate(np.sum)
    assert tga.pformat() == [
        " a   c    d ",
        "--- ---- ---",
        "  0  0.0   4",
        "  1  3.0  13",
        "  2 22.0   6",
    ]

    # Aggregate with a column type that cannot by supplied to the aggregating
    # function.  This raises a warning but still works.
    tg = T1.group_by("a")
    with pytest.warns(AstropyUserWarning, match="Cannot aggregate column"):
        tga = tg.groups.aggregate(np.sum)
    assert tga.pformat() == [
        " a   c    d   q  ",
        "              m  ",
        "--- ---- --- ----",
        "  0  0.0   4  4.0",
        "  1  6.0  18 18.0",
        "  2 22.0   6  6.0",
    ]


def test_table_aggregate_reduceat(T1):
    """
    Aggregate table with functions which have a reduceat method
    """

    # Comparison functions without reduceat
    def np_mean(x):
        return np.mean(x)

    def np_sum(x):
        return np.sum(x)

    def np_add(x):
        return np.add(x)

    # Table with only summable cols
    t1 = T1["a", "c", "d"]
    tg = t1.group_by("a")
    # Comparison
    tga_r = tg.groups.aggregate(np.sum)
    tga_a = tg.groups.aggregate(np.add)
    tga_n = tg.groups.aggregate(np_sum)

    assert np.all(tga_r == tga_n)
    assert np.all(tga_a == tga_n)
    assert tga_n.pformat() == [
        " a   c    d ",
        "--- ---- ---",
        "  0  0.0   4",
        "  1  6.0  18",
        "  2 22.0   6",
    ]

    tga_r = tg.groups.aggregate(np.mean)
    tga_n = tg.groups.aggregate(np_mean)
    assert np.all(tga_r == tga_n)
    assert tga_n.pformat() == [
        " a   c   d ",
        "--- --- ---",
        "  0 0.0 4.0",
        "  1 2.0 6.0",
        "  2 5.5 1.5",
    ]

    # Binary ufunc np_add should raise warning without reduceat
    t2 = T1["a", "c"]
    tg = t2.group_by("a")

    with pytest.warns(AstropyUserWarning, match="Cannot aggregate column"):
        tga = tg.groups.aggregate(np_add)
    assert tga.pformat() == [" a ", "---", "  0", "  1", "  2"]


def test_column_aggregate(T1):
    """
    Aggregate a single table column
    """
    for masked in (False, True):
        tg = QTable(T1, masked=masked).group_by("a")
        tga = tg["c"].groups.aggregate(np.sum)
        assert tga.pformat() == [" c  ", "----", " 0.0", " 6.0", "22.0"]


@pytest.mark.skipif(
    NUMPY_LT_1_22_1, reason="https://github.com/numpy/numpy/issues/20699"
)
def test_column_aggregate_f8():
    """https://github.com/astropy/astropy/issues/12706"""
    # Just want to make sure it does not crash again.
    for masked in (False, True):
        tg = Table({"a": np.arange(2, dtype=">f8")}, masked=masked).group_by("a")
        tga = tg["a"].groups.aggregate(np.sum)
        assert tga.pformat() == [" a ", "---", "0.0", "1.0"]


def test_table_filter():
    """
    Table groups filtering
    """

    def all_positive(table, key_colnames):
        return all(
            np.all(table[colname] >= 0)
            for colname in table.colnames
            if colname not in key_colnames
        )

    # Negative value in 'a' column should not filter because it is a key col
    t = Table.read(
        [
            " a c d",
            " -2 7.0 0",
            " -2 5.0 1",
            " 0 0.0 4",
            " 1 3.0 5",
            " 1 2.0 -6",
            " 1 1.0 7",
            " 3 3.0 5",
            " 3 -2.0 6",
            " 3 1.0 7",
        ],
        format="ascii",
    )
    tg = t.group_by("a")
    t2 = tg.groups.filter(all_positive)
    assert t2.groups[0].pformat() == [
        " a   c   d ",
        "--- --- ---",
        " -2 7.0   0",
        " -2 5.0   1",
    ]
    assert t2.groups[1].pformat() == [" a   c   d ", "--- --- ---", "  0 0.0   4"]


def test_column_filter():
    """
    Table groups filtering
    """

    def all_positive(column):
        if np.any(column < 0):
            return False
        return True

    # Negative value in 'a' column should not filter because it is a key col
    t = Table.read(
        [
            " a c d",
            " -2 7.0 0",
            " -2 5.0 1",
            " 0 0.0 4",
            " 1 3.0 5",
            " 1 2.0 -6",
            " 1 1.0 7",
            " 3 3.0 5",
            " 3 -2.0 6",
            " 3 1.0 7",
        ],
        format="ascii",
    )
    tg = t.group_by("a")
    c2 = tg["c"].groups.filter(all_positive)
    assert len(c2.groups) == 3
    assert c2.groups[0].pformat() == [" c ", "---", "7.0", "5.0"]
    assert c2.groups[1].pformat() == [" c ", "---", "0.0"]
    assert c2.groups[2].pformat() == [" c ", "---", "3.0", "2.0", "1.0"]


def test_group_mixins():
    """
    Test grouping a table with mixin columns
    """
    # Setup mixins
    idx = np.arange(4)
    x = np.array([3.0, 1.0, 2.0, 1.0])
    q = x * u.m
    lon = coordinates.Longitude(x * u.deg)
    lat = coordinates.Latitude(x * u.deg)
    # For Time do J2000.0 + few * 0.1 ns (this requires > 64 bit precision)
    tm = time.Time(2000, format="jyear") + time.TimeDelta(x * 1e-10, format="sec")
    sc = coordinates.SkyCoord(ra=lon, dec=lat)
    aw = table_helpers.ArrayWrapper(x)
    nd = np.array([(3, "c"), (1, "a"), (2, "b"), (1, "a")], dtype="<i4,|S1").view(
        NdarrayMixin
    )

    qt = QTable(
        [idx, x, q, lon, lat, tm, sc, aw, nd],
        names=["idx", "x", "q", "lon", "lat", "tm", "sc", "aw", "nd"],
    )

    # Test group_by with each supported mixin type
    mixin_keys = ["x", "q", "lon", "lat", "tm", "sc", "aw", "nd"]
    for key in mixin_keys:
        qtg = qt.group_by(key)

        # Test that it got the sort order correct
        assert np.all(qtg["idx"] == [1, 3, 2, 0])

        # Test that the groups are right
        # Note: skip testing SkyCoord column because that doesn't have equality
        for name in ["x", "q", "lon", "lat", "tm", "aw", "nd"]:
            assert np.all(qt[name][[1, 3]] == qtg.groups[0][name])
            assert np.all(qt[name][[2]] == qtg.groups[1][name])
            assert np.all(qt[name][[0]] == qtg.groups[2][name])

    # Test that unique also works with mixins since most of the work is
    # done with group_by().  This is using *every* mixin as key.
    uqt = unique(qt, keys=mixin_keys)
    assert len(uqt) == 3
    assert np.all(uqt["idx"] == [1, 2, 0])
    assert np.all(uqt["x"] == [1.0, 2.0, 3.0])

    # Column group_by() with mixins
    idxg = qt["idx"].group_by(qt[mixin_keys])
    assert np.all(idxg == [1, 3, 2, 0])


@pytest.mark.parametrize(
    "col",
    [
        time.TimeDelta([1, 2], format="sec"),
        time.Time([1, 2], format="cxcsec"),
        coordinates.SkyCoord([1, 2], [3, 4], unit="deg,deg"),
    ],
)
def test_group_mixins_unsupported(col):
    """Test that aggregating unsupported mixins produces a warning only"""

    t = Table([[1, 1], [3, 4], col], names=["a", "b", "mix"])
    tg = t.group_by("a")

    with pytest.warns(AstropyUserWarning, match="Cannot aggregate column 'mix'"):
        tg.groups.aggregate(np.sum)


@pytest.mark.parametrize("add_index", [False, True])
def test_group_stable_sort(add_index):
    """Test that group_by preserves the order of the table.

    This table has 5 groups with an average of 200 rows per group, so it is not
    statistically possible that the groups will be in order by chance.

    This tests explicitly the case where grouping is done via the index sort.
    See: https://github.com/astropy/astropy/issues/14882
    """
    a = np.random.randint(0, 5, 1000)
    b = np.arange(len(a))
    t = Table([a, b], names=["a", "b"])
    if add_index:
        t.add_index("a")
    tg = t.group_by("a")
    for grp in tg.groups:
        assert np.all(grp["b"] == np.sort(grp["b"]))
