# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import table
from astropy import units as u
from astropy.table import MaskedColumn, Table
from astropy.time import Time, TimeDelta
from astropy.utils.compat.optional_deps import (
    HAS_DASK,
    HAS_DUCKDB,
    HAS_NARWHALS,
    HAS_PANDAS,
    HAS_POLARS,
    HAS_PYARROW,
)

from .conftest import MIXIN_COLS


@pytest.mark.parametrize(
    "backend,use_legacy",
    [
        pytest.param(
            "pandas",
            True,
            marks=pytest.mark.skipif(not HAS_PANDAS, reason="requires pandas"),
            id="pandas-legacy",
        ),
        pytest.param(
            "pandas",
            False,
            marks=pytest.mark.skipif(
                not HAS_NARWHALS or not HAS_PANDAS,
                reason="requires narwhals and pandas",
            ),
            id="pandas-generic",
        ),
        pytest.param(
            "polars",
            False,
            marks=pytest.mark.skipif(
                not HAS_NARWHALS or not HAS_POLARS,
                reason="requires narwhals and polars",
            ),
            id="polars-generic",
        ),
        pytest.param(
            "pyarrow",
            False,
            marks=pytest.mark.skipif(
                not HAS_NARWHALS or not HAS_PYARROW,
                reason="requires narwhals and pyarrow",
            ),
            id="pyarrow-generic",
        ),
        pytest.param(
            "dask",
            False,
            marks=pytest.mark.skipif(
                not HAS_NARWHALS or not HAS_DASK,
                reason="requires narwhals and dask",
            ),
            id="dask-generic",
        ),
        pytest.param(
            "duckdb",
            False,
            marks=pytest.mark.skipif(
                not HAS_NARWHALS or not HAS_DUCKDB,
                reason="requires narwhals and duckdb",
            ),
            id="duckdb-generic",
        ),
    ],
)
class TestDataFrameConversion:
    """Test DataFrame conversion functionality for both legacy pandas and generic backends."""

    def _to_dataframe(self, table, backend, use_legacy, **kwargs):
        """Convert table to dataframe using appropriate method."""
        if use_legacy:
            if backend != "pandas":
                raise ValueError(
                    "Legacy conversion is only supported for the pandas backend."
                )
            return table.to_pandas(**kwargs)

        # Lazy backends cannot be exported to
        if backend in ("dask", "duckdb"):
            pytest.skip("Lazy backends cannot be converted back to Table")

        return table.to_df(backend, **kwargs)

    def _from_dataframe(self, df, backend, use_legacy, **kwargs):
        """Convert dataframe to table using appropriate method."""
        if use_legacy:
            if backend != "pandas":
                raise ValueError(
                    "Legacy conversion is only supported for the pandas backend."
                )
            return table.Table.from_pandas(df, **kwargs)
        return table.Table.from_df(df, **kwargs)

    def test_simple(self, backend, use_legacy):
        """Test basic endianness and data type handling."""
        t = table.Table()

        for endian in ["<", ">", "="]:
            for kind in ["f", "i"]:
                for byte in ["2", "4", "8"]:
                    dtype = np.dtype(endian + kind + byte)
                    x = np.array([1, 2, 3], dtype=dtype)
                    t[endian + kind + byte] = x.view(x.dtype.newbyteorder(endian))

        t["u"] = ["a", "b", "c"]
        t["s"] = [b"a", b"b", b"c"]

        d = self._to_dataframe(t, backend, use_legacy)

        # Basic round-trip test
        t2 = self._from_dataframe(d, backend, use_legacy)

        for column in t.columns:
            original_col = t[column]
            roundtrip_col = t2[column]

            if column in ("u", "s"):
                assert np.all(original_col == roundtrip_col)
                if backend == "pandas":
                    # Pandas-specific checks
                    # Upstream feature of pandas
                    assert d[column].dtype == np.dtype("O")
            else:
                # Generic comparison with tolerance
                assert_allclose(original_col, roundtrip_col)

                # Compare dtypes by value, not identity (normalize endianness)
                t_dtype = original_col.dtype.newbyteorder("=")
                t2_dtype = roundtrip_col.dtype.newbyteorder("=")

                if backend == "polars" and "f2" in column:
                    # No Polars Float16 support
                    pass
                else:
                    assert t_dtype == t2_dtype

                # Pandas-specific checks
                if backend == "pandas":
                    # Pandas-specific exact comparison
                    assert np.all(t[column] == d[column])
                    if t[column].dtype.isnative:
                        assert d[column].dtype == t[column].dtype
                    else:
                        assert d[column].dtype == t[column].dtype.newbyteorder()

        # Pandas-specific endian tests
        if backend == "pandas":
            # Regression test for astropy/astropy#1156 - the following code gave a
            # ValueError: Big-endian buffer not supported on little-endian
            # compiler. We now automatically swap the endian-ness to native order
            # upon adding the arrays to the data frame.
            # Explicitly testing little/big/native endian separately -
            # regression for a case in astropy/astropy#11286 not caught by #3729.
            d[["<i4", ">i4"]]
            d[["<f4", ">f4"]]

            # Additional round-trip checks
            for column in t.columns:
                if column in ("u", "s"):
                    assert np.all(t[column] == t2[column])
                else:
                    assert_allclose(t[column], t2[column])
                if t[column].dtype.isnative:
                    assert t[column].dtype == t2[column].dtype
                else:
                    assert t[column].dtype.newbyteorder() == t2[column].dtype

    @pytest.mark.parametrize("use_IndexedTable", [False, True])
    def test_to_df_index(self, backend, use_legacy, use_IndexedTable):
        """Test indexing options for both legacy pandas and generic backends."""

        class IndexedTable(table.QTable):
            """Always index the first column"""

            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.add_index(self.colnames[0])

        tm = Time([1998, 2002], format="jyear")
        x = [1, 2]
        table_cls = IndexedTable if use_IndexedTable else table.QTable
        t = table_cls([tm, x], names=["tm", "x"])

        match backend:
            case "pandas":
                import pandas as pd

                row_index = pd.RangeIndex(0, 2, 1)
                tm_index = pd.DatetimeIndex(
                    ["1998-01-01", "2002-01-01"],
                    dtype="datetime64[ns]",
                    name="tm",
                    freq=None,
                )
                tp = self._to_dataframe(t, backend, use_legacy)

            case "polars" | "pyarrow":
                with pytest.raises(
                    ValueError,
                    match="Indexing is only supported for pandas-like backends",
                ):
                    self._to_dataframe(t, backend, use_legacy, index="tm")
                return
            case "dask" | "duckdb":
                # Lazy backends will raise ValueError in _to_dataframe
                self._to_dataframe(t, backend, use_legacy)
                return
            case _:
                raise ValueError(f"Unknown backend: {backend}")

        if not use_IndexedTable:
            assert np.all(tp.index == row_index)
            tp = self._to_dataframe(t, backend, use_legacy, index="tm")
            assert np.all(tp.index == tm_index)
            t.add_index("tm")

        tp = self._to_dataframe(t, backend, use_legacy)
        assert np.all(tp.index == tm_index)
        # Make sure writing to dataframe didn't hack the original table
        assert t["tm"].info.indices

        tp = self._to_dataframe(t, backend, use_legacy, index=True)
        assert np.all(tp.index == tm_index)

        tp = self._to_dataframe(t, backend, use_legacy, index=False)
        assert np.all(tp.index == row_index)

        with pytest.raises(ValueError, match="is not in the table columns"):
            self._to_dataframe(t, backend, use_legacy, index="not a column")

    def test_from_df_index(self, backend, use_legacy):
        """Test index handling in from_dataframe conversion."""
        tm = Time([1998, 2002], format="jyear")
        x = [1, 2]
        t = table.Table([tm, x], names=["tm", "x"])
        match backend:
            case "polars" | "pyarrow":
                # Non-pandas backends don't support indexing
                with pytest.raises(
                    ValueError,
                    match="Indexing is only supported for pandas-like backends",
                ):
                    self._to_dataframe(t, backend, use_legacy, index="tm")
                return
            case "dask" | "duckdb":
                # Lazy backends will raise ValueError in _to_dataframe
                self._to_dataframe(t, backend, use_legacy, index="tm")
                return
            case "pandas":
                tp = self._to_dataframe(t, backend, use_legacy, index="tm")
            case _:
                raise ValueError(f"Unknown backend: {backend}")

        t2 = self._from_dataframe(tp, backend, use_legacy)
        assert t2.colnames == ["x"]

        t2 = self._from_dataframe(tp, backend, use_legacy, index=True)
        assert t2.colnames == ["tm", "x"]
        assert np.allclose(t2["tm"].jyear, tm.jyear)

    def test_units(self, backend, use_legacy):
        """Test handling of units in from_dataframe conversion."""
        data = {"x": [1, 2, 3], "t": [1.3, 1.2, 1.8]}
        match backend:
            case "pandas":
                import pandas as pd

                df = pd.DataFrame(data)
            case "polars":
                import polars as pl

                df = pl.DataFrame(data)
            case "pyarrow":
                import pyarrow as pa

                df = pa.Table.from_pydict(data)
            case "dask":
                import dask.array as da
                import dask.dataframe as dd

                # Need to aggregate data into 2D dask array
                df = dd.from_dask_array(
                    da.stack(list(data.values()), axis=1), columns=list(data.keys())
                )
            case "duckdb":
                import io
                import json

                import duckdb

                # DuckDB can't read directly from dict, need to stream via JSON
                df = duckdb.read_json(io.StringIO(json.dumps(data)))
            case _:
                raise ValueError(f"Unknown backend: {backend}")

        t = self._from_dataframe(df, backend, use_legacy, units={"x": u.m, "t": u.s})

        assert t["x"].unit == u.m
        assert t["t"].unit == u.s

        # test error if not a mapping
        with pytest.raises(TypeError):
            self._from_dataframe(df, backend, use_legacy, units=[u.m, u.s])

        # test warning is raised if additional columns in units dict
        with pytest.warns(UserWarning) as record:
            self._from_dataframe(
                df, backend, use_legacy, units={"x": u.m, "t": u.s, "y": u.m}
            )
        assert len(record) == 1
        assert "{'y'}" in str(record[0].message.args[0])

    @pytest.mark.parametrize("unsigned", ["u", ""])
    @pytest.mark.parametrize("bits", [8, 16, 32, 64])
    def test_nullable_int(self, backend, use_legacy, unsigned, bits):
        """Test nullable integer handling."""
        np_dtype = f"{unsigned}int{bits}"
        c = MaskedColumn([1, 2], mask=[False, True], dtype=np_dtype)
        t = Table([c])
        df = self._to_dataframe(t, backend, use_legacy)

        t2 = self._from_dataframe(df, backend, use_legacy)
        assert str(t2["col0"].dtype) == np_dtype
        assert np.all(t2["col0"].mask == [False, True])
        assert np.all(t2["col0"] == c)

    @pytest.mark.parametrize("ndim", [1, 2, 3])
    def test_nd_columns(self, backend, use_legacy, ndim):
        """Test handling of multidimensional columns."""
        # Add one since we want the dimension of each entry to be ndim
        shape = (3,) * ndim
        colshape = (10,) + shape
        t = table.Table()
        t["a"] = np.arange(np.prod(colshape)).reshape(colshape)

        match backend:
            # Pandas and PyArrow do not support multidimensional columns
            case "pandas" | "pyarrow":
                if ndim > 1:
                    with pytest.raises(
                        ValueError,
                        match="Cannot convert a table with multidimensional columns",
                    ):
                        self._to_dataframe(t, backend, use_legacy)
                    return
            case "dask" | "duckdb":
                # Lazy backends will raise ValueError in _to_dataframe
                self._to_dataframe(t, backend, use_legacy)
                return
            case "polars":
                df = self._to_dataframe(t, backend, use_legacy)
                # Convert and check shape
                assert df["a"].dtype.shape == shape

                # Round-trip conversion
                t2 = self._from_dataframe(df, backend, use_legacy)
                assert t2["a"].shape == colshape
                assert np.all(t["a"] == t2["a"])
            case _:
                raise ValueError(f"Unknown backend: {backend}")

    def test_mixin_columns(self, backend, use_legacy):
        """Test handling of astropy mixin columns."""
        t = table.QTable()
        for name in sorted(MIXIN_COLS):
            if not name.startswith("ndarray"):
                t[name] = MIXIN_COLS[name]

        t["dt"] = TimeDelta([0, 2, 4, 6], format="sec")

        tp = self._to_dataframe(t, backend, use_legacy)
        t2 = self._from_dataframe(tp, backend, use_legacy)

        assert np.allclose(t2["quantity"], [0, 1, 2, 3])
        assert np.allclose(t2["longitude"], [0.0, 1.0, 5.0, 6.0])
        assert np.allclose(t2["latitude"], [5.0, 6.0, 10.0, 11.0])
        assert np.allclose(t2["skycoord.ra"], [0, 1, 2, 3])
        assert np.allclose(t2["skycoord.dec"], [0, 1, 2, 3])
        assert np.allclose(t2["arraywrap"], [0, 1, 2, 3])
        assert np.allclose(t2["arrayswap"], [0, 1, 2, 3])
        assert np.allclose(
            t2["earthlocation.y"], [0, 110708, 547501, 654527], rtol=0, atol=1
        )

        # Time and TimeDelta mixins that round-trip the class
        assert isinstance(t2["time"], Time)
        assert np.allclose(t2["time"].jyear, [2000, 2001, 2002, 2003])
        assert np.all(
            t2["time"].isot
            == [
                "2000-01-01T12:00:00.000",
                "2000-12-31T18:00:00.000",
                "2002-01-01T00:00:00.000",
                "2003-01-01T06:00:00.000",
            ]
        )
        assert t2["time"].format == "isot"

        # TimeDelta
        assert isinstance(t2["dt"], TimeDelta)
        assert np.allclose(t2["dt"].value, [0, 2, 4, 6])
        assert t2["dt"].format == "sec"

    def test_mixin_masked(self, backend, use_legacy):
        """Test handling of masked mixin columns."""
        tm = Time([1, 2, 3], format="cxcsec")
        dt = TimeDelta([1, 2, 3], format="sec")
        tm[1] = np.ma.masked
        dt[1] = np.ma.masked
        t = table.QTable([tm, dt], names=["tm", "dt"])

        tp = self._to_dataframe(t, backend, use_legacy)

        match backend:
            case "pandas":
                tm_nulls = tp["tm"].isnull().to_list()
                dt_nulls = tp["dt"].isnull().to_list()
            case "polars":
                tm_nulls = tp["tm"].is_null().to_list()
                dt_nulls = tp["dt"].is_null().to_list()
            case "pyarrow":
                tm_nulls = tp["tm"].is_null().to_pylist()
                dt_nulls = tp["dt"].is_null().to_pylist()
            case _:
                raise ValueError(f"Unknown backend: {backend}")

        # Common assertion for all backends
        expected_nulls = [False, True, False]
        assert tm_nulls == expected_nulls
        assert dt_nulls == expected_nulls

        # Round-trip conversion
        t2 = self._from_dataframe(tp, backend, use_legacy)

        assert np.all(t2["tm"].mask == tm.mask)
        assert np.ma.allclose(t2["tm"].jd, tm.jd, rtol=1e-14, atol=1e-14)

        assert np.all(t2["dt"].mask == dt.mask)
        assert np.ma.allclose(t2["dt"].jd, dt.jd, rtol=1e-14, atol=1e-14)

    @pytest.mark.parametrize("use_nullable_int", [True, False])
    def test_masking(self, backend, use_legacy, use_nullable_int):
        """Test handling of masked columns."""
        t = table.Table(masked=True)

        t["a"] = [1, 2, 3]
        t["a"].mask = [True, False, True]

        t["b"] = [1.0, 2.0, 3.0]
        t["b"].mask = [False, False, True]

        t["u"] = ["a", "b", "c"]
        t["u"].mask = [False, True, False]

        t["s"] = [b"a", b"b", b"c"]
        t["s"].mask = [False, True, False]

        # https://github.com/astropy/astropy/issues/7741
        t["Source"] = [2584290278794471936, 2584290038276303744, 2584288728310999296]
        t["Source"].mask = [False, False, False]

        if use_nullable_int:  # Default
            df = self._to_dataframe(
                t, backend, use_legacy, use_nullable_int=use_nullable_int
            )
        else:
            with pytest.raises(
                ValueError,
                match="Cannot convert masked integer columns to DataFrame without using nullable integers.",
            ):
                df = self._to_dataframe(
                    t, backend, use_legacy, use_nullable_int=use_nullable_int
                )
            return

        t2 = self._from_dataframe(df, backend, use_legacy)
        for name, column in t.columns.items():
            assert np.all(column.data == t2[name].data)
            if hasattr(t2[name], "mask"):
                assert np.all(column.mask == t2[name].mask)

            if column.dtype.kind == "i":
                if np.any(column.mask) and not use_nullable_int:
                    assert t2[name].dtype.kind == "f"
                else:
                    assert t2[name].dtype.kind == "i"
            else:
                if column.dtype.byteorder in ("=", "|"):
                    assert column.dtype == t2[name].dtype
                else:
                    assert column.dtype.newbyteorder() == t2[name].dtype

    def test_basic_roundtrip(self, backend, use_legacy):
        """Test basic round-trip conversion for different backends."""
        t = table.Table()
        t["a"] = [1, 2, 3]
        t["b"] = [4.0, 5.0, 6.0]
        t["c"] = ["x", "y", "z"]

        # Convert to DataFrame and back
        df = self._to_dataframe(t, backend, use_legacy)
        t2 = self._from_dataframe(df, backend, use_legacy)

        # Check that data is preserved
        assert_allclose(t["a"], t2["a"])
        assert_allclose(t["b"], t2["b"])
        assert np.all(t["c"] == t2["c"])

    def test_units_preservation(self, backend, use_legacy):
        """Test that units are handled correctly through DataFrame conversion."""
        t = table.QTable()
        t["x"] = [1, 2, 3] * u.m
        t["y"] = [4.0, 5.0, 6.0] * u.s

        # Test that units are lost in DataFrame conversion (expected behavior)
        df = self._to_dataframe(t, backend, use_legacy)
        t2 = self._from_dataframe(df, backend, use_legacy)

        # Original table should still have units
        assert t["x"].unit == u.m
        assert t["y"].unit == u.s

        # Round-trip table should not have units
        assert t2["x"].unit is None
        assert t2["y"].unit is None

        # But data should be preserved
        assert_allclose(t["x"].value, t2["x"])
        assert_allclose(t["y"].value, t2["y"])

    def test_masked_int_data(self, backend, use_legacy):
        """Test specific masked integer data handling."""
        data = {"data": [0, 1, 2]}
        t = table.Table(data=data, masked=True)
        t["data"].mask = [1, 1, 0]

        df = self._to_dataframe(t, backend, use_legacy)

        match backend:
            case "pandas":
                val = df["data"].iloc[2]
                nulls_first_two = df["data"].isnull().iloc[:2]
            case "polars":
                val = df["data"][2]
                nulls_first_two = df["data"].is_null()[:2]
            case "pyarrow":
                val = df["data"][2].as_py()
                nulls_first_two = df["data"].is_null()[:2].to_numpy()
            case _:
                raise ValueError(f"Unknown backend: {backend}")

        assert val == 2
        assert nulls_first_two.all()
