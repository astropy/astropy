"""
DataFrame conversion utilities for Astropy Tables.

This module provides utility functions for converting between Astropy Tables
and various DataFrame formats (pandas, polars, pyarrow, etc) via the narwhals
library.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, NewType

import numpy as np

from astropy.utils.compat.optional_deps import HAS_NARWHALS, HAS_PANDAS

from .column import Column, MaskedColumn

__all__ = ["from_df", "from_pandas", "to_df", "to_pandas"]

if TYPE_CHECKING:
    from narwhals.typing import EagerAllowed, IntoBackend
    from units.typing import UnitLike

    from .table import Table

# Sentinel value to indicate pandas-like backend validation
PandasLikeSentinel = NewType("PandasLikeSentinel", object)  # Custom type
PANDAS_LIKE: PandasLikeSentinel = PandasLikeSentinel(object())

# Fixed set of integer dtypes
INTEGER_DTYPE_KINDS = frozenset({"u", "i"})


def _numpy_to_pandas_dtype(dtype: np.dtype) -> str:
    """Convert a numpy dtype to a pandas dtype string, handling nullable integers."""
    dtype_name = dtype.name
    # Special case needed for uint -> UInt
    if dtype_name.startswith("uint"):
        return "UInt" + dtype_name.removeprefix("uint")
    return dtype.name.title()


def _encode_mixins(tbl: Table) -> Table:
    """Encode mixin columns to basic columns for DataFrame compatibility."""
    from astropy.time import TimeBase, TimeDelta

    from . import serialize
    from .column import col_copy

    time_cols = [col for col in tbl.itercols() if isinstance(col, TimeBase)]
    if time_cols:
        new_cols = []
        for col in tbl.itercols():
            new_col = col_copy(col, copy_indices=False) if col.info.indices else col
            new_cols.append(new_col)
        tbl = tbl.__class__(new_cols, copy=False)

        for col in tbl.itercols():
            col.info.indices.clear()

        for col in time_cols:
            if isinstance(col, TimeDelta):
                new_col = (col.sec * 1e9).astype("timedelta64[ns]")
                nat = np.timedelta64("NaT")
            else:
                new_col = col.datetime64.copy()
                nat = np.datetime64("NaT")
            if col.masked:
                new_col[col.mask] = nat
            tbl[col.info.name] = new_col

    # Convert the table to one with no mixins, only Column objects.
    encode_tbl = serialize.represent_mixins_as_columns(tbl)
    return encode_tbl


def _get_backend_impl(backend: str):
    """Get the narwhals backend implementation."""
    if not HAS_NARWHALS:
        raise ModuleNotFoundError(
            "The narwhals library is required for generic DataFrame conversion."
        )
    import narwhals as nw

    return nw.Implementation.from_backend(backend)


def _is_pandas_like(
    backend_impl: IntoBackend[EagerAllowed] | PandasLikeSentinel,
) -> bool:
    return backend_impl is PANDAS_LIKE or (
        not isinstance(backend_impl, str)
        and getattr(backend_impl, "is_pandas_like", lambda: False)()
    )


def _validate_columns_for_backend(
    table: Table,
    *,
    backend_impl: IntoBackend[EagerAllowed] | PandasLikeSentinel = PANDAS_LIKE,
) -> None:
    """Validate that table columns are compatible with the target backend.

    backend_impl may be PANDAS_LIKE to indicate pandas-like validation.

    Raises ValueError if there is a multidimensional column with an unsupported backend.
    """
    # Check for multidimensional columns
    badcols = [name for name, col in table.columns.items() if len(col.shape) > 1]

    if not badcols:
        return

    # Check if pandas-like or pyarrow-like
    if _is_pandas_like(backend_impl) or (
        not isinstance(backend_impl, str)
        and getattr(backend_impl, "is_pyarrow", lambda: False)()
    ):
        raise ValueError(
            f"Cannot convert a table with multidimensional columns to a "
            f"pandas-like or pyarrow DataFrame. Offending columns are: {badcols}\n"
            f"One can filter out such columns using:\n"
            f"names = [name for name in tbl.colnames if len(tbl[name].shape) <= 1]\n"
            f"tbl[names].to_pandas(...)"
        )


def _handle_index_argument(
    table: Table,
    *,
    index: bool | str | None,
    backend_impl: IntoBackend[EagerAllowed] | PandasLikeSentinel,
) -> bool | str:
    """Process the index argument for DataFrame conversion."""
    if index is False:
        return False

    # Check if pandas-like, PANDAS_LIKE is reserved for pandas itself
    # Non-pandas backends don't support indexing
    if not _is_pandas_like(backend_impl):
        if index:
            raise ValueError("Indexing is only supported for pandas-like backends.")
        return False

    has_single_pk = table.primary_key is not None and len(table.primary_key) == 1

    if index in (None, True):
        if has_single_pk:
            return table.primary_key[0]
        elif index is None:
            return False
        else:
            assert index is True
            raise ValueError("index=True requires a single-column primary key.")

    elif isinstance(index, str):
        if index not in table.colnames:
            raise ValueError(f"{index!r} is not in the table columns.")
        return index

    else:
        raise TypeError(
            "index must be None, False, True, or a valid column name. "
            f"Type provided: {type(index).__name__}"
        )


def to_df(
    table: Table,
    *,
    backend: str,
    index: bool | str | None = None,
    use_nullable_int: bool = True,
):
    """Convert an Astropy Table to a DataFrame using the specified backend."""
    if not HAS_NARWHALS:
        raise ModuleNotFoundError(
            "The narwhals library is required for the generic to_df method. "
            "If you want to only convert to pandas, use the `to_pandas` method instead."
        )
    import narwhals as nw

    backend_impl = _get_backend_impl(backend)

    # Using private API while narwhals-dev/narwhals#3150 is not merged
    if not nw._utils.is_eager_allowed(backend_impl):
        raise ValueError("Must export to eager compatible DataFrame")

    # Handle index argument
    index = _handle_index_argument(table, index=index, backend_impl=backend_impl)

    # Encode mixins and validate columns
    tbl = _encode_mixins(table)
    _validate_columns_for_backend(tbl, backend_impl=backend_impl)

    # Convert to narwhals DataFrame
    array = tbl.as_array()
    array_dict = {
        n: array[n].data if tbl.has_masked_columns else array[n] for n in tbl.colnames
    }
    df_nw = nw.from_dict(array_dict, backend=backend_impl)

    # Handle masked columns
    if tbl.has_masked_columns:
        masked_cols = [
            name
            for name, col in tbl.columns.items()
            if isinstance(col, MaskedColumn) and col.mask.any()
        ]
        old_dtypes = {n: df_nw[n].dtype for n in masked_cols}

        df_nw = df_nw.with_columns(
            (nw.when(~nw.new_series(n, array[n].mask, backend=backend_impl)).then(n))
            for n in masked_cols
        )

        for n, old_dtype in old_dtypes.items():
            if old_dtype.is_integer() and not use_nullable_int:
                raise ValueError(
                    "Cannot convert masked integer columns to DataFrame without using nullable integers. "
                    f"Set use_nullable_int=True or remove the offending column: {n}."
                )
            elif old_dtype.is_float():
                df_nw = df_nw.with_columns(nw.col(n).cast(old_dtype).alias(n))

    # Convert to dataframe
    df_native = df_nw.to_native()

    # Fix pandas-like nullable integers
    if backend_impl.is_pandas_like() and tbl.has_masked_columns and use_nullable_int:
        for name in masked_cols:
            dtype = array[name].dtype
            if (dtype := array[name].dtype).kind not in INTEGER_DTYPE_KINDS:
                continue

            df_native[name] = df_native[name].astype(_numpy_to_pandas_dtype(dtype))

    # Pandas-like index
    if index:
        df_native.set_index(index, inplace=True, drop=True)

    return df_native


def from_df(
    cls: type[Table],
    df: Any,
    *,
    index: bool = False,
    units: Mapping[str, UnitLike] | None = None,
) -> Table:
    """Create an instance of ``cls`` from any narwhals-compatible DataFrame."""
    if not HAS_NARWHALS:
        raise ModuleNotFoundError(
            "The narwhals library is required for the generic from_df method. "
            "If you want to only convert from pandas, use the `from_pandas` method instead."
        )
    import narwhals as nw

    # Create output
    out = {}

    # Handle pandas index
    if index:
        if not hasattr(df, "index"):
            raise ValueError(
                "The input dataframe does not have an index. "
                "Are you trying to convert a non-pandas dataframe? "
                "Set `index=False` (or leave it unspecified) to avoid this error."
            )
        index_name = str(df.index.name or "index")
        while index_name in df.columns:
            index_name = "_" + index_name + "_"
        df.reset_index(index_name, inplace=True, drop=False)

    # Narwhals layer, must convert to eager
    df_nw = nw.from_native(df)
    if isinstance(df_nw, nw.LazyFrame):
        df_nw = df_nw.collect()

    # Handle units
    if units is None:
        units = {}
    elif not isinstance(units, Mapping):
        raise TypeError(
            f"Expected a Mapping from column-names to units. Got {units!r} with type {type(units)}"
        )
    not_found = set(units.keys()) - set(df_nw.columns)
    if not_found:
        warnings.warn(f"`units` contains additional columns: {not_found}")

    # Iterate over Narwhals columns
    for column in df_nw.iter_columns():
        # Unpack relevant data
        name = column.name
        dtype = column.dtype
        nw_mask = column.is_null()
        unit = units.get(name)

        # Check first if the dtype is supported
        if isinstance(dtype, nw.Int128):
            raise ValueError(
                "Astropy Table does not support narwhals.Int128, please use a smaller integer type."
            )
        data = column.to_numpy()

        # Handle nullable integers
        if dtype.is_integer() and nw_mask.any():
            nullfilled = column.fill_null(0).to_numpy()
            out[name] = MaskedColumn(
                data=nullfilled, mask=nw_mask.to_numpy(), unit=unit, copy=False
            )
            continue

        # Handle datetime columns
        elif isinstance(dtype, nw.Datetime):
            from astropy.time import Time

            datetime = Time(data, format="datetime64")
            datetime.format = "isot"
            out[name] = datetime
            continue

        # Handle timedelta columns
        elif isinstance(dtype, nw.Duration):
            from astropy.time import TimeDelta

            duration = (
                data.astype("timedelta64[ns]").astype(np.float64, copy=False) / 1e9
            )
            out[name] = TimeDelta(duration, format="sec")
            if nw_mask.any():
                out[name][nw_mask.to_numpy()] = np.ma.masked
            continue

        # Handle string-like columns
        elif isinstance(dtype, (nw.String, nw.Binary, nw.Object)):
            if data.dtype.kind == "O":
                ts = (str, bytes)
                # If all elements of an object array are string-like or None or np.nan
                # then coerce back to a native numpy str/unicode array.
                if all((x is None) or isinstance(x, ts) or math.isnan(x) for x in data):
                    # Force any missing (null) values to b''.  Numpy will
                    # upcast to str/unicode as needed. We go via a list to
                    # avoid replacing objects in a view of the pandas array and
                    # to ensure numpy initializes to string or bytes correctly.
                    data = np.array(
                        [b"" if m else d for (d, m) in zip(data, nw_mask.to_numpy())]
                    )

        if nw_mask.any():
            out[name] = MaskedColumn(
                data=data, mask=nw_mask.to_numpy(), unit=unit, copy=False
            )
        else:
            out[name] = Column(data=data, unit=unit, copy=False)

    return cls(out)


def to_pandas(
    table: Table, *, index: bool | str | None = None, use_nullable_int: bool = True
):
    """Convert an Astropy Table to a pandas DataFrame.

    This mirrors the previous DataFrameConverter.to_pandas method but as a
    module-level function.
    """
    if not HAS_PANDAS:
        raise ModuleNotFoundError(
            "pandas is required for to_pandas conversion. "
            "Install with: pip install pandas"
        )
    from pandas import DataFrame, Series

    # Handle index argument (pandas-specific logic)
    index = _handle_index_argument(
        table, index=index, backend_impl=PANDAS_LIKE
    )  # PANDAS_LIKE for pandas validation

    # Encode mixins and validate columns
    tbl = _encode_mixins(table)
    _validate_columns_for_backend(tbl, backend_impl=PANDAS_LIKE)  # pandas validation

    out = {}

    for name, column in tbl.columns.items():
        if getattr(column.dtype, "isnative", True):
            out[name] = column
        else:
            out[name] = column.data.byteswap().view(column.dtype.newbyteorder("="))

        if isinstance(column, MaskedColumn) and np.any(column.mask):
            if column.dtype.kind in ["i", "u"]:
                pd_dtype = column.dtype.name
                if use_nullable_int:
                    # Convert int64 to Int64, uint32 to UInt32, etc for nullable types
                    pd_dtype = pd_dtype.replace("i", "I").replace("u", "U")
                else:
                    from pandas.errors import IntCastingNaNError

                    raise IntCastingNaNError(
                        "Cannot convert masked integer columns to DataFrame without using nullable integers. "
                        f"Set use_nullable_int=True or remove the offending column: {name}."
                    )
                out[name] = Series(out[name], dtype=pd_dtype)

            elif column.dtype.kind not in ["f", "c"]:
                out[name] = column.astype(object).filled(np.nan)

    kwargs = {}

    if index:
        idx = out.pop(index)
        kwargs["index"] = idx

        # We add the table index to Series inputs (MaskedColumn with int values) to override
        # its default RangeIndex, see #11432
        for v in out.values():
            if isinstance(v, Series):
                v.index = idx

    df = DataFrame(out, **kwargs)
    if index:
        # Explicitly set the pandas DataFrame index to the original table
        # index name.
        df.index.name = idx.info.name

    return df


def from_pandas(
    cls: type[Table],
    dataframe: Any,
    index: bool = False,
    units: Mapping[str, UnitLike] | None = None,
) -> Table:
    """Create an instance of ``cls`` from a pandas DataFrame."""
    out = {}

    names = list(dataframe.columns)
    columns = [dataframe[name] for name in names]
    datas = [np.array(column) for column in columns]
    masks = [np.array(column.isnull()) for column in columns]

    if index:
        index_name = dataframe.index.name or "index"
        while index_name in names:
            index_name = "_" + index_name + "_"
        names.insert(0, index_name)
        columns.insert(0, dataframe.index)
        datas.insert(0, np.array(dataframe.index))
        masks.insert(0, np.zeros(len(dataframe), dtype=bool))

    if units is None:
        units = [None] * len(names)
    else:
        if not isinstance(units, Mapping):
            raise TypeError(
                f"Expected a Mapping from column-names to units. Got {units!r} with type {type(units)}"
            )

        not_found = set(units.keys()) - set(names)
        if not_found:
            warnings.warn(f"`units` contains additional columns: {not_found}")

        units = [units.get(name) for name in names]

    for name, column, data, mask, unit in zip(names, columns, datas, masks, units):
        if column.dtype.kind in ["u", "i", "b"] and np.any(mask):
            # Special-case support for pandas nullable int and bool
            np_dtype = column.dtype.numpy_dtype
            data = np.zeros(shape=column.shape, dtype=np_dtype)
            data[~mask] = column[~mask]
            out[name] = MaskedColumn(
                data=data, name=name, mask=mask, unit=unit, copy=False
            )
            continue

        if data.dtype.kind == "O":
            # If all elements of an object array are string-like or np.nan
            # then coerce back to a native numpy str/unicode array.
            string_types = (str, bytes)
            nan = np.nan
            if all(isinstance(x, string_types) or x is nan for x in data):
                # Force any missing (null) values to b''.  Numpy will
                # upcast to str/unicode as needed. We go via a list to
                # avoid replacing objects in a view of the pandas array and
                # to ensure numpy initializes to string or bytes correctly.
                data = np.array([b"" if m else d for (d, m) in zip(data, mask)])

        # Numpy datetime64
        if data.dtype.kind == "M":
            from astropy.time import Time

            out[name] = Time(data, format="datetime64")
            if np.any(mask):
                out[name][mask] = np.ma.masked
            out[name].format = "isot"

        # Numpy timedelta64
        elif data.dtype.kind == "m":
            from astropy.time import TimeDelta

            data_sec = data.astype("timedelta64[ns]").astype(np.float64) / 1e9
            out[name] = TimeDelta(data_sec, format="sec")
            if np.any(mask):
                out[name][mask] = np.ma.masked

        else:
            if np.any(mask):
                out[name] = MaskedColumn(data=data, name=name, mask=mask, unit=unit)
            else:
                out[name] = Column(data=data, name=name, unit=unit)

    return cls(out)
