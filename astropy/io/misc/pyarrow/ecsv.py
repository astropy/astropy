import collections
import functools
import io
import json
import os
import re
import warnings
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import numpy.typing as npt

from astropy.table import SerializedColumn, Table, meta, serialize

__all__ = ["read_ecsv", "register_pyarrow_ecsv_table", "write_ecsv"]


class ECSVEngine:
    """Base class for ECSV reader engines."""

    engines = {}

    def __init_subclass__(cls, **kwargs):
        """Register the subclass as an ECSV engine."""
        super().__init_subclass__(**kwargs)
        cls.engines[cls.name] = cls


class ECSVEnginePyArrow(ECSVEngine):
    """ECSV reader engine using PyArrow."""

    name = "pyarrow.csv"

    @staticmethod
    def convert_parsetype(parsetype: str) -> str:
        """Convert the parsetype to a converter type."""
        # PyArrow does not support float16, so we need to convert it to float32.
        return "float32" if parsetype == "float16" else parsetype


class ECSVEngineIoAscii(ECSVEngine):
    """ECSV reader engine using io.ascii."""

    name = "ascii.csv"

    @staticmethod
    def convert_parsetype(parsetype: str) -> str:
        """Convert the parsetype to a numpy dtype."""
        # Output compatible with io.ascii `converters` option.
        return np.dtype(parsetype).type


class ECSVEnginePandas(ECSVEngine):
    """ECSV reader engine using pandas."""

    name = "pandas.csv"

    @staticmethod
    def convert_parsetype(parsetype: str) -> np.dtype:
        """Convert the parsetype to a pandas dtype."""
        import pandas as pd

        dtype = np.dtype(parsetype)
        if dtype.kind in ("i", "u"):
            # Convert int64 to Int64, uint32 to UInt32, etc for nullable types
            converter = dtype.name.replace("i", "I").replace("u", "U")
        elif dtype.kind == "b":
            converter = "boolean"
        else:
            converter = parsetype
        return pd.api.types.pandas_dtype(converter)


def is_numpy_dtype(dtype: str) -> bool:
    """Check if the given dtype is a valid numpy dtype."""
    try:
        np.dtype(dtype)
    except Exception:
        return False
    else:
        return True


def get_header_lines(
    input_file: str | os.PathLike | io.BytesIO,
    encoding="utf-8",
) -> tuple[int, list[str], int]:
    """
    Extract header lines from a file or file-like object.

    This function reads a file or file-like object and extracts lines that
    start with a specific header prefix ("# ") while skipping blank lines
    and lines starting with a comment prefix ("##"). The function stops
    reading at the first non-blank, non-comment line that does not match
    the header prefix.

    Parameters
    ----------
    input_file : str | os.PathLike | io.BytesIO
        The input file path or file-like object to read. If a file path is
        provided, the function automatically handles compressed files with
        extensions `.gz` or `.bz2`.
    encoding : str, optional
        The encoding used to decode the file content. Default is "utf-8".

    Returns
    -------
    idx : int
        Index of the last line read.
    lines : list[str]
      List of decoded header lines without the header prefix.
    """
    header_prefix = "# ".encode(encoding)
    comment_prefix = "##".encode(encoding)
    lines = []
    n_empty = 0
    n_comment = 0

    with ExitStack() as stack:
        # Get a file-like object as tmp_input_file
        if isinstance(input_file, (str, os.PathLike)):
            ext = Path(input_file).suffix
            if ext == ".gz":
                import gzip

                opener = gzip.open
            elif ext == ".bz2":
                import bz2

                opener = bz2.open
            else:
                opener = open
            tmp_input_file = stack.enter_context(opener(input_file, "rb"))
        else:
            tmp_input_file = input_file

        for idx, line in enumerate(tmp_input_file):
            line_strip = line.strip()
            if line_strip.startswith(header_prefix):
                lines.append(line_strip[2:].decode(encoding))
            elif not line_strip:
                n_empty += 1
            elif line_strip.startswith(comment_prefix):
                n_comment += 1
            else:
                # Stop iterating on first failed comment match for a non-blank line
                break

    # Need to rewind the input file if it is a file-like object
    if isinstance(input_file, io.BytesIO):
        input_file.seek(0)

    return lines, idx, n_empty, n_comment


@dataclass
class ColumnAttrs:
    name: str
    datatype: str
    subtype: str | None = None
    unit: str | None = None
    description: str | None = None
    format: str | None = None
    meta: dict | None = None

    @functools.cached_property
    def parsetype(self) -> str:
        """Type used to parse the CSV file data"""
        return self.parsetype_dtype_shape[0]

    @functools.cached_property
    def dtype(self) -> str:
        """Numpy dtype in the final column data"""
        return self.parsetype_dtype_shape[1]

    @functools.cached_property
    def shape(self) -> tuple[int, ...]:
        """Shape of the column data"""
        return self.parsetype_dtype_shape[2]

    @functools.cached_property
    def parsetype_dtype_shape(self) -> tuple[str, str, tuple[int, ...]]:
        """Get the parsetype, dtype, and shape of the column from ECSV header."""
        return get_parsetype_dtype_shape(self.datatype, self.subtype, self.name)


def get_parsetype_dtype_shape(
    datatype: str, subtype: str | None, name: str
) -> tuple[str, str, tuple[int, ...]]:
    """Get the parsetype, dtype, and shape of the column from datatype and subtype.

    This function implements most of the complexity of the ECSV data type
    handling. It converts the ECSV ``datatype`` and ``subtype`` to the following:
    - ``parsetype``: the type used to parse the CSV file data. This is a string
      representation of the numpy dtype, e.g. "int64", "float32", "str", "object".
    - ``dtype``: the numpy dtype in the final column data as a string.
    - ``shape``: the shape of the column data (tuple[int]).

    The ECSV standard allows for a wide variety of data types and subtypes, and
    we also need to handle some legacy cases and be liberal in what we accept.
    """
    from astropy.io.ascii.core import InconsistentTableError
    from astropy.io.ascii.ecsv import ECSV_DATATYPES, InvalidEcsvDatatypeWarning

    parsetype = "str" if datatype == "string" else datatype
    dtype = None
    shape = ()

    if datatype not in ECSV_DATATYPES:
        msg = (
            f"unexpected datatype {datatype!r} of column {name!r} "
            f"is not in allowed ECSV datatypes {ECSV_DATATYPES}."
        )
        # Try being liberal on input if the `parsetype` (derived from ECSV
        # `datatype`) looks like a numpy dtype. In this case, parse the column as
        # string and then cast as `parsetype`. This allows for back-compatibility
        # with early versions of io.ascii.ecsv that wrote and read e.g.
        # datatype=datetime64.
        if is_numpy_dtype(parsetype):
            dtype = parsetype
            parsetype = "str"
            warnings.warn(msg, InvalidEcsvDatatypeWarning)
        else:
            # No joy, this is an exception
            raise InconsistentTableError(msg)

    if subtype and parsetype != "str":
        # Note: the "column .. failed to convert" bit is odd here but it is to match
        # the io.ascii.ecsv behavior.
        raise ValueError(
            f"column {name!r} failed to convert: "
            f'datatype of column {name!r} must be "string"'
        )

    if subtype:
        # Subtype can be written like "int64[2,null]" and we want to split this
        # out to "int64" and [2, None].
        if "[" in subtype:
            idx = subtype.index("[")
            dtype = subtype[:idx]
            shape = tuple(json.loads(subtype[idx:]))
        else:
            dtype = subtype

        # Map ECSV types to numpy dtypes
        dtype = {"json": "object", "string": "str"}.get(dtype, dtype)

        # Check if the subtype corresponds to a valid numpy dtype. This is required by
        # the astropy implementation, but not by the ECSV standard. The standard states
        # that an unknown subtype can be ignored, so that is what we do here (but with
        # a warning).
        if not is_numpy_dtype(dtype):
            warnings.warn(
                f"unexpected subtype {subtype!r} set for column "
                f"{name!r}, using dtype={parsetype!r} instead.",
                category=InvalidEcsvDatatypeWarning,
            )
            dtype = parsetype

    if parsetype == "float16":
        # PyArrow does not support float16, so we need to parse it as float32 and then
        # cast it to float16 at the end.
        parsetype = "float32"
        dtype = "float16"

    return parsetype, dtype, shape


def read_header(
    input_file: str | os.PathLike | io.BytesIO,
    encoding: str = "utf-8",
    engine: Literal["pyarrow.csv", "ascii.csv"] = "pyarrow.csv",
) -> tuple[int, list[ColumnAttrs], dict]:
    """
    READ: Initialize the header Column objects from the table ``lines``.
    """
    from astropy.io.ascii.core import InconsistentTableError
    from astropy.io.ascii.ecsv import DELIMITERS

    # Extract non-blank comment (header) lines with comment character stripped
    header_lines, n_header, n_empty, n_comment = get_header_lines(
        input_file, encoding=encoding
    )

    # Validate that this is a ECSV file
    ecsv_header_re = r"""%ECSV [ ]
                            (?P<major> \d+)
                            \. (?P<minor> \d+)
                            \.? (?P<bugfix> \d+)? $"""

    no_header_msg = (
        'ECSV header line like "# %ECSV <version>" not found as first line.'
        "  This is required for a ECSV file."
    )

    if not header_lines:
        raise InconsistentTableError(no_header_msg)

    match = re.match(ecsv_header_re, header_lines[0].strip(), re.VERBOSE)
    if not match:
        raise InconsistentTableError(no_header_msg)

    try:
        header = meta.get_header_from_yaml(header_lines)
    except meta.YamlParseError as e:
        raise InconsistentTableError("unable to parse yaml in meta header") from e

    table_meta = header.get("meta", None)

    delimiter = header.get("delimiter", " ")
    if delimiter not in DELIMITERS:
        raise ValueError(
            "only space and comma are allowed for delimiter in ECSV format"
        )

    # Create list of columns from `header`.
    cols = [ColumnAttrs(**col) for col in header["datatype"]]

    # Start line in the file where the CSV table starts (i.e. line with column names)
    if engine == "ascii.csv":
        csv_table_start = n_header - n_empty
    else:
        csv_table_start = n_header

    return csv_table_start, cols, table_meta, delimiter


def read_data(
    input_file: str | os.PathLike | io.BytesIO,
    csv_table_start: int,
    cols: list[ColumnAttrs],
    table_meta: dict,
    null_values: list[str],
    delimiter: str = ",",
    encoding: str = "utf-8",
    engine: Literal["pyarrow.csv", "ascii.csv", "pandas.csv"] = "pyarrow.csv",
) -> Table:
    """Read the data from the table ``lines``."""
    engine_cls = ECSVEngine.engines[engine]
    converters = {col.name: engine_cls.convert_parsetype(col.parsetype) for col in cols}

    kw = {}
    if engine == "ascii.csv":
        kw["fill_values"] = get_null_values_per_column(cols, table_meta, null_values)
        kw["header_start"] = csv_table_start
        kw["converters"] = converters
        # Fast reader does not support converters (defining types in advance) nor any
        # encoding. Converters are required, e.g. for a string column that looks like
        # floats. Would be nice to fix this, but in mean time use Python CSV reader.
        kw["fast_reader"] = False
    elif engine == "pyarrow.csv":
        kw["null_values"] = null_values
        kw["header_start"] = csv_table_start
        kw["dtypes"] = converters
    elif engine == "pandas.csv":
        fill_values = get_null_values_per_column(cols, table_meta, null_values)
        null_values = collections.defaultdict(list)
        for null_value, _, col_name in fill_values:
            null_values[col_name].append(null_value)

        kw["na_values"] = null_values
        kw["keep_default_na"] = False
        kw["skiprows"] = csv_table_start
        kw["dtype"] = converters
    else:
        raise ValueError(
            f"engine must be 'pyarrow.csv', 'ascii.csv' or 'pandas.csv', not {engine!r}"
        )

    data = Table.read(
        input_file,
        format=engine,
        delimiter=delimiter,
        encoding=encoding,
        **kw,
    )
    return data


def get_str_vals(data):
    """Convert the data to a list of str and get the mask if available."""
    # For masked we need a list because for multidim the data under the mask is set
    # to a compatible value.
    if hasattr(data, "mask"):
        str_vals = data.view(np.ndarray).tolist()
        mask = data.mask
    else:
        str_vals = data
        mask = None
    return str_vals, mask


def convert_column(col: ColumnAttrs, data_in: "npt.NDArray") -> "npt.NDArray":
    """
    Convert the column data from original parsetype to the output numpy dtype.
    """
    try:
        if col.dtype == "object" or col.shape:
            # Handle three distinct column types where each row element is serialized
            # to JSON. First convert the input data array or masked array to a list of
            # str and get the mask where available.
            str_vals, mask = get_str_vals(data_in)

            if col.dtype == "object":
                # Any Python objects serializable to JSON
                process_func = process_1d_Nd_object_data
            elif col.shape[-1] is None:
                # Variable length arrays with shape (n, m, ..., *) for fixed
                # n, m, .. and variable in last axis.
                process_func = process_variable_length_array_data
            else:
                # Multidim columns with consistent shape (n, m, ...).
                process_func = process_fixed_shape_multidim_data

            data_out = process_func(col, str_vals, mask)

        # Regular scalar value column
        else:
            data_out = data_in
            # If we need to cast the data to a different dtype, do it now.
            if col.dtype and data_out.dtype != np.dtype(col.dtype):
                data_out = data_out.astype(col.dtype)

        if data_out.shape[1:] != tuple(col.shape):
            raise ValueError("shape mismatch between value and column specifier")

    except json.JSONDecodeError:
        raise ValueError(
            f"column {col.name!r} failed to convert: column value is not valid JSON"
        )
    except Exception as exc:
        raise ValueError(f"column {col.name!r} failed to convert: {exc}")

    return data_out


def process_1d_Nd_object_data(col_attrs, str_vals, mask):
    if mask is not None:
        for idx in np.nonzero(mask)[0]:
            str_vals[idx] = "0"  # could be "null" but io.ascii uses "0"
    col_vals = [json.loads(val) for val in str_vals]
    np_empty = np.empty if mask is None else np.ma.empty
    data_out = np_empty((len(col_vals),) + tuple(col_attrs.shape), dtype=object)
    data_out[...] = col_vals
    if mask is not None:
        data_out.mask = mask
    return data_out


def process_fixed_shape_multidim_data(col_attrs, str_vals, mask):
    # Change empty (blank) values in original ECSV to something
    # like "[[null, null],[null,null]]" so subsequent JSON
    # decoding works.
    if mask is not None:
        all_none_arr = np.full(shape=col_attrs.shape, fill_value=None, dtype=object)
        fill_value = json.dumps(all_none_arr.tolist())
        for idx in np.nonzero(mask)[0]:
            str_vals[idx] = fill_value

    col_vals = [json.loads(val) for val in str_vals]

    # Make a numpy object array of col_vals to look for None (masked values)
    arr_vals = np.array(col_vals, dtype=object)
    arr_vals_mask = arr_vals == None
    if np.any(arr_vals_mask):
        # Replace all the None with an appropriate fill value
        kind = np.dtype(col_attrs.dtype).kind
        arr_vals[arr_vals_mask] = {"U": "", "S": b""}.get(kind, 0)
        # Finally make a MaskedArray with the filled data + mask
        data_out = np.ma.array(arr_vals.astype(col_attrs.dtype), mask=arr_vals_mask)
    else:
        data_out = arr_vals.astype(col_attrs.dtype)

    return data_out


def process_variable_length_array_data(col_attrs, str_vals, mask):
    """Variable length arrays with shape (n, m, ..., *)

    Shape is fixed for n, m, .. and variable in last axis. The output is a 1-d object
    array with each row element being an ``np.ndarray`` or ``np.ma.masked_array`` of the
    appropriate shape.
    """
    # Empty (blank) values in original ECSV are masked. Instead set the values
    # to "[]" indicating an empty list. This operation also unmasks the values.
    if mask is not None:
        fill_value = "[]"
        for idx in np.nonzero(mask)[0]:
            str_vals[idx] = fill_value

    # Remake as a 1-d object column of numpy ndarrays or
    # MaskedArray using the datatype specified in the ECSV file.
    col_vals = []
    for str_val in str_vals:
        obj_val = json.loads(str_val)  # list or nested lists
        try:
            arr_val = np.array(obj_val, dtype=col_attrs.dtype)
        except TypeError:
            # obj_val has entries that are inconsistent with
            # dtype. For a valid ECSV file the only possibility
            # is None values (indicating missing values).
            vals = np.array(obj_val, dtype=object)
            # Replace all the None with an appropriate fill value
            mask_vals = vals == None
            kind = np.dtype(col_attrs.dtype).kind
            vals[mask_vals] = {"U": "", "S": b""}.get(kind, 0)
            arr_val = np.ma.array(vals.astype(col_attrs.dtype), mask=mask_vals)

        col_vals.append(arr_val)

    col_attrs.shape = ()
    col_attrs.dtype = np.dtype(object)
    np_empty = np.empty if mask is None else np.ma.empty
    data_out = np_empty(len(col_vals), dtype=object)
    data_out[:] = col_vals
    if mask is not None:
        data_out.mask = mask
    return data_out


def get_null_values_per_column(
    cols: list[ColumnAttrs],
    table_meta: dict | None,
    null_values: list[str],
) -> list[tuple[str, str, str]]:
    """Get fill values for individual columns for the ascii.csv engine.

    For ECSV handle the corner case of data that has been serialized using
    the serialize_method='data_mask' option, which writes the full data and
    mask directly, AND where that table includes a string column with zero-length
    string entries ("") which are valid data.

    Normally the super() method will set col.fill_value=('', '0') to replace
    blanks with a '0'.  But for that corner case subset, instead do not do
    any filling.

    Returns
    -------
    fill_values : list[tuple[str, str, str]]
        A list of tuples with (null_value, fill_value, column_name) for each
        column in `cols`. If no fill values are needed, returns an empty list.
    """
    if table_meta is None:
        table_meta = {}

    # Get the serialized columns spec or an empty dict if not present.
    serialized_columns: dict[str, SerializedColumn] = table_meta.get(
        "__serialized_columns__", {}
    )

    # A serialized MaskedColumn column (via `serialize_method="data_mask"`) does not
    # have a fill value, so assemble a set of columns names to skip include the data and
    # the mask columns. For example:
    # - __serialized_columns__:
    #     a:
    #       __class__: astropy.table.column.MaskedColumn
    #       data: !astropy.table.SerializedColumn {name: a}
    #       mask: !astropy.table.SerializedColumn {name: a.mask}
    masked_col_names = set()
    for name, sc in serialized_columns.items():
        if sc["__class__"] == "astropy.table.column.MaskedColumn":
            masked_col_names.add(name)
            masked_col_names.add(name + ".mask")

    fill_values = []
    for col in cols:
        if col.name in masked_col_names:
            continue

        fill_value = "" if col.parsetype == "str" else "0"
        for null_value in null_values:
            fill_values.append((null_value, fill_value, col.name))

        # FIXME: only applies for pandas.csv engine
        # if col.parsetype.startswith("float"):
        #     for nan in ("nan", "NaN"):
        #         fill_values.append((nan, fill_value, col.name))

    return fill_values


def read_ecsv(
    input_file: str | os.PathLike | io.BytesIO,
    encoding: str = "utf-8",
    engine: Literal["pyarrow.csv", "ascii.csv"] = "pyarrow.csv",
    null_values: list[str] | None = None,
) -> Table:
    """
    READ: Read the ECSV file and return a Table object.
    """
    from astropy.io.ascii.core import InconsistentTableError

    if null_values is None:
        null_values = [""]

    if isinstance(input_file, io.StringIO):
        input_file = io.BytesIO(input_file.getvalue().encode(encoding))
    elif isinstance(input_file, str) and "\n" in input_file:
        input_file = io.BytesIO(input_file.encode(encoding))
    elif isinstance(input_file, (list, tuple)):
        # TODO: better way to check for an iterable of str?
        input_file = io.BytesIO("\n".join(input_file).encode(encoding))

    csv_table_start, cols_attrs, table_meta, delimiter = read_header(
        input_file, encoding=encoding, engine=engine
    )

    # Read the data from the CSV file starting at the line after the header. This
    # includes handling that is particular to the engine.
    data_raw = read_data(
        input_file,
        csv_table_start,
        cols_attrs,
        table_meta,
        null_values=null_values,
        delimiter=delimiter,
        encoding=encoding,
        engine=engine,
    )

    # Ensure ECSV header names match the data column names.
    ecsv_header_names = [col_attrs.name for col_attrs in cols_attrs]
    if ecsv_header_names != data_raw.colnames:
        raise InconsistentTableError(
            f"column names from ECSV header {ecsv_header_names} do not "
            f"match names from header line of CSV data {data_raw.colnames}"
        )

    # Convert the column data to the appropriate numpy dtype. This is mostly concerned
    # with JSON-encoded data but also handles cases like pyarrow not supporting float16.
    data = {
        col_attrs.name: convert_column(col_attrs, data_raw[col_attrs.name])
        for col_attrs in cols_attrs
    }

    # Create the Table object
    table = Table(data)

    # Transfer metadata from the ECSV header to the Table columns.
    for col_attrs in cols_attrs:
        col = table[col_attrs.name]
        for attr in ["unit", "description", "format", "meta"]:
            if (val := getattr(col_attrs, attr)) is not None:
                setattr(col.info, attr, val)

    # Add metadata to the table
    if table_meta:
        table.meta.update(table_meta)

    # Construct any mixin columns from the raw columns.
    table = serialize._construct_mixins_from_columns(table)

    return table


def write_ecsv(tbl, output, **kwargs):
    tbl.write(output, format="ascii.ecsv", **kwargs)


def register_pyarrow_ecsv_table():
    """
    Register pyarrow.csv with Unified I/O as a Table reader.
    """
    from astropy.io import registry as io_registry
    from astropy.table import Table

    io_registry.register_reader("ecsv", Table, read_ecsv)
    io_registry.register_writer("ecsv", Table, write_ecsv)
