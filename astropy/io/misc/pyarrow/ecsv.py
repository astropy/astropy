import io
import json
import os
import re
import warnings
from contextlib import ExitStack
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from astropy.table import Table, meta, serialize


def get_header_lines(
    input_file: str | os.PathLike | io.BytesIO,
    encoding="utf-8",
) -> tuple[int, list[str]]:
    comment = "# ".encode(encoding)
    comment_comment = "##".encode(encoding)
    lines = []

    with ExitStack() as stack:
        if isinstance(input_file, (str, os.PathLike)):
            tmp_input_file = stack.enter_context(open(input_file, "rb"))
        else:
            tmp_input_file = input_file

        for idx, line in enumerate(tmp_input_file):
            # Allow blank lines and skip them
            if not (line_strip := line.strip()) or line_strip.startswith(
                comment_comment
            ):
                continue
            if line_strip.startswith(comment):
                lines.append(line_strip[2:].decode(encoding))
            else:
                # Stop iterating on first failed comment match for a non-blank line
                break

    if isinstance(input_file, io.BytesIO):
        input_file.seek(0)

    return idx, lines


@dataclass
class ColumnAttrs:
    name: str
    datatype: str
    subtype: str | None = None
    dtype: str | None = None
    shape: tuple[int, ...] = field(default_factory=tuple)
    unit: str | None = None
    description: str | None = None
    format: str | None = None
    meta: dict | None = None


def read_header(
    input_file: str | os.PathLike | io.BytesIO,
    encoding: str = "utf-8",
) -> tuple[int, list[ColumnAttrs], dict]:
    """
    READ: Initialize the header Column objects from the table ``lines``.
    """
    from astropy.io.ascii.core import InconsistentTableError
    from astropy.io.ascii.ecsv import DELIMITERS, ECSV_DATATYPES

    # Extract non-blank comment (header) lines with comment character stripped
    data_start, header_lines = get_header_lines(input_file, encoding=encoding)

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

    for col in cols:
        # Warn if col datatype is not a valid ECSV datatype, but allow reading for
        # back-compatibility with existing older files that have numpy datatypes
        # like datetime64 or object or python str, which are not in the ECSV standard.
        if col.datatype not in ECSV_DATATYPES:
            msg = (
                f"unexpected datatype {col.datatype!r} of column {col.name!r} "
                f"is not in allowed ECSV datatypes {ECSV_DATATYPES}."
            )
            raise InconsistentTableError(msg)

        # Subtype is written like "int64[2,null]" and we want to split this
        # out to "int64" and [2, None].
        subtype = col.subtype
        if subtype and "[" in subtype:
            idx = subtype.index("[")
            col.subtype = subtype[:idx]
            col.shape = json.loads(subtype[idx:])

        col.dtype = col.datatype

        # Convert ECSV "string" to numpy "str"
        for attr in ("datatype", "dtype", "subtype"):
            if getattr(col, attr) == "string":
                setattr(col, attr, "str")

        # ECSV subtype of 'json' maps to numpy 'object' dtype
        if col.subtype == "json":
            col.subtype = "object"

    return data_start, cols, table_meta, delimiter


def read_data(
    input_file: str | os.PathLike | io.BytesIO,
    data_start: int,
    cols: list[ColumnAttrs],
    delimiter: str = ",",
    encoding: str = "utf-8",
) -> Table:
    """
    READ: Read the data from the table ``lines``.
    """
    # Read the data lines
    dtypes = {col.name: col.datatype for col in cols}
    data = Table.read(
        input_file,
        format="pyarrow.csv",
        delimiter=delimiter,
        encoding=encoding,
        dtypes=dtypes,
        header_start=data_start,
    )
    return data


def _get_str_vals(col, data):
    if col.dtype != "str":
        raise ValueError(f'datatype of column {data.name!r} must be "string"')
    # For masked we need a list because for multidim the data under the mask is set
    # to a compatible value.
    if hasattr(data, "mask"):
        str_vals = data.view(np.ndarray).tolist()
        mask = data.mask
    else:
        str_vals = data
        mask = None
    return str_vals, mask


def convert_column(col_attrs: ColumnAttrs, data_in: "npt.NDArray") -> "npt.NDArray":
    """
    Convert the column data from str to the appropriate numpy dtype.
    """
    from astropy.io.ascii.ecsv import InvalidEcsvDatatypeWarning

    try:
        if col_attrs.subtype == "object" or col_attrs.shape:
            # Handle three distinct column types where each row element is serialized
            # to JSON. First convert the input data array or masked array to a list of
            # str and get the mask where available.
            str_vals, mask = _get_str_vals(col_attrs, data_in)

            if col_attrs.subtype == "object":
                # Any Python objects serializable to JSON
                process_func = process_1d_Nd_object_data
            elif col_attrs.shape[-1] is None:
                # Variable length arrays with shape (n, m, ..., *) for fixed
                # n, m, .. and variable in last axis.
                process_func = process_variable_length_array_data
            else:
                # Multidim columns with consistent shape (n, m, ...).
                process_func = process_fixed_shape_multidim_data

            data_out = process_func(col_attrs, str_vals, mask)

        # Regular scalar value column
        else:
            if col_attrs.subtype:
                warnings.warn(
                    f"unexpected subtype {col_attrs.subtype!r} set for column "
                    f"{col_attrs.name!r}, using dtype={col_attrs.dtype!r} instead.",
                    category=InvalidEcsvDatatypeWarning,
                )
            data_out = data_in

        if data_out.shape[1:] != tuple(col_attrs.shape):
            raise ValueError("shape mismatch between value and column specifier")

    except json.JSONDecodeError:
        raise ValueError(
            f"column {col_attrs.name!r} failed to convert: column value is not valid JSON"
        )
    except Exception as exc:
        raise ValueError(f"column {col_attrs.name!r} failed to convert: {exc}")

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
        kind = np.dtype(col_attrs.subtype).kind
        arr_vals[arr_vals_mask] = {"U": "", "S": b""}.get(kind, 0)
        # Finally make a MaskedArray with the filled data + mask
        data_out = np.ma.array(arr_vals.astype(col_attrs.subtype), mask=arr_vals_mask)
    else:
        data_out = arr_vals.astype(col_attrs.subtype)

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
            arr_val = np.array(obj_val, dtype=col_attrs.subtype)
        except TypeError:
            # obj_val has entries that are inconsistent with
            # dtype. For a valid ECSV file the only possibility
            # is None values (indicating missing values).
            vals = np.array(obj_val, dtype=object)
            # Replace all the None with an appropriate fill value
            mask_vals = vals == None
            kind = np.dtype(col_attrs.subtype).kind
            vals[mask_vals] = {"U": "", "S": b""}.get(kind, 0)
            arr_val = np.ma.array(vals.astype(col_attrs.subtype), mask=mask_vals)

        col_vals.append(arr_val)

    col_attrs.shape = ()
    col_attrs.dtype = np.dtype(object)
    np_empty = np.empty if mask is None else np.ma.empty
    data_out = np_empty(len(col_vals), dtype=object)
    data_out[:] = col_vals
    if mask is not None:
        data_out.mask = mask
    return data_out


def read_ecsv(
    input_file: str | os.PathLike | io.BytesIO,
    encoding: str = "utf-8",
) -> Table:
    """
    READ: Read the ECSV file and return a Table object.
    """
    from astropy.io.ascii.core import InconsistentTableError

    # For testing
    if isinstance(input_file, io.StringIO):
        input_file = io.BytesIO(input_file.getvalue().encode(encoding))
    elif isinstance(input_file, str) and "\n" in input_file:
        input_file = io.BytesIO(input_file.encode(encoding))
    elif isinstance(input_file, (list, tuple)):
        input_file = io.BytesIO("\n".join(input_file).encode(encoding))

    data_start, cols_attrs, table_meta, delimiter = read_header(
        input_file, encoding=encoding
    )
    data_raw = read_data(
        input_file,
        data_start,
        cols_attrs,
        delimiter=delimiter,
        encoding=encoding,
    )

    ecsv_header_names = [col_attrs.name for col_attrs in cols_attrs]
    if ecsv_header_names != data_raw.colnames:
        raise InconsistentTableError(
            f"column names from ECSV header {ecsv_header_names} do not "
            f"match names from header line of CSV data {data_raw.colnames}"
        )

    # Convert the column data to the appropriate numpy dtype
    data = {
        col_attrs.name: convert_column(col_attrs, data_raw[col_attrs.name])
        for col_attrs in cols_attrs
    }

    # Create the Table object
    table = Table(data)

    for col_attrs in cols_attrs:
        col = table[col_attrs.name]
        for attr in ["unit", "description", "format", "meta"]:
            if (val := getattr(col_attrs, attr)) is not None:
                setattr(col.info, attr, val)

    # Add metadata to the table
    if table_meta:
        table.meta.update(table_meta)

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

    io_registry.register_reader("pyarrow.ecsv", Table, read_ecsv)
    io_registry.register_writer("pyarrow.ecsv", Table, write_ecsv)
