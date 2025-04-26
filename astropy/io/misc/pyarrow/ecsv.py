import io
import json
import os
import re
import warnings
from contextlib import ExitStack
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from astropy.io.ascii.core import InconsistentTableError
from astropy.io.ascii.ecsv import (
    DELIMITERS,
    ECSV_DATATYPES,
    InvalidEcsvDatatypeWarning,
    _check_dtype_is_str,
)
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
class EcsvColumn:
    name: str
    datatype: str
    dtype: str | None = None
    shape: tuple[int, ...] = field(default_factory=tuple)
    unit: str | None = None
    description: str | None = None
    format: str | None = None
    meta: dict | None = None
    subtype: str | None = None


def read_header(
    input_file: str | os.PathLike | io.BytesIO,
    encoding: str = "utf-8",
) -> tuple[int, list[dict], dict]:
    """
    READ: Initialize the header Column objects from the table ``lines``.
    """
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
    cols = [EcsvColumn(**col) for col in header["datatype"]]

    for col in cols:
        # Warn if col datatype is not a valid ECSV datatype, but allow reading for
        # back-compatibility with existing older files that have numpy datatypes
        # like datetime64 or object or python str, which are not in the ECSV standard.
        if col.datatype not in ECSV_DATATYPES:
            msg = (
                f"unexpected datatype {col.datatype!r} of column {col.name!r} "
                f"is not in allowed ECSV datatypes {ECSV_DATATYPES}. "
                "Using anyway as a numpy dtype but beware since unexpected "
                "results are possible."
            )
            warnings.warn(msg, category=InvalidEcsvDatatypeWarning)

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
    cols: list[EcsvColumn],
    delimiter: str = ",",
    encoding: str = "utf-8",
) -> tuple[list[tuple], dict]:
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
        # TODO remove this line when the upstream bug is fixed
        null_values=[""],
    )
    return data


def convert_column(col: EcsvColumn, data_in: "npt.NDArray") -> "npt.NDArray":
    """
    Convert the column data from str to the appropriate numpy dtype.
    """
    try:
        # 1-d or N-d object columns are serialized as JSON.
        if col.subtype == "object":
            data_out = process_1d_Nd_object_data(col, data_in)

        # Variable length arrays with shape (n, m, ..., *) for fixed
        # n, m, .. and variable in last axis. Masked values here are
        # not currently supported.
        elif col.shape and col.shape[-1] is None:
            data_out = process_variable_length_array_data(col, data_in)

        # Multidim columns with consistent shape (n, m, ...). These
        # might be masked.
        elif col.shape:
            data_out = process_fixed_shape_multidim_data(col, data_in)

        # Regular scalar value column
        else:
            if col.subtype:
                warnings.warn(
                    f"unexpected subtype {col.subtype!r} set for column "
                    f"{col.name!r}, using dtype={col.dtype!r} instead.",
                    category=InvalidEcsvDatatypeWarning,
                )
            data_out = data_in

        if data_out.shape[1:] != tuple(col.shape):
            raise ValueError("shape mismatch between value and column specifier")

    except json.JSONDecodeError:
        raise ValueError(
            f"column {col.name!r} failed to convert: column value is not valid JSON"
        )
    except Exception as exc:
        raise ValueError(f"column {col.name!r} failed to convert: {exc}")

    return data_out


def process_1d_Nd_object_data(col, data_in):
    _check_dtype_is_str(col)
    col_vals = [json.loads(val) for val in data_in]
    data_out = np.empty([len(col_vals)] + col.shape, dtype=object)
    data_out[...] = col_vals
    return data_out


def process_fixed_shape_multidim_data(col, data_in):
    _check_dtype_is_str(col)

    # Change empty (blank) values in original ECSV to something
    # like "[[null, null],[null,null]]" so subsequent JSON
    # decoding works. Delete `col.mask` so that later code in
    # core TableOutputter.__call__() that deals with col.mask
    # does not run (since handling is done here already).
    if hasattr(data_in, "mask"):
        str_vals = data_in.tolist()
        all_none_arr = np.full(shape=col.shape, fill_value=None, dtype=object)
        all_none_json = json.dumps(all_none_arr.tolist())
        for idx in np.nonzero(data_in.mask)[0]:
            str_vals[idx] = all_none_json
        data_in = np.array(str_vals, dtype=object)

    col_vals = [json.loads(val) for val in data_in]
    # Make a numpy object array of col_vals to look for None
    # (masked values)
    arr_vals = np.array(col_vals, dtype=object)
    mask = arr_vals == None
    if not np.any(mask):
        # No None's, just convert to required dtype
        data_out = data_in.astype(col.subtype)
    else:
        # Replace all the None with an appropriate fill value
        kind = np.dtype(col.subtype).kind
        arr_vals[mask] = {"U": "", "S": b""}.get(kind, 0)
        # Finally make a MaskedArray with the filled data + mask
        data_out = np.ma.array(arr_vals.astype(col.subtype), mask=mask)

    return data_out


def process_variable_length_array_data(col, data_in):
    _check_dtype_is_str(col)

    # Empty (blank) values in original ECSV are masked. Instead set the values
    # to "[]" indicating an empty list. This operation also unmasks the values.
    if hasattr(data_in, "mask"):
        for idx in np.nonzero(data_in.mask)[0]:
            data_in[idx] = "[]"

    # Remake as a 1-d object column of numpy ndarrays or
    # MaskedArray using the datatype specified in the ECSV file.
    col_vals = []
    for str_val in data_in:
        obj_val = json.loads(str_val)  # list or nested lists
        try:
            arr_val = np.array(obj_val, dtype=col.subtype)
        except TypeError:
            # obj_val has entries that are inconsistent with
            # dtype. For a valid ECSV file the only possibility
            # is None values (indicating missing values).
            data = np.array(obj_val, dtype=object)
            # Replace all the None with an appropriate fill value
            mask = data == None
            kind = np.dtype(col.subtype).kind
            data[mask] = {"U": "", "S": b""}.get(kind, 0)
            arr_val = np.ma.array(data.astype(col.subtype), mask=mask)

        col_vals.append(arr_val)

    col.shape = ()
    col.dtype = np.dtype(object)
    # np.array(col_vals_arr, dtype=object) fails ?? so this workaround:
    data_out = np.empty(len(col_vals), dtype=object)
    data_out[:] = col_vals
    return data_out


def read_ecsv(
    input_file: str | os.PathLike | io.BytesIO,
    encoding: str = "utf-8",
) -> tuple[Table, dict]:
    """
    READ: Read the ECSV file and return a Table object.
    """
    data_start, cols, table_meta, delimiter = read_header(input_file, encoding=encoding)
    data_raw = read_data(
        input_file,
        data_start,
        cols,
        delimiter=delimiter,
        encoding=encoding,
    )

    # Convert the column data to the appropriate numpy dtype
    data = {col.name: convert_column(col, data_raw[col.name]) for col in cols}

    # Create the Table object
    table = Table(data)

    # Add metadata to the table
    if table_meta:
        table.meta.update(table_meta)

    table = serialize._construct_mixins_from_columns(table)

    return table, data_raw
