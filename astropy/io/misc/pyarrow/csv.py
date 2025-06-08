# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides functionality to read CSV files into Astropy Tables using PyArrow.
"""

import datetime
import io
import os
from contextlib import ExitStack
from typing import TYPE_CHECKING, BinaryIO, Literal, Union

import numpy as np
import numpy.typing

from astropy.utils.compat.optional_deps import HAS_PYARROW

if TYPE_CHECKING:
    import numpy.typing as npt
    import pyarrow as pa
    import pyarrow.csv

    from astropy.table import Table

__all__ = ["convert_pa_table_to_astropy_table", "read_csv"]


def read_csv(
    input_file: os.PathLike[str] | str | BinaryIO,
    *,
    delimiter: str = ",",
    quotechar: str | Literal[False] = '"',
    doublequote: bool = True,
    escapechar: str | bool = False,
    header_start: int | None = 0,
    data_start: int | None = None,
    names: list[str] | None = None,
    include_names: list[str] | None = None,
    dtypes: dict[str, numpy.typing.DTypeLike] | None = None,
    comment: str | None = None,
    null_values: list[str] | None = None,
    encoding: str = "utf-8",
    newlines_in_values: bool = False,
    timestamp_parsers: list[str] | None = None,
) -> "Table":
    """Read a CSV file into an astropy Table using PyArrow.

    This function allows highly performant reading of text CSV files into an astropy
    ``Table`` using `PyArrow <https://arrow.apache.org/docs/python/csv.html>`_. The
    best performance is achieved for files with only numeric data types, but even for
    files with mixed data types, the performance is still better than the the standard
    ``astropy.io.ascii`` fast CSV reader.

    By default, empty values (zero-length string "") in the CSV file are read as masked
    values in the Table. This can be changed by using the ``null_values`` parameter to
    specify a list of strings to interpret as null (masked) values.

    Entirely empty lines in the CSV file are ignored.

    Columns consisting of only string values ``True`` and ``False`` are parsed as
    boolean data.

    Columns with ISO 8601 date/time strings are parsed as shown below:
    - ``12:13:14.123456``: ``object[datetime.time]``
    - ``2025-01-01``: ``np.datetime64[D]``
    - ``2025-01-01T01:02:03``: ``np.datetime64[s]``
    - ``2025-01-01T01:02:03.123456``: ``np.datetime64[ns]``

    Support for ignoring comment lines in the CSV file is provided by the ``comment``
    parameter. If this is set to a string, any line starting with optional whitespace
    and then this string is ignored. This is done by reading the entire file and
    scanning for comment lines. If the comment lines are all at the beginning of the
    file and both ``header_start`` and ``data_start`` are not specified, then the file
    is read efficiently by setting ``header_start`` to the first line after the
    comments. Otherwise the entire file is read into memory and the comment lines are
    removed before passing to the PyArrow CSV reader. Any values of ``header_start`` and
    ``data_start`` apply to the lines counts *after* the comment lines have been
    removed.

    Parameters
    ----------
    input_file : str, path-like, or file-like
        File path or binary file-like object to read from.
    delimiter : 1-character str, optional (default ",")
        Character delimiting individual cells in the CSV data.
    quotechar : 1-character str or False, optional (default '"')
        Character used optionally for quoting CSV values (`False` if quoting is not
        allowed).
    doublequote : bool, optional (default `True`)
        Whether two quotes in a quoted CSV value denote a single quote in the data.
    escapechar : 1-character str or `False`, optional (default `False`)
        Character used optionally for escaping special characters (`False` if escaping
        is not allowed).
    header_start : int, None, optional (default 0)
        Line index for the header line with column names. If `None` this implies that
        there is no header line and the column names are taken from ``names`` or
        generated automatically ("f0", "f1", ...).
    data_start : int, None, optional (default None)
        Line index for the start of data. If `None`, then data starts one line after the
        header line, or on the first line if there is no header.
    names : list, None, optional (default None)
        List of names for input data columns when there is no header line. If supplied,
        then ``header_start`` must be `None`.
    include_names : list, None, optional (default None)
        List of column names to include in output. If `None`, all columns are included.
    dtypes : dict[str, Any], None, optional (default None)
        If provided, this is a dictionary of data types for output columns. Each key is
        a column name and the value is either a PyArrow data type or a data type
        specifier that is accepted as an argument to `numpy.dtype`. Examples include
        ``pyarrow.Int32()``, ``pyarrow.time32("s")``, ``int``, ``np.float32``,
        ``np.dtype('f4')`` or ``"float32"``. Default is to infer the data types.
    comment : 1-character str or None, optional (default None)
        Character used to indicate the start of a comment. Any line starting with
        optional whitespace and then this character is ignored. Using this option will
        cause the parser to be slower and potentially use more memory as it uses Python
        code to strip comments.

    Other Parameters
    ----------------
    null_values : list, optional (default None)
        List of strings to interpret as null values. By default, only empty strings are
        considered as null values (equivalent to ``null_values=[""]``). Set to ``[]`` to
        disable null value handling.
    encoding: str, optional (default 'utf-8')
        Encoding of the input data.
    newlines_in_values: bool, optional (default False)
        Whether newline characters are allowed in CSV values. Setting this to True
        reduces the performance of multi-threaded CSV reading.
    timestamp_parsers: list, optional
        A sequence of strptime()-compatible format strings, tried in order when
        attempting to infer or convert timestamp values. The default is the special
        value ``pyarrow.csv.ISO8601`` uses the optimized internal ISO8601 parser.

    Returns
    -------
    astropy.table.Table
        An astropy Table containing the data from the CSV file.
    """
    check_has_pyarrow()
    from pyarrow import csv

    parse_options = csv.ParseOptions(
        delimiter=delimiter,
        quote_char=quotechar,
        double_quote=doublequote,
        escape_char=escapechar,
        newlines_in_values=newlines_in_values,
    )

    if comment is not None:
        input_file_stripped, header_start_new = strip_comment_lines(
            input_file, comment, encoding, header_start, data_start
        )
        if header_start_new is not None:
            header_start = header_start_new
        if input_file_stripped is not None:
            input_file = input_file_stripped

    read_options = get_read_options(header_start, data_start, names, encoding)
    convert_options = get_convert_options(
        include_names, dtypes, null_values, timestamp_parsers
    )

    table_pa = csv.read_csv(
        input_file,
        parse_options=parse_options,
        read_options=read_options,
        convert_options=convert_options,
    )
    table_apt = convert_pa_table_to_astropy_table(table_pa)

    return table_apt


def check_has_pyarrow():
    """Check for pyarrow and raise an informative exception if not available.

    Raises
    ------
    ModuleNotFoundError
        If pyarrow is not installed.
    """
    if not HAS_PYARROW:
        raise ModuleNotFoundError(
            "pyarrow is required for astropy.io.misc.pyarrow.csv functions"
        )


def convert_pa_string_array_to_numpy(
    arr: "pa.Array",
) -> "npt.NDArray":
    """
    Convert a PyArrow string array to a NumPy array.

    Parameters
    ----------
    arr : pyarrow.Array
        The PyArrow string array to be converted.

    Returns
    -------
    np.ndarray
        The converted NumPy array.
    """
    check_has_pyarrow()
    import pyarrow as pa
    import pyarrow.compute as pc

    # This implementation is faster than these alternatives:
    # >>> arr.to_numpy().astype(str)
    # >>> np.array(arr.to_pylist(), dtype='str')
    max_length = pa.compute.max(pc.utf8_length(arr))
    out = np.empty(len(arr), dtype=f"U{max_length}")
    out[:] = arr.to_numpy()
    return out


def pyarrow_zero(data_type: "pa.DataType"):
    """
    Return a "zero" value for the given PyArrow data type.

    This function provides a default value corresponding to the specified
    PyArrow data type. The returned value is intended to represent a
    neutral or "zero" equivalent for the type.

    Parameters
    ----------
    data_type : pa.DataType
        The PyArrow data type for which to determine the zero value.

    Returns
    -------
    Any
        A zero-equivalent value for the specified data type. The type of the
        returned value depends on the input data type:
        - Integer or floating types: 0
        - Boolean type: False
        - String or large string types: ""
        - Binary or large binary types: b""
        - Date type: datetime.date(1970, 1, 1)
        - Time type: datetime.time(0, 0, 0)
        - Timestamp type: datetime.datetime(1970, 1, 1)
        - Null type: None
        - Decimal type: decimal.Decimal(0)

    Raises
    ------
    NotImplementedError
        If no zero value is defined for the given data type.
    """
    check_has_pyarrow()
    import pyarrow as pa

    if pa.types.is_integer(data_type) or pa.types.is_floating(data_type):
        return 0
    elif pa.types.is_boolean(data_type):
        return False
    elif pa.types.is_string(data_type) or pa.types.is_large_string(data_type):
        return ""
    elif pa.types.is_binary(data_type) or pa.types.is_large_binary(data_type):
        return b""
    elif pa.types.is_date(data_type):
        return datetime.date(1970, 1, 1)
    elif pa.types.is_time(data_type):
        return datetime.time(0, 0, 0)
    elif pa.types.is_timestamp(data_type):
        return datetime.datetime(1970, 1, 1)
    elif pa.types.is_null(data_type):
        return None
    elif pa.types.is_decimal(data_type):
        import decimal

        return decimal.Decimal(0)
    else:
        raise NotImplementedError(f"No zero value defined for type: {data_type}")


def convert_pa_array_to_numpy(
    arr: Union["pa.Array", "pa.ChunkedArray"],
) -> "npt.NDArray":
    """
    Convert a PyArrow array to a NumPy array.

    Parameters
    ----------
    arr : pyarrow.Array, pyarrow.ChunkedArray
        The PyArrow array to be converted.

    Returns
    -------
    np.ndarray or np.ma.MaskedArray
        The converted NumPy array. If the input array contains null values, a masked
        array is returned with nulls masked out.

    Notes
    -----
    - If the input array does not contain null values, the result is an ndarray.
    - If the input array contains null values, they are replaced with a fill value
      (equivalent to ``np.zeros((), dtype)``)  before conversion, and a masked array is
      returned with the null positions masked.
    - The current implementation supports only int, float, bool, and string types.
    """
    check_has_pyarrow()
    import pyarrow as pa

    # Validate input type and set the fill value
    is_string = pa.types.is_string(arr.type)

    if arr.null_count == 0:
        # No nulls, just return an ndarray view of the pyarray
        out = convert_pa_string_array_to_numpy(arr) if is_string else arr.to_numpy()
    else:
        mask = arr.is_null().to_numpy()
        # Fill nulls in `arr` with zero. We do not know of a zero-copy fill for
        # pyarrow arrays with nulls, so we need to copy the data.
        arr = arr.fill_null(pyarrow_zero(arr.type))
        data = convert_pa_string_array_to_numpy(arr) if is_string else arr.to_numpy()
        out = np.ma.array(data, mask=mask, copy=False)

    return out


def convert_pa_table_to_astropy_table(table_pa) -> "Table":
    """
    Convert a PyArrow Table to an Astropy Table.

    Parameters
    ----------
    table_pa : ``pyarrow.Table``
        The PyArrow Table to be converted.

    Returns
    -------
    astropy.table.Table
        Converted astropy Table.
    """
    from astropy.table import Table

    columns = {
        name: convert_pa_array_to_numpy(col)
        for name, col in zip(table_pa.column_names, table_pa.itercolumns())
    }
    out = Table(columns, copy=False)
    return out


def strip_comment_lines(
    input_file: os.PathLike | str | BinaryIO,
    comment: str,
    encoding: str,
    header_start: int | None,
    data_start: int | None,
) -> tuple[io.BytesIO | None, int | None]:
    """
    Handle specified comment string when reading ``input_file``.

    This function has two modes of operation:

    1. If ``header_start`` is `None` or 0, then read the file line by line, looking for
       lines that start with the specified comment string. If all comment lines are at
       the beginning of the file (a common case), then the lines can be efficiently
       skipped within the pyarrow CSV reader. In this case the return value is (`None`,
       header_start after comments).
    2. Otherwise, the function reads the entire file into memory in a BytesIO object,
       removing all lines that start with the specified comment string along the way. In
       this case the return value is (BytesIO object, `None`).

    Checking for the comment string ignores leading whitespace.

    If ``input_file`` is a file path, it will be opened in binary read mode. The
    ``comment`` string is encoded to bytes using the provided encoding for comparison
    with the file content.

    Parameters
    ----------
    input_file : os.PathLike, str, or BinaryIO
        The input file path or file-like object to read from.
    comment : str
        The comment string that identifies lines to be removed.
    encoding : str
        The encoding to use for the input file.
    header_start : int or None
        The line index where the header starts. See ``read_csv`` for details.
    data_start : int or None
        The line index where the data starts. See ``read_csv`` for details.

    Returns
    -------
    output_file : io.BytesIO | None
        A BytesIO object containing the filtered content with comment lines removed, or
        `None` if all comment lines are at the beginning of the file.
    header_start : int | None
        The line index after the last comment line, or `None` if not all comment lines
        are at the beginning of the file.
    """
    comment_encode = comment.encode(encoding)

    with ExitStack() as stack:
        if isinstance(input_file, (str, os.PathLike)):
            input_file = stack.enter_context(open(input_file, "rb"))

        if header_start in (None, 0) and data_start is None:
            idx_last_comment = -1
            for idx, line in enumerate(input_file):
                if line.lstrip().startswith(comment_encode):
                    if idx - idx_last_comment == 1:
                        idx_last_comment = idx
                    else:
                        # Gap between comment lines, need to reset input file handle and
                        # break out to the logic below.
                        input_file.seek(0)
                        break
            else:
                input_file.seek(0)
                return None, idx_last_comment + 1

        # If we get here, we need to read the whole file and remove comment lines and
        # write into an output BytesIO file.
        output = io.BytesIO()
        for idx, line in enumerate(input_file):
            if not line.lstrip().startswith(comment_encode):
                output.write(line)

        # Set up for the consumer to read it.
        output.seek(0)

    return output, None


def get_convert_options(
    include_names: list[str] | None,
    dtypes: dict[str, "npt.DTypeLike"] | None,
    null_values: list[str] | None,
    timestamp_parsers: list[str] | None,
) -> "pyarrow.csv.ConvertOptions":
    """
    Generate PyArrow CSV convert options.

    Parameters
    ----------
    include_names : list or None
        List of column names to include in the conversion. If None, all columns are
        included.
    dtypes : dict or None
        Dictionary mapping column names to their respective data types. If None, default
        data types are used.
    null_values : list or None, optional (default None)
        List of strings to interpret as null values. By default, only empty strings are
        considered as null values.
    timestamp_parsers : list or None
        A sequence of strptime()-compatible format strings, tried in order when
        attempting to infer or convert timestamp values. The default is the special
        value ``pyarrow.csv.ISO8601`` uses the optimized internal ISO8601 parser.

    Returns
    -------
    pyarrow.csv.ConvertOptions
        PyArrow CSV ConvertOptions object configured with the specified column names and
        data types.
    """
    check_has_pyarrow()
    import pyarrow as pa
    from pyarrow import csv

    convert_options = csv.ConvertOptions()
    convert_options.strings_can_be_null = True
    convert_options.null_values = [""] if null_values is None else null_values

    if include_names is not None:
        convert_options.include_columns = include_names
    if dtypes is not None:
        convert_options.column_types = {
            colname: (
                dtype if isinstance(dtype, pa.DataType) else pa.from_numpy_dtype(dtype)
            )
            for colname, dtype in dtypes.items()
        }
    if timestamp_parsers is not None:
        convert_options.timestamp_parsers = timestamp_parsers
    return convert_options


def get_read_options(
    header_start: int | None,
    data_start: int | None,
    names: list[str] | None,
    encoding: str | None,
) -> "pyarrow.csv.ReadOptions":
    """
    Generate PyArrow CSV read options.

    Parameters
    ----------
    header_start : int or None
        The row index where the header starts. If None, no header is present.
    data_start : int or None
        The row index where the data starts. If None, data starts immediately after the
        header.
    names : list or None
        A list of column names. If None, column names will be autogenerated if no header
        is present.
    encoding : str or None
        The encoding of the CSV file. If None, the default encoding is used.

    Returns
    -------
    pyarrow.csv.ReadOptions
        The configured read options for reading the CSV file.
    """
    check_has_pyarrow()
    from pyarrow import csv

    if header_start is not None and names is not None:
        raise ValueError("cannot specify `names` unless `header_start=None`")

    read_options = csv.ReadOptions()

    # No header: header_start = None, data_start = None or int
    #   Set column_names to an empty list and autogenerate_column_names to True.
    if header_start is None:
        read_options.skip_rows_after_names = 0 if data_start is None else data_start
        if names is None:
            # With these two in combination, column names are generated automatically
            # as "f0", "f1", ...
            read_options.column_names = []
            read_options.autogenerate_column_names = True
        else:
            # Explicitly set column names.
            read_options.column_names = names
            read_options.autogenerate_column_names = False
    else:
        # Header present, header_start is not None
        if data_start is None:
            data_start = header_start + 1
        read_options.skip_rows = header_start
        read_options.skip_rows_after_names = data_start - header_start - 1

    if encoding is not None:
        read_options.encoding = encoding

    return read_options


def register_pyarrow_csv_table():
    """
    Register pyarrow.csv with Unified I/O as a Table reader.
    """
    from astropy.io import registry as io_registry
    from astropy.table import Table

    io_registry.register_reader("pyarrow.csv", Table, read_csv)
