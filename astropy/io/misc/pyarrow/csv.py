# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides functionality to read CSV files into Astropy Tables using PyArrow.
"""

import io
import os
from contextlib import ExitStack
from typing import TYPE_CHECKING, BinaryIO, Literal, Union

import numpy as np

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
    dtypes: dict[str, "npt.DTypeLike"] | None = None,
    comment: str | None = None,
    null_values: list | None = None,
    encoding: str = "utf-8",
    newlines_in_values: bool = False,
) -> "Table":
    """Read a CSV file into an astropy Table using ``pyarrow.csv.read_csv()``.

    Parameters
    ----------
    input_file : str, PathLike, or binary file-like object
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
    dtypes : dict, None, optional (default None)
        If provided, this is a dictionary of data types for output columns. Each key is
        a column name and the value is a data type object that is accepted as an
        argument to `np.dtype`. Examples include ``int``, ``np.float32``,
        ``np.dtype('f4')`` or ``"float32"``. Default is to infer the data types.
    comment : 1-character str or None, optional (default None)
        Character used to indicate the start of a comment. Any line starting with
        optional whitespace and then this character is ignored. Using this option will
        cause the parser to be slower and use more memory as it uses Python code to
        strip comments.

    Other Parameters
    ----------------
    null_values : list, optional (default None)
        List of strings to interpret as null values. By default, only empty strings are
        considered as null values.
    encoding: str, optional (default 'utf-8')
        Encoding of the input data.
    newlines_in_values: bool, optional (default False)
        Whether newline characters are allowed in CSV values. Setting this to True
        reduces the performance of multi-threaded CSV reading.

    Returns
    -------
    Table
        An astropy Table containing the data from the CSV file.
    """
    _, csv = get_pyarrow_csv()

    parse_options = csv.ParseOptions(
        delimiter=delimiter,
        quote_char=quotechar,
        double_quote=doublequote,
        escape_char=escapechar,
        newlines_in_values=newlines_in_values,
    )

    read_options = get_read_options(header_start, data_start, names, encoding)
    convert_options = get_convert_options(include_names, dtypes, null_values)

    if comment is not None:
        input_file = strip_comment_lines(input_file, comment, encoding)

    table_pa = csv.read_csv(
        input_file,
        parse_options=parse_options,
        read_options=read_options,
        convert_options=convert_options,
    )
    table_apt = convert_pa_table_to_astropy_table(table_pa)

    return table_apt


def get_pyarrow_csv():
    """Helper function to import pyarrow and pyarrow.csv."""
    if not HAS_PYARROW:
        raise ModuleNotFoundError(
            "pyarrow is required to read and write pyarrow.csv files"
        )
    import pyarrow as pa
    from pyarrow import csv

    return pa, csv


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
    - If the input array is of string type, it delegates the conversion to
      `convert_pa_string_array_to_numpy`.
    - If the input array contains null values, they are replaced with 0 before
      conversion, and a masked array is returned with the null positions masked.
    - If the input array does not contain null values, the result is an ndarray view of
      the underlying data.
    """
    is_string = arr.type == "string"

    if arr.null_count == 0:
        # No nulls, just return an ndarray view of the pyarray
        out = arr.to_numpy()
        if is_string:
            # `out` is an object array at this point, convert to fixed-length Unicode.
            # TODO: this is taking a lot of memory, there should be a memory-efficient
            # way to do this.
            out = out.astype(str)
    else:
        # Return a Masked Array copy, with nulls filled with zero.
        data = arr.fill_null("" if is_string else 0).to_numpy()
        if is_string:
            data = data.astype(str)
        # pa Bool array may not be zero-copyable, so set zero_copy_only=False to avoid
        # an exception.
        mask = arr.is_null().to_numpy(zero_copy_only=False)
        out = np.ma.array(data, mask=mask, copy=False)

    return out


def convert_pa_table_to_astropy_table(table_pa: "pa.Table") -> "Table":
    """
    Convert a PyArrow Table to an Astropy Table.

    Parameters
    ----------
    table_pa : pa.Table
        The PyArrow Table to be converted.

    Returns
    -------
    Table
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
) -> io.BytesIO:
    """
    Remove lines starting with a specified comment string from a file.

    If ``input_file`` is a file path, it will be opened in binary read mode. The
    ``comment`` string is encoded to bytes using the default encoding for comparison
    with the file content.

    Parameters
    ----------
    input_file : os.PathLike, str, or BinaryIO
        The input file path or file-like object to read from.
    comment : str
        The comment string that identifies lines to be removed.

    Returns
    -------
    io.BytesIO
        A BytesIO object containing the filtered content with comment lines removed.
    """
    comment_encode = comment.encode(encoding)
    output = io.BytesIO()
    with ExitStack() as stack:
        if isinstance(input_file, (str, os.PathLike)):
            input_file = stack.enter_context(open(input_file, "rb"))

        for line in input_file:
            if not line.lstrip().startswith(comment_encode):
                output.write(line)
        output.seek(0)

    return output


def get_convert_options(
    include_names: list | None,
    dtypes: dict[str, "npt.DTypeLike"] | None,
    null_values: list | None,
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

    Returns
    -------
    pyarrow.csv.ConvertOptions
        PyArrow CSV ConvertOptions object configured with the specified column names and
        data types.
    """
    pa, csv = get_pyarrow_csv()

    convert_options = csv.ConvertOptions()
    convert_options.strings_can_be_null = True

    if include_names is not None:
        convert_options.include_columns = include_names
    if null_values is not None:
        convert_options.null_values = null_values
    if dtypes is not None:
        convert_options.column_types = {
            colname: pa.from_numpy_dtype(dtype) for colname, dtype in dtypes.items()
        }
    return convert_options


def get_read_options(
    header_start: int | None,
    data_start: int | None,
    names: list | None,
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
    _, csv = get_pyarrow_csv()

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
