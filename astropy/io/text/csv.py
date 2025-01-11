from pathlib import Path
from typing import BinaryIO

import pyarrow as pa
import pyarrow.csv

import astropy.table as apt

from . import core


def read_csv(
    input_file: Path | str | BinaryIO,
    delimiter: str = ",",
    quotechar: str = '"',
    double_quote: bool = True,
    escape_char: str | bool = False,
    newlines_in_values: bool = False,
    ignore_empty_lines: bool = True,
    header_start: int | None = 0,
    data_start: int | None = None,
    names: list | None = None,
    include_names: list | None = None,
    dtypes: dict | None = None,
    comment: str | None = None,
) -> apt.Table:
    """Read a CSV file into an astropy Table using pyarrow.read_csv.

    Parameters
    ----------
    delimiter: 1-character str, optional (default ",")
        Character delimiting individual cells in the CSV data.
    quote_char: 1-character str or False, optional (default '"')
        Character used optionally for quoting CSV values (`False` if quoting is not
        allowed).
    double_quote: bool, optional (default `True`)
        Whether two quotes in a quoted CSV value denote a single quote in the data.
    escape_char: 1-character str or `False`, optional (default `False`)
        Character used optionally for escaping special characters (`False` if escaping
        is not allowed).
    header_start : int, None, optional (default 0)
        Line index for the header line with column names. If `None`, no header line is
        assumed and the column names are taken from ``names`` or generated
        automatically ("f0", "f1", ...).
    data_start : int, None, optional (default None)
        Line index for the start of data. If `None`, then data starts one line after the
         header, or on the first line if there is no header.
    names : list, None, optional (default None)
        List of names for input data columns. This can be used to override the column
        names inferred from the header line or to specify the column names if there is
        no header line.
    include_names : list, None, optional (default None)
        List of column names to include in output.
    dtypes : dict, None, optional (default None)
        Dictionary of data types for output columns. Each key is a column name and the
        value is a data type object that is accepted as an argument to `np.dtype`.
        Examples include ``int``, ``np.float32``, ``np.dtype('f4')`` or ``"float32"``.
    newlines_in_values: bool, optional (default False)
        Whether newline characters are allowed in CSV values. Setting this to True
        reduces the performance of multi-threaded CSV reading.
    ignore_empty_lines: bool, optional (default True)
        Whether empty lines are ignored in CSV input. If False, an empty line is
        interpreted as containing a single empty value (assuming a one-column CSV file).
    comment: 1-character str or None, optional (default None)
        Character used to indicate the start of a comment. Any line starting with this
        character is ignored. Using this option will cause the parser to be slower and
        use more memory as it uses Python code to strip comments. NOT IMPLEMENTED YET.
    encoding: str, optional (default 'utf-8')
        Encoding of the input data.
    """
    parse_options = pa.csv.ParseOptions(
        delimiter=delimiter,
        quote_char=quotechar,
        double_quote=double_quote,
        escape_char=escape_char,
        newlines_in_values=newlines_in_values,
    )

    # skip_rowsint, optional (default 0) The number of rows to skip before the column names
    # (if any) and the CSV data.
    #
    # skip_rows_after_names int, optional (default 0) The number of rows to skip after the
    # column names. This number can be larger than the number of rows in one block, and
    # empty rows are counted. The order of application is as follows: - skip_rows is applied
    # (if non-zero); - column names are read (unless column_names is set); -
    # skip_rows_after_names is applied (if non-zero).
    #
    # column_names list, optional The column names of the target table. If empty, fall back
    # on autogenerate_column_names.
    #
    # autogenerate_column_names bool, optional (default False) Whether to autogenerate column
    # names if column_names is empty. If true, column names will be of the form “f0”, “f1”…
    # If false, column names will be read from the first CSV row after skip_rows.

    if comment is not None:
        raise NotImplementedError("'comment' parameter is not implemented yet.")

    read_options = get_read_options(header_start, data_start, names)
    convert_options = get_convert_options(include_names, dtypes)

    table_pa = pyarrow.csv.read_csv(
        input_file,
        parse_options=parse_options,
        read_options=read_options,
        convert_options=convert_options,
    )
    table_apt = core.convert_pa_table_to_astropy_table(table_pa)

    return table_apt


def get_convert_options(
    include_names: list | None,
    dtypes: dict | None,
) -> pyarrow.csv.ConvertOptions:
    convert_options = pyarrow.csv.ConvertOptions()
    if include_names is not None:
        convert_options.include_columns = include_names
    if dtypes is not None:
        convert_options.column_types = {
            colname: pa.from_numpy_dtype(dtype) for colname, dtype in dtypes.items()
        }
    return convert_options


def get_read_options(
    header_start: int | None,
    data_start: int | None,
    names: list | None,
) -> pyarrow.csv.ReadOptions:
    read_options = pyarrow.csv.ReadOptions()

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
        # Header present
        if header_start is None:
            header_start = 0
        if data_start is None:
            data_start = header_start + 1
        read_options.skip_rows = header_start
        read_options.skip_rows_after_names = data_start - header_start - 1

    return read_options
