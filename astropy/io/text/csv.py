import io
import os
from typing import TYPE_CHECKING, BinaryIO

import astropy.table as apt

from . import core

if TYPE_CHECKING:
    import pyarrow.csv


def read_csv(
    input_file: os.PathLike | str | BinaryIO,
    delimiter: str = ",",
    quotechar: str = '"',
    double_quote: bool = True,
    escape_char: str | bool = False,
    header_start: int | None = 0,
    data_start: int | None = None,
    names: list | None = None,
    include_names: list | None = None,
    dtypes: dict | None = None,
    comment: str | None = None,
    encoding: str = "utf-8",
    newlines_in_values: bool = False,
    ignore_empty_lines: bool = True,
) -> apt.Table:
    """Read a CSV file into an astropy Table using pyarrow.read_csv.

    Parameters
    ----------
    input_file : str, PathLike, or binary file-like object
        File path or binary file-like object to read from.
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
        assumed and the column names are taken from ``names`` or generated automatically
        ("f0", "f1", ...).
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
        If provided, this is a dictionary of data types for output columns. Each key is
        a column name and the value is a data type object that is accepted as an
        argument to `np.dtype`. Examples include ``int``, ``np.float32``,
        ``np.dtype('f4')`` or ``"float32"``. Default is to infer the data types.
    comment: 1-character str or None, optional (default None)
        Character used to indicate the start of a comment. Any line starting with
        optional whitespace and then this character is ignored. Using this option will
        cause the parser to be slower and use more memory as it uses Python code to
        strip comments.
    encoding: str, optional (default 'utf-8')
        Encoding of the input data.
    newlines_in_values: bool, optional (default False)
        Whether newline characters are allowed in CSV values. Setting this to True
        reduces the performance of multi-threaded CSV reading.
    ignore_empty_lines: bool, optional (default True)
        Whether empty lines are ignored in CSV input. If False, an empty line is
        interpreted as containing a single empty value (assuming a one-column CSV file).
    """
    import pyarrow as pa
    import pyarrow.csv

    parse_options = pa.csv.ParseOptions(
        delimiter=delimiter,
        quote_char=quotechar,
        double_quote=double_quote,
        escape_char=escape_char,
        newlines_in_values=newlines_in_values,
        ignore_empty_lines=ignore_empty_lines,
    )

    if comment is not None:
        input_file = strip_comment_lines(input_file, comment)

    read_options = get_read_options(header_start, data_start, names, encoding)
    convert_options = get_convert_options(include_names, dtypes)

    table_pa = pyarrow.csv.read_csv(
        input_file,
        parse_options=parse_options,
        read_options=read_options,
        convert_options=convert_options,
    )
    table_apt = core.convert_pa_table_to_astropy_table(table_pa)

    return table_apt


def strip_comment_lines(
    input_file: os.PathLike | str | BinaryIO, comment: str
) -> BinaryIO:
    """Strip comment lines from input_file.

    A comment is any line that starts with optional whitespace followed by the comment
    character.

    Parameters
    ----------
    input_file : str, PathLike, or binary file-like object
        File path or binary file-like object to read from.
    comment: 1-character str
        Character used to indicate the start of a comment. Any line starting with this
        character is ignored.

    Returns
    -------
    BinaryIO
        BytesIO object with comments stripped.
    """
    if isinstance(input_file, (str, os.PathLike)):
        with open(input_file, "rb") as f:
            lines = f.readlines()
    else:
        lines = input_file.readlines()

    comment_encode = comment.encode()
    stripped_lines = [
        line for line in lines if not line.lstrip().startswith(comment_encode)
    ]
    return io.BytesIO(b"".join(stripped_lines))


def get_convert_options(
    include_names: list | None,
    dtypes: dict | None,
) -> "pyarrow.csv.ConvertOptions":
    """
    Generate PyArrow CSV conversion options based on included column names and data
    types.

    Parameters
    ----------
    include_names : list or None
        List of column names to include in the conversion. If None, all columns are
        included.
    dtypes : dict or None
        Dictionary mapping column names to their respective data types. If None, default
        data types are used.

    Returns
    -------
    pyarrow.csv.ConvertOptions
        PyArrow CSV ConvertOptions object configured with the specified column names and
        data types.
    """
    import pyarrow as pa
    import pyarrow.csv

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
    encoding: str | None,
) -> "pyarrow.csv.ReadOptions":
    """
    Generate read options for reading a CSV file using pyarrow.

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
    import pyarrow.csv

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

    if encoding is not None:
        read_options.encoding = encoding

    return read_options
