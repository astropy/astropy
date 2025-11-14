"""
ECSV Engine Module
--------------------------
This module provides functionality for reading and writing Enhanced Character Separated
Values (ECSV) files using various backends, including:
- `PyArrow <https://arrow.apache.org/docs/python/>`_
- `pandas <https://pandas.pydata.org/>`_
- `Astropy's ASCII engine <https://docs.astropy.org/en/stable/io/ascii/>`_

ECSV is a human-readable, YAML-encoded table format used in the Astropy
ecosystem for storing tables with metadata, units, and complex data types.

Key Features
------------
- Defines data structures for representing ECSV column and header metadata.
- Implements multiple ECSV reader engines, supporting PyArrow, pandas, and Astropy's
  ASCII CSV readers.
- Handles conversion between ECSV datatypes and numpy/pandas/pyarrow types, including
  support for JSON-encoded columns, multidimensional arrays, and masked data.
- Provides robust parsing of ECSV headers and data, including support for compressed
  files and in-memory file-like objects.
- Ensures compatibility with legacy ECSV files and provides liberal error handling for
  unknown datatypes.
- Integrates with Astropy's Unified I/O registry for seamless reading and writing of
  ECSV files.

Main Classes and Functions
--------------------------
- ``ColumnECSV``: Represents the attributes of a column as described in the ECSV header.
- ``ECSVHeader``: Encapsulates the parsed ECSV header, including column definitions and
  table metadata.
- ``ECSVEngine`` and subclasses: Abstract base class and concrete implementations for
  different CSV parsing engines.
- ``read_ecsv``: Reads an ECSV file and returns an Astropy ``Table`` object, handling all
  necessary conversions and metadata.
- ``write_ecsv``: Writes an Astropy ``Table`` to an ECSV file.
- ``register_pyarrow_ecsv_table``: Registers the PyArrow ECSV reader/writer with
  Astropy's I/O registry.

Usage
-----
This module is intended for internal use within Astropy and for advanced users who need
fine-grained control over ECSV parsing and engine selection. For most users, reading and
writing ECSV files can be accomplished via the high-level ``Table.read`` and
``Table.write`` interfaces, e.g.:
```
Table.read(filename, format="ecsv", engine="pyarrow.csv")
```

Dependencies
------------
- numpy
- astropy.table
- pyarrow (optional, for PyArrow engine)
- pandas (optional, for pandas engine)
"""

import abc
import collections
import functools
import io
import json
import os
import re
import warnings
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Final, NamedTuple

import numpy as np
import numpy.typing as npt

from astropy.utils.data import get_readable_fileobj

if TYPE_CHECKING:
    from astropy.table import SerializedColumn, Table

__all__ = [
    "ColumnECSV",
    "ECSVEngine",
    "ECSVEngineIoAscii",
    "ECSVEnginePandas",
    "ECSVEnginePyArrow",
    "ECSVHeader",
    "read_ecsv",
    "register_ecsv_table",
    "write_ecsv",
]

ECSVEngines: Final[dict[str, "ECSVEngine"]] = {}


class DerivedColumnProperties(NamedTuple):
    """Named tuple for derived properties of a ECSV column specification.

    Attributes
    ----------
    csv_np_type : str
        Numpy type string for the CSV column data, e.g. "int64", "float32", "str".
        This is derived from the ECSV `datatype` and `subtype`.
    dtype : str
        Numpy dtype in the final column data. This may differ from `csv_np_type` in
        some cases, e.g. for JSON-encoded columns.
    shape : tuple[int, ...]
        Shape of the final column data as a tuple of integers. This is derived from the
        ECSV `subtype` if applicable, or an empty tuple for scalar columns.
    """

    csv_np_type: str
    dtype: str
    shape: tuple[int, ...]


@dataclass(frozen=True)
class ColumnECSV:
    """
    Class representing attributes of a column in an ECSV header.

    Attributes
    ----------
    name : str
        The name of the column.
    datatype : str
        The data type of the column as specified in the ECSV header.
    subtype : str or None, optional
        The subtype of the column, if applicable.
    unit : str or None, optional
        The unit of the column values, if specified.
    description : str or None, optional
        A description of the column.
    format : str or None, optional
        The format string for the column values.
    meta : dict or None, optional
        Additional metadata associated with the column.

    Properties
    ----------
    csv_np_type : str
        Numpy type string describing the column CSV data. In practice this is the same
        as the ECSV ``datatype`` except that "string" => "str". This is provided to the
        engine ``convert_np_type()`` method to generate the engine-specific type
        provided to the CSV reader. For instance, for pandas the ``int32`` type gets
        converted to ``Int32`` to read columns as a nullable int32.
    dtype : np.dtype
        Numpy dtype in the final column data. This may be entirely different from
        ``csv_np_type`` in some cases, in particular JSON-encoded fields.
    shape : tuple of int
        Shape of the final column data.
    """

    name: str
    datatype: str
    subtype: str | None = None
    unit: str | None = None
    description: str | None = None
    format: str | None = None
    meta: dict | None = None

    @functools.cached_property
    def csv_np_type(self) -> str:
        """Numpy type string describing the column CSV data."""
        return self._derived_properties.csv_np_type

    @functools.cached_property
    def dtype(self) -> np.dtype:
        """Numpy dtype in the final column data"""
        return np.dtype(self._derived_properties.dtype)

    @functools.cached_property
    def shape(self) -> tuple[int, ...]:
        """Shape of the column data"""
        return self._derived_properties.shape

    @functools.cached_property
    def _derived_properties(self) -> DerivedColumnProperties:
        """Get the csv_np_type, dtype, and shape of the column from ECSV header."""
        return get_csv_np_type_dtype_shape(self.datatype, self.subtype, self.name)


@dataclass(slots=True, frozen=True)
class ECSVHeader:
    """
    Class representing the information in an ECSV header.

    Attributes
    ----------
    n_header : int
        Total number of header lines in the ECSV file (including empty).
    n_empty : int
        Number of empty lines in the header section.
    n_comment : int
        Number of comment lines in the header section.
    cols : list[ColumnECSV]
        List of ``ColumnECSV`` objects describing the attributes of each column in
        header.
    table_meta : dict
        Metadata associated with the table, typically parsed from the ECSV header.
    delimiter : str
        Delimiter character used to separate columns in the ECSV file.
    """

    n_header: int
    n_empty: int
    n_comment: int
    cols: list[ColumnECSV]
    table_meta: dict
    delimiter: str


class ECSVEngine(metaclass=abc.ABCMeta):
    """Base class for ECSV reader engines.

    An engine is responsible for reading the raw CSV data that follows the ECSV header.
    This assumes that the engine has a defined Table Unified I/O interface.

    - `name` and `format` must be defined as class attributes in subclasses.
    - `engines` is a base class-level dictionary that maps engine names to their
      respective engine classes. Subclasses should not modify this directly.

    Properties
    ----------
    name : str
        Name of the engine, used for ``engine`` parameter in a call like:
        ``Table.read(filename, format="ecsv", engine="pyarrow")``.
    format : str
        Format string for the engine CSV reader, e.g. "pyarrow.csv", "ascii.csv", etc.
    engines : dict[str, ECSVEngine]
        Dictionary mapping engine names to their respective engine classes.
    """

    name: str | None = None
    format: str | None = None
    engines: dict[str, "ECSVEngine"] = {}

    def __init_subclass__(cls, **kwargs):
        """Register the subclass as an ECSV engine."""
        super().__init_subclass__(**kwargs)

        # Ensure that the subclass has the required string class attributes.
        for attr in ("name", "format"):
            if not isinstance(val := getattr(cls, attr, None), str):
                raise TypeError(
                    f"Subclasses of ECSVEngine must define a class attribute '{attr}' "
                    f"as a string, got {type(val)}."
                )

        cls.engines[cls.name] = cls

    @abc.abstractmethod
    def convert_np_type(self, np_type: str) -> Any:
        """
        Convert a numpy type string to engine-specific type for parsing.

        For instance, for pandas the ``"int32"`` numpy type gets converted to an
        ``Int32Dtype()`` instance to read columns as a nullable int32.

        Parameters
        ----------
        np_type : str
            The numpy type string to be converted.

        Returns
        -------
        Any
            Corresponding engine-specific type.
        """

    def get_converters(self, header):
        """
        Get a dictionary of converters for the columns in the ECSV header.

        This is used to convert column names to engine-specific types.

        Parameters
        ----------
        header : ECSVHeader
            The ECSV header containing column definitions.

        Returns
        -------
        dict[str, Any]
            Dictionary mapping column names to engine-specific converters.
        """
        return {col.name: self.convert_np_type(col.csv_np_type) for col in header.cols}

    @abc.abstractmethod
    def get_data_kwargs(
        self,
        header: ECSVHeader,
        null_values: list[str],
    ) -> dict[str, Any]:
        """
        Generate a dictionary of keyword arguments for data parsing.

        This accounts for the API variations in each engine CSV reader.

        Parameters
        ----------
        header : ECSVHeader
            ECSVHeader object within header information.
        null_values : list of str
            List of strings with values to be interpreted as null or missing data.

        Returns
        -------
        dict[str, Any]
            Dict of keyword arguments to be passed to engine CSV reader.
        """


class ECSVEnginePyArrow(ECSVEngine):
    """ECSV reader engine using PyArrow."""

    name = "pyarrow"
    format = "pyarrow.csv"

    def convert_np_type(self, np_type: str) -> str:
        # PyArrow does not support float128 and there is no workaround (unlike float16).
        if np_type == "float128":
            raise TypeError(
                "pyarrow engine does not support float128, choose a different engine"
            )

        # PyArrow does not support float16, so we need to convert it to float32.
        # The final output is still cast as float16.
        return "float32" if np_type == "float16" else np_type

    def get_data_kwargs(
        self,
        header: ECSVHeader,
        null_values: list[str],
    ) -> dict[str, Any]:
        # See base method for details.
        kw = {}
        kw["null_values"] = null_values
        kw["header_start"] = header.n_header
        kw["dtypes"] = self.get_converters(header)
        return kw


class ECSVEngineIoAscii(ECSVEngine):
    """ECSV reader engine using astropy.io.ascii Python CSV reader."""

    name = "io.ascii"
    format = "ascii.csv"

    def convert_np_type(self, np_type: str) -> np.generic:
        # Convert the np_type string to a numpy dtype type like np.int32, np.float64,
        # etc. This output is compatible with io.ascii `converters` option where is gets
        # used.
        return np.dtype(np_type).type

    def get_data_kwargs(
        self,
        header: ECSVHeader,
        null_values: list[str],
    ) -> dict[str, Any]:
        kw = {}
        kw["fill_values"] = get_null_values_per_column(
            header.cols, header.table_meta, null_values
        )
        kw["header_start"] = header.n_header - header.n_empty
        kw["converters"] = self.get_converters(header)
        # Fast reader does not support converters (defining types in advance) nor any
        # encoding. Converters are required, e.g. for a string column that looks like
        # floats. Would be nice to fix this, but in mean time use Python CSV reader.
        kw["fast_reader"] = False
        kw["strip_column_names"] = False
        return kw


class ECSVEnginePandas(ECSVEngine):
    """ECSV reader engine using pandas."""

    name = "pandas"
    format = "pandas.csv"

    def convert_np_type(self, np_type: str) -> np.dtype:
        # Convert the np_type to a pandas dtype will support for nullable types.
        import pandas as pd

        dtype = np.dtype(np_type)
        if dtype.kind in ("i", "u"):
            # Convert int64 to Int64, uint32 to UInt32, etc for nullable types
            converter = dtype.name.replace("i", "I").replace("u", "U")
        elif dtype.kind == "b":
            converter = "boolean"
        else:
            converter = np_type
        return pd.api.types.pandas_dtype(converter)

    def get_data_kwargs(
        self,
        header: ECSVHeader,
        null_values: list[str],
    ) -> dict[str, Any]:
        fill_values = get_null_values_per_column(
            header.cols, header.table_meta, null_values
        )
        null_values = collections.defaultdict(list)
        converters = self.get_converters(header)
        for null_value, _, col_name in fill_values:
            null_values[col_name].append(null_value)
            # Pandas parser does not natively parse nan or NaN for floats, so we need
            # to declare this as a null value.
            if converters[col_name].kind == "f":
                for nan in ("nan", "NaN"):
                    null_values[col_name].append(nan)

        kw = {
            "na_values": null_values,
            "keep_default_na": False,
            "comment": "#",
            "dtype": converters,
        }
        # Would prefer setting `"skiprows": header.n_header` above (as in the original
        # implementation prior to #18756) instead of "comment": "#". However there is a
        # bug in pandas.read_csv where skiprows does not work when the line includes a
        # quote character, see https://github.com/pandas-dev/pandas/issues/62739.
        return kw  # noqa: RET504


def is_numpy_dtype(np_type: str) -> bool:
    # Check if the given dtype is a valid numpy dtype.
    try:
        np.dtype(np_type)
    except Exception:
        return False
    else:
        return True


def get_header_lines(
    input_file: str | os.PathLike | io.BytesIO,
    encoding="utf-8",
) -> tuple[list[str], int, int, int]:
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
        provided, the function automatically handles compressed file formats
        that are supported by `~astropy.utils.data.get_readable_fileobj`.
    encoding : str, optional
        The encoding used to decode the file content. Default is "utf-8".

    Returns
    -------
    lines : list[str]
      List of decoded header lines without the header prefix.
    idx : int
        Index of the last line read.
    n_empty : int
        Number of empty lines read.
    n_comment : int
        Number of comment lines read.
    """
    header_prefix = "# ".encode(encoding)
    comment_prefix = "##".encode(encoding)
    lines = []
    n_empty = 0
    n_comment = 0

    with get_readable_fileobj(input_file, encoding="binary") as f:
        for idx, line in enumerate(f):
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


def get_csv_np_type_dtype_shape(
    datatype: str, subtype: str | None, name: str
) -> DerivedColumnProperties:
    """Get the csv_np_type, dtype, and shape of the column from datatype and subtype.

    This function implements most of the complexity of the ECSV data type handling. The
    ECSV standard allows for a wide variety of data types and subtypes, and we also need
    to handle some legacy cases and be liberal in what we accept.

    Parameters
    ----------
    datatype : str
        The data type of the column as specified in the ECSV header.
    subtype : str or None
        The subtype of the column, if applicable. This can include additional
        information like array shape or JSON serialization.
    name : str
        The name of the column, used for error messages.

    Returns
    -------
    CSVNpTypeDtypeShape
        A named tuple containing:
        - `csv_np_type`: Numpy type string for the CSV column data.
        - `dtype`: Numpy dtype in the final column data.
        - `shape`: Shape of the final column data as a tuple of integers.

    Raises
    ------
    ValueError
        If the `datatype` or `subtype` is not recognized or cannot be converted to a
        valid numpy dtype.
    InconsistentTableError
        If the `datatype` is not in the allowed ECSV datatypes and cannot be parsed as a
        numpy dtype.
    """
    from astropy.io.ascii.core import InconsistentTableError
    from astropy.io.ascii.ecsv import ECSV_DATATYPES, InvalidEcsvDatatypeWarning

    csv_np_type = "str" if datatype == "string" else datatype
    dtype = csv_np_type
    shape = ()

    if datatype not in ECSV_DATATYPES:
        msg = (
            f"unexpected datatype {datatype!r} of column {name!r} "
            f"is not in allowed ECSV datatypes {ECSV_DATATYPES}."
        )
        # Try being liberal on input if the `csv_np_type` (derived from ECSV
        # `datatype`) looks like a numpy dtype. In this case, parse the column as
        # string and then cast as `csv_np_type`. This allows for back-compatibility
        # with early versions of io.ascii.ecsv that wrote and read e.g.
        # datatype=datetime64.
        if is_numpy_dtype(csv_np_type):
            dtype = csv_np_type
            csv_np_type = "str"
            warnings.warn(msg, InvalidEcsvDatatypeWarning)
        else:
            # No joy, this is an exception
            raise InconsistentTableError(msg)

    if subtype and csv_np_type != "str":
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
                f"{name!r}, using dtype={csv_np_type!r} instead.",
                category=InvalidEcsvDatatypeWarning,
            )
            dtype = csv_np_type

    return DerivedColumnProperties(csv_np_type, dtype, shape)


def read_header(
    input_file: str | os.PathLike | io.BytesIO,
    encoding: str = "utf-8",
) -> ECSVHeader:
    """
    Read and parse the header of an ECSV (Enhanced Character Separated Values) input.

    This function extracts and validates the ECSV header from the given input file,
    parses the YAML metadata, and constructs the corresponding ECSVHeader object
    containing column definitions and table metadata.

    Parameters
    ----------
    input_file : str, os.PathLike, or io.BytesIO
        The path to the ECSV file or a file-like object containing the ECSV data.
    encoding : str, optional
        The encoding to use when reading the file. Default is 'utf-8'.

    Returns
    -------
    ECSVHeader
        An object containing header information, including the number of header lines,
        number of empty lines, column attributes, table metadata, and delimiter.

    Raises
    ------
    InconsistentTableError
        If the ECSV header is missing, malformed, or the YAML metadata cannot be parsed.
    ValueError
        If the delimiter specified in the header is not supported.

    Notes
    -----
    The function expects the first non-blank comment line to be the ECSV version header,
    and only space and comma are allowed as delimiters in the ECSV format.
    """
    from astropy.io.ascii.core import InconsistentTableError
    from astropy.io.ascii.ecsv import DELIMITERS
    from astropy.table import meta

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
    cols_attrs = [ColumnECSV(**col) for col in header["datatype"]]

    return ECSVHeader(n_header, n_empty, n_comment, cols_attrs, table_meta, delimiter)


def read_data(
    input_file: str | os.PathLike | io.BytesIO,
    header: ECSVHeader,
    null_values: list[str],
    encoding: str = "utf-8",
    engine_name: str = "io.ascii",
) -> "Table":
    """
    Read the data from an ECSV table using the specified engine.

    This function uses an engine-specific class to handle reading and converting
    the data according to the ECSV specification and the selected backend.

    Parameters
    ----------
    input_file : str, os.PathLike, or io.BytesIO
        The path to the input file or a file-like object containing the ECSV data.
    header : ECSVHeader
        The parsed ECSV header containing column definitions and metadata.
    null_values : list of str
        List of string values to interpret as null/missing values in the data.
    encoding : str, optional
        The encoding to use when reading the file. Default is "utf-8".
    engine_name: str, optional
        The backend engine to use for reading the data. Default is "io.ascii".
        Built-in options are "pyarrow", "pandas", and "io.ascii".

    Returns
    -------
    Table
        An Astropy Table containing the data read from the ECSV file.

    Raises
    ------
    InconsistentTableError
        If the column names from the ECSV header do not match the column names
        in the data.
    """
    from astropy.table import Table

    engine = ECSVEngine.engines[engine_name]()

    # Get the engine-specific kwargs for reading the CSV data.
    kwargs = engine.get_data_kwargs(header, null_values)

    data = Table.read(
        input_file,
        format=engine.format,
        delimiter=header.delimiter,
        encoding=encoding,
        **kwargs,
    )

    # Ensure ECSV header names match the data column names.
    ecsv_header_names = [col.name for col in header.cols]
    if ecsv_header_names != data.colnames:
        from astropy.io.ascii.core import InconsistentTableError

        raise InconsistentTableError(
            f"column names from ECSV header {ecsv_header_names} do not "
            f"match names from header line of CSV data {data.colnames}"
        )

    return data


def get_str_vals(
    data: np.ndarray | np.ma.MaskedArray,
) -> tuple[list[str] | npt.NDArray[np.str_], npt.NDArray[np.bool_] | None]:
    """Get the string values and the mask if available.

    This assumes a 1-d input array of strings, possibly masked. This array comes from
    reading the ECSV data, which is always a 1-d array. This function is only called if
    that array is a numpy string array or a masked array of strings.

    For a masked array it converts the data to the equivalent Python representation
    (list of strings) and returns the mask as a separate array.

    A list of strings is required in this case because the subsequent
    ``process_*_data`` functions substitute (in-place) a new string with the appropriate
    JSON for an empty/masked fill value. In particular, if the original input consists
    solely of empty strings (which is legal), the numpy string array will be not be wide
    enough to hold the fill value.

    For regular numpy arrays it simply returns the original data as a numpy array.

    Parameters
    ----------
    data : np.ndarray | np.ma.MaskedArray
        The input data array to extract string values from.

    Returns
    -------
    str_vals : list[str] | npt.NDArray[np.str_]
        A list of strings or a 1D numpy array of strings representing the data.
    mask : npt.NDArray[np.bool_] | None
        A 1D numpy array of booleans indicating the mask, or None if not applicable.
    """
    # For masked we need a list because for multidim the data under the mask is set
    # to a compatible value.
    if hasattr(data, "mask"):
        # TODO: for not NUMPY_LT_2_0, try changing this to:
        # str_vals = data.astype("T")
        str_vals = data.view(np.ndarray).tolist()
        mask = data.mask
    else:
        str_vals = data
        mask = None
    return str_vals, mask


def convert_column(
    col: ColumnECSV,
    data_in: np.ndarray | np.ma.MaskedArray,
) -> np.ndarray | np.ma.MaskedArray:
    """
    Convert column data from original CSV numpy type to specified output numpy dtype.

    This function handles both regular scalar columns and more complex columns such as:

    - Object dtype columns containing arbitrary Python objects serializable to JSON.
    - Variable-length array columns, where the last axis may vary in length.
    - Fixed-shape multidimensional columns.

    Depending on the column's dtype and shape, the function selects the appropriate
    processing routine to convert the data, including deserialization from JSON where
    necessary. For regular scalar columns, it casts the data to the target dtype if
    needed.

    Parameters
    ----------
    col : ColumnECSV
        The column specification, including dtype, shape, and name.
    data_in : np.ndarray | np.ma.MaskedArray
        The input data array to be converted.

    Returns
    -------
    np.ndarray | np.ma.MaskedArray
        The converted data array with the appropriate dtype and shape.

    Raises
    ------
    ValueError
        If the data cannot be converted due to shape mismatch or invalid JSON content.
    """
    try:
        if col.dtype == "object" or col.shape:
            # Handle three distinct column types where each row element is serialized to
            # JSON. In this case ``data_in`` is an ndarray or MaskedArray of
            # fixed-length string which are the JSON-encoded representation of the data.
            # See docstring in `get_str_vals` for explanation of the next step, which
            # has some subtlety.
            str_vals, mask = get_str_vals(data_in)

            if col.dtype == "object":
                # Any Python objects serializable to JSON
                process_func = process_object_data
            elif col.shape[-1] is None:
                # Variable length arrays with shape (n, m, ..., *) for fixed
                # n, m, .. and variable in last axis.
                process_func = process_variable_length_array_data
            else:
                # Multidim columns with consistent shape (n, m, ...).
                process_func = process_fixed_shape_multidim_data

            data_out, col_shape = process_func(col, str_vals, mask)

        # Regular scalar value column
        else:
            data_out = data_in
            # If we need to cast the data to a different dtype, do it now.
            if data_out.dtype != col.dtype:
                data_out = data_out.astype(col.dtype)
            col_shape = col.shape

        if data_out.shape[1:] != tuple(col_shape):
            raise ValueError("shape mismatch between value and column specifier")

    except json.JSONDecodeError:
        raise ValueError(
            f"column {col.name!r} failed to convert: column value is not valid JSON"
        )
    except Exception as exc:
        raise ValueError(f"column {col.name!r} failed to convert: {exc}") from exc

    return data_out


def process_object_data(
    col: ColumnECSV,
    str_vals: list[str] | npt.NDArray[np.str_],
    mask: npt.NDArray[np.bool_] | None,
) -> tuple[np.ndarray | np.ma.MaskedArray, tuple[int, ...]]:
    """
    Handle object columns where each row element is a JSON-encoded object.

    The ECSV format only allows a 1-d column of object type.

    Example::

        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: objects, datatype: string, subtype: json}
        # schema: astropy-2.0
        objects
        "{""a"":1}"
        "{""b"":[2.5,null]}"
        true

    Parameters
    ----------
    col : ColumnECSV
        The column specification, including dtype, shape, and name.
    str_vals : list[str] or 1-D ndarray[str]
        JSON-encoded string representations of the data.
    mask : 1-D ndarray[bool] or None
        Boolean mask array 1-D indicating invalid or missing values. If None, no masking
        is applied.

    Returns
    -------
    data_out : numpy.ndarray or numpy.ma.MaskedArray
        An array of objects reconstructed from `str_vals`, with the same shape as `col`.
        If `mask` is provided, a masked array is returned with the mask applied.
    col_shape : tuple[int, ...]
        Expected shape of data_out, used in final sanity check of reading.
    """
    if mask is not None:
        for idx in np.nonzero(mask)[0]:
            str_vals[idx] = "0"  # could be "null" but io.ascii uses "0"
    col_vals = [json.loads(val) for val in str_vals]
    np_empty = np.empty if mask is None else np.ma.empty
    data_out = np_empty((len(col_vals),) + tuple(col.shape), dtype=object)
    data_out[...] = col_vals
    if mask is not None:
        data_out.mask = mask
    return data_out, col.shape


def process_fixed_shape_multidim_data(
    col: ColumnECSV,
    str_vals: list[str] | npt.NDArray[np.str_],
    mask: npt.NDArray[np.bool_] | None,
) -> tuple[np.ndarray | np.ma.MaskedArray, tuple[int, ...]]:
    """
    Handle fixed-shape multidimensional columns as JSON-encoded strings.

    Example::

        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: array3x2, datatype: string, subtype: 'float64[3,2]'}
        # schema: astropy-2.0
        array3x2
        [[0.0,1.0],[2.0,3.0],[4.0,5.0]]
        [[6.0,7.0],[8.0,null],[10.0,11.0]]

    Parameters
    ----------
    col : ColumnECSV
        The column specification, including dtype, shape, and name.
    str_vals : array-like of str
        Array of string representations of the data, typically JSON-encoded.
    mask : numpy.ndarray or None
        Boolean mask array indicating invalid or missing values. If None, no masking is
        applied.

    Returns
    -------
    data_out : numpy.ndarray or numpy.ma.MaskedArray
        An array of objects reconstructed from `str_vals`, with the same shape as `col`.
        If `mask` is provided, a masked array is returned with the mask applied.
    col_shape : tuple[int, ...]
        Expected shape of data_out, used in final sanity check of reading.
    """
    # Change empty (blank) values in original ECSV to something
    # like "[[null, null],[null,null]]" so subsequent JSON
    # decoding works.
    if mask is not None:
        all_none_arr = np.full(shape=col.shape, fill_value=None, dtype=object)
        fill_value = json.dumps(all_none_arr.tolist())
        for idx in np.nonzero(mask)[0]:
            str_vals[idx] = fill_value

    col_vals = [json.loads(val) for val in str_vals]

    # Make a numpy object array of col_vals to look for None (masked values)
    arr_vals = np.array(col_vals, dtype=object)
    arr_vals_mask = arr_vals == None
    if np.any(arr_vals_mask):
        # Replace all the None with an appropriate fill value
        arr_vals[arr_vals_mask] = {"U": "", "S": b""}.get(col.dtype.kind, 0)
        # Finally make a MaskedArray with the filled data + mask
        data_out = np.ma.array(arr_vals.astype(col.dtype), mask=arr_vals_mask)
    else:
        data_out = arr_vals.astype(col.dtype)

    return data_out, col.shape


def process_variable_length_array_data(
    col: ColumnECSV,
    str_vals: list[str] | npt.NDArray[np.str_],
    mask: npt.NDArray[np.bool_] | None,
) -> tuple[np.ndarray | np.ma.MaskedArray, tuple[int, ...]]:
    """
    Handle variable length arrays with shape (n, m, ..., *) as JSON-encoded strings.

    Shape is fixed for n, m, .. and variable in last axis. The output is a 1-d object
    array with each row element being an ``np.ndarray`` or ``np.ma.masked_array`` of the
    appropriate shape.

    Example::

        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: array_var, datatype: string, subtype: 'int64[null]'}
        # schema: astropy-2.0
        array_var
        [1,2]
        [3,4,5,null,7]
        [8,9,10]

    Parameters
    ----------
    col : ColumnECSV
        The column specification, including dtype, shape, and name.
    str_vals : list[str] or 1-D ndarray[str]
        JSON-encoded string representations of the data.
    mask : 1-D ndarray[bool] or None
        Boolean mask array 1-D indicating invalid or missing values. If None, no masking
        is applied.

    Returns
    -------
    data_out : numpy.ndarray or numpy.ma.MaskedArray
        An array of objects reconstructed from `str_vals`, with the same shape as `col`.
        If `mask` is provided, a masked array is returned with the mask applied.
    col_shape : tuple[int, ...]
        Expected shape of data_out, used in final sanity check of reading.
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
            arr_val = np.array(obj_val, dtype=col.dtype)
        except TypeError:
            # obj_val has entries that are inconsistent with
            # dtype. For a valid ECSV file the only possibility
            # is None values (indicating missing values).
            vals = np.array(obj_val, dtype=object)
            # Replace all the None with an appropriate fill value
            mask_vals = vals == None
            vals[mask_vals] = {"U": "", "S": b""}.get(col.dtype.kind, 0)
            arr_val = np.ma.array(vals.astype(col.dtype), mask=mask_vals)

        col_vals.append(arr_val)

    np_empty = np.empty if mask is None else np.ma.empty
    data_out = np_empty(len(col_vals), dtype=object)
    data_out[:] = col_vals
    if mask is not None:
        data_out.mask = mask
    return data_out, ()


def get_null_values_per_column(
    cols: list[ColumnECSV],
    table_meta: dict | None,
    null_values: list[str],
) -> list[tuple[str, str, str]]:
    """Get null and fill values for individual columns.

    For ECSV to handle the corner case of data that has been serialized using the
    serialize_method='data_mask' option, which writes the full data and mask directly,
    AND where that table includes a string column with zero-length string entries ("")
    which are valid data. Normally the super() method will set col.fill_value=('', '0')
    to replace blanks with a '0'.  But for that corner case subset, instead do not do
    any filling.

    Parameters
    ----------
    cols : list[ColumnECSV]
        List of ColumnECSV objects representing the columns in the ECSV file.
    table_meta : dict or None
        Metadata dictionary from the ECSV header, which may include serialized columns
        and other metadata.
    null_values : list[str]
        List of string values to interpret as null/missing values in every column.
        The upstream default from ``read_ecsv`` is [""] but no default is defined here.

    Returns
    -------
    fill_values : list[tuple[str, str, str]]
        A list of tuples with (null_value, fill_value, column_name) for each column in
        `cols` that is not a MaskedColumn. If no fill values are needed, returns an
        empty list.
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

        fill_value = "" if col.csv_np_type == "str" else "0"
        for null_value in null_values:
            fill_values.append((null_value, fill_value, col.name))

    return fill_values


def read_ecsv(
    input_file: str | os.PathLike | io.BytesIO | io.StringIO | Iterable[str],
    *,
    encoding: str = "utf-8",
    engine: str = "io.ascii",
    null_values: list[str] | None = None,
) -> "Table":
    """
    Read an ECSV (Enhanced Character Separated Values) file and return an Astropy Table.

    Parameters
    ----------
    input_file : str, os.PathLike, io.BytesIO, io.StringIO, Iterable[str]
        The ECSV input to read. This can be a file path, a file-like object, a string
        containing the file contents, or an iterable of strings representing lines of
        the file. Note that providing ``io.StringIO`` or an iterable of strings will be
        less memory efficient, as it will be converted to a bytes stream.
    encoding : str, optional
        The encoding to use when reading the file. Default is "utf-8".
    engine : str, optional
        The engine to use for reading the CSV data. Default is "io.ascii", which uses
        astropy to read the CSV data. Other built-in options are "pyarrow" and "pandas".
        The "pyarrow" engine is optimized for performance and can handle large datasets
        efficiently. The "pandas" engine uses the pandas CSV reader, which is also
        faster than the default "io.ascii" engine.
    null_values : list of str or None, optional
        List of string values to interpret as null/missing values. Default is [""]. The
        ECSV standard requires the null values are represented as empty strings in the
        CSV data, but this allows reading non-compliant ECSV files. A notable example
        are the Gaia source download files which are ECSV but use "null".

    Returns
    -------
    table : astropy.table.Table
        The table read from the ECSV file.

    Raises
    ------
    astropy.io.ascii.core.InconsistentTableError
        If the column names in the ECSV header do not match the column names in the CSV
        data.

    Notes
    -----
    - The function handles various input types, including file paths, file-like objects,
      and in-memory strings or lists of strings.
    - Metadata and column attributes (such as unit, description, format, and meta) are
      transferred from the ECSV header to the resulting Table object.
    - Handles JSON-encoded data and ensures appropriate numpy dtypes for columns.
    """
    from astropy.table import Table, serialize

    if null_values is None:
        null_values = [""]

    # Allow input types that are historically supported by io.ascii. These will not be
    # memory or speed efficient but will still work.
    if isinstance(input_file, io.StringIO):
        input_file = io.BytesIO(input_file.getvalue().encode(encoding))
    elif isinstance(input_file, str) and "\n" in input_file:
        input_file = io.BytesIO(input_file.encode(encoding))
    elif isinstance(input_file, (list, tuple)):
        # TODO: better way to check for an iterable of str?
        input_file = io.BytesIO("\n".join(input_file).encode(encoding))

    # Read the ECSV header from the input.
    header = read_header(input_file, encoding=encoding)

    # Read the CSV data from the input starting at the line after the header. This
    # includes handling that is particular to the engine.
    data_raw = read_data(
        input_file,
        header,
        null_values=null_values,
        encoding=encoding,
        engine_name=engine,
    )

    # Convert the column data to the appropriate numpy dtype. This is mostly concerned
    # with JSON-encoded data but also handles cases like pyarrow not supporting float16.
    data = {col.name: convert_column(col, data_raw[col.name]) for col in header.cols}

    # Create the Table object
    table = Table(data)

    # Transfer metadata from the ECSV header to the Table columns.
    for header_col in header.cols:
        col = table[header_col.name]
        for attr in ["unit", "description", "format", "meta"]:
            if (val := getattr(header_col, attr)) is not None:
                setattr(col.info, attr, val)

    # Add metadata to the table
    if header.table_meta:
        table.meta.update(header.table_meta)

    # Construct any mixin columns from the raw columns.
    table = serialize._construct_mixins_from_columns(table)

    return table  # noqa: RET504


def write_ecsv(tbl, output, engine="io.ascii", **kwargs):
    """Thin wrapper around the ``io.ascii`` ECSV writer to write ECSV files.

    Parameters
    ----------
    tbl : astropy.table.Table
        The table to write to ECSV format.
    output : str or os.PathLike or file-like object
        The output file path or file-like object to write the ECSV data to.
    engine : str, optional
        The engine to use for writing the CSV data. Default is "io.ascii", which uses
        astropy to write the CSV data. Currently this is the only option.
    **kwargs : dict, optional
        Additional keyword arguments passed to the ECSV writer. These can include
        options like ``delimiter``, ``encoding``, and others supported by the
        `astropy.io.ascii.Ecsv` writer.
    """
    if engine != "io.ascii":
        raise ValueError(
            f"{engine=} is not a supported engine for writing, use 'io.ascii'"
        )
    tbl.write(output, format="ascii.ecsv", **kwargs)


def register_ecsv_table():
    """
    Register ECSV reader and writer with Unified I/O as a Table reader.
    """
    from astropy.io import registry as io_registry
    from astropy.table import Table

    io_registry.register_reader("ecsv", Table, read_ecsv)
    io_registry.register_writer("ecsv", Table, write_ecsv)
