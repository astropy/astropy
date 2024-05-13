READ_DOCSTRING = """
    Read the input ``table`` and return the table.  Most of the default behavior for
    various parameters is determined by the ``format`` argument.

    Help on the ``read()`` function arguments is available as shown in this example::

      from astropy.io import ascii
      ascii.read.help()  # Common help for all formats
      ascii.read.help("html")  # Common help plus "html" format-specific args

    See also:

    - https://docs.astropy.org/en/stable/io/ascii/
    - https://docs.astropy.org/en/stable/io/ascii/read.html

    Parameters
    ----------
    table : str, file-like, list, `pathlib.Path` object
        Input table as a file name, file-like object, list of string[s],
        single newline-separated string or `pathlib.Path` object.
    guess : bool
        Try to guess the table format. Defaults to None.
    format : str, `~astropy.io.ascii.BaseReader`
        Input table format
    delimiter : str
        Column delimiter string
    comment : str
        Regular expression defining a comment line in table
    quotechar : str
        One-character string to quote fields containing special characters
    header_start : int
        Line index for the header line not counting comment or blank lines.
        A line with only whitespace is considered blank.
    data_start : int
        Line index for the start of data not counting comment or blank lines.
        A line with only whitespace is considered blank.
    data_end : int
        Line index for the end of data not counting comment or blank lines.
        This value can be negative to count from the end.
    converters : dict
        Dictionary of converters to specify output column dtypes. Each key in
        the dictionary is a column name or else a name matching pattern
        including wildcards. The value is either a data type such as ``int`` or
        ``np.float32``; a list of such types which is tried in order until a
        successful conversion is achieved; or a list of converter tuples (see
        the `~astropy.io.ascii.convert_numpy` function for details).
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output.
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fill_values : tuple, list of tuple
        specification of fill values for bad or missing table values
    fill_include_names : list
        List of names to include in fill_values.
    fill_exclude_names : list
        List of names to exclude from fill_values (applied after ``fill_include_names``)
    fast_reader : bool, str or dict
        Whether to use the C engine, can also be a dict with options which
        defaults to `False`; parameters for options dict:

        use_fast_converter: bool
            enable faster but slightly imprecise floating point conversion method
        exponent_style: str
            One-character string defining the exponent or ``'Fortran'`` to auto-detect
            Fortran-style scientific notation like ``'3.14159D+00'`` (``'E'``, ``'D'``, ``'Q'``),
            all case-insensitive; default ``'E'``, all other imply ``use_fast_converter``
        chunk_size : int
            If supplied with a value > 0 then read the table in chunks of
            approximately ``chunk_size`` bytes. Default is reading table in one pass.
        chunk_generator : bool
            If True and ``chunk_size > 0`` then return an iterator that returns a
            table for each chunk.  The default is to return a single stacked table
            for all the chunks.

    encoding : str
        Allow to specify encoding to read the file (default= ``None``).

    Other Parameters
    ----------------
    inputter_cls : `~astropy.io.ascii.BaseInputter`
        Inputter class
    outputter_cls : `~astropy.io.ascii.BaseOutputter`
        Outputter class
    data_splitter_cls : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split data columns
    header_splitter_cls : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split header columns

    Returns
    -------
    dat : `~astropy.table.Table` or <generator>
        Output table

    """

# Specify allowed types for core read() keyword arguments.  Each entry
# corresponds to the name of an argument and either a type (e.g. int) or a
# list of types.  These get used in io.ascii.ui._validate_read_write_kwargs().
# -  The commented-out kwargs are too flexible for a useful check
# -  'list-list' is a special case for an iterable that is not a string.
READ_KWARG_TYPES = {
    # 'table'
    "guess": bool,
    # 'format'
    # 'reader_cls'
    # 'inputter_cls'
    # 'outputter_cls'
    "delimiter": str,
    "comment": str,
    "quotechar": str,
    "header_start": int,
    "data_start": (int, str),  # CDS allows 'guess'
    "data_end": int,
    "converters": dict,
    # 'data_splitter_cls'
    # 'header_splitter_cls'
    "names": "list-like",
    "include_names": "list-like",
    "exclude_names": "list-like",
    "fill_values": "list-like",
    "fill_include_names": "list-like",
    "fill_exclude_names": "list-like",
    "fast_reader": (bool, str, dict),
    "encoding": str,
}


WRITE_DOCSTRING = """
    Write the input ``table`` to ``filename``.  Most of the default behavior
    for various parameters is determined by the Writer class.

    Help on the ``write()`` function arguments is available as shown in this example::

      from astropy.io import ascii
      ascii.write.help()  # Common help for all formats
      ascii.write.help("html")  # Common help plus "html" format-specific args

    See also:

    - https://docs.astropy.org/en/stable/io/ascii/
    - https://docs.astropy.org/en/stable/io/ascii/write.html

    Parameters
    ----------
    table : `~astropy.io.ascii.BaseReader`, array-like, str, file-like, list
        Input table as a Reader object, Numpy struct array, file name,
        file-like object, list of strings, or single newline-separated string.
    output : str, file-like
        Output [filename, file-like object]. Defaults to``sys.stdout``.
    format : str
        Output table format. Defaults to 'basic'.
    delimiter : str
        Column delimiter string
    comment : str, bool
        String defining a comment line in table.  If `False` then comments
        are not written out.
    quotechar : str
        One-character string to quote fields containing special characters
    formats : dict
        Dictionary of format specifiers or formatting functions
    strip_whitespace : bool
        Strip surrounding whitespace from column values.
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output.
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fast_writer : bool, str
        Whether to use the fast Cython writer.  Can be `True` (use fast writer
        if available), `False` (do not use fast writer), or ``'force'`` (use
        fast writer and fail if not available, mostly for testing).
    overwrite : bool
        If ``overwrite=False`` (default) and the file exists, then an OSError
        is raised. This parameter is ignored when the ``output`` arg is not a
        string (e.g., a file object).

    """
# Specify allowed types for core write() keyword arguments.  Each entry
# corresponds to the name of an argument and either a type (e.g. int) or a
# list of types.  These get used in io.ascii.ui._validate_read_write_kwargs().
# -  The commented-out kwargs are too flexible for a useful check
# -  'list-list' is a special case for an iterable that is not a string.
WRITE_KWARG_TYPES = {
    # 'table'
    # 'output'
    "format": str,
    "delimiter": str,
    "comment": (str, bool),
    "quotechar": str,
    "header_start": int,
    "formats": dict,
    "strip_whitespace": (bool),
    "names": "list-like",
    "include_names": "list-like",
    "exclude_names": "list-like",
    "fast_writer": (bool, str),
    "overwrite": (bool),
}
