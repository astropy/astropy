READ_DOCSTRING = """
    Read the input ``table`` and return the table.  Most of
    the default behavior for various parameters is determined by the Reader
    class.

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
    Inputter : `~astropy.io.ascii.BaseInputter`
        Inputter class
    Outputter : `~astropy.io.ascii.BaseOutputter`
        Outputter class
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
        Dictionary of converters
    data_Splitter : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split data columns
    header_Splitter : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split header columns
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output.
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fill_values : dict
        specification of fill values for bad or missing table values
    fill_include_names : list
        List of names to include in fill_values.
    fill_exclude_names : list
        List of names to exclude from fill_values (applied after ``fill_include_names``)
    fast_reader : bool or dict
        Whether to use the C engine, can also be a dict with options which
        defaults to `False`; parameters for options dict:

        use_fast_converter: bool
            enable faster but slightly imprecise floating point conversion method
        parallel: bool or int
            multiprocessing conversion using ``cpu_count()`` or ``'number'`` processes
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

    Returns
    -------
    dat : `~astropy.table.Table` OR <generator>
        Output table

    """

WRITE_DOCSTRING = """
    Write the input ``table`` to ``filename``.  Most of the default behavior
    for various parameters is determined by the Writer class.

    See also:

    - https://docs.astropy.org/en/stable/io/ascii/
    - https://docs.astropy.org/en/stable/io/ascii/write.html

    Parameters
    ----------
    table : `~astropy.io.ascii.BaseReader`, array_like, str, file_like, list
        Input table as a Reader object, Numpy struct array, file name,
        file-like object, list of strings, or single newline-separated string.
    output : str, file_like
        Output [filename, file-like object]. Defaults to``sys.stdout``.
    format : str
        Output table format. Defaults to 'basic'.
    delimiter : str
        Column delimiter string
    comment : str
        String defining a comment line in table
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
    fast_writer : bool
        Whether to use the fast Cython writer.
    overwrite : bool
        If ``overwrite=None`` (default) and the file exists, then a
        warning will be issued. In a future release this will instead
        generate an exception. If ``overwrite=False`` and the file
        exists, then an exception is raised.
        This parameter is ignored when the ``output`` arg is not a string
        (e.g., a file object).

    """
