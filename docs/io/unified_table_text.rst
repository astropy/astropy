.. _unified_table_text:

Text (CSV, fixed-width, HTML, and specialized)
==============================================

The :meth:`~astropy.table.Table.read` and :meth:`~astropy.table.Table.write` methods can
be used to read and write text-based table data in a wide variety of `supported
formats`_. In addition to common formats like `CSV
<https://en.wikipedia.org/wiki/Comma-separated_values>`__ and :ref:`fixed-width
<fixed_width_gallery>`, the unified interface also supports `specialized formats`_ like
`LaTeX <https://en.wikipedia.org/wiki/LaTeX>`_ tables and the AAS :class:`MRT
<astropy.io.ascii.Mrt>` format.

Most of the formats are provided by :ref:`table_io_ascii`, which is a flexible and
powerful interface for reading and writing text tables. In addition, the interface
provides wrappers around select I/O functions in the `pandas`_ library for additional
flexibility and performance.

.. note::

   For reading large CSV files, the astropy :ref:`PyArrow CSV <table_io_pyarrow_csv>`
   reader is a good option to consider since it can be up to 15 times faster than other
   readers.

Supported Formats
-----------------

Character-delimited Formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^
These formats use a character delimiter to separate columns. This is most commonly a
comma (CSV) or a whitespace character like space or tab.

===========================  =====  ======  ============================================================================================
           Format            Write  Suffix                                          Description
===========================  =====  ======  ============================================================================================
                      ascii    Yes          ASCII table in most supported formats (uses guessing)
                ascii.basic    Yes          :class:`~astropy.io.ascii.Basic`: Basic table with custom delimiters
     ascii.commented_header    Yes          :class:`~astropy.io.ascii.CommentedHeader`: Column names in a commented line
                  ascii.csv    Yes    .csv  :class:`~astropy.io.ascii.Csv`: Basic table with comma-separated values
                 ascii.ecsv    Yes   .ecsv  :class:`~astropy.io.ascii.Ecsv`: Basic table with Enhanced CSV (supporting metadata)
            ascii.no_header    Yes          :class:`~astropy.io.ascii.NoHeader`: Basic table with no headers
                  ascii.rdb    Yes    .rdb  :class:`~astropy.io.ascii.Rdb`: Tab-separated with a type definition header line
                  ascii.tab    Yes          :class:`~astropy.io.ascii.Tab`: Basic table with tab-separated values
                 ascii.tdat    Yes   .tdat  :class:`~astropy.io.ascii.Tdat`: Transportable Database Aggregate Table format
                 pandas.csv    Yes          :func:`pandas.read_csv` and :meth:`pandas.DataFrame.to_csv`
                pyarrow.csv     No          :func:`~astropy.io.misc.pyarrow.csv.read_csv`: Performant CSV reader
===========================  =====  ======  ============================================================================================

Fixed-width Formats
^^^^^^^^^^^^^^^^^^^
These formats use fixed-width columns, where each column has a fixed width in characters.
This can be useful for tables that are intended to also be read by humans.

===========================  =====  ======  ============================================================================================
           Format            Write  Suffix                                          Description
===========================  =====  ======  ============================================================================================
          ascii.fixed_width    Yes          :class:`~astropy.io.ascii.FixedWidth`: Fixed width
ascii.fixed_width_no_header    Yes          :class:`~astropy.io.ascii.FixedWidthNoHeader`: Fixed width with no header
 ascii.fixed_width_two_line    Yes          :class:`~astropy.io.ascii.FixedWidthTwoLine`: Fixed width with second header line
                 pandas.fwf     No          :func:`pandas.read_fwf` (fixed width format)
===========================  =====  ======  ============================================================================================

HTML and JSON Formats
^^^^^^^^^^^^^^^^^^^^^
===========================  =====  ======  ============================================================================================
           Format            Write  Suffix                                          Description
===========================  =====  ======  ============================================================================================
                 ascii.html    Yes   .html  :class:`~astropy.io.ascii.HTML`: HTML table
                   jsviewer    Yes          JavaScript viewer format (write-only)
                pandas.html    Yes          :func:`pandas.read_html` and :meth:`pandas.DataFrame.to_html`
                pandas.json    Yes          :func:`pandas.read_json` and :meth:`pandas.DataFrame.to_json`
===========================  =====  ======  ============================================================================================

Specialized Formats
^^^^^^^^^^^^^^^^^^^^
===========================  =====  ======  ============================================================================================
           Format            Write  Suffix                                          Description
===========================  =====  ======  ============================================================================================
               ascii.aastex    Yes          :class:`~astropy.io.ascii.AASTex`: AASTeX deluxetable used for AAS journals
                  ascii.cds     No          :class:`~astropy.io.ascii.Cds`: CDS format table
              ascii.daophot     No          :class:`~astropy.io.ascii.Daophot`: IRAF DAOphot format table
                 ascii.ipac    Yes          :class:`~astropy.io.ascii.Ipac`: IPAC format table
                ascii.latex    Yes    .tex  :class:`~astropy.io.ascii.Latex`: LaTeX table
                  ascii.mrt    Yes          :class:`~astropy.io.ascii.Mrt`: AAS Machine-Readable Table format
                  ascii.qdp    Yes    .qdp  :class:`~astropy.io.ascii.QDP`: Quick and Dandy Plotter files
                  ascii.rst    Yes    .rst  :class:`~astropy.io.ascii.RST`: reStructuredText simple format table
           ascii.sextractor     No          :class:`~astropy.io.ascii.SExtractor`: SExtractor format table
===========================  =====  ======  ============================================================================================

.. _table_io_ascii:

`astropy.io.ascii`
------------------
The :ref:`astropy.io.ascii <io-ascii>` sub-package provides read and write support for
:ref:`many different formats <supported_formats>`, including astronomy-specific formats
like AAS `Machine-Readable Tables (MRT) <https://journals.aas.org/mrt-standards/>`_.

We **strongly recommend** using the unified interface for reading and writing tables via
the :ref:`astropy.io.ascii <io-ascii>` sub-package. This is done by prefixing the
:ref:`format name <supported_formats>` with the ``ascii.`` prefix. For example to read a
DAOphot table use:

.. doctest-skip::

    >>> from astropy.table import Table
    >>> t = Table.read('photometry.dat', format='ascii.daophot')

Use ``format='ascii'`` in order read a table and guess the table format by successively
trying most of the available formats in a specific order. This can be slow and is not
recommended for large tables.

.. doctest-skip::

  >>> t = Table.read('astropy/io/ascii/tests/t/latex1.tex', format='ascii')
  >>> print(t)
  cola colb colc
  ---- ---- ----
     a    1    2
     b    3    4

When writing a table with ``format='ascii'`` the output is a basic
space-delimited file with a single header line containing the
column names.

All additional arguments are passed to the `astropy.io.ascii`
:func:`~astropy.io.ascii.read` and :func:`~astropy.io.ascii.write`
functions. Further details are available in the sections on
:ref:`io_ascii_read_parameters` and :ref:`io_ascii_write_parameters`. For
example, to change the column delimiter and the output format for the ``colc``
column use:

.. doctest-skip::

  >>> t.write(sys.stdout, format='ascii', delimiter='|', formats={'colc': '%0.2f'})
  cola|colb|colc
  a|1|2.00
  b|3|4.00

.. attention:: **ECSV is recommended**

   For writing and reading tables to text in a way that fully reproduces the table data,
   types, and metadata (i.e., the table will "round-trip"), we highly recommend using
   the :ref:`ecsv_format` with ``format="ascii.ecsv"``. This writes the actual data in a
   space- or comma-delimited format that most text table readers can parse, but also
   includes metadata encoded in a comment block that allows full reconstruction of the
   original columns. This includes support for :ref:`ecsv_format_mixin_columns` (such as
   `~astropy.coordinates.SkyCoord` or `~astropy.time.Time`) and
   :ref:`ecsv_format_masked_columns`.

..
  EXAMPLE END

.. _table_io_pandas:

Pandas
------

.. _pandas: https://pandas.pydata.org/pandas-docs/stable/index.html

``astropy`` `~astropy.table.Table` supports the ability to read or write tables
using some of the `I/O methods <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html>`_
available within pandas_. This interface thus provides convenient wrappers to
the following functions / methods:

.. csv-table::
    :header: "Format name", "Data Description", "Reader", "Writer"
    :widths: 25, 25, 25, 25

    ``pandas.csv``,`CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`__,`read_csv() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-read-csv-table>`_,`to_csv() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-store-in-csv>`_
    ``pandas.json``,`JSON <http://www.json.org/>`__,`read_json() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-json-reader>`_,`to_json() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-json-writer>`_
    ``pandas.html``,`HTML <https://en.wikipedia.org/wiki/HTML>`__,`read_html() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-read-html>`_,`to_html() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-html>`_
    ``pandas.fwf``,Fixed Width,`read_fwf() <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_fwf.html#pandas.read_fwf>`_,

**Notes**:

- This is subject to the limitations discussed in :ref:`astropy-table-and-dataframes`.
- There is no fixed-width writer in pandas_.
- Reading HTML requires `BeautifulSoup4 <https://pypi.org/project/beautifulsoup4/>`_ and
  `html5lib <https://pypi.org/project/html5lib/>`_ to be installed.

When reading or writing a table, any keyword arguments apart from the
``format`` and file name are passed through to pandas, for instance:

.. doctest-skip::

  >>> t.write('data.csv', format='pandas.csv', sep=' ', header=False)
  >>> t2 = Table.read('data.csv', format='pandas.csv', sep=' ', names=['a', 'b', 'c'])

.. _table_io_pyarrow_csv:

PyArrow CSV
-----------

.. _pyarrow: https://arrow.apache.org/docs/python/

The `pyarrow`_ library provides a highly-performant CSV reader that can be used in
Astropy with ``Table.read(input_file, format="pyarrow.csv", ...)``. This can by up to 15
times faster and more memory-efficient than the :ref:`astropy.io.ascii <io-ascii>` fast
reader or the default ``pandas.csv`` reader. The best performance is achieved for files
with only numeric data types, but even for files with mixed data types, the performance
is still better than the standard :ref:`astropy.io.ascii <io-ascii>` fast CSV reader.

This reader uses the :func:`~astropy.io.misc.pyarrow.csv.read_csv` function, which in
turn uses the `PyArrow CSV reader <https://arrow.apache.org/docs/python/csv.html>`__ and
sets the various options to ``pyarrow.csv.read_csv()`` appropriately. The interface is
designed to be similar to the :ref:`io.ascii read interface <io_ascii_read_parameters>`
where possible, but there are differences, most notably:

- Input can only be a string file name, `pathlib.Path`, or a binary file-like object.
- Whitespace in string data fields and header column names is preserved.
- Use ``dtypes`` instead of ``converters`` to specify the column data types.
- Use ``null_values`` instead of ``fill_values`` to specify the null (missing) values.
- No ``guess`` parameter and no guessing of the table format (e.g., ``delimiter``).
- No ``data_end`` parameter.
- No ``exclude_names`` parameter.
- Columns consisting of only string values ``True`` and ``False`` are parsed as
  boolean data.
- Columns with ISO 8601 date/time strings are parsed as shown below:
  - ``12:13:14.123456``: ``object[datetime.time]``
  - ``2025-01-01``: ``np.datetime64[D]``
  - ``2025-01-01T01:02:03``: ``np.datetime64[s]``
  - ``2025-01-01T01:02:03.123456``: ``np.datetime64[ns]``
- Timestamp parsing behavior can be customized with the ``timestamp_parsers``
  parameter.

Using the PyArrow CSV reader directly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :mod:`astropy.io.misc.pyarrow.csv` module also provides the
:func:`~astropy.io.misc.pyarrow.csv.convert_pa_table_to_astropy_table` function to
allow converting a ``pyarrow.Table`` to an `astropy.table.Table`. This allows using
the `PyArrow CSV reader <https://arrow.apache.org/docs/python/csv.html>`__ directly
with custom options that are not available in the astropy interface.
