

.. doctest-skip-all

.. _table_io_parquet:

Parquet
-------

.. _Parquet: https://parquet.apache.org/
.. _pyarrow: https://arrow.apache.org/docs/python/

Reading and writing Parquet_ files is supported with ``format='parquet'``
if the pyarrow_ and `pandas <https://pandas.pydata.org/>`__ packages are installed. For writing, the file extensions ``.parquet`` or
``.parq`` will automatically imply the ``'parquet'`` format. For reading,
Parquet files are automatically identified regardless of the extension
if the first four bytes of the file are ``b'PAR1'``.
In many cases you do not need to explicitly specify ``format='parquet'``,
but it may be a good idea anyway if there is any ambiguity about the
file format.

Multiple-file Parquet datasets are not supported for reading and writing.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading from and Writing to Parquet Files

To read a table from a Parquet file named ``observations.parquet``, you can do::

    >>> t = Table.read('observations.parquet')

To write a table to a new file, simply do::

    >>> t.write('new_file.parquet')

As with other formats, the ``overwrite=True`` argument is supported for
overwriting existing files.

One big advantage of the Parquet files is that each column is stored independently,
and thus reading a subset of columns is fast and efficient.  To find out which
columns are stored in a table without reading the data, use the ``schema_only=True``
as shown below. This returns a zero-length table with the appropriate columns::

    >>> schema = Table.read('observations.parquet', schema_only=True)

To read only a subset of the columns, use the ``include_names`` and/or ``exclude_names`` keywords::

    >>> t_sub = Table.read('observations.parquet', include_names=['mjd', 'airmass'])

..
  EXAMPLE END
