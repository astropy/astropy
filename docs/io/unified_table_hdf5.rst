.. doctest-skip-all

.. _table_io_hdf5:

HDF5
----

Reading/writing from/to |HDF5| files is supported with ``format='hdf5'`` (this
requires |h5py| to be installed). However, the ``.hdf5`` file extension is
automatically recognized when writing files, and HDF5 files are automatically
identified (even with a different extension) when reading in (using the first
few bytes of the file to identify the format), so in most cases you will not
need to explicitly specify ``format='hdf5'``.

Since HDF5 files can contain multiple tables, the full path to the table
should be specified via the ``path=`` argument when reading and writing.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading from and Writing to HDF5 Files

To read a table called ``data`` from an HDF5 file named ``observations.hdf5``,
you can do::

    >>> from astropy.table import QTable
    >>> t = QTable.read('observations.hdf5', path='data')

To read a table nested in a group in the HDF5 file, you can do::

    >>> t = QTable.read('observations.hdf5', path='group/data')

To write a table to a new file, the path should also be specified::

    >>> t.write('new_file.hdf5', path='updated_data')

It is also possible to write a table to an existing file using ``append=True``::

    >>> t.write('observations.hdf5', path='updated_data', append=True)

As with other formats, the ``overwrite=True`` argument is supported for
overwriting existing files. To overwrite only a single table within an HDF5
file that has multiple datasets, use *both* the ``overwrite=True`` and
``append=True`` arguments.

Finally, when writing to HDF5 files, the ``compression=`` argument can be
used to ensure that the data is compressed on disk::

    >>> t.write('new_file.hdf5', path='updated_data', compression=True)

..
  EXAMPLE END

Metadata and Mixin Columns
^^^^^^^^^^^^^^^^^^^^^^^^^^

``astropy`` tables can contain metadata, both in the table ``meta`` attribute
(which is an ordered dictionary of arbitrary key/value pairs), and within the
columns, which each have attributes ``unit``, ``format``, ``description``,
and ``meta``.

By default, when writing a table to HDF5 the code will attempt to store each
key/value pair within the table ``meta`` as HDF5 attributes of the table
dataset. This will fail if the values within ``meta`` are not objects that can
be stored as HDF5 attributes. In addition, if the table columns being stored
have defined values for any of the above-listed column attributes, these
metadata will *not* be stored and a warning will be issued.

serialize_meta
~~~~~~~~~~~~~~

To enable storing all table and column metadata to the HDF5 file, call
the ``write()`` method with ``serialize_meta=True``. This will store metadata
in a separate HDF5 dataset, contained in the same file, which is named
``<path>.__table_column_meta__``. Here ``path`` is the argument provided in
the call to ``write()``::

    >>> t.write('observations.hdf5', path='data', serialize_meta=True)

The table metadata are stored as a dataset of strings by serializing the
metadata in YAML following the `ECSV header format
<https://github.com/astropy/astropy-APEs/blob/main/APE6.rst#header-details>`_
definition. Since there are YAML parsers for most common languages, one can
easily access and use the table metadata if reading the HDF5 in a non-astropy
application.

By specifying ``serialize_meta=True`` one can also store
to HDF5 tables that contain :ref:`mixin_columns` such as `~astropy.time.Time` or
`~astropy.coordinates.SkyCoord` columns.

.. note::
    Certain kind of metadata (e.g., numpy object arrays) cannot be serialized correctly
    using YAML.
