.. doctest-skip-all

.. _read_write_tables:

Reading and writing Table objects
===================================

Astropy provides a unified interface for reading and writing data
in different formats.  For many common cases this will
simplify the process of file I/O and reduce the need to master
the separate details of all the I/O packages within Astropy.  For details and
examples of using this interface see the :ref:`table_io`
section.

Getting started
----------------

The :class:`~astropy.table.Table` class includes two methods,
:meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write`, that make it possible to read from
and write to files. A number of formats are automatically supported (see
:ref:`built_in_readers_writers`) and new file formats and extensions can be
registered with the :class:`~astropy.table.Table` class (see
:ref:`io_registry`).

To use this interface, first import the :class:`~astropy.table.Table` class, then
simply call the :class:`~astropy.table.Table`
:meth:`~astropy.table.Table.read` method with the name of the file and
the file format, for instance ``'ascii.daophot'``::

    >>> from astropy.table import Table
    >>> t = Table.read('photometry.dat', format='ascii.daophot')

It is possible to load tables directly from the Internet using URLs. For example,
download tables from Vizier catalogues in CDS format (``'ascii.cds'``)::

    >>> t = Table.read("ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/snrs.dat",
    ...         readme="ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/ReadMe",
    ...         format="ascii.cds")

For certain file formats, the format can be automatically detected, for
example from the filename extension::

    >>> t = Table.read('table.tex')

Similarly, for writing, the format can be explicitly specified::

    >>> t.write(filename, format='latex')

As for the :meth:`~astropy.table.Table.read` method, the format may
be automatically identified in some cases.

Any additional arguments specified will depend on the format.  For examples of this see the
section :ref:`built_in_readers_writers`.  This section also provides the full list of
choices for the ``format`` argument.

Supported formats
------------------

The  :ref:`table_io` has built-in support for the following data file formats:

* :ref:`table_io_ascii`
* :ref:`table_io_hdf5`
* :ref:`table_io_fits`
* :ref:`table_io_votable`
